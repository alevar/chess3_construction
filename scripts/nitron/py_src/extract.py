#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
	Author: Ales Varabyou
"""

# Extract transcripts based on the predicted introns

import pandas as pd
import numpy as np
import subprocess
import argparse
import random
import pickle
import csv
import sys
import os

from py_src import commons

def separate_known_novel(args,stats_fp):
    # load precomputed known
    known_tids = set()
    with open(args.all_known_assembled_labeled,"r") as inFP:
        for line in inFP.readlines():
            lcs = line.split("\t")
            if len(lcs)<9:
                continue

            if lcs[2]=="transcript":
                tid = lcs[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                known_tids.add(tid)

    # get subset of labels for each tissue
    labels_df = pd.read_csv(args.labels,sep="\t")
    labels_df.rename({"Name":"tid"},axis=1,inplace=True)

    # compare old with new to get mappings
    out_fname = args.output+".tb2merge"

    gffcmp_cmd = "gffcompare -T -r "+args.assembly+" "+args.tb_filtered_tissue_gtf+" -o "+out_fname
    subprocess.call(gffcmp_cmd,shell=True)

    # load mappings
    new2old = dict()
    commons.load_tracking_map(out_fname+".tracking",True,2,new2old)

    # get all ALL_tids with class_code "="
    shared_tids = set(new2old.values())


    # replace tids and reformat attributes in the tissue gtfs

    # need to get the sample count from the original data
    tb_tpms = dict()
    with open(args.tb_filtered_tissue_gtf,"r") as inFP:
        for line in inFP.readlines():
            if line[0]=="#":
                continue
            lineCols = line.strip("\n").split("\t")
            if not lineCols[2]=="transcript":
                continue

            tid = lineCols[8].split("transcript_id \"",1)[1].split("\"",1)[0]
            tpm = float(lineCols[8].split("TPM \"",1)[1].split("\"",1)[0])
            if not tid in new2old: # not "=" code
                continue
            new_tid = new2old[tid]
            if new_tid in tb_tpms:
                tpm = max([tpm,tb_tpms[new_tid]])
            tb_tpms[new_tid] = tpm

    with open(args.output+".tb.gtf","w+") as outFP:
        with open(args.assembly,"r") as inFP:
            for line in inFP:
                if line[0]=="#":
                    outFP.write(line)
                    continue
                cols = line.strip().split("\t")
                new_attrs = ""
                tid = ""
                skip = False # whether to skip the line
                for attr in cols[-1].rstrip(";").split(";"):
                    k,v = [x.strip() for x in attr.rstrip("\"").split(" \"")]
                    if k=="transcript_id":
                        if not v in shared_tids: # not a class_code "="
                            skip=True
                            break
                        tid = v
                    if k=="tpms": # replace for compatibility with "novel_introns"
                        k = "tpms"
                    new_attrs+=k+" \""+v+"\"; "
                if not skip:
                    assert tid in tb_tpms,"tid not found: "+tid
                    if cols[2]=="transcript":
                        new_attrs+="tb_tpms \""+str(tb_tpms[tid])+"\";"
                    cols[-1] = new_attrs.rstrip(" ")
                    outFP.write("\t".join(cols)+"\n")

    # select by num_samples, etc
    all_tids = dict()
    with open(args.output+".tb.gtf","r") as inFP:
        for line in inFP.readlines():
            if line[0]=="#":
                continue
            lineCols = line.strip("\n").split("\t")
            if not lineCols[2]=="transcript":
                continue

            tid = lineCols[8].split("transcript_id \"",1)[1].split("\"",1)[0]
            sample_tpms = lineCols[8].split("; tpms \"",1)[1].split("\"",1)[0]
            sample_tpms = [float(x) for x in sample_tpms.split(",")]
            sample_tpms_mean = sum(sample_tpms)/len(sample_tpms)

            tb_tpm = float(lineCols[8].split("; tb_tpms \"",1)[1].split("\"",1)[0])

            all_tids.setdefault(tid,dict())
            all_tids[tid]["sample_tpms_mean"] = sample_tpms_mean
            all_tids[tid]["sample_count"] = len(sample_tpms)
            all_tids[tid]["tb_tpm"] = tb_tpm

    res_df = pd.DataFrame(all_tids).T.reset_index()
    res_df.replace(np.nan,0,inplace=True)
    cur_cols = list(res_df.columns)
    cur_cols[0] = "tid"
    res_df.columns = cur_cols

    tgtf = pd.read_csv(args.output+".tb.gtf",sep="\t",comment="#",names=commons.gff3cols)
    tgtf["tid"] = tgtf["attributes"].str.split("transcript_id \"",expand=True,n=1)[1].str.split("\"",expand=True,n=1)[0]

    tlabels = labels_df[labels_df["tid"].isin(set(tgtf["tid"]))].reset_index(drop=True)

    tlabels = tlabels.merge(res_df,on="tid",how="outer",indicator=True)
    assert len(set(tlabels["_merge"]))==1 and list(set(tlabels["_merge"]))[0]=="both","tids are not the same"
    tlabels.drop("_merge",axis=1,inplace=True)

    # get all tids for novel_pass (not in known refseq and gencode)
    novel_pass = tlabels[~(tlabels["tid"].isin(known_tids))].reset_index(drop=True)
    novel_pass[~((novel_pass["Noise_comment"].str.contains("repeat"))| \
                 (novel_pass["Noise_comment"].str.contains("readthrough")))].reset_index(drop=True)

    novel_pass_gtf = tgtf[tgtf["tid"].isin(set(novel_pass["tid"]))].reset_index(drop=True)
    novel_pass_gtf[commons.gff3cols].to_csv(args.output+".ordered.novel_pass_labels.gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)
    novel_pass.to_csv(args.output+".novel_pass.labels",index=False)

    # get all tids for known (pass or fail - doesn't matter)
    known_gtf = pd.read_csv(args.assembly,sep="\t",comment="#",names=commons.gff3cols)
    known_gtf["tid"] = known_gtf["attributes"].str.split("transcript_id \"",expand=True,n=1)[1].str.split("\"",expand=True,n=1)[0]
    known = tlabels[tlabels["tid"].isin(known_tids)].reset_index(drop=True)
    known_gtf = known_gtf[known_gtf["tid"].isin(known_tids)].reset_index(drop=True)
    known_gtf[commons.gff3cols].to_csv(args.output+".ordered.known.gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)
    known.to_csv(args.output+".known.labels",index=False)

def run_novel_introns(args,stats_fp):
    # run selection
    out_min_fname = args.output+".selected_min"
    out_allg_fname = args.output+".selected_allg"
    out_all_fname = args.output+".selected_all"

    ni_cmd = [args.novel_introns,
              "-g",args.output+".ordered.novel_pass_labels.gtf",
              "-o",out_min_fname,
              "-i",args.predictions,
              "-t","min"]
    subprocess.call(ni_cmd)

    ni_cmd = [args.novel_introns,
              "-g",args.output+".ordered.novel_pass_labels.gtf",
              "-o",out_allg_fname,
              "-i",args.predictions,
              "-t","allg"]
    subprocess.call(ni_cmd)

    ni_cmd = [args.novel_introns,
              "-g",args.output+".ordered.novel_pass_labels.gtf",
              "-o",out_all_fname,
              "-i",args.predictions,
              "-t","all"]
    subprocess.call(ni_cmd)

    return 0

# load all selected known passing transcripts and their tpms
def load_tid_stats(fname,res):
    with open(fname,"r") as inFP:
        for line in inFP:
            if line[0]=="#":
                continue
            lineCols = line.strip("\n").split("\t")
            if lineCols[2]=="transcript":
                tid = lineCols[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                tpm = np.mean([float(x ) for x in lineCols[8].split("; tpms \"",1)[1].split("\"",1)[0].split(",")])
                num_samples = len(lineCols[8].split("tpms \"",1)[1].split("\"",1)[0].split(","))
                res[tid] = [tpm,num_samples,0] # tpm,num_samples,num_exons
            if lineCols[2]=="exon":
                res[tid][2]+=1

def add_single_exon(args,stats_fp):
    known_tids = dict()
    known_fname = args.output+".ordered.known.gtf"
    load_tid_stats(known_fname,known_tids)

    novel_pass_tids = dict()
    novel_pass_fname = args.output+".ordered.novel_pass_labels.gtf"
    load_tid_stats(novel_pass_fname,novel_pass_tids)

    # find minimum for tpm
    tpm_df = pd.DataFrame.from_dict(known_tids,orient='index').reset_index()
    tpm_df.columns=["tid","tpm","num_samples","num_exons"]
    me_tpm_df = tpm_df[tpm_df["num_exons"]>1].reset_index(drop=True)
    se_tpm_df = tpm_df[tpm_df["num_exons"]==1].reset_index(drop=True)

    me_q25,me_q50,me_q75 = me_tpm_df['tpm'].quantile([0.25,0.5,0.75])
    me_iqr = me_q75-me_q25
    me_high_whisker = min([me_q75+1.5*me_iqr,max(me_tpm_df['tpm'])])
    stats_fp.write("minimum TPM for single exon transcripts based on multi-exon TPM IQR*1.5 is: "+str(me_high_whisker)+"\n")

    se_q25,se_q50,se_q75 = se_tpm_df['tpm'].quantile([0.25,0.5,0.75])
    se_iqr = se_q75-se_q25
    se_high_whisker = min([se_q75+1.5*se_iqr,max(se_tpm_df['tpm'])])
    stats_fp.write("minimum TPM for single exon transcripts based on single-exon TPM IQR*1.5 is: "+str(se_high_whisker)+"\n")


    novel_se_tpm_df = pd.DataFrame.from_dict(novel_pass_tids,orient='index').reset_index()
    novel_se_tpm_df.columns=["tid","tpm","num_samples","num_exons"]
    novel_se_tpm_df = novel_se_tpm_df[novel_se_tpm_df["num_exons"]==1].reset_index(drop=True)

    passing_novel_se_tids = novel_se_tpm_df[(novel_se_tpm_df["num_samples"]>args.min_samples)& \
                                            (novel_se_tpm_df["tpm"]>=me_high_whisker)]["tid"].tolist()
    stats_fp.write("number of novel single_exon with treshold based on multi-exon TPM IQR*1.5: "+str(len(passing_novel_se_tids))+"\n")

    tmp = novel_se_tpm_df[(novel_se_tpm_df["num_samples"]>args.min_samples)& \
                          (novel_se_tpm_df["tpm"]>=se_high_whisker)]["tid"].tolist()
    stats_fp.write("number of novel single_exon with treshold based on multi-exon TPM IQR*1.5: "+str(len(tmp))+"\n")

    tissue_single_exon_df = pd.read_csv(novel_pass_fname,sep="\t",comment="#",names=commons.gff3cols)
    tissue_single_exon_df["tid"] = tissue_single_exon_df["attributes"].str.split("transcript_id \"",expand=True,n=1)[1].str.split("\"",expand=True,n=1)[0]
    tissue_single_exon_df = tissue_single_exon_df[tissue_single_exon_df["tid"].isin(passing_novel_se_tids)].reset_index(drop=True)
    tissue_single_exon_df[commons.gff3cols].to_csv(args.output+".novel_single_exon.pass.gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)

# generate gtf for novel transcripts based on novel_introns1 output
def novel_introns_to_gtf(ni_fname,gtf_fname,out_fname):
    selected_tids = set()
    with open(ni_fname,"r") as inFP:
        for line in inFP:
            tid = line.strip()
            selected_tids.add(tid)

    found_count = 0
    with open(out_fname,"w+") as outFP:
        with open(gtf_fname,"r") as inFP:
            for line in inFP:
                if line[0]=="#":
                    continue
                lineCols = line.strip("\n").split("\t")
                tid = lineCols[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                if tid in selected_tids:
                    outFP.write(line)

def main(args):
    stats_fname = args.output+".nitron.extract.stats"
    if os.path.exists(stats_fname):
        os.remove(stats_fname)
    stats_fp = open(stats_fname,"w+")
    stats_fp.write("assembly: "+args.assembly+"\n")
    stats_fp.write("labels: "+args.labels+"\n")
    stats_fp.write("minimum number of samples: "+str(args.min_samples)+"\n")

    separate_known_novel(args,stats_fp)
    run_novel_introns(args,stats_fp)
    add_single_exon(args,stats_fp)
    for suff in ["_min","_all","_allg"]:
        ni_fname = args.output+".selected"+suff
        novel_fname = args.output+".novel"+suff+".gtf"
        ordered_fname = args.output+".ordered.novel_pass_labels.gtf"
        novel_introns_to_gtf(ni_fname,ordered_fname,novel_fname)


    stats_fp.close()


