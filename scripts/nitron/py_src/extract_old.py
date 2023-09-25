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

def main(args):
    stats_fname = args.output+".nitron.extract.stats"
    if os.path.exists(stats_fname):
        os.remove(stats_fname)
    stats_fp = open(stats_fname,"w+")
    stats_fp.write("assembly: "+args.assembly+"\n")
    stats_fp.write("labels: "+args.labels+"\n")
    stats_fp.write("minimum number of samples: "+str(args.min_samples)+"\n")


    labels_df = pd.read_csv(args.labels,sep="\t")
    labels_df.rename({"Name":"tid"},axis=1,inplace=True)

    all_tids = dict()

    with open(args.assembly,"r") as inFP:
        for line in inFP.readlines():
            if line[0]=="#":
                continue
            lineCols = line.strip("\n").split("\t")
            if not lineCols[2]=="transcript":
                continue

            tid = lineCols[8].split("transcript_id \"",1)[1].split("\"",1)[0]
            tpms = lineCols[8].split("tpms \"",1)[1].split("\"",1)[0]
            tpms = [float(x) for x in tpms.split(",")]
            tpm_mean = sum(tpms)/len(tpms)

            all_tids.setdefault(tid,dict())
            all_tids[tid]["tpm_mean"] = tpm_mean
            all_tids[tid]["sample_count"] = len(tpms)

    res_df = pd.DataFrame(all_tids).T.reset_index()
    res_df.replace(np.nan,0,inplace=True)
    cur_cols = list(res_df.columns)
    cur_cols[0] = "tid"
    res_df.columns = cur_cols

    tgtf = pd.read_csv(args.assembly,sep="\t",comment="#",names=commons.gff3cols)
    tgtf["tid"] = tgtf["attributes"].str.split("transcript_id \"",expand=True,n=1)[1].str.split("\"",expand=True,n=1)[0]

    tlabels = labels_df[labels_df["tid"].isin(set(tgtf["tid"]))].reset_index(drop=True)

    tlabels = tlabels.merge(res_df,on="tid",how="outer",indicator=True)
    assert len(set(tlabels["_merge"]))==1 and list(set(tlabels["_merge"]))[0]=="both","tids are not the same"
    tlabels.drop("_merge",axis=1,inplace=True)

    # now select only those labels which we are interested in
    # Those are:
    # repeat
    # readthrough
    # tpm_mean<1
    # sample_count<min([10,len(tissue2samples[tissue])])
    known = tlabels[((tlabels["Noise_comment"].str.contains("refseq"))| \
                     (tlabels["Noise_comment"].str.contains("gencode")))].reset_index(drop=True)
    known_pass = known[((known["tpm_mean"]>=1)& \
                        (known["sample_count"]>=args.min_samples))].reset_index(drop=True)

    novel_pass = tlabels[~(tlabels["Noise_comment"].str.contains("refseq"))&
                         ~(tlabels["Noise_comment"].str.contains("gencode"))].reset_index(drop=True)
    novel_pass = novel_pass[~((novel_pass["Noise_comment"].str.contains("repeat"))| \
                              (novel_pass["Noise_comment"].str.contains("readthrough"))| \
                              (novel_pass["tpm_mean"]<1)| \
                              (novel_pass["sample_count"]<args.min_samples))].reset_index(drop=True)
    all_df = pd.concat([known,novel_pass],axis=0).reset_index(drop=True)
    stats_fp.write("total number of assembled transcripts: "+str(len(tlabels))+"\n")
    stats_fp.write("number of known assembled transcripts: "+str(len(known))+"\n")
    stats_fp.write("number of known passing assembled transcripts: "+str(len(known_pass))+"\n")
    stats_fp.write("number of novel assembled transcripts: "+str(len(novel_pass))+"\n")
    stats_fp.write("number of passing assembled transcripts: "+str(len(all_df))+"\n")

    # get subset of the tissue gtf based on the labels
    known_gtf = tgtf[tgtf["tid"].isin(set(known["tid"]))].reset_index(drop=True)
    known_gtf[commons.gff3cols].to_csv(args.output+".known.all.gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)
    known.to_csv(args.output+".known.all.labels",index=False)

    known_pass_gtf = tgtf[tgtf["tid"].isin(set(known_pass["tid"]))].reset_index(drop=True)
    known_pass_gtf[commons.gff3cols].to_csv(args.output+".known.pass.gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)
    known_pass.to_csv(args.output+".known.pass.labels",index=False)

    novel_pass_gtf = tgtf[tgtf["tid"].isin(set(novel_pass["tid"]))].reset_index(drop=True)
    novel_pass_gtf[commons.gff3cols].to_csv(args.output+".novel.pass.gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)
    novel_pass.to_csv(args.output+".novel.pass.labels",index=False)

    tgtf = tgtf[tgtf["tid"].isin(set(all_df["tid"]))].reset_index(drop=True)
    tgtf[commons.gff3cols].to_csv(args.output+".all.pass.gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)
    all_df.to_csv(args.output+".all.pass.labels",index=False)

    # run selection
    ni_cmd = [args.novel_introns,
              "-g",args.output+".all.pass.gtf",
              "-o",args.output+".selected",
              "-i",args.predictions]
    subprocess.call(ni_cmd)