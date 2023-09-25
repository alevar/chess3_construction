#!/usr/bin/env python

#====================================
# preprocess and organize the data
#====================================

import os
import csv
import sys
import time
import glob
import shutil
import argparse
import subprocess
import numpy as np
import pandas as pd

gff3cols = ["seqid","source","type","start","end","score","strand","phase","attributes"]

def run_get_intron_transcript_counts(args,t2s):
    outdir = args.outdir.rstrip("/")+"/"

    for tissue,paths in t2s.items():
        print("run_get_intron_transcript_counts: "+tissue)
        tissue_dir=outdir+tissue+"/"
        gtf_fname = outdir+"ALL."+tissue+".gtf"

        full_gtf = pd.read_csv(gtf_fname,sep="\t",comment="#",names=gff3cols)
        full_gtf["tid"] = full_gtf["attributes"].str.split("transcript_id \"",expand=True,n=1)[1].str.split("\"",expand=True,n=1)[0]
        full_gtf["gid"] = full_gtf["attributes"].str.split("gene_id \"",expand=True,n=1)[1].str.split("\"",expand=True,n=1)[0]
        full_gtf["type"]=pd.Categorical(full_gtf["type"],categories=["transcript","exon"],ordered=True)
        full_gtf.sort_values(by=["gid","tid","type","start"],ascending=True,inplace=True)
        full_gtf[gff3cols].to_csv(tissue_dir+tissue+".ordered.gtf",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)
        
        full_gtf = full_gtf[full_gtf["type"]=="exon"].reset_index(drop=True)
        full_gtf.sort_values(by=["tid","seqid","strand","start","end"],ascending=True,inplace=True)
        full_gtf.reset_index(drop=True,inplace=True)
        # get introns
        full_gtf[["next_start","next_tid"]] = full_gtf[["start","tid"]].shift(-1)
        full_gtf = full_gtf[full_gtf["tid"]==full_gtf["next_tid"]].reset_index(drop=True)
        full_gtf.dropna(inplace=True)
        full_gtf["start"]=full_gtf["end"]
        full_gtf["end"]=full_gtf["next_start"].astype(int)
        full_gtf = full_gtf[["seqid","strand","start","end"]]
        full_gtf["tx_count"] = 1
        full_gtf = full_gtf.groupby(by=["seqid","strand","start","end"]).count().reset_index()
        full_gtf.to_csv(tissue_dir+tissue+".assembled.introns",index=False,header=False)

        os.replace(tissue_dir+tissue+".ordered.gtf",gtf_fname)

def run_tissue_tiebrush(args,t2s):
    outdir = args.outdir.rstrip("/")+"/"
    tb_fname = outdir+"tb.txt"
    with open(tb_fname,"w+") as tbFP:
        for tissue,paths in t2s.items():
            tissue_dir=outdir+tissue+"/"
            if not os.path.exists(tissue_dir):
                os.makedirs(tissue_dir)
            
            lst_fname = outdir+tissue+".lst"
            with open(lst_fname,"w+") as lstFP:
                lstFP.write("\n".join(paths))
            
            tbFP.write(args.tiewrap+" --batch-size "+str(args.num_samples_per_call)+" --output "+tissue_dir+tissue+".tb.bam"+" "+" "+lst_fname+"\n")
            
            
    # now can run the file in parallel
    parallel_cmd = "parallel -j "+str(args.threads)+" < "+tb_fname
    subprocess.call(parallel_cmd,shell=True)
    
    print("indexing tissue merges")
    for tissue,paths in t2s.items():
        tissue_dir=outdir+tissue+"/"
        idx_cmd = ["samtools","index",tissue_dir+tissue+".tb.bam"]
        subprocess.call(idx_cmd)    

    print("merging all tmps")

    # lastly run the final tiebrush to merge them all together
    tiewrap_cmd = [args.tiewrap,
                   "--batch-size",str(args.num_samples_per_call),
                   "--output",outdir+"all.tb.bam"]
    for tissue,paths in t2s.items():
        tissue_dir=outdir+tissue+"/"
        tiewrap_cmd.append(tissue_dir+tissue+".tb.bam")
    print(" ".join(tiewrap_cmd))
    subprocess.call(tiewrap_cmd)

    idx_cmd = ["samtools","index",outdir+"all.tb.bam"]

def run_tissue_tiecov(args,t2s):
    outdir = args.outdir.rstrip("/")+"/"
    with open(outdir+"tiecov.parallel","w+") as outFP:
        for tissue,paths in t2s.items():
            tissue_dir = outdir+tissue+"/"
            outFP.write(args.tiecov+" -l 75 -c "+tissue_dir+tissue+".cov.bed -j "+tissue_dir+tissue+".def.sjs -s "+tissue_dir+tissue+".sjs.stats "+tissue_dir+tissue+".tb.bam\n")

    # now can run the file in parallel
    parallel_cmd = "parallel -j "+str(args.threads)+" < "+outdir+"tiecov.parallel"
    subprocess.call(parallel_cmd,shell=True)

def run_tissue_tiecov_default(args,t2s):
    outdir = args.outdir.rstrip("/")+"/"
    with open(outdir+"tiecov_default.parallel","w+") as outFP:
        for tissue,paths in t2s.items():
            print("run_tissue_tiecov_default: "+tissue)
            tissue_dir = outdir+tissue+"/"
            outFP.write(args.tiecov_default+" -s "+tissue_dir+tissue+".def.sample -j "+tissue_dir+tissue+".def.junctions -c "+tissue_dir+tissue+".def.coverage "+tissue_dir+tissue+".tb.bam\n")

    # now can run the file in parallel
    parallel_cmd = "parallel -j "+str(args.threads)+" < "+outdir+"tiecov_default.parallel"
    subprocess.call(parallel_cmd,shell=True)

def run_convert_to_bigwig(args,t2s):
    outdir = args.outdir.rstrip("/")+"/"
    with open(outdir+"bigwig.parallel","w+") as outFP:
        for tissue,paths in t2s.items():
            print("run_convert_to_bigwig: "+tissue)
            tissue_dir = outdir+tissue+"/"
            outFP.write("bedtools sort -i "+tissue_dir+tissue+".def.coverage.bedgraph > "+tissue_dir+tissue+".def.coverage.sorted.bedgraph && bedGraphToBigWig "+tissue_dir+tissue+".def.coverage.sorted.bedgraph ~/genomes/human/hg38/hg38_p12_ucsc.no_alts.no_fixs.fa.fai "+tissue_dir+tissue+".def.coverage.bigwig\n")

    # now can run the file in parallel
    parallel_cmd = "parallel -j "+str(args.threads)+" < "+outdir+"bigwig.parallel"
    
    subprocess.call(parallel_cmd,shell=True)

def run_tiebrush_tiecov_bigwig_all(args,t2s):
    print("run_tiebrush_tiecov_bigwig_all")
    outdir = args.outdir.rstrip("/")+"/"

    tiecov_cmd = [args.tiecov_default,
                "-s",outdir+"all.def.sample",
                "-j",outdir+"all.def.junctions",
                "-c",outdir+"all.def.coverage",outdir+"all.tb.bam"]
    print(" ".join(tiecov_cmd))
    subprocess.call(tiecov_cmd)

    sort_cmd =  "bedtools sort -i "+outdir+"all.def.coverage.bedgraph > "+outdir+"all.def.coverage.sorted.bedgraph"
    print(sort_cmd)
    subprocess.call(sort_cmd,shell=True)

    bw_cmd =  ["bedGraphToBigWig",
                 outdir+"all.def.coverage.sorted.bedgraph",args.reference,outdir+"all.def.coverage.bigwig"]
    print(" ".join(bw_cmd))
    subprocess.call(bw_cmd)


def run_assembly_stats(args,t2s):
    print("run_assembly_stats")
    outdir = args.outdir.rstrip("/")+"/"

    # run assembly stats to get additional data
    stats_default_cmd = [args.assembly_stats,
                       "-g",args.ALL_gtf,
                       "-o",outdir+"ALL",
                       "-a",args.ALL_tracking,
                       "-t",args.tissue_trackings,
                       "-c",args.tissue_gtfs]
    print(" ".join(stats_default_cmd))
    subprocess.call(stats_default_cmd)


def run_step1(args):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    assert os.path.exists(args.input),"input file does not exist"

    # first need to form a dictionary of tissues to samples
    tissue2samples = dict()
    with open(args.input,"r") as inFP:
        for line in inFP.readlines():
            line = line.strip()
            tissue,cram_fp = line.split(",")
            tissue2samples.setdefault(tissue,[]).append(cram_fp)

#     run_tissue_tiebrush(args,tissue2samples)
#     run_tissue_tiecov(args,tissue2samples)
    run_assembly_stats(args,tissue2samples)
    run_get_intron_transcript_counts(args,tissue2samples)
    run_tissue_tiecov_default(args,tissue2samples)
    run_convert_to_bigwig(args,tissue2samples)
    run_tiebrush_tiecov_bigwig_all(args,tissue2samples)

def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument('--input',
                        required=True,
                        type=str,
                        help="Input file in CSV format where column #1 is the path to a cram file and column #2 is the tissue name")
    parser.add_argument('--outdir',
                        required=True,
                        type=str,
                        help="Output directory in which all output and temporary data will be stored")
    parser.add_argument("--threads",
                        required=False,
                        type=int,
                        default=1,
                        help="number of threads to be used by GNU parallel")
    parser.add_argument("--num_samples_per_call",
                        required=False,
                        type=int,
                        default=20,
                        help="number of samples to process with tiebrush within a single batch")
    parser.add_argument("--ALL_tracking",
                        required=True,
                        type=str,
                        help="path to the tracking file for the ALL level of assembly as generated by gffcompare")
    parser.add_argument("--ALL_gtf",
                        required=True,
                        type=str,
                        help="path to the GTF file for the ALL level of assembly as generated by gffcompare")
    parser.add_argument("--tissue_trackings",
                        required=True,
                        type=str,
                        help="path to a file containing a list of tissue tracking files")
    parser.add_argument("--tissue_gtfs",
                        required=True,
                        type=str,
                        help="path to a TSV file containing paths to the gtf names (1st column) for each tissue and names of tissues (2nd column)")
    parser.add_argument("--tiebrush",
                        required=False,
                        type=str,
                        default="tiebrush",
                        help="path to the tiebrush executable")
    parser.add_argument("--tiewrap",
                        required=False,
                        type=str,
                        default="tiewrap.py",
                        help="path to the tiewrap executable")
    parser.add_argument("--tiecov",
                        required=True,
                        type=str,
                        help="path to the tiecov_sjs executable")
    parser.add_argument("--tiecov_default",
                        required=False,
                        type=str,
                        default="tiecov",
                        help="path to the standard tiecov executable")
    parser.add_argument("--assembly_stats",
                        required=True,
                        type=str,
                        help="path to the assembly_stats executable")
    parser.add_argument("--reference",
                        required=True,
                        type=str,
                        help="path to the reference genome")
    

    parser.set_defaults(func=run_step1)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main(sys.argv[1:])
