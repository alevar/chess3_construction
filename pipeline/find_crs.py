#!/usr/bin/env python3

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
from pybedtools import BedTool

gff3cols = ["seqid","source","type","start","end","score","strand","phase","attributes"]

def run(args):
    stats_fname = args.output+".find_crs.stats"
    if os.path.exists(stats_fname):
        os.remove(stats_fname)
    stats_fp = open(stats_fname,"w+")


    crfs_df = pd.read_csv(args.crs)
    chess_df = pd.read_csv(args.input,sep="\t",comment="#",names=gff3cols)
    chess_df["tid"] = chess_df["attributes"].str.split("transcript_id \"",expand=True,n=1)[1].str.split("\"",expand=True,n=1)[0]
    chess_bed = chess_df[chess_df["type"] == "exon"][["seqid","start","end","tid","score","strand"]].reset_index(drop=True)

    crfs_bed = BedTool.from_dataframe(crfs_df[["chr_hg38","start_hg38","end_hg38","CRS","pscore","strand_hg38"]])
    chess_bed = BedTool.from_dataframe(chess_bed)

    intersect=chess_bed.intersect(crfs_bed, wao=True)
    int_df=pd.read_table(intersect.fn,names=["seqid",
                                      "start",
                                      "end",
                                      "tid",
                                      "score",
                                      "strand",
                                      "cseqid",
                                      "cstart",
                                      "cend",
                                      "cname",
                                      "cscore",
                                      "cstrand",
                                      "distance"])
    print("total number of CRSs with partial and full overlaps: "+str(len(set(int_df[int_df["distance"]>0]["cname"]))))
    print("total number of transcripts with partial and full overlaps: "+str(len(set(int_df[int_df["distance"]>0]["tid"]))))

    int_df["full_contain"] = np.where((int_df["cstart"]>=int_df["start"])&(int_df["cend"]<=int_df["end"]),True,False)

    print("number of CRSs fully contained within an exon: "+str(len(set(int_df[int_df["full_contain"]]["cname"]))))
    print("number of transcripts with a fully contained CRS: "+str(len(set(int_df[int_df["full_contain"]]["tid"]))))

    tdf = chess_df[chess_df["type"]=="transcript"].reset_index(drop=True)
    tdf["status"] = tdf["attributes"].str.split("type \"",expand=True,n=1)[1].str.split("\"",expand=True,n=1)[0]

    full_cont_df = int_df[int_df["full_contain"]].reset_index(drop=True)
    full_cont_df = full_cont_df.merge(tdf[["tid","status"]],how="left")

    print("number of novel transcripts with a fully contained CRS: "+str(len(set(full_cont_df[full_cont_df["status"]=="novel"]["tid"]))))
    print("number of known transcripts with a fully contained CRS: "+str(len(set(full_cont_df[full_cont_df["status"]=="known"]["tid"]))))

    # also want to know the number of CRS in novel only
    cname_in_novel = set(full_cont_df[full_cont_df["status"]=="novel"]["cname"])
    cname_in_known = set(full_cont_df[full_cont_df["status"]=="known"]["cname"])

    cname_novel_only = cname_in_novel-cname_in_known
    print("number of CRSs in novel transcripts only: "+str(len(cname_novel_only)))

    # get known transcript ids
    known_tids = set(full_cont_df[full_cont_df["status"]=="known"]["tid"])

    tid2crs = dict()
    for name, grp in full_cont_df[["tid","cname"]].groupby(by="tid"):
        tid2crs[name] = ",".join(grp["cname"].tolist())

    # add CRSs to the attributes of transcripts

    max_num = 0

    with open(args.output,"w+") as outFP:
        with open(args.input,"r") as inFP:
            for line in inFP:
                if line[0]=="#":
                    outFP.write(line)
                    continue
                cols = line.split("\t")
                if cols[2]=="transcript":
                    line = line.rstrip("\n").rstrip(";")+";"
                    tid = cols[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                    
                    crs_status = tid in tid2crs
                    line+=" CRS_status \""+str(crs_status)+"\";"
                    if crs_status:
                        line+=" CRS_names \""+tid2crs[tid]+"\";"
                    
                    line+="\n";
                    outFP.write(line)
                else:
                    outFP.write(line)

    stats_fp.close()

def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument("--input",
                        required=True,
                        type=str,
                        help="Input file to be analyzed")
    parser.add_argument("--output",
                        required=True,
                        type=str,
                        help="Output filename")
    parser.add_argument("--crs",
                        required=True,
                        type=str,
                        help="file with a list of CRS to be intersected with the input annotation")
    

    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main(sys.argv[1:])