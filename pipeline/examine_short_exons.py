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
    stats_fname = args.output+".examine_short_exons.stats"
    if os.path.exists(stats_fname):
        os.remove(stats_fname)
    stats_fp = open(stats_fname,"w+")


    # get known tids
    known_tids = set()
    tx_num_exons = dict()

    with open(args.input,"r") as inFP:
        prev_tid = None
        for line in inFP:
            if line[0]=="#":
                continue
            
            lcs = line.rstrip().split("\t")
            attrs = lcs[-1]
            
            if lcs[2]=="transcript":
                status = lcs[8].split("type \"",1)[1].split("\"",1)[0]
                tid = lcs[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                if status=="known":
                    known_tids.add(tid)
                    
                if prev_tid is not None:
                    tx_num_exons[prev_tid]=exon_no
                exon_no = 0
                prev_tid = tid
                    
            if lcs[2]=="exon":
                exon_no+=1
                    
    tx_num_exons[list(tx_num_exons)[0]]

    # are any of the exons novel or are they part of the known set?
    known_exons = dict()

    for name,fname in [(x.split(":")[0],x.split(":")[1]) for x in args.annotations.split(",")]:
        with open(fname,"r") as inFP:
            for line in inFP:
                if line[0]=="#":
                    continue
                
                lcs = line.rstrip().split("\t")
                if lcs[2]=="exon":
                    exon = lcs[0]+lcs[6]+lcs[3]+"-"+lcs[4]
                    known_exons.setdefault(exon,set())
                    known_exons[exon].add(name)

    # find short exons
    num_short_exons = 0
    short_exon_tids = set()
    short_exon_order = list()

    unique_short_known_exons = set()
    unique_short_novel_exons = set()

    short_exons = dict()

    novel_tids_with_short_exons = set()

    with open(args.input,"r") as inFP:
        exon_no = 0
        for line in inFP:
            if line[0]=="#":
                continue
            
            lcs = line.rstrip().split("\t")
            attrs = lcs[-1]
            
            if lcs[2]=="exon":
                exon_no+=1
                if int(lcs[4])-int(lcs[3])<3:
                    tid = lcs[8].split("transcript_id \"",1)[1].split("\"",1)[0]
    #                 if not tid in known_tids:
                    short_exon_tids.add(tid)
                    known_exon_status = "novel"
                    exon = lcs[0]+lcs[6]+lcs[3]+"-"+lcs[4]
                    exon_igv = lcs[0]+":"+lcs[3]+"-"+lcs[4]
                    if exon in known_exons:
                        known_exon_status = ",".join(list(known_exons[exon]))
                        unique_short_known_exons.add(exon)
                    else:
                        unique_short_novel_exons.add(exon)
                        novel_tids_with_short_exons.add(tid)
                    print("novel short exon ("+str(int(lcs[4])-int(lcs[3]))+" - "+exon+"    "+exon_igv+"): "+tid+" at position: "+str(exon_no)+" out of "+str(tx_num_exons[tid])+" : is in: "+known_exon_status)
                    short_exon_order.append(exon_no)
                    short_exons.setdefault((int(lcs[4])-int(lcs[3]))+1,set())
                    short_exons[(int(lcs[4])-int(lcs[3]))+1].add(exon)
                        
            else:
                exon_no = 0

    # add tag for short exons
    # if novel with short exon and the short exon is 1st or last - remove the exon and adjust coordinates
    with open(args.output,"w+") as outFP:
        with open(args.input,"r") as inFP:
            for line in inFP:
                lcs = line.rstrip().split("\t")
                if len(lcs)<9:
                    outFP.write(line)

                if lcs[2]=="transcript":
                    tid = lcs[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                    if tid in novel_tids_with_short_exons:
                        line = line.rstrip("\n").rstrip()
                        if not line[-1]==";":
                            line+=";"
                        line+=" notes \"short_exon\";\n"
                        
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
    parser.add_argument("--annotations",
                        required=True,
                        type=str,
                        help="Comma separated list of reference annotations to use in the following format: <name1>:<path1>,<name2>:<path2>")
    

    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main(sys.argv[1:])