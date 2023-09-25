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

def get_atts(att_string):
    atts_tmp = lcs[8].split(";")
    atts = dict()
    for kv in atts_tmp:
        k,v = atts_tmp.split(" \"")
        assert len(k)>0,"empty key"
        assert len(v)>0,"empty value"    
        v = v.rtrim("\"")
        atts[k]=v
    return atts

def run(args):
    stats_fname = args.output+".reader.stats"
    if os.path.exists(stats_fname):
        os.remove(stats_fname)
    stats_fp = open(stats_fname,"w+")

    # run reader on the reference annotations
    reader_res = dict()
    for name,fname in [(x.split(":")[0],x.split(":")[1]) for x in args.annotations.split(",")]:
        ref_out_gtf_fname = args.output+".reader_"+name+".gtf"
        reader_res[name]={"infname":fname,"outfname":ref_out_gtf_fname,"status":dict()}

        reader_cmd = [args.reader,
                     "-i",args.input,
                     "-r",fname,
                     "-o",ref_out_gtf_fname]

        stats_fp.write(" ".join(reader_cmd))
        subprocess.call(reader_cmd)

        with open(ref_out_gtf_fname,"r") as inFP:
            for line in inFP:
                lcs=line.split("\t")
                if lcs[2]=="transcript":
                    tid = lcs[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                    status = lcs[8].split("reader_status \"",1)[1].split("\"",1)[0]
                    reader_res[name]["status"][tid]=status


    # consolidate the two versions by using the union of two sets of identified readthroughs
    # this way we are still being strict for prot-coding (even if present in refgen as readthrough) but also allow for possible readthroughs in reagions not covered by mane

    # if the transcript has status pass in mane - include reagardless - otherwise decide based on refgen
    with open(args.output,"w+") as outFP:
        with open(args.input,"r") as inFP:
            for line in inFP:
                lcs=line.rstrip("\n").rstrip().rstrip(";").split("\t")
                if lcs[2]=="transcript":
                    tid = lcs[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                    status = "pass"
                    for refname in set(reader_res):
                        ref_status = reader_res[refname]["status"][tid]
#                         if ref_status=="pass" or status=="pass":
#                             status="pass"
#                         elif ref_status=="fail" and status=="unknown":
#                             status="fail"
#                         else:
#                             continue
                        if ref_status=="fail" or status=="fail":
                            status="fail"
                        lcs[8] = lcs[8]+"; reader_"+refname+"_status \""+ref_status+"\""

                    lcs[8] = lcs[8]+"; reader_status \""+status+"\";"
                
                outFP.write("\t".join(lcs)+"\n")

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
    parser.add_argument("--reader",
                        required=True,
                        type=str,
                        help="Path to the reader executable")
    parser.add_argument("--annotations",
                        required=True,
                        type=str,
                        help="Comma separated list of reference annotations to use in the following format: <name1>:<path1>,<name2>:<path2>")
    

    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main(sys.argv[1:])