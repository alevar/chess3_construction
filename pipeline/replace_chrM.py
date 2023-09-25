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
    stats_fname = args.output+".replace_chrM.stats"
    if os.path.exists(stats_fname):
        os.remove(stats_fname)
    stats_fp = open(stats_fname,"w+")

    # replace mitochondrial annotation with the one from GENCODE/RefSeq (should agree)
    with open(args.output,"w+") as outFP:
        with open(args.input,"r") as inFP: # write everything from chess except for the chrM
            for line in inFP:
                if line[0]=="#":
                    outFP.write(line)
                lcs = line.split("\t")
                if not lcs[0]=="chrM":
                    outFP.write(line)
                    
        with open(args.annotation,"r") as inFP:
            for line in inFP:
                if line[0]=="#":
                    continue
                lcs = line.split("\t")
                if lcs[0]=="chrM":
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
    parser.add_argument("--annotation",
                        required=True,
                        type=str,
                        help="Annotation to use for the mitochondrial sequence")
    

    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main(sys.argv[1:])