#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
	Author: Ales Varabyou
"""

# train and classify introns from transcriptomes

import argparse
import sys
import os

from py_src import model
from py_src import extract

def nitron(argv):

    parser = argparse.ArgumentParser(description='''Help Page''')
    subparsers = parser.add_subparsers(help='sub-command help')

    #===========================================
    #===================MODEL===================
    #===========================================

    parser_model=subparsers.add_parser('model',help='build help')
    parser_model.add_argument('-r',
                              '--reference',
                              required=True,
                              type=str,
                              help="Reference annotation in GTF to use for selecting positive junctions for training")
    parser_model.add_argument('-s',
                              '--junctions',
                              required=True,
                              type=str,
                              help="File containing junction data")
    parser_model.add_argument('-i',
                                '--introns',
                                required=True,
                                type=str,
                                help="File containing introns and trancsript counts")
    parser_model.add_argument('-o',
                                '--output',
                                required=True,
                                type=str,
                                help="output file")
    parser_model.add_argument('-t',
                              '--threads',
                              required=False,
                              type=int,
                              default=1,
                              help="number of threads")
    parser_model.add_argument('-c',
                              '--coverage',
                              required=False,
                              type=int,
                              default=100,
                              help="Coverage threshold for known splice junctions")
    parser_model.add_argument('-f',
                              '--fraction',
                              required=False,
                              type=int,
                              default=80,
                              help="Fraction (%) of the dataset to allocate for training. The rest to be used in testing")

    parser_model.set_defaults(func=model.main)

    #===========================================
    #==================EXTRACT==================
    #===========================================

    parser_extract=subparsers.add_parser('extract',help='build help')
    parser_extract.add_argument('-l',
                              '--labels',
                              required=True,
                              type=str,
                              help="File containing labels for the transcripts to be classified")
    parser_extract.add_argument('-p',
                                '--predictions',
                                required=True,
                                type=str,
                                help="File containing intron scores and predictions from nitron model")
    parser_extract.add_argument('-a',
                                '--assembly',
                                required=True,
                                type=str,
                                help="File containing assembled transcripts to be selected from")
    parser_extract.add_argument('-n',
                                '--novel_introns',
                                required=True,
                                type=str,
                                help="Path to the executable of novel_introns to run selection of transcripts")
    parser_extract.add_argument('-o',
                              '--output',
                              required=True,
                              type=str,
                              help="output file")
    parser_extract.add_argument('-t',
                              '--threads',
                              required=False,
                              type=int,
                                default=1,
                              help="number of threads")
    parser_extract.add_argument('-s',
                                '--min_samples',
                                required=True,
                                type=int,
                                default=10,
                                help="minimum number of samples in which transcript was assembled")
    parser_extract.add_argument('--all_known_assembled_labeled',
                                required=True,
                                type=str,
                                help="Path to a file containing all assembled transcripts labeled with the reference ids")
    parser_extract.add_argument('--tb_filtered_tissue_gtf',
                                required=True,
                                type=str,
                                help="Path to a file containing transcripts assembled from filtered tiebrushed alignments")

    parser_extract.set_defaults(func=extract.main)

    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    nitron(sys.argv[1:])
