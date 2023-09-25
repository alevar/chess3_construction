#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
	Author: Ales Varabyou
"""

# train and classify introns from transcriptomes

import pandas as pd
import numpy as np
import argparse
import random
import pickle
import csv
import sys
import os

gff3cols = ["seqid","source","type","start","end","score","strand","phase","attributes"]

sjs_cols = ["seqid","strand","start","end","med_num_samples","mean_qual","total_cov","unique_cov", \
            "unique_cov_ratio","strand_cov_ratio","left_cov","right_cov","left_cov_ratio","right_cov_ratio", \
            "max_start","max_end","entropy_starts","entropy_ends","combined_entropy","unique_starts", \
            "unique_ends","frac_span"]
use_cols = ["med_num_samples","unique_cov","unique_cov_ratio","overhang_cov_ratio", \
            "max_start","max_end","entropy_starts","entropy_ends","tx_count"]
y_cols = ["type","seqid","start","end","strand"]


def load_data(gff_fname,sjs_fname,introns_fname,min_known_cov,stats_fp=None):
    global gff3cols
    global sjs_cols
    global use_cols

    if not stats_fp is None:
        stats_fp.write("load_data\n")

    gtf_df = pd.read_csv(gff_fname,sep="\t",comment="#",names=gff3cols)
    gtf_df = gtf_df[gtf_df["type"]=="exon"][["seqid","start","end","strand","attributes"]].reset_index(drop=True)
    gtf_df["tid"] = gtf_df["attributes"].str.split("Parent=",expand=True,n=1)[1].str.split(";",expand=True,n=1)[0]
    gtf_df.sort_values(by=["tid","seqid","strand","start","end"],ascending=True,inplace=True)
    # get introns
    gtf_df[["next_start","next_tid"]] = gtf_df[["start","tid"]].shift(-1)
    gtf_df = gtf_df[gtf_df["tid"]==gtf_df["next_tid"]].reset_index(drop=True)
    gtf_df.dropna(inplace=True)
    gtf_df["start"]=gtf_df["end"]
    gtf_df["end"]=gtf_df["next_start"].astype(int)
    gtf_df.drop(["next_start","next_tid","attributes"],axis=1,inplace=True)
    gtf_df.drop_duplicates(["seqid","strand","start","end"],inplace=True)

    # sjs removes any multimappers and adds frac_ends information
    df = pd.read_csv(sjs_fname,sep="\t",skiprows=1,names=sjs_cols)
    df["start"] = df["start"]
    df["end"] = df["end"]+1
    if not stats_fp is None:
        stats_fp.write("Total number of introns: "+str(len(df))+"\n")
    df = df.merge(gtf_df,on=["seqid","strand","start","end"],how="left")
    df.rename(columns={"tid":"type"},inplace=True)
    df["type"] = np.where(df["type"].isna(),0,1)

    # lastly we need to remove anything that is not in the assembled transcript set
    # and add fraction of transcription explained by each intron
    introns = pd.read_csv(introns_fname,names=["seqid","strand","start","end","tx_count"])
    df = df.merge(introns,on=["seqid","strand","start","end"],how="inner")
    if not stats_fp is None:
        stats_fp.write("Number of assembled introns: "+str(len(df))+"\n")
        stats_fp.write("Number of assembled introns not in alignments: "+str(len(introns)-len(df))+"\n")

    df.reset_index(drop=True)
    if not stats_fp is None:
        stats_fp.write("Total number of annotated introns: "+str(len(df[df["type"]==1]))+"\n")
        stats_fp.write("Number of annotated introns below minimum coverage threshold: "+str(len(df[(df["type"]==1)&(df["unique_cov"]<min_known_cov)]))+"\n")
    df["type"] = np.where(df["unique_cov"]<min_known_cov,0,df["type"])

    # aggregate some of the datapoints for better accuracy
    df["overhang_cov_ratio"] = df["left_cov_ratio"]*df["right_cov_ratio"]

    if not stats_fp is None:
        stats_fp.write("\n========\n")
    return df

def get_train_test(df,num_train_known,num_train_novel,stats_fp=None):
    global gff3cols
    global use_cols
    global y_cols

    if stats_fp is not None:
        stats_fp.write("get_train_test\n")

    known_df = df[df["type"]==1].reset_index(drop=True)
    novel_df = df[df["type"]==0].reset_index(drop=True)

    assert num_train_known>0 and \
           num_train_novel>0 and \
           num_train_known<=len(df[df["type"]==1]) and \
           num_train_novel<=len(df[df["type"]==0]),"incorrect number of taining given"

    # get training and test data
    known_train_idxs = random.sample(known_df.index.to_list(),num_train_known)
    novel_train_idxs = random.sample(novel_df.index.to_list(),num_train_novel)

    known_train_df = known_df[known_df.index.isin(known_train_idxs)].reset_index(drop=True)
    novel_train_df = novel_df[novel_df.index.isin(novel_train_idxs)].reset_index(drop=True)

    train_df = pd.concat([known_train_df,novel_train_df],axis=0)
    train_df = train_df.sample(frac=1).reset_index(drop=True) # shuffle

    train_x = train_df.drop(y_cols,axis=1)
    train_y = train_df[y_cols]

    # the remainder can then be testing dataset

    known_test_df = known_df[~(known_df.index.isin(known_train_idxs))].reset_index(drop=True)
    novel_test_df = novel_df[~(novel_df.index.isin(novel_train_idxs))].reset_index(drop=True)

    test_df = pd.concat([known_test_df,novel_test_df],axis=0)
    test_df = test_df.sample(frac=1).reset_index(drop=True) # shuffle

    test_x = test_df.drop(y_cols,axis=1)
    test_y = test_df[y_cols]

    if stats_fp is not None:
        # how many known and novel are allocated into train/test
        stats_fp.write("number of training samples: "+str(len(train_x))+"\n")
        stats_fp.write("number of testing samples: "+str(len(test_x))+"\n")
        stats_fp.write("number of known items in training data: "+str(len(known_train_df))+"\n")
        stats_fp.write("number of known items in testing data: "+str(len(known_test_df))+"\n")
        stats_fp.write("number of novel items in training data: "+str(len(novel_train_df))+"\n")
        stats_fp.write("number of novel items in testing data: "+str(len(novel_test_df))+"\n")

        stats_fp.write("\n========\n")
    return train_x,train_y,test_x,test_y

def split_train(df):
    global use_cols
    global y_cols

    train_x = df[use_cols]
    train_y = df[y_cols]
    return train_x,train_y

def analyze(clf,df,all_x,all_y,out_fname):
    pred = clf.predict_proba(all_x)
    pred_df = pd.concat([pd.DataFrame(pred,columns=["0","1"]),pd.DataFrame(all_y)],axis=1)
    introns_df = pred_df.merge(df[["seqid","strand","start","end"]+use_cols],on=["seqid","strand","start","end"],how="inner")
    assert len(introns_df)==len(df),"wrong lengths"
    introns_df.to_csv(out_fname,index=False)

# create a map to link new transcripts with old
def load_tracking_map(tracking_fname,eq,ref_col_no,res): # eq - only class_code "="
    with open(tracking_fname,"r") as inFP:
        for line in inFP:
            cols = line.strip().split("\t")
            class_code = cols[3]
            if eq and not class_code=="=":
                continue

            new_tid = cols[ref_col_no].split("|")[-1]
            for tmp in cols[4:]:
                if tmp == "-":
                    continue
                orig_tid = tmp.split(":",1)[1].split("|",2)[1]
                res[orig_tid]=new_tid