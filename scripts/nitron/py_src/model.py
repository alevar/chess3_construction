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

from lightgbm import LGBMClassifier

from py_src import commons

def main(args):
    out_model_fname = args.output+".nitron.model"
    out_preds_fname = args.output+".nitron.preds"

    stats_fname = args.output+".nitron.model.stats"
    if os.path.exists(stats_fname):
        os.remove(stats_fname)
    stats_fp = open(stats_fname,"w+")

    stats_fp.write("reference: "+args.reference+"\n")
    stats_fp.write("introns: "+args.introns+"\n")
    stats_fp.write("junctions: "+args.junctions+"\n")
    stats_fp.write("fraction: "+str(args.fraction)+"\n")
    stats_fp.write("coverage: "+str(args.coverage)+"\n")
    stats_fp.write("model: "+out_model_fname+"\n")

    # load data
    ann_df = commons.load_data(args.reference,args.junctions,args.introns,args.coverage,stats_fp)

    # split data
    num_train_known = int(len(ann_df[ann_df["type"]==1])*(args.fraction/100))
    num_train_novel = int(len(ann_df[ann_df["type"]==0])*(args.fraction/100))
    train_x,train_y,test_x,test_y = commons.get_train_test(ann_df,num_train_known,num_train_novel,stats_fp)

    # classify
    lgb = LGBMClassifier(
        application='binary',
        objective='binary',
        metric='auc',
        is_unbalance='true',
        boosting='gbdt',
        num_leaves=10,
        max_depth=3,
        min_data_in_leaf=100,
        feature_fraction=0.5,
        bagging_fraction=0.3,
        bagging_freq=1,
        learning_rate=0.005,
        num_boost_round=50000,
        early_stopping_rounds=500,
        n_jobs=args.threads,
        verbose=-1
    )
    lgb.fit(train_x[commons.use_cols],train_y["type"],eval_set=[(test_x[commons.use_cols],test_y["type"])])

    lgbPickle = open(out_model_fname,'wb')
    pickle.dump(lgb,lgbPickle)
    lgbPickle.close()

    # run prediction
    all_x,all_y = commons.split_train(ann_df)
    commons.analyze(lgb,ann_df,all_x,all_y,out_preds_fname)

    stats_fp.close()