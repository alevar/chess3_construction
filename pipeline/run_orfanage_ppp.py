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

gff3Cols = ["seqid","source","type","start","end","score","strand","phase","attributes"]

def correct_ppp(ifn,ofn): # add transcript_id to the CDS records of phylocsf
    with open(ofn,"w+") as outFP:
        with open(ifn,"r") as inFP:
            cur_tid = ""
            for line in inFP.readlines():
                if line[0]=="#":
                    outFP.write(line)
                    continue
                lineCols = line.rstrip("\n").split("\t")
                if lineCols[2]=="transcript":
                    cur_tid = lineCols[-1].split("transcript_id \"")[1].split("\";")[0]
                if lineCols[2]=="CDS":
                    line = "\t".join(lineCols[:8])+"\ttranscript_id \""+cur_tid+"\"; "+lineCols[8]+"\n"
                outFP.write(line)

def load_fasta(fname):
    res = dict()
    with open(fname,"r") as inFP:
        cur_tid = ""
        cur_seq = ""
        first = True
        for line in inFP.readlines():
            if line[0]==">":
                if not first:
                    res[cur_tid]=cur_seq
                first=False
                cur_seq = ""
                cur_tid = line[1:].strip()
            else:
                cur_seq+=line.strip()
        # add the last one
        if not first:
            res[cur_tid]=cur_seq
    return res

def trim_to_len(cds_chain,strand,new_len):
    num_bp_left = new_len
    new_chain = []
    if strand=='+':
        for cs,ce in cds_chain:
            assert ce>=cs,"incorrect pos-strand coords"
            cur_len = (ce-cs)+1
            if cur_len>=num_bp_left:
                ne = cs+(num_bp_left-1)
                new_chain.append((cs,ne))
                break
            else:
                new_chain.append((cs,ce))
                num_bp_left-=cur_len
    else:
        for cs,ce in cds_chain:
            assert cs>=ce,"incorrect neg-strand coords"
            cur_len = (cs-ce)+1
            if(cur_len>=num_bp_left):
                ne = cs-(num_bp_left-1)
                new_chain.append((cs,ne))
                break
            else:
                new_chain.append((cs,ce))
                num_bp_left-=cur_len
                
    return new_chain

def get_introns(chain):
    introns = []
    first=True
    pe = -1
    for cs,ce in chain:
        if not first:
            introns.append((pe,cs))
        first=False
        pe=ce    
    return introns

def chain_len(chain):
    return sum((max(x[1],x[0])-min(x[1],x[0]))+1 for x in chain)
    
def load_chains(fname,aas,remove_stop=False):
    res = list()
    with open(fname,"r") as inFP:
        chain = ""
        for line in inFP.readlines():
            if line[0]=="#":
                continue
            lineCols = line.strip().split("\t")
            if lineCols[2]=="transcript":
                tid = lineCols[-1].split("transcript_id \"")[1].split("\";")[0]
                aa = ""
                if tid in aas:
                    aa = aas[tid]
                strand = lineCols[6]
                res.append([tid,[],[],-1,-1,0,0,aa,strand]) # tid, cds chain, cds introns, start, end, length, num_segments,amino sequence, strand
            elif lineCols[2]=="CDS":
                if(tid in aas):
                    strand = lineCols[6]
                    cs = int(lineCols[3])
                    ce = int(lineCols[4])
                    res[-1][1].append((cs,ce)) # add segment to the CDS chain
            else:
                continue
           
        for r in res:
            if len(r[1])>0: # not empty
                if r[8]=='-': # reverse the chain to reflect the transcription direction
                    r[1] = [(x[1],x[0]) for x in r[1][::-1]]
                if remove_stop:
                    r[1] = trim_to_len(r[1],r[8],len(r[7].rstrip(".")*3))

                r[2] = get_introns(r[1])
                r[3] = r[1][0][0]
                r[4] = r[1][-1][1]
                r[5] = chain_len(r[1])
                r[6] = len(r[1])
                
    return res

def load_gtf(fname):
    df = pd.read_csv(fname,sep="\t",names=gff3Cols,comment="#")
    df = df[df["type"]=="transcript"].reset_index(drop=True)
    df["tid"]=df.attributes.str.split("transcript_id \"",expand=True,n=1)[1].str.split("\";",expand=True,n=1)[0]
    return df

def run(args):
    stats_fname = args.output+".run_orfanage_ppp.stats"
    if os.path.exists(stats_fname):
        os.remove(stats_fname)
    stats_fp = open(stats_fname,"w+")

    # run orfanage
    orfanage_cmd = [args.orfanage,
                    "-i",args.input,
                    "-c",args.annotations,
                    "-r",args.reference,
                    "-o",args.output+".reflist",
                    "-l"]
    stats_fp.write(" ".join(orfanage_cmd))
    subprocess.call(orfanage_cmd)


    # merge gencode and refseq to get a reference set for comparison
    gffcmp_cmd = ["gffcompare",
                  "-o",args.output+".reflist_gffcmp_merge", # refseq_merge_gencode
                  "-T","--no-merge"]
    gffcmp_cmd.extend(args.annotations.split(","))
    stats_fp.write(" ".join(gffcmp_cmd))
    subprocess.call(gffcmp_cmd)

    # run gffcompare against both refseq and gencode to get the number of truly novel transcripts in the perfect_fit set
    refgen_gtf_fname = args.output+".reflist_gffcmp_merge.combined.gtf"
    gffcmp_cmd = ["gffcompare",
                  "-o",args.output+".reflist.perfect",
                  "-r",refgen_gtf_fname,
                  args.output+".reflist.perfect.gtf"]
    stats_fp.write(" ".join(gffcmp_cmd))
    subprocess.call(gffcmp_cmd)

    # find orfs on imperfect transcripts that overlap protein-coding genes using phylocsf
    ppp_cmd = [args.phylocsfpp,"find-cds",
               "--output",args.output+".reflist.imperfect.phylocsf",
               args.reference,
               args.ppp_idx,
               args.output+".reflist.imperfect.gtf"]
    stats_fp.write(" ".join(ppp_cmd))
    subprocess.call(ppp_cmd)


    # run on fitted CDS for comparisons
    ppp_cmd = [args.phylocsfpp,"find-cds",
               "--output",args.output+".reflist.perfect.phylocsf",
               args.reference,
               args.ppp_idx,
               args.output+".reflist.perfect.gtf"]
    stats_fp.write(" ".join(ppp_cmd))
    subprocess.call(ppp_cmd)

    out_base = args.output.split("/")[-1]
    ppp_fname = args.output+".reflist.imperfect.phylocsf/"+out_base+".reflist.imperfect.PhyloCSF++.gtf"
    ppp_out_fname = args.output+".reflist.imperfect.phylocsf/"+out_base+".reflist.imperfect.PhyloCSF++.cor.gtf"
    correct_ppp(ppp_fname,ppp_out_fname)

    # extract aa sequences from ppp
    ppp_out_aa_fname = args.output+".reflist.imperfect.phylocsf/"+out_base+".reflist.imperfect.PhyloCSF++.cor.aa.fa"
    gffread_cmd = ["gffread",
                   "-y",ppp_out_aa_fname,
                   "-g",args.reference,
                   ppp_out_fname]
    stats_fp.write(" ".join(gffread_cmd))
    subprocess.call(gffread_cmd)

    # extract nt sequences from ppp
    ppp_out_nt_fname = args.output+".reflist.imperfect.phylocsf/"+out_base+".reflist.imperfect.PhyloCSF++.cor.nt.fa"
    gffread_cmd = ["gffread",
                   "-x",ppp_out_nt_fname,
                   "-g",args.reference,
                   ppp_out_fname]
    stats_fp.write(" ".join(gffread_cmd))
    subprocess.call(gffread_cmd)

    # extract aa sequences from orfanage
    orf_out_aa_fname = args.output+".reflist.imperfect.aa.fa"
    gffread_cmd = ["gffread",
                   "-y",orf_out_aa_fname,
                   "-g",args.reference,
                   args.output+".reflist.imperfect.gtf"]
    stats_fp.write(" ".join(gffread_cmd))
    subprocess.call(gffread_cmd)

    # extract nt sequences from orfanage
    orf_out_nt_fname = args.output+".reflist.imperfect.nt.fa"
    gffread_cmd = ["gffread",
                   "-x",orf_out_nt_fname,
                   "-g",args.reference,
                   args.output+".reflist.imperfect.gtf"]
    stats_fp.write(" ".join(gffread_cmd))
    subprocess.call(gffread_cmd)


    # load tids
    chess3_tids = set(load_gtf(args.input)["tid"])

    # load perfect and imperfect orfanage results
    orf_perfect_gtf_fname = args.output+".reflist.perfect.gtf"
    orf_imperfect_gtf_fname = args.output+".reflist.imperfect.gtf"

    orf_perfect_df = load_gtf(orf_perfect_gtf_fname)
    orf_imperfect_df = load_gtf(orf_imperfect_gtf_fname)
    ppp_imperfect_df = load_gtf(ppp_out_fname)

    stats_fp.write("number of transcripts in Chess3: "+str(len(chess3_tids)))

    # check for duplicates in the outputs
    perfect_tids = orf_perfect_df[orf_perfect_df["type"]=="transcript"]["tid"].to_list()
    assert len(perfect_tids)==len(set(perfect_tids)),"found duplicate tids"

    perfect_tids = set(perfect_tids)

    perfect_chess3_tids = perfect_tids.intersection(chess3_tids)

    stats_fp.write("number of Chess3 transcripts in perfect: "+str(len(perfect_chess3_tids)))

    imperfect_tids = orf_imperfect_df[orf_imperfect_df["type"]=="transcript"]["tid"].to_list()
    assert len(imperfect_tids)==len(set(imperfect_tids)),"found duplicate tids"

    imperfect_tids = set(imperfect_tids)

    imperfect_chess3_tids = imperfect_tids.intersection(chess3_tids)

    stats_fp.write("number of Chess3 transcripts in imperfect: "+str(len(imperfect_chess3_tids)))

    # load fasta for each
    imperfect_orf_aas = load_fasta(orf_out_aa_fname)
    imperfect_ppp_aas = load_fasta(ppp_out_aa_fname)


    imperfect_chess3_chains = load_chains(orf_imperfect_gtf_fname,imperfect_orf_aas,False)
    ppp_chess3_chains = load_chains(ppp_out_fname,imperfect_ppp_aas,True)
    assert len(ppp_chess3_chains)==len(imperfect_chess3_chains),"lengths do not match"

    imperfect_ppp_chain_df = pd.DataFrame(imperfect_chess3_chains,columns=["tid",
                                                                           "cds_chain_orf",
                                                                           "intron_chain_orf",
                                                                           "cds_start_orf",
                                                                           "cds_end_orf",
                                                                           "cds_len_orf",
                                                                           "num_seg_orf",
                                                                           "aa_orf",
                                                                           "strand_orf"]).merge(pd.DataFrame(ppp_chess3_chains,columns=["tid",
                                                                                                                                        "cds_chain_ppp",
                                                                                                                                        "intron_chain_ppp",
                                                                                                                                        "cds_start_ppp",
                                                                                                                                        "cds_end_ppp",
                                                                                                                                        "cds_len_ppp",
                                                                                                                                        "num_seg_ppp",
                                                                                                                                        "aa_ppp",
                                                                                                                                        "strand_ppp"]),on="tid")
    assert len(imperfect_ppp_chain_df)==len(ppp_chess3_chains),"lengths do not match"

    # need to know whether the stop codon is present
    missing_stops = dict()
    with open(orf_imperfect_gtf_fname,"r") as inFP:
        for line in inFP:
            if line[0]=="#":
                continue
            lineCols = line.strip().split("\t")
            if lineCols[2]=="transcript":
                tid = lineCols[-1].split("transcript_id \"")[1].split("\";")[0]
                last_codon = lineCols[-1].split("last_codon \"")[1].split("\";")[0]
                missing_stops[tid] = last_codon

    # add missing_stops to the dataframe
    tmp = pd.DataFrame.from_dict(missing_stops,orient='index').reset_index()
    tmp.columns = ["tid","last_codon_orf"]
    imperfect_ppp_chain_df = imperfect_ppp_chain_df.merge(tmp,on="tid",how="left")

    # add missing start and stop for orfanage
    imperfect_ppp_chain_df["missing_start_orf"] = np.where(imperfect_ppp_chain_df["aa_orf"].str.startswith("M"),0,1)
    imperfect_ppp_chain_df["missing_stop_orf"] = np.where(imperfect_ppp_chain_df["last_codon_orf"]==".",0,1)

    stats_fp.write("total imperfect: "+str(len(imperfect_ppp_chain_df)))
    left_df = imperfect_ppp_chain_df.copy(deep=True)

    # remove everything where there is no start or stop codon
    no_start_stop_df = left_df[((left_df["missing_start_orf"]==1)|\
                                (left_df["missing_stop_orf"]==1))&\
                               (left_df["cds_start_ppp"]==-1)].reset_index(drop=True)
    left_df = left_df[~(((left_df["missing_start_orf"]==1)|\
                                (left_df["missing_stop_orf"]==1))&\
                               (left_df["cds_start_ppp"]==-1))].reset_index(drop=True)
    stats_fp.write("orfanage no start/stop and no phylocsf: "+str(len(no_start_stop_df)))
    stats_fp.write("left: "+str(len(left_df)))

    orf_ppp_match_df = left_df[(left_df["cds_start_orf"]>=0)&\
                               (left_df["cds_chain_orf"]==left_df["cds_chain_ppp"])].reset_index(drop=True)
    left_df = left_df[~((left_df["cds_start_orf"]>=0)&\
                        (left_df["cds_chain_orf"]==left_df["cds_chain_ppp"]))].reset_index(drop=True)
    stats_fp.write("orfanage+phylocsf match: "+str(len(orf_ppp_match_df)))
    stats_fp.write("left: "+str(len(left_df)))

    orf_but_no_ppp = left_df[(left_df["cds_start_orf"]>=0)&\
                             (left_df["cds_start_ppp"]==-1)].reset_index(drop=True)
    left_df = left_df[~((left_df["cds_start_orf"]>=0)&\
                        (left_df["cds_start_ppp"]==-1))].reset_index(drop=True)
    stats_fp.write("orfanage but no phylocsf: "+str(len(orf_but_no_ppp)))
    stats_fp.write("left: "+str(len(left_df)))

    no_orf_but_ppp = left_df[(left_df["cds_start_orf"]==-1)&\
                             (left_df["cds_start_ppp"]>=0)].reset_index(drop=True)
    left_df = left_df[~((left_df["cds_start_orf"]==-1)&\
                        (left_df["cds_start_ppp"]>=0))].reset_index(drop=True)
    stats_fp.write("no orfanage but phylocsf: "+str(len(no_orf_but_ppp)))
    stats_fp.write("left: "+str(len(left_df)))

    # create dictionary with CDS and attributes for each transcript
    orfs = dict()
    for index, row in no_start_stop_df.iterrows():
        orfs[row["tid"]] = {"CDS_TYPE":"no_start_stop","CDS_ID":"","CDS_text":""}

    for index, row in orf_ppp_match_df.iterrows():
        orfs[row["tid"]] = {"CDS_TYPE":"orf_and_ppp","CDS_ID":"","CDS_text":""}

    for index, row in orf_but_no_ppp.iterrows():
        orfs[row["tid"]] = {"CDS_TYPE":"orf_no_ppp","CDS_ID":"","CDS_text":""}
        
    for index, row in no_orf_but_ppp.iterrows():
        orfs[row["tid"]] = {"CDS_TYPE":"ppp_no_orf","CDS_ID":"","CDS_text":""}

    # label everything else for PPP only
    for index, row in left_df.iterrows():
        orfs[row["tid"]] = {"CDS_TYPE":"select_ppp","CDS_ID":"","CDS_text":""}

    # add perfect matches
    with open(args.output+".reflist.perfect.gtf","r") as inFP:
        for line in inFP:
            if line[0]=="#":
                continue
            cols = line.split("\t")
            if cols[2]=="transcript":
                tid = cols[8].split("transcript_id \"",1)[1].split("\"",1)[0] 
                cid = cols[8].split("CDS_tid \"",1)[1].split("\"",1)[0]
                orfs[tid] = {"CDS_TYPE":"cds_match","CDS_ID":cid,"CDS_text":""}

    # add non-overlapping 
    with open(args.output+".reflist.nonoverlapping.gtf","r") as inFP:
        for line in inFP:
            if line[0]=="#":
                continue
            cols = line.split("\t")
            if cols[2]=="transcript":
                tid = cols[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                orfs[tid] = {"CDS_TYPE":"noncoding","CDS_ID":"","CDS_text":""}
                
    # add cds_tid for imperfect found by orfanage
    with open(args.output+".reflist.imperfect.gtf","r") as inFP:
        for line in inFP:
            if line[0]=="#":
                continue
            cols = line.split("\t")
            if cols[2]=="transcript":
                tid = cols[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                if orfs[tid]["CDS_TYPE"]=="orf_and_ppp" or orfs[tid]["CDS_TYPE"]=="orf_no_ppp":
                    cid = cols[8].split("CDS_tid \"",1)[1].split("\"",1)[0]
                    orfs[tid]["CDS_ID"] = cid

    # add CDS records to the orfs dict
    with open(args.output+".reflist.perfect.gtf","r") as inFP:
        for line in inFP:
            if line[0]=="#":
                continue
            cols = line.split("\t")
            if cols[2]=="CDS":
                tid = cols[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                orfs[tid]["CDS_text"] += line
                
    with open(args.output+".reflist.imperfect.gtf","r") as inFP:
        for line in inFP:
            if line[0]=="#":
                continue
            cols = line.split("\t")
            if cols[2]=="CDS":
                tid = cols[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                if orfs[tid]["CDS_TYPE"]=="orf_and_ppp" or orfs[tid]["CDS_TYPE"]=="orf_no_ppp":
                    orfs[tid]["CDS_text"] += line
                    
    with open(args.output+".reflist.imperfect.phylocsf/"+out_base+".reflist.imperfect.PhyloCSF++.cor.gtf","r") as inFP:
        for line in inFP:
            if line[0]=="#":
                continue
            cols = line.split("\t")
            if cols[2]=="CDS":
                tid = cols[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                if orfs[tid]["CDS_TYPE"]=="ppp_no_orf" or orfs[tid]["CDS_TYPE"]=="select_ppp":
                    cols[8] = "transcript_id \""+tid+"\";\n"
                    line = "\t".join(cols)
                    orfs[tid]["CDS_text"] += line


    # create the final GTF with CDS and all attributes
    assert not os.path.exists(args.output),"output file already exists"
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
                    orf_record = orfs[tid]
                    line+=" CDS_TYPE \""+orf_record["CDS_TYPE"]+"\";"
                    if len(orf_record["CDS_ID"])>0:
                        line+=" CDS_ID \""+orf_record["CDS_ID"]+"\";"
                    line+="\n"
                    if len(orf_record["CDS_text"])>0:
                        line+=orf_record["CDS_text"]
                    outFP.write(line)
                elif cols[2]=="exon":
                    # remove any attributes that are not "CDS"
                    tid = cols[8].split("transcript_id \"",1)[1].split("\"",1)[0]
                    cols[8] = "transcript_id \""+tid+"\";\n"
                    line = "\t".join(cols)
                    outFP.write(line)
                elif cols[2]=="CDS":
                    continue
                else:
                    assert False,"unknown type of record: "+line

    stats_fp.close()

def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument('--orfanage',
                        required=True,
                        type=str,
                        help="Path to the orfanage executable")
    parser.add_argument('--phylocsfpp',
                        required=True,
                        type=str,
                        help="Path to the phylocsf++ executable")
    parser.add_argument("--ppp_idx",
                        required=True,
                        type=str,
                        help="Path to the index for phylocsf")
    parser.add_argument("--reference",
                        required=True,
                        type=str,
                        help="Path to the reference genome")
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
                        help="Comma-separated list of reference annotations")
    parser.add_argument("--threads",
                        required=False,
                        type=int,
                        help="number of threads")
    

    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main(sys.argv[1:])