# contains reusable funcitons used throughout the experiments


# 1s cds chains of two transcripts and returns a set of statistics:
# 1. number match
#    - inframe
#    - outframe
# 2. number mismatch
# 3. match start
# 4. match stop

import os
import json
import random
import pyBigWig
import upsetplot
import subprocess
import numpy as np
import pandas as pd

gff3cols=["seqid","source","type","start","end","score","strand","phase","attributes"]

def load_fasta_dict(fa_fname,rev=False,upper=False):
    res = dict()
    with open(fa_fname,"r") as inFP:
        cur_nm = None
        
        for line in inFP:
            if line[0]==">":
                cur_nm = line.strip()[1:].split()[0]
                assert cur_nm not in res,"duplicate record name: "+nm
                res[cur_nm] = ""
            else:
                assert cur_nm is not None,"empty record name"
                res[cur_nm]+=line.strip().upper()
                
    if rev:
        im = dict()
        for k, v in res.items():
            im[v] = im.get(v, []) + [k]
        
        res = im
    return res

def get_gff3cols():
    return gff3cols

def test_defs(): # simple test to check whether the file has been loaded
    print("test passed")
    
def merge_intervals(intervals):
    merged = []
    for interval in intervals:
        if not merged or merged[-1][1] < interval[0]:
            merged.append(list(interval))
        else:
            merged[-1][1] = max(merged[-1][1], interval[1])
    return [tuple(x) for x in merged]

def union_of_intervals(intervals1, intervals2):
    merged_intervals = sorted(intervals1 + intervals2, key=lambda x: x[0])
    return merge_intervals(merged_intervals)

def intersect(s1,s2):
    res = [0,-1,0]
    tis = max(s1[0],s2[0])
    tie = min(s1[1],s2[1])
    if(tis<=tie):
        res[0] = tis
        res[1] = tie
        return (tie-tis)+1,res
    return 0,res

def split(s1,s2):
    left  = [0,-1,-1]
    right = [0,-1,-1]
    
    il,inter = intersect(s1,s2)
    if il>0:
        if inter[0]>s1[0]:
            left[0] = s1[0]
            left[1] = inter[0]-1
            left[2] = s1[2]
        if inter[1]<s1[1]:
            right[0] = inter[1]+1
            right[1] = s1[1]
            right[2] = s1[2]
    else:
        if s1[0]<s2[0]:
            left = s1
        else:
            right = s1
        
    return left,inter,right

def slen(s):
    return (s[1]-s[0])+1

def clen(chain):
    res = 0
    for c in chain:
        res+=slen(c)
    return res

def compare(i1,i2):
    intervals = []
    for i in i1:
        intervals.append([i[0],i[1]])
        intervals[-1].append(-1)
    for i in i2:
        intervals.append([i[0],i[1]])
        intervals[-1].append(1)
    intervals.sort()
    
    if len(i1)==0 and len(i2)==0:
        return []
    
    stack = []
    stack.append(intervals[0])
    for i in intervals[1:]:
        
        left,inter,right = split(stack[-1],i)
        if slen(right)>0:
            assert slen(inter)==slen(i) # must be intirely contained within
        else:
            tmp,inter2,right = split(i,stack[-1])
            if(slen(tmp)>0):
                t2 = stack[-1]
                stack[-1]=tmp
                stack.append(t2)
                
            else:
                assert slen(tmp)<=0,str(tmp)+","+str(inter2)+","+str(right)
            assert inter==inter2
        
        stack.pop()
        
        if slen(left)>0:
            stack.append(left)
        if slen(inter)>0:
            inter[2] = 0
            stack.append(inter)
        if slen(right)>0:
            stack.append(right)
        
    return stack

def subset_gtf(in_gtf_fname,out_gtf_fname,tids):
    with open(out_gtf_fname,"w+") as outFP:
        with open(in_gtf_fname,"r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                if not len(lcs) == 9:
                    continue

                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if tid in tids:
                    outFP.write(line)
                    
def load_gtf(gtf_fname):
    assert os.path.exists(gtf_fname),"input file does not exist: "+gtf_fname
    res = pd.read_csv(gtf_fname,sep="\t",names=gff3cols,comment="#")
    return res

def extract_gtf_attribute(gtf_df,k,column="attributes"):
    assert column in set(gtf_df.columns), "please speify the correct column to extract attributes from"
    return gtf_df[column].str.split(k+" \"",expand=True)[1].str.split("\"",expand=True)[0]

def load_map(gtf_fname, qname, tname, pass_codes):
    assert os.path.exists(gtf_fname),"input file does not exist: "+gtf_fname
    
    res = dict()
    
    with open(gtf_fname,"r") as inFP:
        for line in inFP.readlines():
            lcs = line.split("\t")
            if not len(lcs) == 9:
                continue

            if not lcs[2]=="transcript":
                continue

            code = lcs[8].split("class_code \"",1)[1].split("\"",1)[0]
            if code not in pass_codes:
                continue

            qtid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
            ttid = lcs[8].split("cmp_ref \"", 1)[1].split("\"", 1)[0]

            assert qtid not in res,"duplicate qtid: "+qtid

            res[qtid] = [code,ttid,lcs[0],lcs[6]]
    res = pd.DataFrame.from_dict(res,orient="index").reset_index()
    res.columns = [qname,"code",tname,"seqid","strand"]
    return res

def get_attribute(gtf_fname,attrs,cols=None,feature="transcript",gff=False):
    clean_attrs = []
    if type(attrs)==list:
        clean_attrs = attrs
    else: # is string
        clean_attrs = [attrs]
        
    tid2name = dict()
    with open(gtf_fname,"r") as inFP:
        for line in inFP:
            lcs = line.split("\t")
            if not len(lcs) == 9:
                continue

            if lcs[2]==feature:
                tid = ""
                if not gff:
                    tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                else:
                    print(lcs)
                    tid = lcs[8].split("ID=",1)[1].split(";",1)[0]
                
                tid2name[tid]=[]
                if cols is not None:
                    for c in cols:
                        tid2name[tid].append(lcs[c])
                for a in clean_attrs:
                    gn = "-"
                    if a in lcs[8]:
                        if not gff:
                            gn = lcs[8].split(a+" \"", 1)[1].split("\"", 1)[0]
                        else:
                            gn = lcs[8].split(a+"=",1)[1].split(";",1)[0]
                    tid2name[tid].append(gn)
    res = pd.DataFrame.from_dict(tid2name,orient="index").reset_index()
    tmp_cols = ["tid"]
    if cols is not None:
        for c in cols:
            tmp_cols.append(c)
    for a in clean_attrs:
        tmp_cols.append(a)
    res.columns = tmp_cols
    return res

def get_chains(gtf_fname,feature_type,coords,phase=False):
    res = dict()
    with open(gtf_fname,"r") as inFP:
        for line in inFP:
            lcs = line.split("\t")
            if not len(lcs) == 9:
                continue

            if lcs[2]=="transcript":
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]

                if(coords):
                    res.setdefault(tid,[0,lcs[0],lcs[6],lcs[0]+":"+str(lcs[3])+"-"+str(lcs[4]),[]])
                else:
                    res.setdefault(tid,[0,[]])

            if not lcs[2]==feature_type:
                continue
                
            tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]

            res[tid][-1].append((int(lcs[3]),int(lcs[4])))
            if phase:
                res[tid][-1][-1].append(lcs[7])
            res[tid][0]=1
            
    for k,v in res.items():
        res[k] = v[:-1]
        res[k].append(sorted(v[-1]))
    res = pd.DataFrame.from_dict(res,orient="index").reset_index()
    if(coords):
        res.columns = ["tid","has_cds","seqid","strand","coords","chain"]
    else:
        res.columns = ["tid","has_cds","chain"]
    return res

def chain_inv(chain):
    # inverts a chain
    if len(chain)<=1:
        return list()
    
    res = []
    for i in range(1,len(chain)):
        res.append((chain[i-1][1],chain[i][0]))
        
    return res

# runs compare() funciton and labels all matches as in and out of frame accordingly
def compare_label_frame(chain1,chain2,strand):
    if chain2 is np.nan or len(chain2)==0:
        [[x[0],x[1],-1] for x in chain1]
    if chain1 is np.nan or len(chain1)==0:
        [[x[0],x[1],1] for x in chain2]

    mod_chain = compare(chain1,chain2)

    if strand=="-":
        mod_chain.reverse()

    t_frame = 0
    q_frame = 0
    
    for mc in mod_chain:
        if (mc[2] == -1): # extra positions in the query
            q_frame += slen(mc)
        elif (mc[2] == 1): # template positions missing from the query
            t_frame += slen(mc)
        elif (mc[2] == 0): # matching positions between query and template
            if (q_frame % 3 == t_frame % 3):
                mc[2] = 100 # inframe
            else:
                mc[2] = -100 # outframe
        else:
            print("wrong code")
            return

    return mod_chain

def compare_and_extract(chain1,chain2,strand):
    if chain2 is np.nan or len(chain2)==0:
        return pd.Series([[[x[0],x[1],-1] for x in chain1],-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])
    if chain1 is np.nan or len(chain1)==0:
        return pd.Series([[[x[0],x[1],1] for x in chain2],-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])

    # 1. compute the total number of matching positions between query and template
    # 2. compute the number of matching positions in frame between query and template
    mod_chain = compare(chain1,chain2)
    
    c1len = clen(chain1)
    c2len = clen(chain2)
    
    ulen = clen(union_of_intervals(chain1,chain2))

    if strand=="-":
        mod_chain.reverse()

    num_bp_extra = 0
    num_bp_missing = 0
    num_bp_inframe = 0
    num_bp_match = 0
    num_bp_outframe = 0

    t_frame = 0
    q_frame = 0
    
    for mc in mod_chain:
        if (mc[2] == -1): # extra positions in the query
            num_bp_extra += slen(mc)
            q_frame += slen(mc)
        elif (mc[2] == 1): # template positions missing from the query
            num_bp_missing += slen(mc)
            t_frame += slen(mc)
        elif (mc[2] == 0): # matching positions between query and template
            num_bp_match += slen(mc)
            if (q_frame % 3 == t_frame % 3):
                num_bp_inframe += slen(mc) # TODO: shouldn't this be stranded?
            else:
                num_bp_outframe += slen(mc)
        else:
            print("wrong code")
            return

    # compute lpd, ilpd, mlpd, etc
    lpi = int((100.0 * (float(c1len) / float(ulen))))
    ilpi = int((100.0 * (float(num_bp_inframe) / float(ulen))))
    mlpi = int((100.0 * (float(num_bp_match) / float(ulen))))

    match_start = chain1[0][0]==chain2[0][0] if strand=='+' else chain1[-1][1]==chain2[-1][1]
    match_end = chain1[-1][1]==chain2[-1][1] if strand=='+' else chain1[0][0]==chain2[0][0]

    return pd.Series([mod_chain,c1len,c2len,match_start,match_end,num_bp_extra,num_bp_missing,num_bp_inframe,num_bp_match,num_bp_outframe,lpi,ilpi,mlpi])

def load_tid2aa(fname):
    tid2aa = dict()

    with open(fname,"r") as inFP:
        cur_tid = ""
        cur_aa = ""
        for line in inFP:
            if line[0]==">":

                if not len(cur_tid)==0:
                    tid2aa[cur_tid]=cur_aa

                cur_tid = line[1:].rstrip()
                cur_aa = ""
            else:
                cur_aa += line.rstrip()

        if not len(cur_tid)==0:
            tid2aa[cur_tid]=cur_aa

    res = pd.DataFrame.from_dict(tid2aa,orient="index").reset_index()
    res.columns = ["tid","aa"]
    return res

def merge(segs):
    segs.sort()
    res = [[segs[0][0],segs[0][1]]]
    for s in segs:
        prev = res[-1]
        if s[0] <= prev[1]:
            prev[1] = max(prev[1], s[1])
        else:
            res.append([s[0],s[1]])
            
    return res

def load_segments(fname,feature_type,strandless):
    
    res = dict({"+":dict(),
                "-":dict()})
    if strandless:
        res = dict()
    
    with open(fname,"r") as inFP:
        for line in inFP:
            lcs = line.split("\t")
            if not len(lcs) == 9:
                continue

            if not lcs[2]==feature_type:
                continue
                
            if strandless:
                res.setdefault(lcs[0],set())
                res[lcs[0]].add((int(lcs[3]),int(lcs[4])))
            else:
                res[lcs[6]].setdefault(lcs[0],set())
                res[lcs[6]][lcs[0]].add((int(lcs[3]),int(lcs[4])))
            
    for k,v in res.items():
        if strandless:
            res[k] = merge(list(v))
        else:
            for k2,v2 in v.items():
                res[k][k2] = merge(list(v2))
            
    return res

def extract_from_comp(segs): # separated "shared,left,right" into separate objects
    left=[]
    shared=[]
    right=[]
    for s in segs:
        if s[2]==-1:
            left.append(s[:2])
        if s[2]==1:
            right.append(s[:2])
        if s[2]==0:
            shared.append(s[:2])
            
    return left,shared,right

# extract three-way comparisons
def upset_data(segs):
    res = {111:dict(),
           100:dict(),
             1:dict(),
            10:dict(),
           101:dict(),
           110:dict(),
            11:dict()}
    
    for seqid in segs[0]:
        for k,v in res.items():
            v.setdefault(seqid,list())
            
        s0 = segs[0][seqid]
        s1 = segs[1][seqid]
        s2 = segs[2][seqid]
        
        s0s1 = compare(s0,s1)
        s2s1 = compare(s2,s1)
        s0s2 = compare(s0,s2)
        
        s0_not_s1,s0_and_s1,s1_not_s0 = extract_from_comp(s0s1)
        s2_not_s1,s2_and_s1,s1_not_s2 = extract_from_comp(s2s1)
        s0_not_s2,s0_and_s2,s2_not_s0 = extract_from_comp(s0s2)
        
        # shared between all three
        tmp = compare(s0_and_s1,s2)
        x,s0_and_s1_and_s2,y = extract_from_comp(tmp)
        # in s0 only
        tmp = compare(s0_not_s1,s2)
        s0_not_s1_not_s2,x,y = extract_from_comp(tmp)
        # in s2 only
        tmp = compare(s2_not_s1,s0)
        s2_not_s1_not_s0,x,y = extract_from_comp(tmp)
        # in s1 only
        tmp = compare(s1_not_s0,s2)
        s1_not_s0_not_s2,x,y = extract_from_comp(tmp)
        # s0 and s2
        tmp = compare(s0_and_s2,s1)
        s0_and_s2_not_s1,x,y = extract_from_comp(tmp)
        # s0 and s1
        tmp = compare(s0_and_s1,s2)
        s0_and_s1_not_s2,x,y = extract_from_comp(tmp)
        # s1 and s2
        tmp = compare(s2_and_s1,s0)
        s2_and_s1_not_s0,x,y = extract_from_comp(tmp)

        res[111][seqid] = s0_and_s1_and_s2
        res[100][seqid] = s0_not_s1_not_s2
        res[1  ][seqid] = s2_not_s1_not_s0
        res[10 ][seqid] = s1_not_s0_not_s2
        res[101][seqid] = s0_and_s2_not_s1
        res[110][seqid] = s0_and_s1_not_s2
        res[11 ][seqid] = s2_and_s1_not_s0
        
    return res

# extract scores
def get_scores(scores_fname,ud,num_random=None):
    res = dict()
    min_score = 0
    max_score = 0
    
    bw = pyBigWig.open(scores_fname)   
    
    for k,v in ud.items():
        res.setdefault(k,[])
        for seqid,v2 in v.items():
            if not seqid in bw.chroms():
                continue
            for r in v2:
                res[k].extend(bw.values(seqid,r[0],r[1]+1))
        if len(res[k])>0:
            min_score = min(min_score,min(res[k]))
            max_score = max(max_score,max(res[k]))
        
    if num_random is not None:
        for k,v in res.items():
            if len(v)>num_random:
                res[k] = random.sample(v,num_random)
                
    return res,min_score,max_score

# get position count matrix
def pos_mat(ud,labels,scores=None):
    res =  upsetplot.from_memberships( [[],
                                        [labels[0]],
                                        [labels[1]],
                                        [labels[1], labels[0]],
                                        [labels[2]],
                                        [labels[2], labels[0]],
                                        [labels[2], labels[1]],
                                        [labels[2], labels[1], labels[0]],],data=[0,
                                                                                sum([clen(v) for k,v in ud[100].items()]),
                                                                                sum([clen(v) for k,v in ud[10].items()]),
                                                                                sum([clen(v) for k,v in ud[110].items()]),
                                                                                sum([clen(v) for k,v in ud[1].items()]),
                                                                                sum([clen(v) for k,v in ud[101].items()]),
                                                                                sum([clen(v) for k,v in ud[11].items()]),
                                                                                sum([clen(v) for k,v in ud[111].items()])])
    
    res_scores = [(sum([clen(v) for k,v in ud[100].items()]),[]),
                (sum([clen(v) for k,v in ud[10].items()]),[]),
                (sum([clen(v) for k,v in ud[110].items()]),[]),
                (sum([clen(v) for k,v in ud[1].items()]),[]),
                (sum([clen(v) for k,v in ud[101].items()]),[]),
                (sum([clen(v) for k,v in ud[11].items()]),[]),
                (sum([clen(v) for k,v in ud[111].items()]),[])]
    if scores is not None:
        res_scores = [(sum([clen(v) for k,v in ud[100].items()]),scores[100]),
                    (sum([clen(v) for k,v in ud[10].items()]),scores[10]),
                    (sum([clen(v) for k,v in ud[110].items()]),scores[110]),
                    (sum([clen(v) for k,v in ud[1].items()]),scores[1]),
                    (sum([clen(v) for k,v in ud[101].items()]),scores[101]),
                    (sum([clen(v) for k,v in ud[11].items()]),scores[11]),
                    (sum([clen(v) for k,v in ud[111].items()]),scores[111])]
    
    return [res,res_scores]

# get_mean_score
def mean_score(score_fname):
    means = []
    bw = pyBigWig.open(score_fname)
    for seqid,l in bw.chroms().items():
        means.append(bw.stats(seqid,0,l-1))
    return np.mean(means)

# extract sashimi and gtf for a specified set of transcripts based on several annotations
def extract_sashimi(sbin,cmp_gtf_fname,ref_gtf_fname,q_gtf_fname,out_base_fname,cmp_tid,qtids,title_str):
    out_gtf_fname = out_base_fname+".gtf"
    out_svg_fname = out_base_fname+".svg"
    
    with open(out_gtf_fname,"w+") as outFP:
        # first write out MANE
        with open(cmp_gtf_fname,"r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                if not len(lcs) == 9:
                    continue
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if tid == cmp_tid:
                    outFP.write(line)
        # next write out ORFanage
        with open(q_gtf_fname,"r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                if not len(lcs) == 9:
                    continue
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if tid in qtids or tid == qtids:
                    lcs[8] = "transcript_id \"ORFanage:"+tid+"\""
                    line = "\t".join(lcs)+"\n"
                    outFP.write(line)
        # lastly write out regular RefSeq
        with open(ref_gtf_fname,"r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                if not len(lcs) == 9:
                    continue
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if tid in qtids or tid == qtids:
                    outFP.write(line)
    
    sashimi_cmd = [sbin,
                   "--compare",cmp_tid,
                   "--title",title_str,
                   "--gtf",out_gtf_fname,
                   "-o",out_svg_fname]
    print(" ".join(sashimi_cmd))
    subprocess.call(sashimi_cmd)
    
def get_poly_gids(gtf_fname):
    df = get_chains(gtf_fname,"CDS",True)
    df = df[df["has_cds"]==1].reset_index(drop=True)
    df["seqid"]=df["coords"].str.split(":",n=1,expand=True)[0]
    df["start"] = df["chain"].apply(lambda row: row[0][0])
    df["end"] = df["chain"].apply(lambda row: row[-1][1])
    # add gene ids
    gid=pd.read_csv(gtf_fname,sep="\t",names=gff3cols,comment="#")
    gid=gid[gid["type"]=="transcript"].reset_index(drop=True)
    gid["tid"]=gid["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    gid["gid"]=gid["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    gid = gid[["gid","tid"]]

    df = df.merge(gid,on="tid",how="left",indicator=False)

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    df.sort_values(by=["seqid","strand","start","end"],ascending=True,inplace=True)

    df = df.groupby(by=["seqid","strand","gid"]).agg({"start":min,"end":max}).reset_index()
    df.sort_values(by=["seqid","strand","start","end"],ascending=True,inplace=True)
    df["nc"]=df.seqid.shift(-1)
    df["nt"]=df.strand.shift(-1)
    df["ns"]=df.start.shift(-1)
    df["nid"]=df.gid.shift(-1)
    df.fillna(0,inplace=True)
    df["od"] = np.where((df["seqid"]==df["nc"]) & 
                               (df["strand"]==df["nt"]) & 
                               (df["end"]>df["ns"]),1,0)
    pids = set(df[df["od"]==1]["gid"]).union(set(df[df["od"]==1]["nid"]))
    
    return pids

def clean(in_gtf_fname,out_gtf_fname,gffread_bin,fa_fname,remove_single_exon,keep_seqids = None,keep_gids=None):
    # deal with stops
    tmp1_fname = out_gtf_fname+".tmp1.gtf"
    cmd = [gffread_bin,
           "-g",fa_fname,
           "--adj-stop","-V","-T","-F","-J",
           "-o",tmp1_fname,
           in_gtf_fname]

    print(" ".join(cmd))
    subprocess.call(cmd)
    
    discard_gids = set()
    if keep_gids is not None and len(keep_gids)>0:
        t2g = get_attribute(in_gtf_fname,"gene_id")
        discard_gids = set(t2g["gene_id"])-keep_gids
        print("discarding these many genes that are not in the kept_genes: "+str(len(discard_gids)))
    
    seids = set()
    if remove_single_exon:
        df=pd.read_csv(tmp1_fname,sep="\t",names=gff3cols,comment="#")
        df=df[df["type"]=="exon"].reset_index(drop=True)
        df["tid"]=df["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
        df["gid"]=df["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
        gdf = df[["tid","seqid"]].groupby(by="tid").count().reset_index()
        gdf.columns = ["tid","num_exons"]
        gdf = gdf.merge(df[["tid","gid"]],on="tid",how="left")
        ggdf = gdf[["gid","num_exons"]].groupby(by="gid").max().reset_index()
        ggdf.columns = ["gid","max_num_exons"]
        seids = set(ggdf[ggdf["max_num_exons"]==1]["gid"])
        print("number of single exon genes discarded: "+str(len(seids)))
        
    # deal with poly
    pids = get_poly_gids(tmp1_fname)
    print("number of polycistronic genes: "+str(len(pids)))
    
    # deal with exceptions
    df=pd.read_csv(tmp1_fname,sep="\t",names=gff3cols,comment="#")
    df=df[df["type"]=="transcript"].reset_index(drop=True)
    df["tid"]=df["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    df["gid"]=df["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]

    # retain chromosomes
    cids = set()
    if not keep_seqids is None:
        cids = df[~(df["seqid"].isin(keep_seqids))]["gid"]
        print("number of transcripts not on selected seqids: "+str(len(cids)))
    
    try:
        df["attr"]=df["attributes"].str.split("exception \"",expand=True)[1].str.split("\"",expand=True)[0]
    except:
        df["attr"]=None
    df["seleno"] = df["attributes"].str.lower().str.contains("selen")

    sids = set(df[df["seleno"]]["gid"])
    print("number of seleno: "+str(len(sids)))

    print("exceptions: "+", ".join(list(set(df[~(df["attr"].isna())]["attr"].tolist()))))
    eids = set(df[~(df["attr"].isna())]["gid"])
    print("number of exceptions: "+str(len(eids)))
    
    dirty_gids = seids.union(pids).union(sids).union(eids).union(cids).union(discard_gids)
    print("number of genes to discard: "+str(len(dirty_gids)))
    
    with open(out_gtf_fname,"w+") as outFP:
        with open(tmp1_fname,"r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                gid = lcs[8].split("gene_id \"",1)[1].split("\"",1)[0]
                if not gid in dirty_gids:
                    outFP.write(line)
                    
def get_coding_gids(gtf_fname):
    df=pd.read_csv(gtf_fname,sep="\t",names=gff3cols,comment="#")
    df = df[df["type"]=="CDS"].reset_index(drop=True)
    df["gid"] = df["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    return set(df["gid"])

def subset_gtf_by_seqid(in_gtf_fname,out_gtf_fname,seqids):
    with open(out_gtf_fname,"w+") as outFP:
        with open(in_gtf_fname,"r") as inFP:
            for line in inFP:
                lcs = line.rstrip().split("\t")
                if not seqids == False and lcs[0] not in seqids:
                    continue
                    
                outFP.write(line)
    

def subset_gtf(in_gtf_fname,out_gtf_fname,gids,tids):
    writing_tid = ""
    with open(out_gtf_fname,"w+") as outFP:
        with open(in_gtf_fname,"r") as inFP:
            for line in inFP:
                lcs = line.rstrip().split("\t")
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if lcs[2]=="transcript":
                    gid = lcs[8].split("gene_id \"", 1)[1].split("\"", 1)[0]
                    if not gids == False and gid in gids:
                        outFP.write(line)
                        writing_tid = tid
                        continue
                
                if not tids == False and tid in tids:
                    outFP.write(line)
                    continue
                    
                # handle non transcript ffeatures without gene_id for whihc a transcript was found based on gene-id
                if writing_tid == tid:
                    outFP.write(line)
                    continue
                
                
def get_coding_stats(gtf_fname):
    # average number of isoforms per gene

    cgids = get_coding_gids(gtf_fname)
    t2g = get_attribute(gtf_fname,"gene_id")
    t2g.columns = ["tid","gid"]
    chains = get_chains(gtf_fname,"CDS",False)
    chains = chains[chains["has_cds"]==1].reset_index(drop=True)
    chains["chain_str"] = chains.apply(lambda row: ",".join(str(x[0])+"-"+str(x[1]) for x in row["chain"]),axis=1)

    t2g_c = t2g[t2g["gid"].isin(cgids)].reset_index(drop=True)
    t2g_c_g = t2g_c.groupby(by="gid").count().reset_index()
    print("mean number of transcripts per coding gene: "+str(t2g_c_g["tid"].mean()))

    chains = chains.merge(t2g_c,on="tid",how="inner")
    chains = chains[["chain_str","gid"]]
    chains.drop_duplicates(inplace=True)
    chains_g = chains.groupby(by="gid").count().reset_index()
    print("mean number of unique proteins per coding gene: "+str(chains_g["chain_str"].mean()))
    
# extracts num_pos coordinates from the chain
# if reverse - will extract from the end
def get_coords(chain,num_pos,reverse):
    res_coords = []
    
    tmp_chain = chain
    inc = 1
    if reverse:
        inc = -1
        tmp_chain = [[x[1],x[0]] for x in tmp_chain[::-1]]
    
    for c in tmp_chain:
        for i in range(c[0],c[1]+1,inc):
            res_coords.append(i)
            if len(res_coords)>=num_pos:
                return res_coords
            
            
def contained_intervals(is1,is2,inverse=False): # if inverse - return intervals of is2 hich contain intervals in is1
    res = []
    for i1 in is1:
        for i2 in is2:
            if i1[0] >= i2[0] and i1[1] <= i2[1]:
                if inverse:
                    res.append(i2)
                else:
                    res.append(i1)
                break
    return res

def extract_attributes(attribute_str:str,gff=False)->dict: # extract attribute key values into dictionary
    attrs = attribute_str.rstrip().rstrip(";").split(";")
    attrs = [x.strip() for x in attrs]
    attrs = [x.strip("\"") for x in attrs]
    attrs_dict = dict()
    sep = " \""
    if gff:
        sep = "="
    for at in attrs:
        k,v = at.split(sep)
        attrs_dict.setdefault(k,v)
        
    return attrs_dict

def rename_attributes(attrs:dict,rename_dict=dict)->dict: # renames attributes
    res_dict = {}
    for k,v in attrs.items():
        if k in rename_dict:
            res_dict[rename_dict[k]]=v
        else:
            res_dict[k]=v
    return res_dict

def to_attribute_string(attrs:dict,gff=False,feature_type=None)->str: # converts attribute key values back into string
    order = ["ID","Parent","transcript_id","gene_id","gene_name","gene_type","db_xref","description","max_TPM","sample_count","assembly_id","tag"]
    res = ""
    sep = " "
    quote = "\""
    end = "; "
    if gff:
        assert feature_type in ["gene","transcript","exon","CDS"],"wrong type: "+str(feature_type)
        sep = "="
        quote = ""
        end = ";"
        
    for k in order:
        if k in attrs:
            if gff:
                assert ";" not in attrs[k],"invalid character in attribute: "+attrs[k]
            
            if gff and feature_type=="gene" and k=="transcript_id":
                continue
            elif gff and feature_type=="gene" and k=="gene_id":
                res+="ID="+quote+attrs[k]+quote+end
            elif gff and feature_type=="transcript" and k=="transcript_id":
                res+="ID="+quote+attrs[k]+quote+end
            elif gff and feature_type=="transcript" and k=="gene_id":
                res+="Parent="+quote+attrs[k]+quote+end
            elif gff and feature_type in ["exon","CDS"] and k=="transcript_id":
                res+="Parent="+quote+attrs[k]+quote+end
            elif gff and feature_type in ["exon","CDS"] and k=="gene_id":
                continue
            else:        
                res+=k+sep+quote+attrs[k]+quote+end
    
    # add any other attributes in sorted order
    for k in sorted(list(attrs)):
        if k not in order:
            if gff:
                assert ";" not in attrs[k],"invalid character in attribute: "+attrs[k]
            res+=k+sep+quote+attrs[k]+quote+end
    
    if not gff:
        res = res.rstrip()
    if gff:
        res = res.rstrip(";")
    return res

def get_intervals(gtf_fname,feature = "exon", invert=False):
    res_intervals = dict() # seqids and strand as keys, lists of introns as values1
    
    intervals = {}
    with open(gtf_fname, 'r') as inFP:
        for line in inFP:
            if line[0] == "#":
                continue
            lcs = line.strip().split('\t')
            if lcs[2] == feature:
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if tid not in intervals:
                    intervals[tid] = {"seqname": lcs[0],
                                  "strand": lcs[6],
                                  "intervals": []}
                intervals[tid]["intervals"].append((int(lcs[3]), int(lcs[4])))

    for tid, idata in intervals.items():
        if invert: # get introns
            for ii in range(1,len(idata["intervals"][1:]),1):
                key = (idata["seqname"],idata["strand"])
                res_intervals.setdefault(key,dict())
                rit = (idata["intervals"][ii-1][1]+1,idata["intervals"][ii][0]-1)
                res_intervals[key].setdefault(rit,set())
                res_intervals[key][rit].add(tid)
        else:
            for ii in range(len(idata["intervals"])):
                key = (idata["seqname"],idata["strand"])
                res_intervals.setdefault(key,dict())
                rit = (idata["intervals"][ii][0],idata["intervals"][ii][1])
                res_intervals[key].setdefault(rit,set())
                res_intervals[key][rit].add(tid)

    return res_intervals

def quant_AS(qry_gtf_fname: str, tmpl_gtf_fname: str) -> pd.DataFrame:
    # load tmpl orfs and extract start coordinates and check for match
    tmpl_orfs = get_chains(tmpl_gtf_fname,"CDS",True)
    tmpl_orfs["frv"] = tmpl_orfs.apply(lambda row: get_coords(row["chain"],3,reverse=False),axis=1)
    tmpl_orfs["rev"] = tmpl_orfs.apply(lambda row: get_coords(row["chain"],3,reverse=True),axis=1)
    tmpl_orfs["start"] = np.where(tmpl_orfs["strand"]=="+",tmpl_orfs["frv"],tmpl_orfs["rev"])
    tmpl_orfs["end"] = np.where(tmpl_orfs["strand"]=="+",tmpl_orfs["rev"],tmpl_orfs["frv"])
    tmpl_orfs.drop(["frv","rev"],axis=1,inplace=True)
    
    # load tmpl orfs and extract start coordinates and check for match
    qry_orfs = get_chains(qry_gtf_fname,"CDS",True)
    qry_orfs = qry_orfs[qry_orfs["has_cds"]==1].reset_index(drop=True)
    qry_orfs["frv"] = qry_orfs.apply(lambda row: get_coords(row["chain"],3,reverse=False),axis=1)
    qry_orfs["rev"] = qry_orfs.apply(lambda row: get_coords(row["chain"],3,reverse=True),axis=1)
    qry_orfs["start"] = np.where(qry_orfs["strand"]=="+",qry_orfs["frv"],qry_orfs["rev"])
    qry_orfs["end"] = np.where(qry_orfs["strand"]=="+",qry_orfs["rev"],qry_orfs["frv"])
    qry_orfs.drop(["frv","rev"],axis=1,inplace=True)
    qry_orfs = qry_orfs[~(qry_orfs["start"].isna())].reset_index(drop=True)
    qry_orfs = qry_orfs[~(qry_orfs["end"].isna())].reset_index(drop=True)
    
    tmpl_orf_starts = set([",".join([str(x[0]),str(x[1]),str(x[2])]) for x in tmpl_orfs["start"].to_list()])
    tmpl_orf_ends = set([",".join([str(x[0]),str(x[1]),str(x[2])]) for x in tmpl_orfs["end"].to_list()])
    qry_orf_starts = set([",".join([str(x[0]),str(x[1]),str(x[2])]) for x in qry_orfs["start"].to_list()])
    qry_orf_ends = set([",".join([str(x[0]),str(x[1]),str(x[2])]) for x in qry_orfs["end"].to_list()])
    
    num_alt_starts = qry_orf_starts - tmpl_orf_starts
    num_alt_ends = qry_orf_ends - tmpl_orf_ends
    
    # next we need a script to extract exon skipping events
    # exon skipping: tmpl exon entirely contained within an intron of query
    qry_introns = get_intervals(qry_gtf_fname,feature="exon",invert=True)
    tmpl_exons = get_intervals(tmpl_gtf_fname,feature="exon",invert=False)
    
    # total number of tmpl exons
    tmpl_exon_count = 0
    for k,v in tmpl_exons.items():
        tmpl_exon_count+=len(v)
    
    skipped_tmpl_exon_count = 0
    tmpl_exon_skipping_event_tids = set() # chess transcripts hich have one or more exon skipping events
    
    for (seqid,strand),tmpl_intervals in tmpl_exons.items():
        res_contained = contained_intervals(list(tmpl_intervals),list(qry_introns[(seqid,strand)]),False)
        res_contained_in = contained_intervals(list(tmpl_intervals),list(qry_introns[(seqid,strand)]),True)
        for r in res_contained:
            skipped_tmpl_exon_count+=1
        for r in res_contained_in:
            tmpl_exon_skipping_event_tids.update(qry_introns[(seqid,strand)][r])
            
    
    # subset of the same but using only tmpl CDS pieces and see how many of those are gone
    
    # cds skipping: tmpl cds entirely contained within an intron of query
    qry_introns = get_intervals(qry_gtf_fname,feature="exon",invert=True)
    tmpl_cds = get_intervals(tmpl_gtf_fname,feature="CDS",invert=False)
    
    # total number of tmpl cds
    tmpl_cds_count = 0
    for k,v in tmpl_cds.items():
        tmpl_cds_count+=len(v)
    
    skipped_tmpl_cds_count = 0
    tmpl_cds_skipping_event_tids = set() # query transcripts hich have one or more cds skipping events
    
    for (seqid,strand),tmpl_intervals in tmpl_cds.items():
        res_contained = contained_intervals(list(tmpl_intervals),list(qry_introns[(seqid,strand)]),False)
        res_contained_in = contained_intervals(list(tmpl_intervals),list(qry_introns[(seqid,strand)]),True)
        for r in res_contained:
            skipped_tmpl_cds_count+=1
        for r in res_contained_in:
            tmpl_cds_skipping_event_tids.update(qry_introns[(seqid,strand)][r])
    
    # intron retention - how many of the tmpl introns are entirely contained within query exons
    
    qry_exons = get_intervals(qry_gtf_fname,feature="exon",invert=False)
    tmpl_introns = get_intervals(tmpl_gtf_fname,feature="exon",invert=True)
    
    # total number of tmpl introns
    tmpl_intron_count = 0
    for k,v in tmpl_introns.items():
        tmpl_intron_count+=len(v)
    
    retained_tmpl_intron_count = 0
    tmpl_intron_retention_event_tids = set()
    
    for (seqid,strand),tmpl_intervals in tmpl_introns.items():
        res_contained = contained_intervals(list(tmpl_intervals),list(qry_exons[(seqid,strand)]),False)
        res_contained_in = contained_intervals(list(tmpl_intervals),list(qry_exons[(seqid,strand)]),True)
        for r in res_contained:
            retained_tmpl_intron_count+=1
        for r in res_contained_in:
            tmpl_intron_retention_event_tids.update(qry_exons[(seqid,strand)][r])
    
    # intron retention - using cds this time
    
    qry_exons = get_intervals(qry_gtf_fname,feature="exon",invert=False)
    tmpl_introns = get_intervals(tmpl_gtf_fname,feature="CDS",invert=True)
    
    # total number of tmpl introns
    tmpl_cds_intron_count = 0
    for k,v in tmpl_introns.items():
        tmpl_cds_intron_count+=len(v)
    
    retained_tmpl_cds_intron_count = 0
    tmpl_cds_intron_retention_event_tids = set()
    
    for (seqid,strand),tmpl_intervals in tmpl_introns.items():
        res_contained = contained_intervals(list(tmpl_intervals),list(qry_exons[(seqid,strand)]),False)
        res_contained_in = contained_intervals(list(tmpl_intervals),list(qry_exons[(seqid,strand)]),True)
        for r in res_contained:
            retained_tmpl_cds_intron_count+=1
        for r in res_contained_in:
            # print(r)
            # print(qry_introns[(seqid,strand)][r])
            tmpl_cds_intron_retention_event_tids.update(qry_exons[(seqid,strand)][r])
    
    # report stop/start, intron retention, exon-skipping
    df = pd.DataFrame([],columns=["count","type"])
    df.loc[len(df.index)]=[len(num_alt_starts),"Alternative Starts"]
    df.loc[len(df.index)]=[len(num_alt_ends),"Alternative Ends"]
    df.loc[len(df.index)]=[tmpl_exon_count,"Total number of Exons in Template"]
    df.loc[len(df.index)]=[skipped_tmpl_exon_count,"Template Exons Skipped"]
    df.loc[len(df.index)]=[len(tmpl_exon_skipping_event_tids),"Query Transcripts Skipping tmpl Exons"]
    df.loc[len(df.index)]=[tmpl_cds_count,"Total number of Coding Exons in Template"]
    df.loc[len(df.index)]=[skipped_tmpl_cds_count,"Template Coding Exons Skipped"]
    df.loc[len(df.index)]=[len(tmpl_cds_skipping_event_tids),"Query Transcripts Skipping tmpl Coding Exons"]
    df.loc[len(df.index)]=[tmpl_intron_count,"Total number of introns in Template"]
    df.loc[len(df.index)]=[retained_tmpl_intron_count,"Template Introns Retained"]
    df.loc[len(df.index)]=[len(tmpl_intron_retention_event_tids)," Query Transcripts Retaining tmpl Introns"]
    df.loc[len(df.index)]=[tmpl_cds_intron_count,"Total number of Coding Introns in Template"]
    df.loc[len(df.index)]=[retained_tmpl_cds_intron_count,"Template Coding Introns Retained"]
    df.loc[len(df.index)]=[len(tmpl_cds_intron_retention_event_tids),"Query Transcripts Retaining tmpl Coding Introns"]
    
    return df


# converts gtf to bed format. for every gene id get the minimum start and maximum end as a single bed interval
def gtf_to_gene_bed(gtf_fname: str,bed_fname: str, offset:int=1):
    df = load_gtf(gtf_fname)
    df = df[df["type"]=="transcript"].reset_index()
    df["gid"] = df["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    df = df[["gid","seqid","start","end"]].groupby(by="gid").agg({"seqid":"min","start":"min","end":"max"}).reset_index()
    df["start"] = np.where((df["start"]-offset)<1,1,df["start"]-offset)
    df["end"] = df["end"]+offset
    df.sort_values(by=["seqid","start","end"],ascending=True)
    df[["seqid","start","end"]].to_csv(bed_fname,sep="\t",index=False,header=False)
    
    
# counts the number of transcripts per gene
def count_gene_transcripts(gtf_fname: str):
    df = load_gtf(gtf_fname)
    df = df[df["type"]=="transcript"].reset_index()
    df["gid"] = df["attributes"].str.split("gene_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    df["tid"] = df["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    
    df = df[["gid","tid"]].groupby(by="gid").count().reset_index()
    df.columns = ["gid","num_txs"]
    return df


# find longest ORF in each transcript
# what do we do if there is multiple ORFs of the same length? - just count for now. could also just skip those genes alltogether
def find_longest_orfs(seq):
    longest = []
    max_len = 0

    matches = re.finditer(r'(?=(ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)))', seq)
    for match in matches:
        result = match.group(1)
        coords = [match.start(),match.start()+len(result)-1]

        if max_len<len(result):
            longest = [coords]
            max_len = len(result)
            continue
        if max_len==len(result):
            longest.append(coords)
            continue

    return longest

def trans2genome(chain, strand, zero_pos):
    chain_pos = -1
    left_to_stop = zero_pos
    found_pos = False
    if strand=='+':
        for i in range(len(chain)):
            clen = slen(chain[i])
            if left_to_stop<clen: # found the segment with the stop codon
                chain_pos = chain[i][0]+left_to_stop
                found_pos = True
                break
            
            left_to_stop-=clen
        
        if not found_pos: # return the last position
            chain_pos = chain[-1][1]
        
    else:
        for i in range(len(chain)-1,-1,-1):
            clen = slen(chain[i])
            if left_to_stop<clen: # found the cds segment with the stop codon
                chain_pos = chain[i][1]-left_to_stop
                found_pos = True
                break
            
            left_to_stop-=clen
            
        if not found_pos: # return the last position
            chain_pos = chain[0][0]
        
    assert chain_pos>=0,"unexpected chain_pos<0"
    return chain_pos


def cut_chain(chain, start, end):
    res = []
    for cs,ce in chain:
        new_cs = cs
        new_ce = ce
        if new_cs<=start and new_ce>=start:
            new_cs = start
        if new_ce>=end:
            new_ce = end
            res.append((new_cs,new_ce))
            break
        if new_ce<start or new_cs>end:
            continue
        res.append((new_cs,new_ce))
    return res


# extract number of exons for each transcript in gtf
def num_exons(gtf_fname):
    df=pd.read_csv(gtf_fname,sep="\t",names=gff3cols,comment="#")
    df=df[df["type"]=="exon"].reset_index(drop=True)
    df["tid"]=df["attributes"].str.split("transcript_id \"",expand=True)[1].str.split("\"",expand=True)[0]
    gdf = df[["tid","seqid"]].groupby(by="tid").count().reset_index()
    gdf.columns = ["tid","num_exons"]
    return gdf

# check if the quant.sf file has all the transcripts
def check_slmn(slmn_dir:str)->bool:
    meta_fname = slmn_dir.rstrip("/")+"/aux_info/meta_info.json"
    qsf_fname = slmn_dir.rstrip("/")+"/quant.sf"
    if not os.path.exists(meta_fname):
        return False
    if not os.path.exists(qsf_fname):
        print(2)
        return False
    
    with open(meta_fname, 'r') as metaFP:
        meta = json.load(metaFP)
        if not "end_time" in meta:
            return False
        
        return meta["end_time"] is not None