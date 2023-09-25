//
// Created by sparrow on 3/4/20.
//

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <iomanip>

#include "TrackingStats.h"
#include "gff.h"
#include "GFaSeqGet.h"

#ifndef ASSEMBLY_STATS_TRACKINGTREE_H
#define ASSEMBLY_STATS_TRACKINGTREE_H

// Types of transcripts
enum TYPE {KNOWN_TX       = 1,
           NOVEL_TX       = 2,
           UNDEFINED_TX   = -1};

// Minimum Transcript
struct MT{
    float cov; // coverage
    float tpm; // tpm
    float fpkm; // fpkm
    std::string tid; // transcript id
    std::string locus; // locus id
    std::string tissue;
    std::string sample;
    int type = TYPE::UNDEFINED_TX;
    MT(std::string tid,std::string locus,float fpkm,float tpm,float cov,std::string tissue,std::string sample):tid(tid), locus(locus), fpkm(fpkm), tpm(tpm), cov(cov), tissue(tissue),sample(sample) { }
    void set_cov(int cov){cov=cov;}
    void set_tpm(int tpm){tpm=tpm;}
    std::string get_strg(){
        return tissue+"\t"+sample+"\t"+tid+":"+locus+"|"+std::to_string(fpkm)+"|"+std::to_string(tpm)+"|"+std::to_string(cov);
    }
    void set_type(int type){
//        if(std::strcmp(this->tid.c_str(),"128-cell__STRG.2.1")==0){
//            std::cout<<"found"<<std::endl;
//        }
        if(this->type!=TYPE::UNDEFINED_TX && this->type!=type){
            std::cerr<<"repeating types detected on sample level"<<std::endl;
            exit(-1);
        }
        this->type=type;
    }
};

typedef std::map<std::string,MT> MTM; // map of transcripts per sample
typedef std::pair<MTM::iterator,bool> MTM_IT;
typedef std::map<std::string,std::vector<MTM_IT>> Loci_Sample;
typedef std::pair<Loci_Sample::iterator,bool> Loci_Sample_IT;

// Minimum transcript to tissue
struct MTT{
    std::string tid;
    std::string all_tid;
    std::string locus;
    std::string all_loc;
    std::string tissue;
    std::set<std::string> samples;
    std::vector<MTM_IT> txs;
    int type = TYPE::UNDEFINED_TX;
    MTT(std::string tid):tid(tid){};
    void set_tissue(std::string tissue){this->tissue = tissue;}
    void set_locus(std::string locus){this->locus = locus;}
    void set_all_tid(std::string new_tid){this->all_tid=new_tid;}
    void set_all_loc(std::string new_loc){this->all_loc=new_loc;}
    MTT(std::string tid,std::string locus,std::string tissue):tid(tid),locus(locus),tissue(tissue){}
    void add_tx(MTM_IT mtm_it){
        this->txs.push_back(mtm_it);
        this->samples.insert(mtm_it.first->second.sample);
    }

    std::string get_strg(){
        std::string res = "";
        for(auto& v : txs){
            res+=tid+"\t"+locus+"\t"+v.first->second.get_strg()+"\n";
        }
        return res;
    }

    int get_num_txs(){
        return txs.size();
    }

    float get_sum_tpms(){
        float res=0;
        for(auto& t : txs){
            res+=t.first->second.tpm;
        }

        return res;
    }

    void set_type(int type){
        if(this->type!=TYPE::UNDEFINED_TX && this->type!=type){
            std::cerr<<"repeating types detected on tissue level"<<std::endl;
            exit(-1);
        }
        this->type=type;
        // propagate down
        for(auto& t : txs){
            t.first->second.set_type(type);
        }
    }
};

typedef std::map<std::string,MTT> MTTM; // map of tissue merged transcripts to sample transcripts
typedef std::pair<MTTM::iterator,bool> MTTM_IT;
typedef std::set<std::string> Tissues;
typedef std::pair<std::set<std::string>::iterator,bool> Tissues_IT;
typedef std::map<std::string,std::vector<MTTM_IT>> Loci_Tissue;
typedef std::pair<Loci_Tissue::iterator,bool> Loci_Tissue_IT;

struct MATT{
    std::string tid;
    std::string locus;
    int num_exons,elen;
    std::vector<MTTM_IT> txs;
    std::set<std::string> tissues;
    int type = TYPE::UNDEFINED_TX;
    MATT(std::string tid,std::string loc):tid(tid),locus(loc){}
    void add_tx(MTTM_IT mttm_it){
        this->txs.push_back(mttm_it);
        this->tissues.insert(mttm_it.first->second.tissue);
    }
    int num_tissues(){
        return this->tissues.size();
    }
    std::string get_strg(){
        std::string res = "";
        for(auto& v : txs){
            res+=tid+"\t"+locus+"\t"+v.first->second.get_strg()+"\n";
        }

        return res;
    }

    int get_num_sample_txs(){
        int res=0;
        for(auto& v : txs){
            res += v.first->second.get_num_txs();
        }

        return res;
    }

    float get_sum_sample_tpms(){
        float res=0;
        for(auto& t : txs){
            res += t.first->second.get_sum_tpms();
        }

        return res;
    }

    int set_type(int type){
        if(this->type!=TYPE::UNDEFINED_TX){
            std::cerr<<"Type was already set for: "<<tid<<std::endl;
            exit(-1);
        }
        this->type=type;
        // now propagate down to the tissue and sample levels
        for(auto& v : txs){
            v.first->second.set_type(type);
        }
        return 1;
    }

    void set_elen(int el){
        this->elen = el;
    }
    void set_num_exons(int ne){
        this->num_exons = ne;
    }
};

typedef std::map<std::string,MATT> MATTM; // map of all merged transcripts to tissue merged transcripts
typedef std::pair<MATTM::iterator,bool> MATTM_IT;
typedef std::map<std::string,std::pair<std::vector<MATTM_IT>,int>> Loci; // value contains 1. iterators to transcripts; 2. type of locus (real,noise)
typedef std::pair<Loci::iterator,bool> Loci_IT;

struct TR{
    std::string new_tid,new_loc,old_tid,old_loc,qj,fpkm,cov;
    float tpm;
    TR() = default;
    TR(std::string new_tid,std::string new_loc,std::string old_tid,std::string old_loc,std::string qj,std::string fpkm,float tpm,std::string cov):new_tid(new_tid),new_loc(new_loc),old_tid(old_tid),old_loc(old_loc),qj(qj),fpkm(fpkm),tpm(tpm),cov(cov){}
};

class TrackingTree {
public:
    TrackingTree(std::string ttf,std::string atf,std::string agf,float tpm_thresh):tissue_tracking_fname(ttf),all_tracking_fname(atf),all_gtf_fname(agf),tpm_thresh(tpm_thresh){}
    void set_tmap_fname(std::string tmapf){
        this->tmap_fname = tmapf;
        this->tmap_fname_set = true;
    }
    ~TrackingTree() = default;

    void load();

    void convert(std::string gtf_list_fname,std::string base_out_fname);

    void get_stats(TrackingStats& stats,std::string base_out_fname){
        get_tx_tpms(base_out_fname);
        get_loc_tpms(base_out_fname);
        get_num_tissue_sample_per_tx(stats);
        get_num_tissue_sample_per_loc(stats);
        std::cout<<"<<<computing stats"<<std::endl;
    }

private:
    float tpm_thresh;
    MTM mtm;
    MTM_IT mtm_it;
    MTTM mttm;
    MTTM_IT mttm_it,mttm_it_tmp;
    MATTM mattm;
    MATTM_IT mattm_it;
    Loci loci;
    Loci_IT loci_it;
    Tissues tissues;
    Tissues_IT tissues_it;
    Loci_Tissue loci_tissue;
    Loci_Tissue_IT loci_tissue_it;
    Loci_Sample loci_sample;
    Loci_Sample_IT loci_sample_it;

    std::map<std::string,std::set<std::string>> tts; // tissue to sample map; stores the sample names for each tissue
    std::pair<std::map<std::string,std::set<std::string>>::iterator,bool> tts_it;

    std::string get_tissue_name(std::string fname);

    void break_tracking(std::vector<TR>& trs,std::string& tline);

    void add_tracking(std::string fname);

    void load_tt();

    void load_at();

    void load_types_gtf();
    void load_types_tmap();

    std::string tissue_tracking_fname,all_tracking_fname,all_gtf_fname,tmap_fname;
    bool tmap_fname_set = false;

    void _check_dups(std::map<std::string,std::pair<std::vector<std::string>,bool>>& dups,std::string& gtf_fname);
    void _convert_gtf(std::string gtf_fname,std::string tissue,std::string base_out_fname,std::map<std::string,std::pair<std::vector<std::string>,bool>>& dups);

    // GETTERS

    void get_tx_tpms(std::string base_out_fname); // get all tpms for each transcript across the dataset - separated by tissue (/) then by sample (;)
    void get_loc_tpms(std::string base_out_fname); // get all tpms for each locus across the dataset - groupped by tissue (/) then by sample (;)
    void get_num_tissue_sample_per_tx(TrackingStats& stats); // get number of tissues and samples and samples per tissue each transcript occurs in
    void get_num_tissue_sample_per_loc(TrackingStats& stats); // get number of tissues and samples and samples per tissue each locus occurs in
};


#endif //ASSEMBLY_STATS_TRACKINGTREE_H
