//
// Created by sparrow on 3/4/20.
//

# include <cmath>
#include <numeric>

#include "TrackingStats.h"

void TrackingStats::save_num_tissue_sample_per_tx(std::string base_out_fname){
    std::string tx_counts_fname = base_out_fname+".tx_counts";
    std::ofstream tx_counts_ss(tx_counts_fname.c_str());
    tx_counts_ss<<"tid,code,num_tissues,num_samples_per_tissue_pass,num_samples_pass,num_samples_per_tissue_fail,num_samples_fail"<<std::endl;

    for(auto &atx : this->num_tissue_sample_per_tx){ // iterate over ALL transcripts
        tx_counts_ss<<std::get<0>(atx)<<","<<char(std::get<1>(atx))<<","<<std::get<4>(atx).size()<<",";
        int num_samples_pass = 0;
        for(auto &val : std::get<4>(atx)){
            num_samples_pass = num_samples_pass+val.size();
            tx_counts_ss<<val.size()<<";";
        }
        tx_counts_ss.seekp(-1, std::ios_base::end);
        tx_counts_ss<<","<<num_samples_pass<<",";

        int num_samples_fail = 0;
        for(auto &val : std::get<5>(atx)){
            num_samples_fail = num_samples_fail+val.size();
            tx_counts_ss<<val.size()<<";";
        }
        tx_counts_ss.seekp(-1, std::ios_base::end);
        tx_counts_ss<<","<<num_samples_fail<<std::endl;
    }

    tx_counts_ss.close();
}

void TrackingStats::save_num_tissue_sample_per_loc(std::string base_out_fname){
    std::string loc_counts_fname = base_out_fname+".loc_counts";
    std::ofstream loc_counts_ss(loc_counts_fname.c_str());
    loc_counts_ss<<"lid,code,num_tissues,num_samples_per_tissue,num_samples"<<std::endl;

    for(auto &atx : this->num_tissue_sample_per_loc){ // iterate over ALL loci
        loc_counts_ss<<atx.first<<","<<std::get<0>(atx.second)<<","<<std::get<1>(atx.second).size()<<",";
        int num_samples = 0;
        for(auto &val : std::get<1>(atx.second)){
            num_samples = num_samples+val.second.size();
            loc_counts_ss<<val.second.size()<<";";
        }
        loc_counts_ss.seekp(-1, std::ios_base::end);
        loc_counts_ss<<","<<num_samples<<std::endl;
    }

    loc_counts_ss.close();
}

void get_tx_stats(std::vector<std::vector<float> >& samples_tpms,int& ns,int& med_ns,float& mean_ns,float& sd_ns,float& med_tpm,float& mean_tpm,float& sd_tpm){
    if(samples_tpms.size()==0){
        return;
    }

    std::sort(samples_tpms.begin(),samples_tpms.end(),[](const std::vector<float> & a, const std::vector<float> & b){ return a.size() < b.size();});

    // get median number of samples per tissue
    if (samples_tpms.size() % 2 != 0) { // even
        med_ns = samples_tpms[samples_tpms.size() / 2].size();
    }
    else{ // odd
        med_ns = (samples_tpms[(samples_tpms.size()-1)/2].size() + samples_tpms[samples_tpms.size()/2].size())/2.0;
    }

    // get total number of samples, mean and sd number of samples per tissue
    for(auto &val : samples_tpms){
        ns = ns+val.size();
    }
    mean_ns = ((float)ns)/((float)samples_tpms.size());

    // get sd number of samples per tissue
    float total = 0;
    for(auto& t : samples_tpms){
        total += (t.size()-mean_ns)*(t.size()-mean_ns);
    }
    sd_ns = sqrt(total/ns);

    // get all tpms into a single vector
    std::vector<float> all_tx_tpms;
    for(auto& t : samples_tpms){
        for(auto tpm : t){
            all_tx_tpms.push_back(tpm);
        }
    }
    std::sort(all_tx_tpms.begin(),all_tx_tpms.end());
    med_tpm = all_tx_tpms[all_tx_tpms.size() / 2];

    // get total number of samples, mean and sd number of samples per tissue
    float sum_tpms = std::accumulate(all_tx_tpms.begin(),all_tx_tpms.end(), 0.0);
    mean_tpm = ((float)sum_tpms)/((float)all_tx_tpms.size());

    // get sd number of samples per tissue
    float total_tpm = 0;
    for(auto& t : all_tx_tpms){
        total_tpm += (t-mean_tpm)*(t-mean_tpm);
    }
    sd_tpm = sqrt(total_tpm/sum_tpms);
}

void TrackingStats::save_merged_features(std::string base_out_fname){
    std::string agg_fname = base_out_fname+".agg";
    std::ofstream agg_ss(agg_fname.c_str());
    agg_ss<<"tid,code,elen,num_exons,"
          <<"num_tissues_tid_pass,num_samples_tid_pass,"
          <<"median_num_samples_per_tissue_tid_pass,mean_num_samples_per_tissue_tid_pass,sd_num_samples_per_tissue_tid_pass,"
          <<"median_tpm_tid_pass,mean_tpm_tid,sd_tpm_tid_pass,"
          <<"num_tissues_tid_fail,num_samples_tid_fail,"
          <<"median_num_samples_per_tissue_tid_fail,mean_num_samples_per_tissue_tid_fail,sd_num_samples_per_tissue_tid_fail,"
          <<"median_tpm_tid_fail,mean_tpm_tid,sd_tpm_tid_fail,"
          <<"lid,num_tissues_lid,num_samples_lid"<<std::endl;
    for(auto &atx : this->num_tissue_sample_per_tx){ // iterate over ALL transcripts
        if(std::get<4>(atx).empty()){ // if empty array - do not report transcript since none pass thresholds
            continue;
        }
        agg_ss<<std::get<0>(atx)<<","<<char(std::get<1>(atx))<<","<<std::get<2>(atx)<<","<<std::get<3>(atx)<<","<<std::get<4>(atx).size()<<",";
        // sort the num_samples vector
        int ns_pass=0,med_ns_pass=0;
        float mean_ns_pass=0,sd_ns_pass=0,med_tpm_pass=0,mean_tpm_pass=0,sd_tpm_pass=0;
        get_tx_stats(std::get<4>(atx),ns_pass,med_ns_pass,mean_ns_pass,sd_ns_pass,med_tpm_pass,mean_tpm_pass,sd_tpm_pass);
        agg_ss<<ns_pass<<","<<med_ns_pass<<","<<mean_ns_pass<<","<<sd_ns_pass<<","<<med_tpm_pass<<","<<mean_tpm_pass<<","<<sd_tpm_pass<<",";

        int ns_fail=0,med_ns_fail=0;
        float mean_ns_fail=0,sd_ns_fail=0,med_tpm_fail=0,mean_tpm_fail=0,sd_tpm_fail=0;
        get_tx_stats(std::get<5>(atx),ns_fail,med_ns_fail,mean_ns_fail,sd_ns_fail,med_tpm_fail,mean_tpm_fail,sd_tpm_fail);
        agg_ss<<ns_fail<<","<<med_ns_fail<<","<<mean_ns_fail<<","<<sd_ns_fail<<","<<med_tpm_fail<<","<<mean_tpm_fail<<","<<sd_tpm_fail<<",";


        // Gene information
        agg_ss<<std::get<6>(atx)<<",";
        // find corresponding entry in the locus specific map
        this->ntspl_it.first = this->num_tissue_sample_per_loc.find(std::get<6>(atx));
        if(this->ntspl_it.first == this->num_tissue_sample_per_loc.end()){
            std::cerr<<"locus not found: "<<std::get<6>(atx)<<std::endl;
            exit(-1);
        }

        int num_loc_samples = 0;
        for(auto &val : std::get<1>(this->ntspl_it.first->second)){
            num_loc_samples += val.second.size();
        }

        float mean_num_loc_samples = ((float)num_loc_samples)/((float)std::get<1>(this->ntspl_it.first->second).size());

        // TODO: need to sort values and get remaining stuff

        agg_ss<<std::get<1>(this->ntspl_it.first->second).size()<<","<<num_loc_samples<<","<<0<<","<<mean_num_loc_samples<<std::endl;
    }

    agg_ss.close();
}