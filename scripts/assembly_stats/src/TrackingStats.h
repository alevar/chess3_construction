//
// Created by Ales Varabyou on 3/4/20.
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

#ifndef ASSEMBLY_STATS_TRACKINGSTATS_H
#define ASSEMBLY_STATS_TRACKINGSTATS_H

#define FIXED_FLOAT(x)

class TrackingStats{
public:
    TrackingStats() = default;
    ~TrackingStats() = default;

    std::vector<std::tuple<std::string,int,int,int,std::vector<std::vector<float> >,std::vector<std::vector<float> >,std::string > > num_tissue_sample_per_tx; // tid, class_code, vector in which each element is a tissue and the value of the element is the code,elen,number of exons and a vector of tpms for each sample within that tissue, lid. Two vectors are one for tpms above the threshold and the second for tpms below threshold

    std::map<std::string,std::tuple<int,std::map<std::string,std::set<std::string> > > > num_tissue_sample_per_loc; // lid, class_code, vector in which each element is a map from tissue to samples
    std::pair<std::map<std::string,std::tuple<int,std::map<std::string,std::set<std::string> > > >::iterator,bool> ntspl_it;

//    std::map<std::string,std::vector<std::tuple<std::vector<int> > > > tx_tpms;
//    std::pair<std::map<std::string,std::vector<std::tuple<std::vector<int> > > >::iterator,bool> tt_it;

    void save_stats(std::string base_out_fname){
        std::cout<<"saving stats"<<std::endl;
        save_num_tissue_sample_per_tx(base_out_fname);
        save_num_tissue_sample_per_loc(base_out_fname);
        save_merged_features(base_out_fname);
        std::cout<<"done saving stats"<<std::endl;
    }

private:
    void save_num_tissue_sample_per_tx(std::string base_out_fname);
    void save_num_tissue_sample_per_loc(std::string base_out_fname);

    void save_merged_features(std::string base_out_fname);
};


#endif //ASSEMBLY_STATS_TRACKINGSTATS_H
