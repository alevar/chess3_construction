#include "GArgs.h"
#include "GStr.h"
#include "GVec.hh"
#include "htslib/sam.h"
#include <set>
#include <functional>
#include <iostream>
#include <map>
#include <tuple>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>

#define VERSION "0.0.2"

const char* USAGE="TieCov v" VERSION " usage:\n"
                  " tiecov [-o out.bedgraph] in.bam\n"
                  " Other options: \n"
                  "  -C   : use alternate (faster) coverage calculation\n"
                  "  -N   : maximum NH score (if available) to include when reporting coverage\n"
                  "  -Q   : minimum mapping quality to include when reporting coverage\n"
                  "  --all,-a  : report all positions and do not group into intervals\n"
                  "  --no_cov,-n    : aggregate coverage without taking coverage into account. Useful for getting a list of positions that pass filters\n"
                  "  -s       : compute coverage across splice junctions and write to file\n";

struct Filters{
    int max_nh = MAX_INT;
    int min_qual = -1;
} filters;

struct OutputFMT{
    bool report_all = false;
    bool no_cov     = false;
    bool report_sjs = false;
    std::string sjs_fname  = "";
    std::string sjs_combined = "";
} outFmt;

GStr outfname, infname;
FILE* outf=NULL;
bool debugMode=false;
bool verbose=false;
bool alesCov=true;

void processOptions(int argc, char* argv[]);

std::map<std::pair<uint16_t,uint32_t>,uint16_t,std::less<std::pair<uint16_t,uint32_t>> > cov_pq; // since RB+ Tree is sorted by std::greater - smallest position (pair of seqid and pos) will always be at the top (begin())
std::pair<std::map<std::pair<uint16_t,uint32_t>,uint16_t,std::less<std::pair<uint16_t,uint32_t>> >::iterator,bool>  cpq_it;

std::map<std::pair<int,int>,int> donor_ends,acceptor_starts; // holds coverage data about the last/first exonic bases for the donor and acceptor sites.
std::pair<std::map<std::pair<int,int>,int>::iterator,bool> da_it;

struct CovInterval{
    int tid=-1;
    int start=-1; // start position
    int end=-1; // end position
    int cov=-1; // current coverage for the interval

    CovInterval() = default;
    ~CovInterval(){
        if(cov_ss.is_open()){
            this->cov_ss.close();
        }
    }

    void set_outfname(std::string out_fname){
        if(out_fname=="-" || out_fname.empty()){
            return;
        }
        this->cov_ss.open(out_fname,std::ios::out);
    }

    void extend(bam_hdr_t *al_hdr,int new_tid,int pos,int new_cov){
        if(tid == -1){
            tid=new_tid;
            start=pos;
            end=pos;
            cov=new_cov;
        }
        else if(tid==new_tid && pos-end==1 && (cov==new_cov || outFmt.no_cov)){ // position can extend the interval
            end++;
        }
        else{
            write(al_hdr);
            tid=new_tid;
            start=pos;
            end=pos;
            cov=new_cov;
        }
    }

    void write(sam_hdr_t *al_hdr){
        if(outFmt.report_sjs){
            for(int cur_pos=start;cur_pos<=end;cur_pos++){
                da_it.first = donor_ends.find(std::make_pair(tid,cur_pos));
                if(da_it.first != donor_ends.end()){
                    da_it.first->second = cov;
                }

                da_it.first = acceptor_starts.find(std::make_pair(tid,cur_pos));
                if(da_it.first != acceptor_starts.end()){
                    da_it.first->second = cov;
                }
            }
        }
        if(outFmt.report_all){
            for(int cur_pos=start;cur_pos<=end;cur_pos++){
                if(cov_ss.is_open()){
                    cov_ss<<al_hdr->target_name[tid]<<"\t"<<cur_pos<<"\t"<<cov<<std::endl;
                }
                else{
                    fprintf(stdout,"%s\t%d\t%d\n",al_hdr->target_name[tid],cur_pos,cov);
                }
            }
            return;
        }
        if(cov_ss.is_open()){
            cov_ss<<al_hdr->target_name[tid]<<"\t"<<start-1<<"\t"<<end<<"\t"<<cov<<std::endl;
        }
        else{
            fprintf(stdout,"%s\t%d\t%d\t%d\n",al_hdr->target_name[tid],start-1,end,cov);
        }
    }
private:
    std::fstream cov_ss;
} covInterval;

void cleanPriorityQueue(sam_hdr_t* al_hdr){ // cleans all entries
    for(auto& pq : cov_pq){
        covInterval.extend(al_hdr,pq.first.first,pq.first.second,pq.second);
    }
    cov_pq.clear();
}

void cleanPriorityQueue(sam_hdr_t* al_hdr,int tid,int pos){ // process and remove all positions upto but not including the given position
    if(cov_pq.empty()){
        return;
    }
    while(!cov_pq.empty() && cov_pq.begin()->first < std::pair<uint16_t,uint32_t>(tid,pos)){
        covInterval.extend(al_hdr,cov_pq.begin()->first.first,cov_pq.begin()->first.second,cov_pq.begin()->second);
        cov_pq.erase(cov_pq.begin());
    }
}

class Intron{
public:
    Intron(int tid,int start,int end):tid(tid),start(start),end(end){
        this->stats = std::map<int,Stats>{
                {-1,Stats()},
                {0,Stats()},
                {1,Stats()}
        };
    };
    ~Intron()=default;

    int get_total_cov(int strand){
        return this->stats[strand].cov + this->stats[strand].cov_multi;
    }
    int get_total_cov(){
        int total_cov = 0;
        for(auto& stat : this->stats){
            total_cov += stat.second.cov + stat.second.cov_multi;
        }
        return total_cov;
    }

    int get_unique_cov(int strand){
        return this->stats[strand].cov;
    }
    int get_unique_cov(){
        int total_cov = 0;
        for(auto& stat : this->stats){
            total_cov += stat.second.cov;
        }
        return total_cov;
    }

    int get_multi_cov(int strand){
        return this->stats[strand].cov_multi;
    }
    int get_multi_cov(){
        int total_cov = 0;
        for(auto& stat : this->stats){
            total_cov += stat.second.cov_multi;
        }
        return total_cov;
    }

    float get_unique_cov_ratio(int strand){
        return (float)get_unique_cov(strand)/(float)(get_total_cov(strand));
    }
    float get_unique_cov_ratio(){
        return (float)get_unique_cov()/(float)(get_total_cov());
    }

    float get_cov_ratio(int strand){
        return (float)get_multi_cov(strand)/(float)(get_total_cov(strand));
    }
    float get_cov_ratio(){
        return (float)get_multi_cov()/(float)(get_total_cov());
    }

    float get_left_cov(){
        da_it.first = donor_ends.find(std::make_pair(this->tid,this->start));
        if(da_it.first == donor_ends.end()){
            std::cerr<<"left cov does not exist...."<<std::endl;
            exit(-1);
        }
        return da_it.first->second;
    }
    float get_right_cov(){
        da_it.first = acceptor_starts.find(std::make_pair(this->tid,this->end));
        if(da_it.first == acceptor_starts.end()){
            std::cerr<<"right cov does not exist...."<<std::endl;
            exit(-1);
        }
        return da_it.first->second;
    }

    float get_left_cov_ratio(){
//        if((float)get_total_cov()/(float)(get_left_cov()) > 1.0){
//            std::cout<<"found"<<std::endl;
//        }
        return (float)get_total_cov()/(float)(get_left_cov());
    }

    float get_right_cov_ratio(){
        return (float)get_total_cov()/(float)(get_right_cov());
    }

    int get_unique_starts(int strand){
        return std::set<int>(this->stats[strand].l_lens.begin(),this->stats[strand].l_lens.end() ).size();
    }
    int get_unique_starts(){
        std::set<int> starts;
        for(auto& stat : stats){
            std::copy(stat.second.l_lens.begin(),stat.second.l_lens.end(),std::inserter(starts,starts.end()));
        }
        return starts.size();
    }

    int get_unique_ends(int strand){
        return std::set<int>(this->stats[strand].r_lens.begin(),this->stats[strand].r_lens.end() ).size();
    }
    int get_unique_ends(){
        std::set<int> ends;
        for(auto& stat : stats){
            std::copy(stat.second.r_lens.begin(),stat.second.r_lens.end(),std::inserter(ends,ends.end()));
        }
        return ends.size();
    }

    int get_max_start_dist(int strand){
        return *std::max_element(this->stats[strand].l_lens.begin(),this->stats[strand].l_lens.end());
    }
    int get_max_start_dist(){
        std::set<int> starts;
        for(auto& stat : stats){
            std::copy(stat.second.l_lens.begin(),stat.second.l_lens.end(),std::inserter(starts,starts.end()));
        }
        return *std::max_element(starts.begin(),starts.end());
    }

    int get_max_end_dist(int strand){
        return *std::max_element(this->stats[strand].r_lens.begin(),this->stats[strand].r_lens.end());
    }
    int get_max_end_dist(){
        std::set<int> ends;
        for(auto& stat : stats){
            std::copy(stat.second.r_lens.begin(),stat.second.r_lens.end(),std::inserter(ends,ends.end()));
        }
        return *std::max_element(ends.begin(),ends.end());
    }

    int get_expected_unique_starts(int strand){
//        int total_cov = get_total_cov(strand);
        int max_start = get_max_start_dist(strand);
//        if(total_cov > max_start){
            return max_start;
//        }
//        return total_cov;
    }
    int get_expected_unique_starts(){
//        int total_cov = get_total_cov();
        int max_start = get_max_start_dist();
//        if(total_cov > max_start){
            return max_start;
//        }
//        return total_cov;
    }

    int get_expected_unique_ends(int strand){
//        int total_cov = get_total_cov(strand);
        int max_end = get_max_end_dist(strand);
//        if(total_cov > max_end){
            return max_end;
//        }
//        return total_cov;
    }
    int get_expected_unique_ends(){
//        int total_cov = get_total_cov();
        int max_end = get_max_end_dist();
//        if(total_cov > max_end){
            return max_end;
//        }
//        return total_cov;
    }

    float get_frac_startsNends(int strand){ // compute fraction of expected starts and ends which are present
        int expected_starts = get_expected_unique_starts(strand);
        int expected_ends = get_expected_unique_ends(strand);
        int starts = get_unique_starts(strand);
        int ends = get_unique_ends(strand);

        return (float)(starts+ends)/(float)(expected_starts+expected_ends);
    }
    float get_frac_startsNends(){ // compute fraction of expected starts and ends which are present
        int expected_starts = get_expected_unique_starts();
        int expected_ends = get_expected_unique_ends();
        int starts = get_unique_starts();
        int ends = get_unique_ends();

        return (float)(starts+ends)/(float)(expected_starts+expected_ends);
    }

    float get_mean_qual(int strand){
        return std::accumulate(this->stats[strand].quals.begin(),this->stats[strand].quals.end(),0.0)/(float)this->stats[strand].quals.size();
    }

    float get_mean_qual(){
        float total_sum = 0.0;
        float total_size = 0.0;
        for(auto& stat : stats){
            total_sum += std::accumulate(stat.second.quals.begin(),stat.second.quals.end(),0.0);
            total_size += stat.second.quals.size();
        }
        return total_sum/total_size;
    }

    int get_median_num_samples(int strand){
        const auto median_it = this->stats[strand].num_samples.begin() + (int)(this->stats[strand].num_samples.size() / 2);
        std::nth_element(this->stats[strand].num_samples.begin(),median_it,this->stats[strand].num_samples.end());
        return *median_it;
    }

    int get_median_num_samples(){
        std::vector<int> num_samples_all;
        for(auto& stat : stats){
            num_samples_all.insert(num_samples_all.end(),stat.second.num_samples.begin(),stat.second.num_samples.end());
        }
        const auto median_it = num_samples_all.begin() + (int)(num_samples_all.size() / 2);
        std::nth_element(num_samples_all.begin(),median_it,num_samples_all.end());
        return *median_it;
    }

    std::string get_header(){
        std::string res = "seqid,"
                          "strand,"
                          "start,"
                          "end,"
                          "med_num_samples,"
                          "mean_qual,"
                          "total_cov,"
                          "unique_cov,"
                          "unique_cov_ratio,"
                          "cov_ratio,"
                          "left_cov_ratio,"
                          "right_cov_ratio,"
                          "max_start,"
                          "max_end,"
                          "unique_starts,"
                          "unique_ends,"
                          "frac_span";
        return res;
    }

    std::string get_tab_str_sep_strands(sam_hdr_t* al_hdr){
        std::string res = "";
        for(auto& stat : this->stats){ // for each strand
            if(stat.first==-1){
                continue;
            }
            if(stat.second.cov<=0 && stat.second.cov_multi<=0){
                continue;
            }
            res += al_hdr->target_name[this->tid];
            res += "\t" + std::string(1,(char)stat.first) +
                   "\t" + std::to_string(this->start) +
                   "\t" + std::to_string(this->end) +
                   "\t" + std::to_string(get_median_num_samples(stat.first)) +
                   "\t" + std::to_string(get_mean_qual(stat.first)) +
                   "\t" + std::to_string(get_total_cov(stat.first)) +
                   "\t" + std::to_string(get_unique_cov(stat.first)) +
                   "\t" + std::to_string(get_unique_cov_ratio(stat.first)) +
                   "\t" + std::to_string(get_cov_ratio(stat.first)) +
                   "\t" + std::to_string(get_left_cov_ratio()) +
                   "\t" + std::to_string(get_right_cov_ratio()) +
                   "\t" + std::to_string(get_max_start_dist(stat.first)) +
                   "\t" + std::to_string(get_max_end_dist(stat.first)) +
                   "\t" + std::to_string(get_unique_starts(stat.first)) +
                   "\t" + std::to_string(get_unique_ends(stat.first)) +
                   "\t" + std::to_string(get_frac_startsNends(stat.first))+"\n";
        }
        if(!res.empty()){
            res.pop_back();
        }
        return res;
    }

    std::string get_comb_header(){
        std::string res = "seqid,"
                          "start,"
                          "end,"
                          "med_num_samples,"
                          "mean_qual,"
                          "total_cov,"
                          "unique_cov,"
                          "unique_cov_ratio,"
                          "cov_ratio,"
                          "left_cov_ratio,"
                          "right_cov_ratio,"
                          "max_start,"
                          "max_end,"
                          "unique_starts,"
                          "unique_ends,"
                          "frac_span";
        return res;
    }

    std::string get_tab_str_combined(sam_hdr_t* al_hdr){
        std::string res = "";
        res += al_hdr->target_name[this->tid];
        res += "\t" + std::to_string(this->start) +
                "\t" + std::to_string(this->end) +
                "\t" + std::to_string(get_median_num_samples()) +
                "\t" + std::to_string(get_mean_qual()) +
                "\t" + std::to_string(get_total_cov()) +
                "\t" + std::to_string(get_unique_cov()) +
                "\t" + std::to_string(get_unique_cov_ratio()) +
                "\t" + std::to_string(get_cov_ratio()) +
                "\t" + std::to_string(get_left_cov_ratio()) +
                "\t" + std::to_string(get_right_cov_ratio()) +
                "\t" + std::to_string(get_max_start_dist()) +
                "\t" + std::to_string(get_max_end_dist()) +
                "\t" + std::to_string(get_unique_starts()) +
                "\t" + std::to_string(get_unique_ends()) +
                "\t" + std::to_string(get_frac_startsNends());
        return res;
    }

    void add_read(bam1_t* al){ // process read to populate the stats information
        int strand = this->getXS(al);
        int multi = this->getNH(al);
        int factor = this->getYC(al);
        if(multi==1){
            this->stats[strand].cov+=(1*factor);
        }
        else{ // multimapper
            this->stats[strand].cov_multi+=(1*factor);
        }

        // look at the read and add data to the stat vectors
        this->_add_read(al);
    }
private:
    int tid=-1;
    int start=-1;
    int end=-1;

    struct Stats{
        int cov=0;
        int cov_multi=0;
        int multi_weight=0; // Coverage of multimappers weighted by the number of multimappers per read
        std::vector<int> l_lens; // length of the segment on the left side
        std::vector<int> r_lens; // length of the segment on the right side
        std::vector<int> mismatches; // positions of mismatches relative to the intron position
        std::vector<int> quals; // quality values of reads spanning an intron
        std::vector<int> edits; // edit distances of reads
        std::vector<int> num_samples; // number of samples each read is observed in
    }; // separate statistics for two strands in case the intron is observed on both

    std::map<int,Stats> stats;

    void parse_md(std::string& md_str,std::vector<int> md_pairs){
        int pos = 0; // position on the read
    }

    void _add_read(bam1_t* al){ // populate vectors by cycling through the read positions
        int pos=al->core.pos+1; // 1-based
        int query_pos = 0; // position on the query/read

        bool first_match = true; // indicates whether the first match position has been encountered - triggers cleanup of the priority queue
        bool soft_end = false; // soft-clipped end
        int first_al_base = 0;
        int last_al_base = 0;

        int intron_coord_query = 0; // poition of the intron on the query

        for (uint8_t c=0;c<al->core.n_cigar;++c){
            uint32_t *cigar_full=bam_get_cigar(al);
            int opcode=bam_cigar_op(cigar_full[c]);
            int oplen=bam_cigar_oplen(cigar_full[c]);

            switch(opcode){
                case BAM_CINS: // no change in coverage and position
                    query_pos+=oplen;
                    break;
                case BAM_CDEL: // skip to the next position - no change in coverage
                    pos+=oplen;
                    break;
                case BAM_CREF_SKIP: // skip to the next position - no change in coverage
                    // check if matches the current intron - then get left and right lengths (only aligned bases of the query - no clipping)
                    if(pos-1==start && pos+oplen==end){ // current intron
                        intron_coord_query = query_pos;
                    }
                    pos+=oplen;
                    break;
                case BAM_CSOFT_CLIP: // skip to the next position - no change in coverage
                    if(!first_match){ // means this clipping is at the end of the read since a match was already found
                        soft_end = true;
                        last_al_base = query_pos;
                    }
                    query_pos+=oplen;
                    break;
                case BAM_CMATCH: // base match - increment coverage
                    if(first_match){
                        first_match=false;
                        first_al_base = query_pos;
                    }
                    for(int i=0;i<oplen;i++){ // use to look get positions of the mismatches
                        pos++;
                        query_pos++;
                    }
                    break;
                default:
                    fprintf(stderr,"ERROR: unknown opcode: %c from read: %s",bam_cigar_opchr(opcode),bam_get_qname(al));
                    exit(-1);
            }
        }
        if(!soft_end){
            last_al_base = query_pos;
        }
        int factor = getYC(al);
        int strand = getXS(al);
        int num_samples = getYX(al);
        for(int i=0;i<factor;i++){
            this->stats[strand].r_lens.push_back(last_al_base-intron_coord_query);
            this->stats[strand].l_lens.push_back(intron_coord_query-first_al_base);
            this->stats[strand].num_samples.push_back(num_samples);
            this->stats[strand].quals.push_back(al->core.qual);
        }
    }

    int getXS(bam1_t *al){
        uint8_t* ptr_xs=bam_aux_get(al,"XS");
        if(ptr_xs){
            return bam_aux2A(ptr_xs);
        }
        return -1; // -1 indicates that strand information is not available
    }

    int getNH(bam1_t *al){
        uint8_t* ptr_nh=bam_aux_get(al,"NH");
        if(ptr_nh){
            return bam_aux2i(ptr_nh);
        }
        return 1;
    }

    int getYC(bam1_t *in_rec){
        uint8_t* ptr_yc=bam_aux_get(in_rec,"YC");
        if(ptr_yc){
            return bam_aux2i(ptr_yc);
        }
        return 1;
    }

    int getYX(bam1_t *in_rec){
        uint8_t* ptr_yx=bam_aux_get(in_rec,"YX");
        if(ptr_yx){
            return bam_aux2i(ptr_yx);
        }
        return 1;
    }

    char* getMD(bam1_t *in_rec){
        uint8_t* ptr_md=bam_aux_get(in_rec,"MD");
        if(ptr_md){
            return bam_aux2Z(ptr_md);
        }
        return NULL;
    }
};

std::map<std::tuple<int,int,int>,Intron> sjs;
std::pair<std::map<std::tuple<int,int,int>,Intron>::iterator,bool> sj_it;

void write_sjs(sam_hdr_t *al_hdr){
    std::fstream int_ss;
    int_ss.open(outFmt.sjs_fname,std::ios::out);

    std::fstream int_comb_ss;
    int_comb_ss.open(outFmt.sjs_combined,std::ios::out);

    int_ss<<sjs.begin()->second.get_header()<<std::endl;
    int_comb_ss<<sjs.begin()->second.get_comb_header()<<std::endl;

    std::string res_strand_sep = "";
    for(auto& sj : sjs){
        res_strand_sep = sj.second.get_tab_str_sep_strands(al_hdr);
        if(!res_strand_sep.empty()){
            int_ss<<res_strand_sep<<std::endl;
        }
        int_comb_ss<<sj.second.get_tab_str_combined(al_hdr)<<std::endl;
    }

    int_ss.close();
    int_comb_ss.close();
}

void incCov(sam_hdr_t *al_hdr,bam1_t *in_rec,int dupCount){
    // remove all elements from the priority queue which indicate positions smaller than current
    // iterate over positions following cigar and add to the priority queue
    int pos=in_rec->core.pos+1; // 1-based

    bool first_match = true; // indicates whether the first match position has been encountered - triggers cleanup of the priority queue

    for (uint8_t c=0;c<in_rec->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(in_rec);
        int opcode=bam_cigar_op(cigar_full[c]);
        int oplen=bam_cigar_oplen(cigar_full[c]);

        switch(opcode){
            case BAM_CINS: // no change in coverage and position
                break;
            case BAM_CDEL: // skip to the next position - no change in coverage
                pos+=oplen;
                break;
            case BAM_CREF_SKIP: // skip to the next position - no change in coverage
                if(outFmt.report_sjs){
                    sj_it = sjs.insert(std::make_pair(std::make_tuple(in_rec->core.tid,pos-1,pos+oplen),Intron(in_rec->core.tid,pos-1,pos+oplen)));
                    sj_it.first->second.add_read(in_rec);

                    // add to donor/acceptor maps
                    da_it = donor_ends.insert(std::make_pair(std::make_pair(in_rec->core.tid,pos-1),0));
                    da_it = acceptor_starts.insert(std::make_pair(std::make_pair(in_rec->core.tid,pos+oplen),0));
                }
                pos+=oplen;
                break;
            case BAM_CSOFT_CLIP: // skip to the next position - no change in coverage
                //pos+=oplen; //WRONG
                break;
            case BAM_CMATCH: // base match - increment coverage
                if(first_match){
                    cleanPriorityQueue(al_hdr,in_rec->core.tid,pos);
                    first_match=false;
                }
                for(int i=0;i<oplen;i++){ // add coverage for each position
                    cpq_it = cov_pq.insert(std::make_pair(std::make_pair(in_rec->core.tid,pos),0));
                    cpq_it.first->second+=dupCount;
                    pos+=1;
                }
                break;
            default:
                fprintf(stderr,"ERROR: unknown opcode: %c from read: %s",bam_cigar_opchr(opcode),bam_get_qname(in_rec));
                exit(-1);
        }
    }
}

void addCov(bam1_t *in_rec, int dupCount, GVec<uint16_t>& bcov, int b_start) {
    int pos=in_rec->core.pos; // 0-based
    for (uint8_t c=0;c<in_rec->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(in_rec);
        int opcode=bam_cigar_op(cigar_full[c]);
        int oplen=bam_cigar_oplen(cigar_full[c]);
        switch(opcode){
            case BAM_CINS: // no change in coverage and position
                break;
            case BAM_CDEL: // skip to the next position - no change in coverage
                pos+=oplen;
                break;
            case BAM_CREF_SKIP: // skip to the next position - no change in coverage
                pos+=oplen;
                break;
            case BAM_CSOFT_CLIP:
                break;
            case BAM_CMATCH: // base match - add coverage
                for(int i=0;i<oplen;i++) {
                    bcov[pos-b_start]+=dupCount;
                    pos++;
                }
                break;
            default:
                GError("ERROR: unknown opcode: %c from read: %s",bam_cigar_opchr(opcode),bam_get_qname(in_rec));
        }
    }
}


int getYC(bam1_t *in_rec){
    uint8_t* ptr_yc_1=bam_aux_get(in_rec,"YC");
    if(ptr_yc_1){
        return bam_aux2i(ptr_yc_1);
    }
    else{
        return 1;
    }
}

int getNH(bam1_t *in_rec){
    uint8_t* ptr_nh_1=bam_aux_get(in_rec,"NH");
    if(ptr_nh_1){
        return bam_aux2i(ptr_nh_1);
    }
    else{
        return -1;
    }
}

void flushCoverage(sam_hdr_t* hdr, GVec<uint16_t>& bcov,  int tid, int b_start) {
    if (tid<0 || b_start<0) return;
    int i=0;
    if(outFmt.report_all){
        while (i<bcov.Count()) {
            uint16_t icov=bcov[i];
            if (icov!=0){
                fprintf(outf, "%s\t%d\t%d\n", hdr->target_name[tid], b_start+i, icov);
            }
        }
        return;
    }
    while (i<bcov.Count()) {
        uint16_t icov=bcov[i];
        int j=i+1;
        while (j<bcov.Count() && icov==bcov[j]) {
            j++;
        }
        if (icov!=0)
            fprintf(outf, "%s\t%d\t%d\t%d\n", hdr->target_name[tid], b_start+i, b_start+j, icov);
        i=j;
    }
}

// >------------------ main() start -----
int main(int argc, char *argv[])  {
    processOptions(argc, argv);
    htsFile* hts_file=hts_open(infname.chars(), "r");
    if (hts_file==NULL)
        GError("Error: could not open alignment file %s \n",infname.chars());
    sam_hdr_t* hdr = sam_hdr_read(hts_file);
    int prev_tid=-1;
    bam1_t* b = bam_init1();
    if (alesCov) {
        covInterval.set_outfname(outfname.chars());
        while (sam_read1(hts_file, hdr, b) >= 0) {
            int nh = getNH(b);
            if(nh>filters.max_nh){continue;}
            if (b->core.qual<filters.min_qual) { continue; }
            int tid=b->core.tid;
            if (tid!=prev_tid) {
                cleanPriorityQueue(hdr);
                prev_tid=tid;
            }
            int accYC = getYC(b);
            incCov(hdr, b, accYC);
        }
        cleanPriorityQueue(hdr);
        covInterval.write(hdr);

        if(outFmt.report_sjs){
            write_sjs(hdr);
        }
    }
    else {
        
    }
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
}// <------------------ main() end -----

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;all;no_cov;anCDVho:s:N:Q:");
    args.printError(USAGE, true);
    if (args.getOpt('h') || args.getOpt("help") || args.startNonOpt()==0) {
        GMessage(USAGE);
        exit(1);
    }

    if (args.getOpt('h') || args.getOpt("help")) {
        fprintf(stdout,"%s",USAGE);
        exit(0);
    }

    outFmt.report_all=(args.getOpt("all")!=NULL || args.getOpt('a')!=NULL);
    outFmt.no_cov=(args.getOpt("no_cov")!=NULL || args.getOpt('n')!=NULL);

    GStr max_nh_str=args.getOpt('N');
    if (!max_nh_str.is_empty()) {
        filters.max_nh=max_nh_str.asInt();
    }
    GStr min_qual_str=args.getOpt('Q');
    if (!min_qual_str.is_empty()) {
        filters.min_qual=min_qual_str.asInt();
    }

    debugMode=(args.getOpt("debug")!=NULL || args.getOpt('D')!=NULL);
    verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);
    if (args.getOpt('C'))
        alesCov=false ;
    if (args.getOpt("version")) {
        fprintf(stdout,"%s\n", VERSION);
        exit(0);
    }
    //verbose=(args.getOpt('v')!=NULL);
    if (verbose) {
        fprintf(stderr, "Running TieCov " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }
    outfname=args.getOpt('o');
    if (args.startNonOpt()==0) GError("Error: no alignment file given!\n");
    infname=args.nextNonOpt();

    if(args.getOpt('s')){
        outFmt.report_sjs = true;
        outFmt.sjs_fname = args.getOpt('s');
        std::string sjs_filename = args.getOpt('s');
        outFmt.sjs_combined = sjs_filename+".combined";
    }
}
