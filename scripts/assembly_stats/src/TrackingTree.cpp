//
// Created by sparrow on 3/4/20.
//

#include <cmath>
#include "TrackingTree.h"

void TrackingTree::load(){
    // load tissue_tracking
    std::cout<<"<<<loading tissue tracking info"<<std::endl;
    load_tt();

    // add all tracking
    std::cout<<"<<<loading all tracking info"<<std::endl;
    load_at();

    // load transcript IDs for ALL that are true
    std::cout<<"<<<setting types for transcripts"<<std::endl;
    if(this->tmap_fname_set){
        load_types_tmap();
    }
    else{
        load_types_gtf();
    }
}

std::string TrackingTree::get_tissue_name(std::string fname){
    std::size_t found_sample = fname.rfind("/");
    if (found_sample!=std::string::npos){

        std::size_t found_ext = fname.rfind(".");
        if (found_ext!=std::string::npos){
            return fname.substr(found_sample+1,found_ext-found_sample-1);
        }
        else{
            std::cerr<<"no extension separator found in the path: "<<fname<<std::endl;
            exit(-1);
        }
    }
    else{
        std::cerr<<"no sample separator found in the path: "<<fname<<std::endl;
        exit(-1);
    }
}

void TrackingTree::break_tracking(std::vector<TR>& trs,std::string& tline){
    std::string new_tid,new_loc;
    // get tid by finding first tab
    std::size_t found_tid = tline.find("\t");
    if (found_tid!=std::string::npos){
        new_tid = tline.substr(0,found_tid);
    }
    else{
        std::cerr<<"ERROR #1"<<std::endl;
        exit(-1);
    }
    // get locus
    std::size_t found_gid = tline.find("\t",found_tid+1);
    if(found_gid!=std::string::npos){
        new_loc = tline.substr(found_tid+1,found_gid-found_tid - 1);
    }
    // get individual transcripts
    std::size_t found_name = tline.find("\t",found_gid+1);
    std::size_t found_code = tline.find("\t",found_name+1);

    std::stringstream *txs_stream = new std::stringstream(tline.substr(found_code+1,tline.size()-found_code));
    std::string txs;
    while(std::getline(*txs_stream,txs,'\t')) {
        if (std::strcmp(txs.c_str(), "-") == 0) { // transcript does not exist in a sample
            continue;
        }

        // get qj
        std::string qj = "";
        std::size_t found_qj = txs.find(":");
        if(found_qj!=std::string::npos){
            qj = txs.substr(0,found_qj);
        }
        else {
            std::cerr << "ERROR #2" << std::endl;
            exit(-1);
        }

        std::stringstream *tx_stream = new std::stringstream(txs.substr(found_qj+1,txs.size()-found_qj));
        std::string tx;
        while(std::getline(*tx_stream,tx,',')) {
            TR tr;
            tr.new_loc = new_loc;
            tr.new_tid = new_tid;

            tr.qj = qj;

            // get old_gid
            std::size_t found_old_gid = tx.find("|");
            if (found_old_gid != std::string::npos) {
                tr.old_loc = tx.substr(0, found_old_gid);
            } else {
                std::cerr << "ERROR #3" << std::endl;
                exit(-1);
            }
            // get old_tid
            std::size_t found_old_tid = tx.find("|", found_old_gid + 1);
            if (found_old_tid != std::string::npos) {
                tr.old_tid = tx.substr(found_old_gid + 1, found_old_tid - found_old_gid - 1);
            } else {
                std::cerr << "ERROR #4" << std::endl;
                exit(-1);
            }
            // skip num_exons
            std::size_t found_ne = tx.find("|", found_old_tid + 1);
            // get fpkm
            std::size_t found_fpkm = tx.find("|", found_ne + 1);
            if (found_fpkm != std::string::npos) {
                if ((found_fpkm - found_ne - 1) == 0) {
                    tr.fpkm = "0";
                } else {
                    tr.fpkm = tx.substr(found_ne + 1, found_fpkm - found_ne - 1);
                }
            } else {
                std::cerr << "ERROR #5" << std::endl;
                exit(-1);
            }
            // get tpm
            std::size_t found_tpm = tx.find("|", found_fpkm + 1);
            if (found_tpm != std::string::npos) {
                if ((found_fpkm - found_ne - 1) == 0) {
                    tr.tpm = 0;
                } else {
                    tr.tpm = std::stof(tx.substr(found_fpkm + 1, found_tpm - found_fpkm - 1));
                }
            } else {
                std::cerr << "ERROR #6" << std::endl;
                exit(-1);
            }
            // get cov
            std::size_t found_cov = tx.find("|", found_tpm + 1);
            if (found_cov != std::string::npos) {
                if ((found_fpkm - found_ne - 1) == 0) {
                    tr.cov = "0";
                } else {
                    tr.cov = tx.substr(found_tpm + 1, found_cov - found_tpm - 1);
                }
            } else {
                std::cerr << "ERROR #7" << std::endl;
                exit(-1);
            }
            trs.push_back(tr);
        }
        delete tx_stream;
    }
    delete txs_stream;
}

void TrackingTree::add_tracking(std::string fname){
    std::string tissue = get_tissue_name(fname);

    // parse lines
    std::ifstream tissue_track_stream;
    tissue_track_stream.open(fname.c_str(),std::ios::in);
    if (!tissue_track_stream.good()){
        std::cerr<<"@ERROR::Couldn't open tracking file: "<<fname<<std::endl;
        exit(1);
    }
    std::ios::sync_with_stdio(false);

    std::string tline,uid,ulid,sid;
    std::vector<TR> trs;
    int line_n = 0;
    while (std::getline(tissue_track_stream,tline)) {
        break_tracking(trs,tline);
        line_n+=1;
        for(auto tr : trs){
            this->tts_it = this->tts.insert(std::make_pair(tissue,std::set<std::string>{})); // add to tissue_sample map
            this->tts_it.first->second.insert(tr.qj);

            uid = tissue+"_"+tr.qj+"_"+tr.old_tid;
            ulid = tissue+"_"+tr.qj+"_"+tr.old_loc;
            sid = tissue+"_"+tr.qj; // sample id
            // now that we have this information - we need to store it in a realtionship
            this->mtm_it = this->mtm.insert(std::make_pair(uid,MT(uid,tr.old_loc,std::stof(tr.fpkm),tr.tpm,std::stof(tr.cov),tissue,sid)));
            // create an entry about the locus information here
            this->loci_sample_it = loci_sample.insert(std::make_pair(ulid,std::vector<MTM_IT>{}));
            this->loci_sample_it.first->second.push_back(this->mtm_it);
            // now to create/update entry in the tissue to sample map
            this->mttm_it = this->mttm.insert(std::make_pair(tr.new_tid,MTT(tr.new_tid,tr.new_loc,tissue)));
//            this->mttm_it.first = this->mttm.find(tr.new_tid);
//            if(this->mttm_it.first == this->mttm.end()){
//                std::cerr<<"tracking transcript "<<tr.new_tid<<" not found"<<std::endl;
//                exit(-1);
//            }
            this->mttm_it.first->second.set_locus(tr.new_loc);
            this->mttm_it.first->second.set_tissue(tissue);
            this->mttm_it.first->second.add_tx(this->mtm_it);
            // add to tissue locus map for downstream stats computation
            this->loci_tissue_it = loci_tissue.insert(std::make_pair(tissue+"_"+tr.new_loc,std::vector<MTTM_IT>{}));
            this->loci_tissue_it.first->second.push_back(this->mttm_it);
            this->tissues.insert(tissue);
        }
        trs.clear();
    }
    tissue_track_stream.close();
}

void TrackingTree::load_tt(){
    std::ifstream tissue_track_stream;
    tissue_track_stream.open(tissue_tracking_fname.c_str(),std::ios::in);
    if (!tissue_track_stream.good()){
        std::cerr<<"@ERROR::Couldn't open file with the list of sample paths: "<<tissue_tracking_fname<<std::endl;
        exit(1);
    }
    std::ios::sync_with_stdio(false);

    std::string aline,tissue;
    while (std::getline(tissue_track_stream,aline)) {
        std::cout<<aline<<std::endl;
        add_tracking(aline);
    }
    tissue_track_stream.close();
}

void TrackingTree::load_at(){
    std::ifstream all_track_stream;
    all_track_stream.open(all_tracking_fname.c_str(),std::ios::in);
    if (!all_track_stream.good()){
        std::cerr<<"@ERROR::Couldn't open all tracking file: "<<all_tracking_fname<<std::endl;
        exit(1);
    }
    std::ios::sync_with_stdio(false);

    std::string aline;
    std::vector<TR> trs;
    while (std::getline(all_track_stream,aline)) {
        break_tracking(trs,aline);
        for(auto tr : trs){
//            if(std::strcmp(tr.new_tid.c_str(),"ALL_00000747")==0){
//                std::cout<<"found"<<std::endl;
//            }
//            this->mttm_it = this->mttm.insert(std::make_pair(tr.old_tid,MTT(tr.old_tid)));
            this->mttm_it.first = this->mttm.find(tr.old_tid);
            if(this->mttm_it.first == this->mttm.end()){
                std::cerr<<"tracking transcript "<<tr.old_tid<<" not found"<<std::endl;
                exit(-1);
            }

            this->mttm_it.first->second.set_all_tid(tr.new_tid);
            this->mttm_it.first->second.set_all_loc(tr.new_loc);

            // add entry to the map linking all transcripts and
            this->mattm_it = mattm.insert(std::make_pair(tr.new_tid,MATT(tr.new_tid,tr.new_loc)));
            this->mattm_it.first->second.add_tx(this->mttm_it);

            // also add to locus map linking all loci to all transcripts
            this->loci_it = loci.insert(std::make_pair(tr.new_loc,std::make_pair(std::vector<MATTM_IT>{},-1)));
            this->loci_it.first->second.first.push_back(this->mattm_it);
        }
        trs.clear();
    }
    all_track_stream.close();
}

void TrackingTree::load_types_gtf(){
    FILE* gff_file = fopen(this->all_gtf_fname.c_str(), "r");
    if (gff_file == nullptr){
        std::cerr << "@ERROR::Couldn't open annotation: " << this->all_gtf_fname<< std::endl;
        exit(1);
    }

    GffReader gffReader;
    gffReader.init(gff_file,true);
    gffReader.readAll(true);

    GffObj *p_gffObj;

    std::set<std::string> conflicts; // conflicting loci = those having both known and unknown codes

    for (int i = 0; i < gffReader.gflst.Count(); ++i){
        p_gffObj = gffReader.gflst.Get(i);
        if (p_gffObj->isDiscarded() || p_gffObj->exons.Count()==0){
            continue;
        }
        p_gffObj->store_elen();

        // find the transcript in the index
        this->mattm_it.first = this->mattm.find(p_gffObj->getID());
        if(this->mattm_it.first==this->mattm.end()){
//            std::cerr<<"transcript not found: "<<p_gffObj->getID()<<" in: "<<all_gtf_fname<<std::endl;
//            exit(-1);
            continue;
        }

        // now need to get classification code
        int class_code_attid = p_gffObj->names->attrs.getId("class_code");
        if(class_code_attid==-1){ // atribute not found
            std::cerr<<"no attribute class_code found in file"<<std::endl;
            exit(-1);
        }
        else{
            if(!p_gffObj->attrs->getAttr(class_code_attid)){ // atribute not found - something is wrong...
                std::cerr<<"attribute class_code not found in current transcript: "<<p_gffObj->getID()<<std::endl;
                exit(-1);
            }

            std::string class_code = p_gffObj->attrs->getAttr(class_code_attid);

            this->mattm_it.first->second.set_type(class_code[0]);
            this->mattm_it.first->second.set_elen(p_gffObj->elen);
            this->mattm_it.first->second.set_num_exons(p_gffObj->exons.Count());

            // now also find the corresponding locus in the index and set its type according to what's been precomputed
            this->loci_it.first = loci.find(p_gffObj->getGeneID());
            if(this->loci_it.first != loci.end()){
//            if(class_code == "u" || class_code == "y" || class_code == "s" || class_code == "x" || class_code == "r" || class_code == "p"){ // intergenic locus
                if(class_code=="u"){
                    this->loci_it.first->second.second = 0;
                }
                else{
                    this->loci_it.first->second.second = 1;
                }
            }
            else{ // danger zone
                std::cerr<<"locus: "<<p_gffObj->getGeneID()<<" from type is not found in index"<<std::endl;
                exit(-1);
            }
        }
    }
    std::cout<<"number of conflicting loci: "<<conflicts.size()<<std::endl;
}

void TrackingTree::load_types_tmap(){
    std::ifstream tmap_stream;
    tmap_stream.open(tmap_fname.c_str(),std::ios::in);
    if (!tmap_stream.good()){
        std::cerr<<"@ERROR::Couldn't open the tmap file: "<<tmap_fname<<std::endl;
        exit(1);
    }
    std::ios::sync_with_stdio(false);

    std::string aline;
    std::getline(tmap_stream,aline); // skip first line with the header
    while (std::getline(tmap_stream,aline)) {
        int col_no = 0;
        std::stringstream *col_stream = new std::stringstream(aline);
        std::string col;
        std::string tid,gid,class_code;
        while(std::getline(*col_stream,col,'\t')) {
            if(col_no == 2){
                class_code = col;
            }
            if(col_no == 3){
                gid = col;
            }
            if(col_no == 4){
                tid = col;
            }
            col_no++;
        }
        delete col_stream;


        this->mattm_it.first = this->mattm.find(tid);
        if(this->mattm_it.first==this->mattm.end()){
//            std::cerr<<"transcript not found: "<<tid<<" in: "<<all_gtf_fname<<std::endl;
//            exit(-1);
            continue;
        }

        this->mattm_it.first->second.set_type(class_code[0]);

        // now also find the corresponding locus in the index and set its type according to what's been precomputed
        this->loci_it.first = loci.find(gid);
        if(this->loci_it.first != loci.end()){
//            if(class_code == "u" || class_code == "y" || class_code == "s" || class_code == "x" || class_code == "r" || class_code == "p"){ // intergenic locus
            if(class_code=="u"){
                this->loci_it.first->second.second = 0;
            }
            else{
                this->loci_it.first->second.second = 1;
            }
        }
        else{ // danger zone
            std::cerr<<"locus: "<<gid<<" from type is not found in index"<<std::endl;
            exit(-1);
        }
    }
    tmap_stream.close();
}

void TrackingTree::_convert_gtf(std::string gtf_fname,std::string tissue,std::string base_out_fname,std::map<std::string,std::pair<std::vector<std::string>,bool>>& dups){
    std::string gtf_out_fname = base_out_fname+"."+tissue+".gtf";
    std::ofstream gtf_out_ss(gtf_out_fname.c_str());

    FILE* gff_file = fopen(gtf_fname.c_str(), "r");
    if (gff_file == nullptr){
        std::cerr << "@ERROR::Couldn't open annotation: " <<gtf_fname<<" for conversion"<<std::endl;
        exit(1);
    }

    GffReader gffReader;
    gffReader.init(gff_file,true);
    gffReader.readAll(true);

    GffObj *p_gffObj;
    GffExon *exon;

    std::set<std::string> conflicts; // conflicting loci = those having both known and unknown codes

    std::map<std::string,std::pair<std::vector<std::string>,bool>>::iterator dit;

    for (int i = 0; i < gffReader.gflst.Count(); ++i){
        p_gffObj = gffReader.gflst.Get(i);
        if (p_gffObj->isDiscarded() || p_gffObj->exons.Count()==0){
            continue;
        }
        p_gffObj->store_elen();

        // find the transcript in the index
        this->mttm_it.first = this->mttm.find(p_gffObj->getID());
        if(this->mttm_it.first==this->mttm.end()){
            std::cerr<<"transcript not found: "<<p_gffObj->getID()<<" in: "<<gtf_fname<<std::endl;
            exit(-1);
        }

        // set new tid, novel vs known, tpms, etc
        if(this->mttm_it.first->second.all_tid.empty()){
//            std::cerr<<"all tid not set"<<std::endl;
//            exit(-1);
            continue;
        }

        // check if duplicates exist
        dit = dups.find(this->mttm_it.first->second.all_tid);
        if(dit != dups.end()){
            if(dit->second.second){ // duplicates have been previously observed and written - nothing to do
                continue;
            }
            dit->second.second = true; // set that duplicate merged
        }

        gtf_out_ss<<p_gffObj->getGSeqName()<<"\t"
                  <<"astats"<<"\t"
                  <<"transcript"<<"\t"
                  <<p_gffObj->start<<"\t"
                  <<p_gffObj->end<<"\t"
                  <<"."<<"\t"
                  <<p_gffObj->strand<<"\t"
                  <<"."<<"\t"
                  <<"transcript_id \""<<this->mttm_it.first->second.all_tid<<"\";"
                  <<" gene_id \""<<this->mttm_it.first->second.all_loc<<"\";"
                  <<" type \""<<(char)this->mttm_it.first->second.type<<"\";";
        // add all tpms for the transcript
        gtf_out_ss<<" tpms \"";
        if(dit != dups.end()){
            for(auto& dup : dit->second.first){
                this->mttm_it_tmp.first = this->mttm.find(dup);
                if(this->mttm_it_tmp.first==this->mttm.end()){
                    std::cerr<<"transcript not found: "<<dup<<" in duplicates"<<std::endl;
                    exit(-1);
                }
                if(this->mttm_it_tmp.first->second.all_tid.empty()){
                    continue;
                }
                for(auto& tx : this->mttm_it_tmp.first->second.txs){
                    gtf_out_ss<<tx.first->second.tpm<<",";
                }
            }
            gtf_out_ss.seekp(-1, std::ios_base::end);
            gtf_out_ss<<"\";";

            gtf_out_ss<<" original_tids \"";
            for(auto& dup : dit->second.first){
                this->mttm_it_tmp.first = this->mttm.find(dup);
                if(this->mttm_it_tmp.first==this->mttm.end()){
                    std::cerr<<"transcript not found: "<<dup<<" in duplicates"<<std::endl;
                    exit(-1);
                }
                if(this->mttm_it_tmp.first->second.all_tid.empty()){
                    continue;
                }
                gtf_out_ss<<dup<<",";
            }
            gtf_out_ss.seekp(-1, std::ios_base::end);
            gtf_out_ss<<"\";";
        }
        else{
            for(auto& tx : this->mttm_it.first->second.txs){
                gtf_out_ss<<tx.first->second.tpm<<",";
            }
            gtf_out_ss.seekp(-1, std::ios_base::end);
            gtf_out_ss<<"\";";
            gtf_out_ss<<" original_tids \""<<p_gffObj->getID()<<"\";";
        }

        gtf_out_ss<<std::endl;

        // write exons
        for(int eidx=0;eidx<p_gffObj->exons.Count();eidx++){
            exon = p_gffObj->exons.Get(eidx);
            gtf_out_ss<<p_gffObj->getGSeqName()<<"\t"
                      <<"astats"<<"\t"
                      <<"exon"<<"\t"
                      <<exon->start<<"\t"
                      <<exon->end<<"\t"
                      <<"."<<"\t"
                      <<p_gffObj->strand<<"\t"
                      <<"."<<"\t"
                      <<"transcript_id \""<<this->mttm_it.first->second.all_tid<<"\";"
                      <<"gene_id \""<<this->mttm_it.first->second.all_loc<<"\";"
                      <<" type \""<<(char)this->mttm_it.first->second.type<<"\";"
                      <<std::endl;
        }
    }

    gtf_out_ss.close();
}

void TrackingTree::_check_dups(std::map<std::string,std::pair<std::vector<std::string>,bool>>& dups,std::string& gtf_fname){
    FILE* gff_file = fopen(gtf_fname.c_str(), "r");
    if (gff_file == nullptr){
        std::cerr << "@ERROR::Couldn't open annotation: " <<gtf_fname<<" for conversion"<<std::endl;
        exit(1);
    }

    GffReader gffReader;
    gffReader.init(gff_file,true);
    gffReader.readAll(true);

    GffObj *p_gffObj;

    std::pair<std::map<std::string,std::pair<std::vector<std::string>,bool>>::iterator,bool> dit;

    for (int i = 0; i < gffReader.gflst.Count(); ++i) {
        p_gffObj = gffReader.gflst.Get(i);
        if (p_gffObj->isDiscarded() || p_gffObj->exons.Count() == 0) {
            continue;
        }
        p_gffObj->store_elen();

        // find the transcript in the index
        this->mttm_it.first = this->mttm.find(p_gffObj->getID());
        if (this->mttm_it.first == this->mttm.end()) {
            std::cerr << "transcript not found: " << p_gffObj->getID() << " in: " << gtf_fname << std::endl;
            exit(-1);
        }

        // set new tid, novel vs known, tpms, etc
        if (this->mttm_it.first->second.all_tid.empty()) {
            continue;
        }

        dit = dups.insert(std::make_pair(this->mttm_it.first->second.all_tid,std::make_pair(std::vector<std::string>{},false)));
        dit.first->second.first.push_back(this->mttm_it.first->second.tid);
    }

    dit.first = dups.begin();
    while (dit.first != dups.end()) {
        if (dit.first->second.first.size()<2) {
            dit.first = dups.erase(dit.first);
        } else {
            ++dit.first;
        }
    }
    std::cout<<dups.size()<<std::endl;
}

void TrackingTree::convert(std::string gtf_list_fname,std::string base_out_fname){
    std::ifstream gtf_list_stream;
    gtf_list_stream.open(gtf_list_fname.c_str(),std::ios::in);
    if (!gtf_list_stream.good()){
        std::cerr<<"@ERROR::Couldn't open file with the list of tissue GTF paths: "<<gtf_list_fname<<std::endl;
        exit(1);
    }
    std::ios::sync_with_stdio(false);

    std::string aline;
    while (std::getline(gtf_list_stream,aline)) {
        std::string gtf_fname = aline.substr(0, aline.find("\t"));
        aline.erase(0,aline.find("\t")+1);
        std::string tissue = aline.substr(0,aline.find("\n"));

        if(tissue.empty()){
            std::cerr<<"tissue not present"<<std::endl;
            exit(-1);
        }
        // before converting - check if any transcripts result in the same ALL transcript - if found - not so that the records can be merged
        std::map<std::string,std::pair<std::vector<std::string>,bool>> dups;
        _check_dups(dups,gtf_fname);
        _convert_gtf(gtf_fname,tissue,base_out_fname,dups);
    }
    gtf_list_stream.close();
}

// get all tpms for each ALL transcript across the dataset - separated by tissue (/) then by sample (;)
void TrackingTree::get_tx_tpms(std::string base_out_fname){
    std::cout<<"aggregating all TPMs for each transcript"<<std::endl;

    std::string tx_tpms_fname = base_out_fname+".tx_tpms";
    std::ofstream tx_tpms_ss(tx_tpms_fname.c_str());
    tx_tpms_ss<<"tid,code,tpms"<<std::endl;

    for(auto &atx : this->mattm){ // iterate over ALL transcripts
        int tx_type = atx.second.type;
        tx_tpms_ss<<atx.first<<","<<char(tx_type)<<",";
        std::string tx_locus = atx.second.locus;
        for(auto& ttx : atx.second.txs){ // tissue level
            for(auto& stx : ttx.first->second.txs){ // sample level
                tx_tpms_ss<<stx.first->second.tpm<<";";
            }
            tx_tpms_ss.seekp(-1, std::ios_base::end);
            tx_tpms_ss<<"|";
        }
        tx_tpms_ss.seekp(-1, std::ios_base::end);
        tx_tpms_ss<<std::endl;
    }

    tx_tpms_ss.close();
}

// get all tpms for each ALL locus across the dataset - separated by transcript (:) then by tissue (/) then by sample (;)
void TrackingTree::get_loc_tpms(std::string base_out_fname){
    std::cout<<"aggregating all TPMs for each locus"<<std::endl;

    std::string loc_tpms_fname = base_out_fname+".loc_tpms";
    std::ofstream loc_tpms_ss(loc_tpms_fname.c_str());
    loc_tpms_ss<<"lid,real,tpms"<<std::endl;

    for(auto& loc : this->loci){ // begin iterating over all loci
        loc_tpms_ss<<loc.first<<","<<loc.second.second<<",";
        for(auto& atx : loc.second.first){
            for(auto& ttx : atx.first->second.txs){ // tissue level
                for(auto& stx : ttx.first->second.txs){ // sample level
                    loc_tpms_ss<<stx.first->second.tpm<<";";
                }
                loc_tpms_ss.seekp(-1, std::ios_base::end);
                loc_tpms_ss<<"|";
            }
            loc_tpms_ss.seekp(-1, std::ios_base::end);
            loc_tpms_ss<<":";
        }
        loc_tpms_ss.seekp(-1, std::ios_base::end);
        loc_tpms_ss<<std::endl;
    }

    loc_tpms_ss.close();
}

// get number of tissues and samples and samples per tissue each transcript occurs in
void TrackingTree::get_num_tissue_sample_per_tx(TrackingStats& stats){
    std::cout<<"Computing number of transcript occurrences"<<std::endl;

    for(auto &atx : this->mattm){ // iterate over ALL transcripts
        int tx_type = atx.second.type;
        stats.num_tissue_sample_per_tx.push_back(std::make_tuple(atx.first,tx_type,atx.second.elen,atx.second.num_exons,std::vector<std::vector<float> >{},std::vector<std::vector<float> >{},atx.second.locus));
        for(auto& ttx : atx.second.txs){ // tissue level
            std::get<4>(stats.num_tissue_sample_per_tx.back()).push_back(std::vector<float>{});
            std::get<5>(stats.num_tissue_sample_per_tx.back()).push_back(std::vector<float>{});
            for(auto& stx : ttx.first->second.txs){ // sample level
                if(stx.first->second.tpm>this->tpm_thresh){
                    std::get<4>(stats.num_tissue_sample_per_tx.back()).back().push_back(stx.first->second.tpm);
                }
                else{
                    std::get<5>(stats.num_tissue_sample_per_tx.back()).back().push_back(stx.first->second.tpm);
                }
            }
            // now check if any passed the threshold - if empty - remove tissue
            if(std::get<4>(stats.num_tissue_sample_per_tx.back()).back().empty()){
                std::get<4>(stats.num_tissue_sample_per_tx.back()).pop_back();
            }
            if(std::get<5>(stats.num_tissue_sample_per_tx.back()).back().empty()){
                std::get<5>(stats.num_tissue_sample_per_tx.back()).pop_back();
            }
        }
    }
}

// get number of tissues and samples and samples per tissue each locus occurs in
void TrackingTree::get_num_tissue_sample_per_loc(TrackingStats& stats){
    std::cout<<"Computing number of locus occurrences"<<std::endl;

    for(auto& loc : this->loci){ // begin iterating over all loci
        stats.ntspl_it = stats.num_tissue_sample_per_loc.insert(std::make_pair(loc.first,std::make_tuple(loc.second.second,std::map<std::string,std::set<std::string>>{})));
        if(!stats.ntspl_it.second){
            std::cerr<<"locus already exists: "<<loc.first<<std::endl;
            exit(-1);
        }
        std::pair<std::map<std::string,std::set<std::string>>::iterator,bool> mit;
        for(auto& atx : loc.second.first){
            for(auto& ttx : atx.first->second.txs){ // tissue level
                mit = std::get<1>(stats.ntspl_it.first->second).insert(std::make_pair(ttx.first->second.tissue,std::set<std::string>{}));
                for(auto& stx : ttx.first->second.txs){ // sample level
                    mit.first->second.insert(stx.first->second.sample);
                }
            }
        }
    }
}