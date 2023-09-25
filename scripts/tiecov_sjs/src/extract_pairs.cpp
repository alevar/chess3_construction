#include "GArgs.h"
#include "GStr.h"
#include "GVec.hh"
#include "htslib/sam.h"
#include "GSam.h"
#include <set>
#include <functional>
#include <iostream>
#include <map>
#include <tuple>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cstring>

#define VERSION "0.0.3"

const char* USAGE="ExtractPairs v" VERSION " usage:\n"
                  " extract_pairs [-b out.flt.bam] in.bam\n";

GStr bfname, infname;

bool debugMode=false;
bool verbose=false;

void processOptions(int argc, char* argv[]);

void write_idx(uint32_t next_pos,std::fstream& fp){
    char bytes_next[4];
    bytes_next[0] = (next_pos >> 24) & 0xFF;
    bytes_next[1] = (next_pos >> 16) & 0xFF;
    bytes_next[2] = (next_pos >> 8) & 0xFF;
    bytes_next[3] = next_pos & 0xFF;
    fp.write(bytes_next,4);
}

struct PairedMates{
public:
    PairedMates() = default;
    ~PairedMates() = default;

    uint32_t add_read(bam1_t* rec,uint32_t pos){ // return -1 if no mate was found. If mate found - returns the position of the first mate
        this->fit = this->firsts.insert(std::make_pair(bam_get_qname(rec),pos));
        if(!this->fit.second){ // second mate found
            uint32_t return_pos = this->fit.first->second;
            this->firsts.erase(this->fit.first);
            return return_pos;
        }
        return -1;
    }
private:
    std::map<std::string,uint32_t> firsts; // readname to position of the first mate
    std::pair<std::map<std::string,uint32_t>::iterator,bool> fit;
} mates;

std::map<uint32_t,uint32_t> sorted_offsets;

// >------------------ main() start -----
int main(int argc, char *argv[])  {
    processOptions(argc, argv);

    GSamReader samreader(infname.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

    GSamFileType oftype=(bfname=="-") ?
                        GSamFile_SAM : GSamFile_BAM;

    std::fstream out_stream(bfname.chars(),std::ios::out | std::ios::binary);

    GSamRecord brec;
    int rec_pos = 0;
    while (samreader.next(brec)) {
        if(!brec.get_b()->core.flag & 0x1){ // if not paired - just write offset of 0
            sorted_offsets.insert(std::make_pair(rec_pos,0));
            rec_pos++;
            continue;
        }
        if(brec.get_b()->core.flag & 0x100){ // not a primary alignment - can continue
            continue;
        }
        else{
            // primary alignment - need to find the mate
            uint32_t first_mate_pos = mates.add_read(brec.get_b(),rec_pos);
            if(first_mate_pos != -1){ // read can be written
                sorted_offsets.insert(std::make_pair(first_mate_pos,rec_pos-first_mate_pos));
            }
            else{
                rec_pos++;
            }
        }

    } //while GSamRecord emitted
    for(auto& offset : sorted_offsets){
        write_idx(offset.second,out_stream);
    }

    out_stream.close();


}// <------------------ main() end -----

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;DVhb:");
    args.printError(USAGE, true);
    if (args.getOpt('h') || args.getOpt("help") || args.startNonOpt()==0) {
        GMessage(USAGE);
        exit(1);
    }

    if (args.getOpt('h') || args.getOpt("help")) {
        fprintf(stdout,"%s",USAGE);
        exit(0);
    }

    debugMode=(args.getOpt("debug")!=NULL || args.getOpt('D')!=NULL);
    verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);

    if (args.getOpt("version")) {
        fprintf(stdout,"%s\n", VERSION);
        exit(0);
    }
    //verbose=(args.getOpt('v')!=NULL);
    if (verbose) {
        fprintf(stderr, "Running UnTie " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }
    bfname=args.getOpt('b');
    if (args.startNonOpt()!=1) GError("Error: no alignment file given!\n");
    infname=args.nextNonOpt();
}
