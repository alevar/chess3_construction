#include "GArgs.h"
#include "GStr.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>

#define VERSION "0.0.5"

// Merge two indices together
const char* USAGE="TieMe v" VERSION " usage:\n"
                  " tieme -o output.tbd index1.tbd index2.tbd\n";

GStr outfname,infname1,infname2;

bool debugMode=false;
bool verbose=false;

void processOptions(int argc, char* argv[]);

class Index{
public:
    Index() = default;
    Index(GStr& ifname):index_fname(ifname.chars()){
        this->index_ss.open(ifname.chars(),std::ios::in | std::ios::binary);
        this->index_ss.unsetf(std::ios::skipws);
    };
    ~Index(){
        index_ss.close();
    };

    void set_infname(GStr& ifname){
        this->index_fname = ifname.chars();
        this->index_ss.open(ifname.chars(),std::ios::in | std::ios::binary);
        this->index_ss.unsetf(std::ios::skipws);
    }

    int next(uint32_t &val) {
        if (!index_ss.is_open() || !index_ss.good())
            GError("Warning: Index::next() called with no open file.\n");
        char buffer[4];
        index_ss.read(buffer,4);
        if(!index_ss) {
            return 0;
        }
        val = ((uint8_t)buffer[0] << 24) | ((uint8_t)buffer[1] << 16) | ((uint8_t)buffer[2] << 8) | (uint8_t)buffer[3];
        return 1;
    }

    void merge_in(Index& ix2,std::fstream* out_fp){

    }

    void add(uint64_t dupcount,std::fstream* out_fp){
        buffer[nr*4]   = (dupcount >> 24) & 0xFF;
        buffer[(nr*4)+1] = (dupcount >> 16) & 0xFF;
        buffer[(nr*4)+2] = (dupcount >> 8) & 0xFF;
        buffer[(nr*4)+3] = dupcount & 0xFF;

        nr += 1;

        if(nr*4 == 1024*4096){
            write(out_fp);
        }
    }
    void clear(std::fstream* out_fp){
        write(out_fp);
    }

    std::ifstream::pos_type file_size(){
        return index_ss.tellg();
    }

private:
    std::string index_fname = "";
    std::fstream index_ss;

    char buffer[1024*4096];
    int nr = 0;

    void write(std::fstream* out_fp){
        out_fp->write(buffer,nr*4);
        nr=0;
    }
};

bool fp1_larger_fp2(GStr& fname1,GStr& fname2){
    std::ifstream fp1(fname1,std::ifstream::ate | std::ifstream::binary);
    std::ifstream fp2(fname2,std::ifstream::ate | std::ifstream::binary);

    std::ifstream::pos_type fp1_size = fp1.tellg();
    std::ifstream::pos_type fp2_size = fp2.tellg();

    fp1.close();
    fp2.close();

    return fp1_size>fp2_size;
}

// >------------------ main() start -----
int main(int argc, char *argv[])  {
    processOptions(argc, argv);

    Index more_collapsed_tbd;
    Index less_collapsed_tbd;
    if(fp1_larger_fp2(infname1,infname2)){
        more_collapsed_tbd.set_infname(infname1);
        less_collapsed_tbd.set_infname(infname2);
    }

    Index out_tbd;
    std::fstream* out_fp = new std::fstream();
    out_fp->open(outfname,std::ios::out | std::ios::binary);

    uint32_t more_collapsed;
    uint32_t less_collapsed;
    while(more_collapsed_tbd.next(more_collapsed)){
        if(more_collapsed==0){
            out_tbd.add(0,out_fp);
        }
        else{
            while(less_collapsed_tbd.next(less_collapsed)){
                if(less_collapsed>0){
                    out_tbd.add(less_collapsed,out_fp);
                    break;
                }
            }
        }
    }
    out_fp->close();
    delete out_fp;
}// <------------------ main() end -----

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;DVhb:i:p:N:Q:");
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
    outfname=args.getOpt('o');
    if(outfname==NULL){
        GError("No duplicity index file provided\n");
    }

    if (args.startNonOpt()!=1) GError("Error: no index fils given (must be two)!\n");
    infname1=args.nextNonOpt();
    infname2=args.nextNonOpt();
}