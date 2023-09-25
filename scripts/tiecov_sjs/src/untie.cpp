#include "GArgs.h"
#include "GStr.h"
#include "htslib/sam.h"
#include "GSam.h"
#include <set>
#include <tuple>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <unordered_map>

#define VERSION "0.0.5"

const char* USAGE="UnTie v" VERSION " usage:\n"
                  " untie [-b out.flt.bam] [-i index.tbd] [-p index.tbp] in.bam\n";

struct Filters{
    int max_nh = MAX_INT;
    int min_qual = -1;
} filters;

GStr ifname, pfname, bfname, infname;
FILE* boutf=NULL;

GSamWriter* outfile=NULL;

bool debugMode=false;
bool verbose=false;

void processOptions(int argc, char* argv[]);

class Index{
public:
    Index(GStr& ifname):index_fname(ifname.chars()){
        this->index_ss.open(ifname.chars(),std::ios::in | std::ios::binary);
        this->index_ss.unsetf(std::ios::skipws);
    };
    ~Index(){
        index_ss.close();
    };

    void next(uint32_t &val) {
        if (!index_ss.is_open() || !index_ss.good())
            GError("Warning: Index::next() called with no open file.\n");
        char buffer[4];
        index_ss.read(buffer,4);
        if(!index_ss) {
            GError("Error: only %d bytes could be loaded\n",index_ss.gcount());
        }
        val = ((uint8_t)buffer[0] << 24) | ((uint8_t)buffer[1] << 16) | ((uint8_t)buffer[2] << 8) | (uint8_t)buffer[3];
    }

private:
    std::string index_fname = "";
    std::fstream index_ss;
};

struct Mates{
public:
    Mates() = default;
    ~Mates() = default;

    void add_pos(uint32_t name,uint32_t pos,uint32_t offset){
        this->mit = this->mates.insert(std::make_pair(pos+offset,std::vector<uint32_t>{}));
        this->mit.first->second.push_back(name);
    }

    int pop_mate(uint32_t pos){ // given current line number - checks if the mate for this has already been added and returns the name to be assigned
        this->mit.first = this->mates.find(pos);
        if(this->mit.first != this->mates.end()){
            uint32_t name = this->mit.first->second.front();
            this->mit.first->second.erase(this->mit.first->second.begin());
            if(this->mit.first->second.empty()){
                this->mates.erase(this->mit.first);
            }
            return name;
        }
        return -1;
    }
private:
    std::unordered_map<uint32_t,std::vector<uint32_t>> mates; // key is the line_number of the mate and value is a vector of readnames to be assigned for pairing. readnames are sorted by the order in which they will appear in the output
    std::pair<std::unordered_map<uint32_t,std::vector<uint32_t>>::iterator,bool> mit;
} mates;

uint32_t readid = 0;
uint32_t first_mates_count = 0;
uint32_t both_mates_count = 0;

// >------------------ main() start -----
int main(int argc, char *argv[])  {
    processOptions(argc, argv);
    //htsFile* hts_file=hts_open(infname.chars(), "r");
    //if (hts_file==NULL)
    //   GError("Error: could not open alignment file %s \n",infname.chars());
    Index tie_idx(ifname);
    Index pair_idx(pfname);

    GSamReader samreader(infname.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

    GSamFileType oftype=(bfname=="-") ?
                        GSamFile_SAM : GSamFile_BAM;
    outfile=new GSamWriter(bfname,samreader.header(),oftype);

    GSamRecord brec;
    while (samreader.next(brec)) {
        if(brec.get_b()->core.flag & 0x100){
            GError("Cannot untie files with supplementary alignments\n");
        }
        // TODO: if(passes_options()) continue;
        // get the next position of the index
        uint32_t dupcount;
        tie_idx.next(dupcount);

        for(int i=0;i<dupcount;i++){
            int bm = both_mates_count;
            int fm = first_mates_count;
            if(brec.get_b()->core.flag & 0x1){
                // check if the current position already exists in the stored mates
                int name = mates.pop_mate(both_mates_count);
                if(name==-1){ // mate not found
                    // get next pair information
                    uint32_t mate_offset;
                    pair_idx.next(mate_offset);
                    name = first_mates_count;
                    mates.add_pos(first_mates_count,both_mates_count,mate_offset);
                    first_mates_count++;
                }
                both_mates_count++;
                name = std::stoi("2" + std::to_string(name)); // adds 2 to tell it's a paired read
                brec.replace_qname(name);
                outfile->write(&brec);
            }
            else{
                int name = std::stoi("1" + std::to_string(readid)); // adds 1 to tell it's a single-end read
                brec.replace_qname(name);
                readid++;
                outfile->write(&brec);
            }
        }
    } //while GSamRecord emitted
    delete outfile;
}// <------------------ main() end -----

// TODO: need to read in the options from the BAM header to guide the re-pairing of the file
//       (make sure the same policies are applied to count the line numbers)
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
    ifname=args.getOpt('i');
    if(ifname==NULL){
        GError("No duplicity index file provided\n");
    }
    pfname=args.getOpt('p');
    if(pfname==NULL){
        GError("No pair index file provided\n");
    }
    if (args.startNonOpt()!=1) GError("Error: no alignment file given!\n");
    infname=args.nextNonOpt();
}
