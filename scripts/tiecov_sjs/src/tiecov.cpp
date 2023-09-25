#include <vector>
#include <math.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>


#include "GArgs.h"
#include "GStr.h"
#include "GVec.hh"
#include "GSam.h"

#define VERSION "0.0.5"

const char* USAGE="TieCov v" VERSION " usage:\n"
                  " tiecov [-b out.flt.bam] [-s out.sample.bed] [-c out.coverage.bedgraph] [-j out.junctions.bed] in.bam\n"
                  " Otions: \n"
                  "  -b   : bam file after applying filters (-N/-Q)\n"
                  "  -s   : BED file with number of samples which contain alignments for each interval. Binning of intervals within bundles is controlled by -B\n"
                  "  -c   : BedGraph file with coverage for all mapped bases. Binning of intervals within bundles is controlled by -B\n"
                  "  -j   : BED file with coverage of all splice-junctions in the input file.\n"
                  "  -B   : default 1. If set to N - average coverage will be reported for each N bases . If set to 0 - will output average value for each complete bundle\n"
                  "  -I   : index file or a file containing a list of paths to indices. If provided, only samples corresponding to the indices will be used in the calculations\n"
                  "  -N   : maximum NH score (if available) to include when reporting coverage\n"
                  "  -Q   : minimum mapping quality to include when reporting coverage\n";

// TODO: make sure filetypes are appended correctly (like bedgraph)

struct Filters{
    int max_nh = MAX_INT;
    int min_qual = -1;
} filters;

GStr covfname, jfname, bfname, infname, sfname;
FILE* boutf=NULL;
FILE* coutf=NULL;
FILE* joutf=NULL;
FILE* soutf=NULL;

bool debugMode=false;
bool verbose=false;
int juncCount=0;

struct CJunc {
	int start, end;
	char strand;
	uint64_t dupcount;
	CJunc(int vs=0, int ve=0, char vstrand='+', uint64_t dcount=1):
	  start(vs), end(ve), strand(vstrand), dupcount(dcount) { }

	bool operator==(const CJunc& a) {
		return (strand==a.strand && start==a.start && end==a.end);
	}
	bool operator<(const CJunc& a) {
		if (start==a.start) return (end<a.end);
		else return (start<a.start);
	}

	void add(CJunc& j) {
       dupcount+=j.dupcount;
	}

	void write(FILE* f, const char* chr) {
		juncCount++;
		fprintf(f, "%s\t%d\t%d\tJUNC%08d\t%ld\t%c\n",
				chr, start-1, end, juncCount, (long)dupcount, strand);
	}
};

GArray<CJunc> junctions(64, true);

void addJunction(GSamRecord& r, int dupcount) {
	char strand = r.spliceStrand();
	if (strand!='+' && strand!='-') return;
	for (int i=1;i<r.exons.Count();i++) {
		CJunc j(r.exons[i-1].end+1, r.exons[i].start-1, strand,
				dupcount);
		int ei;
		int r=junctions.AddIfNew(j, &ei);
		if (r==-1) {//existing junction, update
			junctions[ei].add(j);
		}
	}
}

void flushJuncs(FILE* f, const char* chr) {
    for (int i=0;i<junctions.Count();i++) {
    	junctions[i].write(f,chr);
    }
    junctions.Clear();
    junctions.setCapacity(128);
}

void processOptions(int argc, char* argv[]);

void addMean(GSamRecord& r, int val, std::vector<std::pair<float,uint64_t>>& bvec, int b_start){ // for YX (number of samples) we are not interested in the sum but rather the average number of smaples that describe the position. Giving a heatmap
    bam1_t* in_rec=r.get_b();
    int pos=in_rec->core.pos; // 0-based
    b_start--; //to make it 0-based
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
                    bvec[pos-b_start].first+=(val-bvec[pos-b_start].first)/bvec[pos-b_start].second; // dynamically compute average
                    bvec[pos-b_start].second++;
                    pos++;
                }
                break;
            default:
                GError("ERROR: unknown opcode: %c from read: %s",bam_cigar_opchr(opcode),bam_get_qname(in_rec));
        }
    }
}

//b_start MUST be passed 1-based
void addCov(GSamRecord& r, int val, GVec<uint64_t>& bvec, int b_start) {
	bam1_t* in_rec=r.get_b();
    int pos=in_rec->core.pos; // 0-based
    b_start--; //to make it 0-based
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
                    bvec[pos-b_start]+=val;
                    pos++;
                }
                break;
            default:
                GError("ERROR: unknown opcode: %c from read: %s",bam_cigar_opchr(opcode),bam_get_qname(in_rec));
        }
    }
}

//b_start MUST be passed 1-based
void flushCoverage(FILE* outf,sam_hdr_t* hdr, GVec<uint64_t>& bvec,  int tid, int b_start) {
  if (tid<0 || b_start<=0) return;
  int i=0;
  b_start--; //to make it 0-based;
  while (i<bvec.Count()) {
     uint64_t ival=bvec[i];
     int j=i+1;
     while (j<bvec.Count() && ival==bvec[j]) {
    	 j++;
     }
     if (ival!=0)
       fprintf(outf, "%s\t%d\t%d\t%ld\n", hdr->target_name[tid], b_start+i, b_start+j, (long)ival);
     i=j;
  }
}

void flushCoverage(FILE* outf,sam_hdr_t* hdr, std::vector<std::pair<float,uint64_t>>& bvec,  int tid, int b_start) {
    if (tid<0 || b_start<=0) return;
    int i=0;
    b_start--; //to make it 0-based;
    while (i<bvec.size()) {
        uint64_t ival=bvec[i].second;
        float hval = bvec[i].first;
        int j=i+1;
        while (j<bvec.size() && ival==bvec[j].second) {
            j++;
        }
        if (ival!=0)
            fprintf(outf, "%s\t%d\t%d\t%ld\t%f\n", hdr->target_name[tid], b_start+i, b_start+j, (long)ival,hval);
        i=j;
    }
}

void discretize(std::vector<std::pair<float,uint64_t>>& bvec1){
    for(auto& val : bvec1){
        val.second = std::ceil(val.first);
        val.first = 0;
    }
}

void average(std::vector<uint64_t>& bvec,float thresh){
    // iterate
    // find min and max of the range of values
    // find percentage by which values can be similar to group together
    // find difference between two points
    // if within difference - keep average and compare next observation with an average
    // if not in range - write all previous values with the average and start again
    uint64_t max_val = *std::max_element(std::begin(bvec), std::end(bvec));
    uint64_t min_val = *std::min_element(std::begin(bvec), std::end(bvec));

    for(auto& val : bvec){

    }
}

void normalize(std::vector<std::pair<float,uint64_t>>& bvec,float mint, float maxt, int num_samples){ // normalizes values to a specified range
    float denom = num_samples-1;
    float mult = (maxt-mint);

    for (auto& val : bvec){
        val.first = ((val.second-1)/denom)*mult+mint;
    }
}

bool parse_pg_sample_line(std::string& line){ // returns true if is sample pg line
    std::stringstream *line_stream = new std::stringstream(line);
    std::string col;

    // make sure it's PG
    std::getline(*line_stream,col,'\t');
    if(std::strcmp(col.c_str(),"@PG")!=0){
        delete line_stream;
        return false;
    }

    // check if ID == SAMPLE
    std::getline(*line_stream,col,'\t');
    if(std::strcmp(col.c_str(),"ID:SAMPLE")!=0){
        delete line_stream;
        return false;
    }

    std::getline(*line_stream,col,'\t');
    std::stringstream *col_stream = new std::stringstream(col);
    std::string kv;
    std::getline(*col_stream,kv,':');
    if(std::strcmp(kv.c_str(),"SP")!=0){
        delete line_stream;
        delete col_stream;
        return false;
    }
    std::getline(*col_stream,kv,'\t');
    line = kv;
    delete line_stream;
    delete col_stream;
    return true;
}

void load_sample_info(sam_hdr_t* hdr,std::vector<std::string>& info){
    bool found_sample_line = false;
    int line_pos = 0;
    std::string line;
    while(true){
        kstring_t str = KS_INITIALIZE;
        if (sam_hdr_find_line_pos(hdr, "PG", line_pos, &str)!=0) {
            if(!found_sample_line){
                GError("Error: no sample lines found in header");
            }
            break;
        }
        else{
            // parse line to check if indeed SAMPLE PG
            line = std::string(str.s);
            bool ret = parse_pg_sample_line(line);
            if(ret){
                found_sample_line=true;
                info.push_back(line);
            }
            line_pos++;
        }
        ks_free(&str);
    }
}

// >------------------ main() start -----
int main(int argc, char *argv[])  {
	processOptions(argc, argv);
    //htsFile* hts_file=hts_open(infname.chars(), "r");
    //if (hts_file==NULL)
    //   GError("Error: could not open alignment file %s \n",infname.chars());
	GSamReader samreader(infname.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
    if (!covfname.is_empty()) {
       if (covfname=="-" || covfname=="stdout")
    	   coutf=stdout;
       else {
          coutf=fopen(covfname.chars(), "w");
          if (coutf==NULL) GError("Error creating file %s\n",
        		  covfname.chars());
          fprintf(coutf, "track type=bedGraph\n");
       }
    }
    if (!jfname.is_empty()) {
       joutf=fopen(jfname.chars(), "w");
       if (joutf==NULL) GError("Error creating file %s\n",
        		  jfname.chars());
       fprintf(joutf, "track name=junctions\n");
    }
    if (!sfname.is_empty()) {
        soutf=fopen(sfname.chars(), "w");
        if (soutf==NULL) GError("Error creating file %s\n",
                                sfname.chars());
        fprintf(soutf, "track type=bedGraph name=\"Sample Count Heatmap\" description=\"Sample Count Heatmap\" visibility=full graphType=\"heatmap\" color=200,100,0 altColor=0,100,200\n");
    }

    // load sample info
    std::vector<std::string> sample_info;
    load_sample_info(samreader.header(),sample_info);

    int prev_tid=-1;
    GVec<uint64_t> bcov(2048*1024);
    std::vector<std::pair<float,uint64_t>> bsam(2048*1024,{0,1}); // number of samples. 1st - current average; 2nd - total number of values
    int b_end=0; //bundle start, end (1-based)
    int b_start=0; //1 based
    GSamRecord brec;
	while (samreader.next(brec)) {
		    int nh = brec.tag_int("NH");
		    if(nh>filters.max_nh)  continue;
		    if (brec.mapq()<filters.min_qual) continue;
		    int endpos=brec.end;
		    if (brec.refId()!=prev_tid || (int)brec.start>b_end) {
		    	if (coutf) {
                    flushCoverage(coutf,samreader.header(), bcov, prev_tid, b_start);
                }
                if (soutf) {
                    discretize(bsam);
                    normalize(bsam,0.1,1.5,sample_info.size());
                    flushCoverage(soutf,samreader.header(),bsam,prev_tid,b_start);
                }
			    if (joutf) {
                    flushJuncs(joutf, samreader.refName(prev_tid));
                }
			    b_start=brec.start;
			    b_end=endpos;
			    if (coutf) {
			        bcov.setCount(0);
			        bcov.setCount(b_end-b_start+1);
			    }
			    if (soutf) {
			        bsam.clear();
                    bsam.resize(b_end-b_start+1,{0,1});
			    }
			    prev_tid=brec.refId();
		    } else { //extending current bundle
			    if (b_end<endpos) {
				    b_end=endpos;
				    bcov.setCount(b_end-b_start+1, (int)0);
				    if (soutf){
                        bsam.resize(b_end-b_start+1,{0,1});
				    }
			    }
		    }
		    int accYC = brec.tag_int("YC", 1);
		    float accYX = (float)brec.tag_int("YX", 1);
		    if(coutf){
		        addCov(brec, accYC, bcov, b_start);
		    }
		    if(soutf){
		        addMean(brec, accYX, bsam, b_start);
		    }
		    if (joutf && brec.exons.Count()>1) {
		    	addJunction(brec, accYC);
		    }
	} //while GSamRecord emitted
	if (coutf) {
       flushCoverage(coutf,samreader.header(), bcov, prev_tid, b_start);
       if (coutf!=stdout) fclose(coutf);
	}
	if (soutf) {
        discretize(bsam);
        normalize(bsam,0.1,1.5,sample_info.size());
	    flushCoverage(soutf,samreader.header(),bsam,prev_tid,b_start);
	    if (soutf!=stdout) fclose(soutf);
	}
	if (joutf) {
		flushJuncs(joutf, samreader.refName(prev_tid));
		fclose(joutf);
	}

}// <------------------ main() end -----

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;DVhc:s:j:b:N:Q:");
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
        fprintf(stderr, "Running TieCov " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }
    covfname=args.getOpt('c');
    jfname=args.getOpt('j');
    bfname=args.getOpt('b');
    sfname=args.getOpt('s');
    if (args.startNonOpt()!=1) GError("Error: no alignment file given!\n");
    infname=args.nextNonOpt();
}
