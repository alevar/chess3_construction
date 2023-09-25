#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "arg_parse.h"
#include "TrackingStats.h"
#include "TrackingTree.h"

enum Opt {TISSUES   = 't',
          GTF       = 'g',
          ALL       = 'a',
          OUTPUT    = 'o',
          TMAP      = 'm',
          THRESH    = 'r',
          CONVERT   = 'c'};

int main(int argc, char** argv) {

    ArgParse args("assembly_stats");
    args.add_string(Opt::TISSUES, "tissues", "", "File containing a list of paths to tracking files per tissue", true);
    args.add_string(Opt::GTF, "gtf", "", "File containing the merged gtf of all tissues. The file is expected to have classification codes", true);
    args.add_string(Opt::ALL, "all", "", "Path to the tracking for all samples across all tissues", true);
    args.add_string(Opt::OUTPUT, "output", "", "Basename for the output files", true);
    args.add_string(Opt::TMAP,"tmap","","File with the classification codes in the tmap format from gffcompare for the top level of assembly merging",false);
    args.add_double(Opt::THRESH,"thresh",0,"Minimum TPM of observations to count",false);
    args.add_string(Opt::CONVERT,"convert","","file with a list of path to gtfs of tissues to convert to IDs from ALL",true);

    // we probably only the GTFs to get effective lengths of the transcripts, since we can get coverage and TPM information from the tissue tracking files

    if (argc <= 1 || strcmp(argv[1], "--help") == 0) {
        std::cerr << args.get_help() << std::endl;
        exit(1);
    }

    args.parse_args(argc, argv);

    // first create the execution string
    std::string cl = "assembly_stats ";
    for (int i = 0; i < argc; i++) {
        if (i == 0) {
            cl += argv[i];
        } else {
            cl += " ";
            cl += argv[i];
        }
    }

    // now can load the tissue to sample transcript relationship table
    TrackingTree tt(args.get_string(Opt::TISSUES), args.get_string(Opt::ALL), args.get_string(Opt::GTF),args.get_double(Opt::THRESH));
    if(args.is_set(Opt::TMAP)){
        tt.set_tmap_fname(args.get_string(Opt::TMAP));
    }
    tt.load();
    std::cerr<<"converting"<<std::endl;
    tt.convert(args.get_string(Opt::CONVERT),args.get_string(Opt::OUTPUT));

//    TrackingStats stats;
//    tt.get_stats(stats,args.get_string(Opt::OUTPUT));
//    stats.save_stats(args.get_string(Opt::OUTPUT));

    return 0;
}