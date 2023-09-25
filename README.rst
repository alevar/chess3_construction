# Collection of scripts and tools used in the processing of raw data and partial construction of the CHESS Catalogue of the Human Genes and Transcripts.


The following collection is an extends the description of the methods provided in the corresponding publication. However, some parts of the protocol are described exclusively in the manuscript (eg. alignment, assembly, reference selection, etc). Furthermore, while this protocol aims to provide an account of the automated steps of the CHESS catalogue construction, multiple steps were manually curated and are not covered by the present repository, rather are reflected in the manuscript itself. Please do not hesistate to reach out to us regarding any parts of the protocol and we will gladly assist with any questions to the best of our ability.

Contents:
1. pipeline
    - step1_process_assembly.ipynb and process_assembly.py - This notebook executes the post-assembly scripts and prepares all data for the sucessive evaluation and analysis. The notebook also fetches all required referece data, and performs comparisons between assembled and known transcripts. Sequencing and assembly data are used to propagate names of transcripts across the assembly hierarchy (Dataset -> tissues -> samples) and create summaries with tierus and tiecov. The notebook contains details of the execution for each step.
    - step2_select_transcripts_ml.ipynb - This notebook runs the standalone nitron utility which trains a LightGBM model at each tissue level, assigns splice junction confidence and selects transcripts. The notebook terminates after creating tissue-level draft catalogues along with a merged catalogue for the entire dataset.
    - step3_refine.ipynb - This notebook executes various scripts, cleaning up the draft catalogue. Namely, the notebook performs a) replacement of the mitochondrial genes with established models; b) ORF prediction and CDS assignment using ORFanage; c) single-exon extraction; d) 3' and 5' re-assignment based on available evidence; e) re-introducing missing validated known isoforms which were not observed in the data but are believed to be valid; f) readthrough transcript detection; and g) preliminary gene and transcript ID assignment based on the previous versions of the dataset. Multiple steps in the notebook do not add or remove features from the dataset, rather flagging them in the attributes for further deliberation downstream in the analysis. 
    - step4_remove_tags.ipynb - This notebook performs cleanup of the features in the draft dataset based on the attributes computed up to this point in the analysis.
    - step5_namer.ipynb - Performs preliminary gene and transcript ID assignment..
    - step6_add_missing.ipynb - Performs recovery of the missing genes and transcripts from other available datasets.
    - sashimi.ipynb - Builds sashimi figures for selected genes.
2. Reviews - Collection of notebooks adressing reviewer comments and suggestions.
3. Scripts - Various scripts mentioned and used in the notebooks.
