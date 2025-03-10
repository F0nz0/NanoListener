# NanoListener

<p align="center">
<img src="https://github.com/F0nz0/NanoListener/blob/master/NatMet_NanoListener_NanoSpeech_Fig1_HD_GitHub.png" align="center">
</p>

<p align="justify">
NanoListener is a small suite of Python scripts to create custom datasets for training direct-RNA modification-aware basecaller models. Raw signals stored in fast5, or converted pod5 files, are basecalled via modification-unaware basecalling models and mapped against reference sequences. Alignments are filtered collecting reads falling on regions of interest and respecting strict criteria, using NanoListener accessory scripts. 
Ionic currents needs to be re-squiggled onto the reference via f5c eventalign. NanoListener takes advantage of these re-squiggled signals to retrieve context information and extract random chunks of electric measurements trying to avoid problematic regions (in red) due to modified nucleotides. It randomly exploits as anchors, low-noise “unmodified” flanking regions to extract whole raw signals from fast5, filling in turn, re-squiggling gaps. For every extracted chunk, the deduced k-mer is marked for modified nucleotides in accordance with tables of per-read modification positions. Finally, a training dataset is returned, consisting in pairs of chunks/annotated k-mers where a padding is added to currents measurements making these, uniform in length. NanoListener needs positional information, at a per-read level, to annotate all the output k-mers. These can be obtained either using synthetic in-vitro transcribed molecules or using other orthogonal methods. 
</p>

## **Required Software**:
NanoListener uses internally (and not) some software that should be installed preferably into a new conda enviroment. \
After the activation of the conda enviroment install the following softwares:
1) Python >= 3.9
2) Samtools >= 1.21
3) Minimap2 == 2.24
4) f5c >= 1.5

## **Installation**:

1) Download the source code from GitHub repository at the url:
        
    https://github.com/F0nz0/NanoListener

2) Create a new virtual environment (it's suggested to create, use and activate a base conda environment with all the required software):

		# create a new conda environment
		conda create --name NanoListener python=3.9

		# activate the conda env
		conda activate NanoListener

		# install samtools
		conda install -c bioconda samtools >= 1.21

		# install minimap2
		conda install -c bioconda minimap2 == 2.24

		# install f5c
		conda install -c bioconda f5c == 1.5

		# create virtual environment inside Conda NanoListener environment
		python3 -m venv NanoListener_venv

4) Activate the venv:
    
		source NanoListener_venv/bin/activate

5) Upgrade pip version:
    	
		python3 -m pip install --upgrade pip
    
6) Install wheel package via pip:
    	
		pip install wheel
    
7) Install required Python packages using the requirements.txt file:
    
		python -m pip install -r requirements.txt

## **Basic Usage**:
NanoListener has 1 main scripts that is **NanoListener.py** and 4 accessory scripts used for preliminary preprocessing of input data and post-processing of output datasets. 

The first step is to obtain fast5 files (or converted pod5 via pod5 library at https://pod5-file-format.readthedocs.io/en/0.1.21/docs/tools.html#pod5-convert-to-fast5). Then these raw data have to be basecalled into fastQ files by a basecaller model (either modification -unaware or -aware, for recursive trainings) and mapped against reference sequences: we stronly suggest to use reference transcriptomes to avoid introns in interrupted genes. Alignments in BAM format then need to be filtered out from supplemenary and secondary alignments retaining only reads mapped on the forward strand. These BAM files then can be further filtered using NanoListener accessory scripts depending on the experimental design.


The first (facultative) accessory script which can be used is *create_transcript_white_list.py*. It allows to retrieve reference regions respecting a set of given criteria and takes as input the following parameters:

	python3 create_transcript_white_list.py -h

	optional arguments:
	  -h, --help            show this help message and exit
	  -b BAM_FILEPATH, --bam_filepath BAM_FILEPATH
	                        --bam_filepath: a <str> for the BAM file path to be filtered. The filtered outputfile will be at the path
	                        <BAM_basename>.filtered.bam
	  -l MIN_LEN, --min_len MIN_LEN
	                        minimum length [1000]
	  -c MIN_COV, --min_cov MIN_COV
	                        minimum coverage percentage [90]
	  -r MIN_NUMREADS, --min_numreads MIN_NUMREADS
	                        minimum number of reads per transcript [20]
	  -d MIN_MEANDEPTH, --min_meandepth MIN_MEANDEPTH
	                        minimum average depth for transcript [15]
	  -q MIN_MEANBASEQ, --min_meanbaseq MIN_MEANBASEQ
	                        minimum average base quality per transcript [15]
	  -mp MIN_MEANMAPQ, --min_meanmapq MIN_MEANMAPQ
	                        minimum average mapping quality per transcript [40]

The output of the *create_transcript_white_list.py* script is a white-list of reference regions/transcripts where NanoListener will focuse its computations. This list can be feed to the second NanoListerer accessory script, *filter_bam_file_for_training_dataset.py* which will produce a filtered BAM file containing only reads covering a minimum percentage of the corresponding reference sequence. It can be run as following:

	python3 filter_bam_file_for_training_dataset.py -h
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -f FASTA_FILEPATH, --fasta_filepath FASTA_FILEPATH
	                        --fasta_filepath: a <str> for the fasta file path.
	  -b BAM_FILEPATH, --bam_filepath BAM_FILEPATH
	                        --bam_filepath: a <str> for the BAM file path to be filtered. The filtered outputfile will be at the path
	                        <BAM_basename>.filtered.bam
	  -t PERC_THRESHOLD, --perc_threshold PERC_THRESHOLD
	                        --perc_threshold: a <float> indicating the maximum amount of unligned reads with respect reference transcript to be
	                        filtered. [0.3]
	  -r REGION_WHITELIST, --region_whitelist REGION_WHITELIST
	                        --region_whitelist: a <str> optional indicating the path for a list of regions to be extracted from BAM file. [None]
	  -o FILT_BAM_FILE_SUFFIX, --filt_bam_file_suffix FILT_BAM_FILE_SUFFIX
	                        --filt_bam_file_suffix: a <str> optional indicating the suffix to be used to write the filtered BAM output file on
	                        disk.. [None]

After that, all reads falling on the regions in the white-list and covering a given minimum percentage of the reference sequence will be retrieved and collected into the filtered output BAM file. The subsequent step is to launch the f5c eventalign program (Gamaarachchi et al., 2020, GitHub: https://github.com/hasindu2008/f5c) on this BAM file to perform the re-squiggling of raw currents against the selected reference. The f5c eventalign program (version >= 1.4 for RNA004) has to be run with the following configurations and steps: 
A) FAST5/FASTQ files indexing:

	f5c index -t 20 --iop 25 -d $FAST5 $FASTQ

B) Eventalign program:

 	f5c eventalign --iop $N_IOP -t $N_THREADS -r $FASTQ \
		-b $BAM \
		-g $GENOME \
		--rna \
		--scale-events \
		--print-read-names \
		--samples \
		--signal-index \
		--min-mapq 0 \ ##### --pore rna004 \ ADD THIS LINE FOR THE NEW ONT RNA PORE
		--summary $SUMMOUT > $EVOUT

At the end of f5c computation, the eventalign table and its summary table will be provided to **NanoListener.py** main script. It will produce an extended training datasets using re-squiggled events as anchors to extract and annotate, pairs of raw signals chunks directly from fast5 files along with their annotated k-mer. NanoListeners receives several inputs and configurations such as the required interval of chunks lengths and the pA-scale limit values to clip extracted currents chunks. Here an explaination of its functioning:


	python3 NanoListener.py -h
	usage: NanoListener.py [-h] -e EVENTALIGN_FILEPATH [-s SUMMARY_TABLE_FILEPATH] [-t THREADS_N] [-ml MIN_LENGTH_CHUNKS] [-l MAX_LENGTH_CHUNKS]
	                       [-m MOD_POS_TABLES_DIR] [-n N_READS_TO_PROCESS] [-c CLIP_OUTLIERS] -r REF_FILEPATH [-o OUTDIR_PREFIX]
	
	NanoListener. A features extractor from eventalign and fast5 files: chunks of a given length will be extracted along with related sequence
	and other useful informations. These chunks can be used to train a model from ionic currents measurements.
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -e EVENTALIGN_FILEPATH, --eventalign_filepath EVENTALIGN_FILEPATH
	                        --eventalign_filepath: a <str> indicating the fullpath for the eventalign table to be used to extract ionic currents
	                        features/chunks.
	  -s SUMMARY_TABLE_FILEPATH, --summary_table_filepath SUMMARY_TABLE_FILEPATH
	                        --summary_table_filepath: a <str> with the full_path for the summary file produced by the eventalign sofware. It is
	                        needed if you want to go back to pA values. [None]
	  -t THREADS_N, --threads_n THREADS_N
	                        --threads_n: a <int> indicating the number of threads/consumers to generate in parallel (do not use more threads than
	                        available CPUs). [1]
	  -ml MIN_LENGTH_CHUNKS, --min_length_chunks MIN_LENGTH_CHUNKS
	                        --min_length_chunks: a <int> indicating the minimum length of the chunks to extract. [300]
	  -l MAX_LENGTH_CHUNKS, --max_length_chunks MAX_LENGTH_CHUNKS
	                        --max_length_chunks: a <int> indicating the maximum length of the chunks to extract. [300]
	  -m MOD_POS_TABLES_DIR, --mod_pos_tables_dir MOD_POS_TABLES_DIR
	                        -mod_pos: a <str> indicating the fullpath for the directory containing tsv tables with the expected modified
	                        positions for each read (region pos0based alt_base). [None]
	  -n N_READS_TO_PROCESS, --n_reads_to_process N_READS_TO_PROCESS
	                        -n_reads_to_process: a <int> indicating the total number of reads to be processed. [100]
	  -c CLIP_OUTLIERS, --clip_outliers CLIP_OUTLIERS
	                        -clip_outliers: a <str> indicating the min and max pA values to clip with mean. [40-165]
	  -r REF_FILEPATH, --ref_filepath REF_FILEPATH
	                        --ref_filepath: a <str> indicating the reference fasta file used to align reads and to perform the f5c eventalign
	                        procedure.
	  -o OUTDIR_PREFIX, --outdir_prefix OUTDIR_PREFIX
	                        -outdir_prefix: a <str> indicating the ouput directory prefix to use. By default the eventalign basename will be used
	                        as prefix. [None]

To annotate k-mers of each chunks for modified nucleotides, NanoListener expects in input the full-path for a directory containing tsv files with all the known positions for alternative bases, and their symbols (1 or more modifications can be used) at a per-read level. Counterwise, it will work in un-modifided mode. The tsv files (one for every read) have to be named as *<read_id>.tsv* and respect the following content where the first 2 columns indicate the mapping coordinates (0-based) and the third one a given char for each expected modification:

	Curlcake1	17	Y
	Curlcake1	21	I
	Curlcake1	23	Y
	Curlcake1	31	Y
	Curlcake1	37	M
	Curlcake1	39	Y
	Curlcake1	45	D
	Curlcake1	49	Y

During its multi-thread computation, for every k-mer, NanoListeners will look inside these tables: if the k-mer contains any modification, it will be annotated with the selected char identifing the modified nucleotide.
The final output of NanoListener is a directory named as <BAM_FILE_NAME>.<MIN_LENGTH_CHUNKS>to<MAX_LENGTH_CHUNKS>.chunks and containg one tsv file for each thread/worker. These tsv files are composed of pairs (one line per pair) of padded chunks (to MAX_LENGTH_CHUNKS value) and additional columns containg chunks metadata and the extracted (and annotated if required) output k-mer.

These tsv files can be considered as intermediate training datasets which require additional filtering with the *nanolistener_cons_out_filt.py* accessory script. This will eliminate possible duplicates or abnormal extracted chunks/k-mers. This scripts will also consolidate all the intermediated tsv files into training and test datasets, containing in turn, all the examples derived from the analyzed sample. 

	python3 nanolistener_cons_out_filt.py -h
	usage: nanolistener_cons_out_filt.py [-h] -nld NANOLISTENER_CHUNKS_DIR -tml TARGET_MAX_LEN [-mb MOD_BASE] [-kli KMER_LEN_INTERVAL]
	                                     [-bal BALANCING_STRATEGY] [-ttr TRAIN_TEST_RATIO] [-n SAMPLE_N] [-c CONS_OUT_N] [-o OUTSUFFIX]
	
	nanolistener_cons_out_filt.py. A program to filter and balance consumers outputs produced by NanoListener before NanoSpeech model training.
	It merges all the selected consumer outputs into X and y_meta dataframes and it then splits them into X_train, X_test, y_train_meta and
	y_test_meta datasets saving these on disk at the same folder-path of the NanoListener chunks directory.
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -nld NANOLISTENER_CHUNKS_DIR, --nanolistener_chunks_dir NANOLISTENER_CHUNKS_DIR
	                        --nanolistener_chunks_dir: a <str> indicating the fullpath for the NanoListener output folder.
	  -tml TARGET_MAX_LEN, --target_max_len TARGET_MAX_LEN
	                        --target_max_len: a <int> with target max length of the chunks. It will try to fit to this value the input chunks.
	  -mb MOD_BASE, --mod_base MOD_BASE
	                        --mod_base: a <str> indicating the modified base to be searched into output kmer (at least one mod_base needed to
	                        mantain the chunks). [None]
	  -kli KMER_LEN_INTERVAL, --kmer_len_interval KMER_LEN_INTERVAL
	                        --kmer_len_interval: a <str> indicating the min and maximum allowed length of output kmers. [5-100]
	  -bal BALANCING_STRATEGY, --balancing_strategy BALANCING_STRATEGY
	                        --balancing_strategy: a <int> indicating (0: no balancing, 1: balancing on regions). [1]
	  -ttr TRAIN_TEST_RATIO, --train_test_ratio TRAIN_TEST_RATIO
	                        -train_test_ratio: a <float> indicating the ratio of chunks to be splitted into test set. [0.03]
	  -n SAMPLE_N, --sample_n SAMPLE_N
	                        -sample_n: a <int> indicating the number of chunks to sample randomly. [None]
	  -c CONS_OUT_N, --cons_out_n CONS_OUT_N
	                        -cons_out_n: a <int> indicating how many consumer outputs to considers. [None]
	  -o OUTSUFFIX, --outsuffix OUTSUFFIX
	                        -outsuffix: a <str> with a suffix to be added to the output X and y merged and filtered datasets. [None]

The *nanolistener_cons_out_filt.py* script saves output training and test datasets into the NanoListener main program output directory separating data (X_train or X_test files) from metada + k-mers (y_train, y_test):

	1) <BAM_FILE_NAME>.<MIN_LENGTH_CHUNKS>to<MAX_LENGTH_CHUNKS>.chunks/X_train_<OUTSUFFIX>.tsv
	2) <BAM_FILE_NAME>.<MIN_LENGTH_CHUNKS>to<MAX_LENGTH_CHUNKS>.chunks/X_test_<OUTSUFFIX>.tsv
	3) <BAM_FILE_NAME>.<MIN_LENGTH_CHUNKS>to<MAX_LENGTH_CHUNKS>.chunks/y_train_<OUTSUFFIX>.tsv
	4) <BAM_FILE_NAME>.<MIN_LENGTH_CHUNKS>to<MAX_LENGTH_CHUNKS>.chunks/y_test_<OUTSUFFIX>.tsv

It's very likely that the final training dataset will be composed of a combination of different organisms and runs so, because of that, an additional accessory script is provided to unify filtered training/test datasets into a single tabular file, ready to be used for training and testing purposes. This is the *make_global_dataset.py* accessory script, which takes in input this mandatory arguments:

  	 python3 make_global_dataset.py $1 $2 $3 $4 $5

    	 where the 4 arguments are:

  	 1) single_datasets_list_fp ===> a csv with one row for every dataset to be merged and 2 columns containing:
    					 A) the sample/run name, B) the full-path to the NanoListener directory with X_* / y_* filtered files;
	 2) outdir ====================> the output directory where the merged datasets will be saved;
	 3) dataset_group =============> the dataset group (train or test, alternatively);
	 4) dataset_suffix ============> if a suffix has been used provide it here;
	 5) overwrite =================> whether the script can overwrite the output if it does already exist. (False or True, ### DANGER ###: if True, be very careful!)

In the output folder <outdir> the final merged, shuffled and consolidated training (or test, depending on the used <dataset_group>) dataset will be saved as: 

	<outdir>/concat.x_y_meta.<dataset_group>.shuffled.tsv

The final datasets is ready to be used to train modification-aware (or unaware if NanoListener has been used in unmodified mode without per-read tables) basecalling models for dRNA experiments (R9 or RNA004).
