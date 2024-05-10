import os,pysam
from datetime import datetime
import pandas as pd
import argparse

def create_transcripts_white_list(bam_filepath, 
                                  min_length = 1000,
                                  min_coverage = 90,
                                  min_numreads = 20,
                                  min_meandepth = 15,
                                  min_meanbaseq = 15,
                                  min_meanmapq = 40):
    # open bam file and start computation
    with pysam.AlignmentFile(bam_filepath) as bam:
        # create coverage file for each region Transcript and load in memory
        bam_coverage_filepath = bam_filepath+".coverage"
        print(f"[{datetime.now()}] Generating BAM coverage table using samtools.", flush=True)
        os.system(f"samtools coverage {bam_filepath} > {bam_coverage_filepath}")
        bam_coverage = pd.read_table(bam_coverage_filepath)
        print(f"[{datetime.now()}] Loading generated coverage table via samtools coverage command.", flush=True)
        print(f"[{datetime.now()}] Coverage Table shape:", flush=True)
        print(bam_coverage.shape)
        print(f"\n[{datetime.now()}] Printing coverage table statistics:", flush=True)
        print(bam_coverage.describe())
        print("", flush=True)
        # select only transcripts with:
        print(f"[{datetime.now()}] Starting generation of white list for transcripts and regions from the aligment BAM file...")
        bam_cov_filt = bam_coverage[(bam_coverage["endpos"]>=min_length)&
                                    (bam_coverage["coverage"]>=min_coverage)&
                                    (bam_coverage["numreads"]>=min_numreads)&
                                    (bam_coverage["meandepth"]>=min_meandepth)&
                                    (bam_coverage["meanbaseq"]>=min_meanbaseq)&
                                    (bam_coverage["meanmapq"]>=min_meanmapq)].copy()
        print(f"[{datetime.now()}] White list of transcripts length:", flush=True)
        print(bam_cov_filt.shape, flush=True)
        print("", flush=True)
        print(f"[{datetime.now()}] Save on disk the transcript white list to {bam_filepath+'.white_list'}", flush=True)
        bam_cov_filt["#rname"].to_csv(bam_filepath+".white_list", header=None, index=None)


if __name__ == "__main__":
    # bam_filepath, min_length = 1000, min_coverage = 90, 
    # min_numreads = 20, min_meandepth = 15, min_meanbaseq = 15, min_meanmapq = 40

    parser = argparse.ArgumentParser(description=f"create_transcript_white_list_KO")

    parser.add_argument("-b",
                        "--bam_filepath",
                        required=True,
                        type=str,
                        help="--bam_filepath: \t a <str> for the BAM file path to be filtered. The filtered outputfile will be at the path <BAM_basename>.filtered.bam")
    
    parser.add_argument("-l",
                        "--min_len",
                        required=False,
                        type=int,
                        default=1000,
                        help="minimum length [1000]")
    
    parser.add_argument("-c",
                        "--min_cov",
                        required=False,
                        type=int,
                        default=90,
                        help="minimum coverage percentage [90]")
    
    parser.add_argument("-r",
                        "--min_numreads",
                        required=False,
                        type=int,
                        default=20,
                        help="minimum number of reads per transcript [20]")
    
    parser.add_argument("-d",
                        "--min_meandepth",
                        required=False,
                        type=int,
                        default=15,
                        help="minimum average depth for transcript [15]")

    parser.add_argument("-q",
                        "--min_meanbaseq",
                        required=False,
                        type=int,
                        default=15,
                        help="minimum average base quality per transcript [15]")
    
    parser.add_argument("-mp",
                        "--min_meanmapq",
                        required=False,
                        type=int,
                        default=40,
                        help="minimum average mapping quality per transcript [40]")

    # bam_filepath, min_length, min_coverage, min_numreads, 
    # min_meandepth, min_meanbaseq, min_meanmapq
    args = parser.parse_args()
    bam_filepath = args.bam_filepath
    min_length = args.min_len
    min_coverage = args.min_cov
    min_numreads = args.min_numreads
    min_meandepth = args.min_meandepth
    min_meanbaseq = args.min_meanbaseq
    min_meanmapq = args.min_meanmapq

    # print some starting info related to version, used program and to the input arguments
    print(f"[{datetime.now()}] create_transcript_white_list_KO.py program.", flush=True)
    print(f"[{datetime.now()}] Input arguments:", flush=True)
    for argument in args.__dict__.keys():
        print(f"\t- {argument} --> {args.__dict__[argument]}", flush=True)

    # launch main program
    # filter starting bam file
    create_transcripts_white_list(bam_filepath, 
                                  min_length = min_length,
                                  min_coverage = min_coverage,
                                  min_numreads = min_numreads,
                                  min_meandepth = min_meandepth,
                                  min_meanbaseq = min_meanbaseq,
                                  min_meanmapq = min_meanmapq)