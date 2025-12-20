# import basic modules
from __init__ import __version__, __author__
import os, sys, pysam
from shutil import rmtree
from datetime import datetime
from tqdm import tqdm
from math import ceil


##### Author: Adriano Fonzino, Ph.D.
##### Description: A script to sequentially split BAM files (only primary forward reads) into N parts, to produce smaller eventalign file and
##### use a map-reduce like strategy when producing training datasets by NanoListener.

# main function
def split_bam_file(bam_fp, n_splits=10, overwrite=False):
    print("#############################", flush=True)
    print(f"[{datetime.now()}] NanoListener {__version__}.", flush=True)
    print(f"[{datetime.now()}] split_bam_file.py program.", flush=True)
    print(f"- Input BAM file: {bam_fp}", flush=True)
    print(f"- N° of splits: {n_splits}", flush=True)
    out_dir = bam_fp + f".{n_splits}_split_bams"
    print(f"[{datetime.now()}] Making output directory where split BAM file will be stored: {out_dir}", flush=True)
    if os.path.exists(out_dir):
        if overwrite:
            rmtree(out_dir)
        else:
            sys.exit(f"[{datetime.now()}] ERROR! Output directory {out_dir} at already exists. Please set overwrite to True. EXITING.")
    os.mkdir(out_dir)
    # start computations
    bam = pysam.AlignmentFile(bam_fp, "rb")
    n_reads = bam.mapped + bam.unmapped
    n_reads_per_split = ceil(n_reads / n_splits)
    pbar_update_interval = int(n_reads / 100)
    print(f"[{datetime.now()}] n° of reads into BAM file: {n_reads}", flush=True)
    print(f"[{datetime.now()}] n° reads per split: {n_reads_per_split}", flush=True)
    total_reads_written = 0
    n_reads_written_currents_split = 0
    current_split = 1
    current_bam_fp = os.path.join(out_dir, os.path.basename(bam_fp) + f".split{current_split}.bam")
    current_bam = pysam.AlignmentFile(current_bam_fp, "wb", header=bam.header)
    print(f"[{datetime.now()}] Working on split {current_split} and writing on: {current_bam_fp}", flush=True)
    with tqdm(total=n_reads, desc="Processed reads") as pbar:
        for read in bam:
            if n_reads_written_currents_split == n_reads_per_split:
                print(f"\n[{datetime.now()}] BAM file split: {current_split} completed. Total number of reads written {n_reads_written_currents_split}. Closing, indexing and opening the next one", flush=True)
                current_bam.close()
                pysam.index(current_bam_fp)
                current_split += 1
                n_reads_written_currents_split = 0
                current_bam_fp = os.path.join(out_dir, os.path.basename(bam_fp) + f".split{current_split}.bam")
                current_bam = pysam.AlignmentFile(current_bam_fp, "wb", header=bam.header)
                print(f"[{datetime.now()}] Working on split {current_split} and writing on: {current_bam_fp}", flush=True)
            current_bam.write(read)
            total_reads_written += 1
            n_reads_written_currents_split += 1
            if total_reads_written % pbar_update_interval == 0: pbar.update(pbar_update_interval)
    print(f"\n[{datetime.now()}] BAM file split: {current_split} (last one) completed. Total number of reads written {n_reads_written_currents_split}. Closing, indexing and opening the next one", flush=True)
    current_bam.close()
    pysam.index(current_bam_fp)
    bam.close()
    print(f"[{datetime.now()}] Computation Finished.", flush=True)

# define Inputs
bam_fp = sys.argv[1]
n_splits = int(sys.argv[2])
overwrite = bool(sys.argv[3])

# launch program
split_bam_file(bam_fp, n_splits, True)