from __init__ import __version__
from collections import Counter
import numpy as np
import pandas as pd
import os, sys, argparse
from tqdm import tqdm
from datetime import datetime

def GetNucleotidesDistributionsFromNanoListenerFinalDS(ds_fp, mods_string='imfyezdxrkwb'):
    '''
    NanoListener: GetNucleotidesDistributionsFromNanoListenerFinalDS.py program
    Description: A script to count modified and unmodified chunks and nucleotides from a final NanoListener dataset.", flush=True)
    '''
    print(f"[{datetime.now()}] Input NanoListener Final DS: {ds_fp}", flush=True)
    # start Computations
    print(f"[{datetime.now()}] Start Computations...", flush=True)
    start_time = datetime.now()
    mods = list(mods_string)
    with tqdm() as pbar:
        with open(ds_fp, "r") as ds:
            total_chunks = 0
            mod_chunks = 0
            unmod_chunks = 0
            nt_counts = {}
            for c, r in enumerate(ds):
                kmer_status = "UNMOD"
                kmer = r.strip().split("\t")[-1].lower()
                kmer_dict = dict(Counter(kmer))
                for nt in kmer_dict.keys():
                    if not nt in nt_counts.keys():
                        nt_counts[nt] = kmer_dict[nt]
                    else:
                        nt_counts[nt] += kmer_dict[nt]
                for b in kmer_dict.keys():
                    if b in mods:
                        kmer_status = "MOD"
                        mod_chunks += 1
                        break
                if not kmer_status == "MOD":
                    unmod_chunks += 1
                if c % 10000 == 0: pbar.update(10000)
                total_chunks += 1

    # print logs
    print("#############")
    print("Computation Finished. Elapsed Time:", (datetime.now() - start_time), flush=True)
    print(f"\nTotal Processed Chunks [mod] [unmod]: {total_chunks} [{mod_chunks}] [{unmod_chunks}]", flush=True)
    print("\nNucleotides distribution:", flush=True)
    print(pd.Series(nt_counts).sort_index(), flush=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"NanoListener {__version__}: GetNucleotidesDistributionsFromNanoListenerFinalDS.py program.\nDescription: A script to count modified and unmodified chunks and nucleotides from a final NanoListener dataset.")

    parser.add_argument("-ds",
                        "--dataset_filepath",
                        required=True,
                        type=str,
                        help="--dataset_filepath: \t a <str> with the full_path for the NanoListener final dataset to be analyzed.")

    parser.add_argument("-m",
                        "--mods_string",
                        required=False,
                        type=str,
                        default='imfyezdxrkwb',
                        help="--mods_string: \t a <str> containing all the modification chars to be looked for. [imfyezdxrkwb]")

    args = parser.parse_args()
    dataset_filepath = args.dataset_filepath
    mods_string = args.mods_string

    # print some starting info related to version, used program and to the input arguments
    print(f"[{datetime.now()}] NanoListener {__version__}", flush=True)
    print(f"[{datetime.now()}] GetNucleotidesDistributionsFromNanoListenerFinalDS.py program.", flush=True)
    for argument in args.__dict__.keys():
        print(f"\t- {argument} --> {args.__dict__[argument]}", flush=True)

    # launch main program
    GetNucleotidesDistributionsFromNanoListenerFinalDS(ds_fp=dataset_filepath, mods_string=mods_string)