# import required modules/libraries
from __init__ import __version__
import os, pysam, shutil, argparse, sys
import numpy as np
import pandas as pd
from datetime import datetime
from tqdm import tqdm

##### Author: Adriano Fonzino, Ph.D.
##### Description: an utility script to collect Dorado modifications calls produced via ModKit Extract Calls.
##### These modifications are then used to convert dorado call_codes into nanospeech call codes and
##### to produce, eventually, per-read tsv tables with modified positions (contig,pos0based,mod_char)
##### used by NanoListener to mark modified bases/nucleotides.
##### V2 --> in this version a user can decide to feed NanoListener with ModKit information about high quality SNV base quality

# define main function
def CreateModPosTablesFromDoradoCalls(bam_fp, dorado_calls_fp, mods_code_table_fp, min_call_p=0., min_base_q=0., retrieve_snvs=False, discard_mods=False, output_dir=None, overwrite=False):
    start_time = datetime.now()

    # print some starting logs
    print(f"[{datetime.now()}] BAM file: {bam_fp}", flush=True)
    print(f"[{datetime.now()}] Dorado ModKit Extract Calls file: {dorado_calls_fp}", flush=True)
    print(f"[{datetime.now()}] MODS_CODE.csv file: {mods_code_table_fp}", flush=True)
    if min_call_p > 0:
        print(f"[{datetime.now()}] Minimal call probability filtering activated: {min_call_p}", flush=True)
    else:
        print(f"[{datetime.now()}] Minimal call probability filtering deactivated.", flush=True)
    if min_base_q > 0:
        print(f"[{datetime.now()}] Minimal base quality filtering activated: {min_base_q}", flush=True)
    else:
        print(f"[{datetime.now()}] Minimal base quality filtering deactivated.", flush=True)
    print(f"[{datetime.now()}] Retrieve SNVs info: {retrieve_snvs}", flush=True)
    print(f"[{datetime.now()}] Discard all Modifications (useful for IVTs runs to collect only SNVs): {discard_mods}", flush=True)

    if not retrieve_snvs and discard_mods:
        sys.exit(f"[{datetime.now()}] ERROR! You are discarding both SNVs and modification...the ModKit table would be not used at all. Exiting...")

    if not output_dir:
        print(f"[{datetime.now()}] Output directory not provided. Producing one using basic BAM file radix path...",
              flush=True)
        output_dir = bam_fp + ".mod_pos_tables"
    print(f"[{datetime.now()}] Output Directory: {output_dir}", flush=True)
    if os.path.exists(output_dir):
        if overwrite:
            print(f"[{datetime.now()}] Output Directory already exists. Overwriting...", flush=True)
            shutil.rmtree(output_dir)
        else:
            sys.exit(
                f"[{datetime.now()}] ERROR! Output directory {output_dir} at already exists. Please set overwrite to True. EXITING.")
    os.mkdir(output_dir)

    # open input files and perform some checks
    print(f"[{datetime.now()}] Loading MODS_CODE.csv containing modifications code:", flush=True)
    mods_code_table = pd.read_csv(mods_code_table_fp)
    print(mods_code_table, flush=True)

    bam = pysam.AlignmentFile(bam_fp)

    print(f"[{datetime.now()}] Loading Dorado-ModKit Extract Calls table in-memory and printing some statistics...",
          flush=True)
    dorado_calls = pd.read_table(dorado_calls_fp,
                                 names="read_id,forward_read_position,ref_position,chrom,mod_strand,ref_strand,ref_mod_strand,fw_soft_clipped_start,fw_soft_clipped_end,alignment_start,alignment_end,read_length,call_prob,call_code,base_qual,ref_kmer,query_kmer,canonical_base,modified_primary_base,fail,inferred,within_alignment,flag,motifs".split(
                                     ","))
    print("\n##### Modifications Table Calls Stats:", flush=True)
    print(dorado_calls.describe(), flush=True)
    print("\n##### Modifications Table Calls Codes Distribution:", flush=True)
    print(dorado_calls[["call_code", "canonical_base", "modified_primary_base"]].value_counts(), flush=True)

    print("\n##### Modifications Table Distributions of Regions/Contigs:", flush=True)
    print(dorado_calls["chrom"].value_counts(), flush=True)

    print("\n##### Strand check:", flush=True)
    print(dorado_calls[["mod_strand", "ref_strand", "ref_mod_strand"]].value_counts(), flush=True)

    # discarding unmodified calls rows if required
    if not retrieve_snvs:
        dorado_calls = dorado_calls.query("call_code != '-'").copy()
        print(f"[{datetime.now()}] Discarding unmodified (only modified nucleotides - no SNVs) positions...", flush=True)
        print("\n##### Modifications Table Calls Codes Distribution after unmodified calls filtering:", flush=True)
        print(dorado_calls[["call_code", "canonical_base", "modified_primary_base"]].value_counts(), flush=True)

    # discarding modifications rows if required
    if discard_mods:
        dorado_calls = dorado_calls.query("call_code == '-'").copy()
        print(f"[{datetime.now()}] Discarding modified rows (only unmodified nucleotides - only SNVs)...", flush=True)
        print("\n##### Modifications Table Calls Codes Distribution after modified calls filtering:", flush=True)
        print(dorado_calls[["call_code", "canonical_base", "modified_primary_base"]].value_counts(), flush=True)

    # filtering for call probability and base quality, if requested...
    if min_call_p > 0 or min_base_q > 0:
        print(f"[{datetime.now()}] Filtering for call probability and/or base quality:", flush=True)
        dorado_calls = dorado_calls.query(f"call_prob >= {min_call_p}").query(f"base_qual >= {min_base_q}").copy()
        print("\n##### Modifications Table Calls Codes Distribution after call probability and/or base quality filtering:", flush=True)
        print(dorado_calls[["call_code", "canonical_base", "modified_primary_base"]].value_counts(), flush=True)

    print(f"\n[{datetime.now()}] Start Computations...", flush=True)
    total_number_of_reads = bam.mapped + bam.unmapped  # retrieve total number of reads into bam file
    pbar_update = int(total_number_of_reads / 100)
    # define some counters
    processed_reads_unmapped = 0
    processed_reads_w_mods = 0
    processed_reads_no_mods = 0
    processed_reads_total = 0
    # starting interation over bam file
    with tqdm(total=total_number_of_reads, desc="Processed Reads") as pbar:
        for c, read in enumerate(bam):
            if not read.is_unmapped:
                read_id = read.qname
                # retrieve unmodified calls for current read-id and sorting per-coordinate
                read_calls = dorado_calls.query(f"read_id == '{read_id}'").copy()
                read_calls.sort_values(by=["ref_position"], inplace=True)
                if not read_calls.empty:
                    if read_calls.chrom.value_counts().shape[0] == 1:  # check if multiple alignment for the same reads are absent
                        mod_pos_table = []
                        for base in read_calls.itertuples():
                            # detect call_code and convert to NanoSpeech single char encoding
                            if not base.call_code == "-": # it's a putative modification
                                dorado_mod_code = base.call_code
                                nanospeech_mod_code = mods_code_table[mods_code_table["Dorado_code (pseudo ChEBI)"] == dorado_mod_code]["ModificationCharNanoSpeech"].iloc[0]
                                mod_pos_table.append([base.chrom, base.ref_position, nanospeech_mod_code])
                            else: # it is a putative snv (if not requested this part of code will be not executed at all)
                                # take ref_kmer[2] != query_kmer[2] (for 5-mer extracted kmers by modkit)
                                if base.ref_kmer[2] != base.query_kmer[2]:
                                    mod_pos_table.append([base.chrom, base.ref_position, base.query_kmer[2].upper()]) # add the predicted SNV base
                        mod_pos_table = pd.DataFrame(mod_pos_table, columns=["contig", "pos0based", "mod_code_nanospeech"])
                        # save to disk
                        mod_pos_table.to_csv(os.path.join(output_dir, f"{read_id}.tsv"), header=None, index=None, sep="\t")
                        processed_reads_w_mods += 1 # or snv detected
                    else:
                        sys.exit(f"ERROR!!! Multiple regions for the same read-id {read_calls}")
                else:
                    processed_reads_no_mods += 1
            else:
                processed_reads_unmapped += 1
            processed_reads_total += 1
            if c % pbar_update == 0:
                pbar.update(pbar_update)
    bam.close()

    print(
        f"[{datetime.now()}] Iteration Finished. Total processed reads (unmapped|mapped_w_mods|mapped_no_mods): {processed_reads_total} ({processed_reads_unmapped}|{processed_reads_w_mods}|{processed_reads_no_mods})",
        flush=True)
    print(f"[{datetime.now()}] Computation Finished. Elapsed time: {datetime.now() - start_time}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="NanoListener: CreateModPosTablesFromDoradoCalls_V2.py program. A tool to produce mod_pos tables for NanoListener/NanoSpeech starting from Dorado modification calls and ModKit Extract Calls toolkit.")
    # bam_fp, dorado_calls_fp, mods_code_table_fp, output_dir = None, overwrite = False
    parser.add_argument("-b",
                        "--bam_fp",
                        required=True,
                        type=str,
                        help="--bam_fp: \t a <str> indicating the fullpath for the filtered (mod)BAM file (only strand +).")

    parser.add_argument("-d",
                        "--dorado_calls_fp",
                        required=True,
                        type=str,
                        help="--dorado_calls_fp: \t a <str> indicating the fullpath for the Dorado/ModKitExtractCalls file produced from the input (mod) BAM.")

    parser.add_argument("-m",
                        "--mods_code_table_fp",
                        required=True,
                        type=str,
                        help="--mods_code_table_fp: \t a <str> indicating the fullpath for the mods_code_table file used for binding Dorado and NanoSpeech call_codes.")

    parser.add_argument("-p",
                        "--min_call_p",
                        required=False,
                        default=0.,
                        type=float,
                        help="--min_call_p: \t a <float> indicating the minimum call probability to retain a call (from 0. to 1.). [0. --> deactivated]")

    parser.add_argument("-q",
                        "--min_base_q",
                        required=False,
                        default=0.,
                        type=float,
                        help="--min_base_q: \t a <float> indicating the minimum base quality to retain a call. [0. --> deactivated]")

    parser.add_argument("-snv",
                        "--retrieve_snvs",
                        required=False,
                        type=str,
                        default=False,
                        help="--retrieve_snvs: \t a <bool> indicating if you want to retrieve also SNVs information from the Dorado/ModKitExtractCalls file. [False]")

    parser.add_argument("-D",
                        "--discard_mods",
                        required=False,
                        type=str,
                        default=False,
                        help="--discard_mods: \t a <bool> indicating if you want to discard modifications rows (retrieve only SNVs information) from the Dorado/ModKitExtractCalls file. [False]")

    parser.add_argument("-o",
                        "--output_dir",
                        required=False,
                        type=str,
                        default=None,
                        help="--output_dir: \t a <str> with the full_path for the output directory. If not provided by default will use the BAM file radix + .mod_pos_tables pattern [None]")

    parser.add_argument("-f",
                        "--overwrite",
                        required=False,
                        type=str,
                        default=False,
                        help="--overwrite: \t a <bool> indicating if you want to delete the output directory if it already exists. [False]")

    args = parser.parse_args()
    bam_fp = args.bam_fp
    dorado_calls_fp = args.dorado_calls_fp
    mods_code_table_fp = args.mods_code_table_fp
    min_call_p = args.min_call_p
    min_base_q = args.min_base_q
    retrieve_snvs = args.retrieve_snvs
    if retrieve_snvs == "None" or retrieve_snvs == "False":
        retrieve_snvs = False
    elif retrieve_snvs == "True":
        retrieve_snvs = True
    discard_mods = args.discard_mods
    if discard_mods == "None" or discard_mods == "False":
        discard_mods = False
    elif discard_mods == "True":
        discard_mods = True
    output_dir = args.output_dir
    if output_dir == "None": output_dir = None
    overwrite = args.overwrite
    if overwrite == "None" or overwrite == "False":
        overwrite = False
    elif overwrite == "True":
        overwrite = True

        # print some starting info related to version, used program and to the input arguments
        print(f"[{datetime.now()}] NanoListener version {__version__}.", flush=True)
        print(f"[{datetime.now()}] CreateModPosTablesFromDoradoCalls_V2.py program.", flush=True)
        print(f"[{datetime.now()}] Input arguments:", flush=True)
        for argument in args.__dict__.keys():
            print(f"\t- {argument} --> {args.__dict__[argument]}", flush=True)

    # launch main program
    CreateModPosTablesFromDoradoCalls(bam_fp=bam_fp,
                                      dorado_calls_fp=dorado_calls_fp,
                                      mods_code_table_fp=mods_code_table_fp,
                                      min_call_p=min_call_p,
                                      min_base_q=min_base_q,
                                      retrieve_snvs=retrieve_snvs,
                                      discard_mods=discard_mods,
                                      output_dir=output_dir,
                                      overwrite=overwrite)


