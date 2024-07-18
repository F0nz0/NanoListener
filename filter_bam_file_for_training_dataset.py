# import basic modules
import pysam, os, argparse
import pandas as pd
from datetime import datetime


def filter_bam_file_by_perc_map(fasta_filepath, bam_filepath, perc_threshold=0.3, region_list=None, filt_bam_file_suffix=None):
    reads_to_exclude_filepath = bam_filepath+".reads_to_exclude"
    fasta = pysam.FastaFile(fasta_filepath)
    bam = pysam.AlignmentFile(bam_filepath)
    reads_to_exclude = open(reads_to_exclude_filepath, "w")
    # load the list of regions to be mantained only (the remaining part will be excluded)
    # this is to focus on a given set of transcripts
    if region_list:
        print(f"[{datetime.now()}] Loading regions whitelist from {region_list}.")
        rl = set()
        with open(region_list, "r") as l:
            for r in l:
                rl.add(r.rstrip())
        diffs_perc = []
        reads_excluded = 0
        for total,r in enumerate(bam):
            if not r.is_unmapped:
                if r.reference_name in rl:
                    # it is in the white list...evaluate diffs into mapped read and len of 
                    # transcript for the curre read
                    fasta_ref = fasta.fetch(r.reference_name)
                    diff = len(fasta_ref)-r.reference_length
                    diff_perc = diff/len(fasta_ref)
                    diffs_perc.append(diff_perc)
                    if diff_perc >= perc_threshold:
                        reads_to_exclude.write(r.qname+"\n")
                        reads_excluded += 1
                else:
                    # not in the whitelist
                    reads_to_exclude.write(r.qname+"\n")
                    reads_excluded += 1
            else:
                reads_to_exclude.write(r.qname+"\n")
                reads_excluded += 1

    # no reference white list provided filtering only of perc_mapped 
    # reads with respect the whole ref seq
    else:
        print(f"[{datetime.now()}] Regions whitelist not provided: working only on perc. mapped read.")
        diffs_perc = []
        reads_excluded = 0
        for total,r in enumerate(bam):
            if not r.is_unmapped:
                fasta_ref = fasta.fetch(r.reference_name)
                diff = len(fasta_ref)-r.reference_length
                diff_perc = diff/len(fasta_ref)
                diffs_perc.append(diff_perc)
                if diff_perc >= perc_threshold:
                    reads_to_exclude.write(r.qname+"\n")
                    reads_excluded += 1
            else:
                reads_to_exclude.write(r.qname+"\n")
                reads_excluded += 1

    print(f"[{datetime.now()}] Reads to be excluded are {reads_excluded} on a total of {total}.", flush=True)
    bam.close()
    fasta.close()
    reads_to_exclude.close()

    # print some statistics of percentage
    df_diffs_perc = pd.DataFrame(diffs_perc, columns=["difference_perc"])
    print(f"[{datetime.now()}] Statistics about percentages of covered reference sequence from alignments of the input BAM file:", flush=True)
    print(df_diffs_perc.describe(), flush=True)
    print("", flush=True)

    # filtering
    if filt_bam_file_suffix == None:
        output_bam = os.path.splitext(bam_filepath)[0]+".filtered.bam"
    else:
        output_bam = os.path.splitext(bam_filepath)[0]+"."+filt_bam_file_suffix+".bam"
    print(f"[{datetime.now()}] Filtering via samtools. Output at: {output_bam}", flush=True)
    os.system(f"samtools view -N {reads_to_exclude_filepath} -U {output_bam} -o {output_bam+'.excluded'} {bam_filepath}")
    os.system(f"rm {reads_to_exclude_filepath}")
    os.system(f"rm {output_bam+'.excluded'}")
    os.system(f"samtools index {output_bam}")
    print(f"[{datetime.now()}] Computation finished.", flush=True)


if __name__ == "__main__":
    #fasta_filepath, bam_filepath, perc_threshold=0.3, region_list=None, filt_bam_file_suffix=None
    parser = argparse.ArgumentParser(description=f"filter_bam_file_by_perc_map")

    parser.add_argument("-f",
                        "--fasta_filepath",
                        required=True,
                        type=str,
                        help="--fasta_filepath: \t a <str> for the fasta file path.")

    parser.add_argument("-b",
                        "--bam_filepath",
                        required=True,
                        type=str,
                        help="--bam_filepath: \t a <str> for the BAM file path to be filtered. The filtered outputfile will be at the path <BAM_basename>.filtered.bam")
    
    parser.add_argument("-t",
                        "--perc_threshold",
                        required=False,
                        type=float,
                        default=0.3,
                        help="--perc_threshold: \t a <float> indicating the maximum amount of unligned reads with respect reference transcript to be filtered. [0.3]")

    parser.add_argument("-r",
                        "--region_whitelist",
                        required=False,
                        type=str,
                        default=None,
                        help="--region_whitelist: \t a <str> optional indicating the path for a list of regions to be extracted from BAM file. [None]")

    parser.add_argument("-o",
                        "--filt_bam_file_suffix",
                        required=False,
                        type=str,
                        default=None,
                        help="--filt_bam_file_suffix: \t a <str> optional indicating the suffix to be used to write the filtered BAM output file on disk.. [None]")

    args = parser.parse_args()
    fasta_filepath = args.fasta_filepath
    bam_filepath = args.bam_filepath
    perc_threshold = args.perc_threshold
    region_whitelist = args.region_whitelist
    filt_bam_file_suffix = args.filt_bam_file_suffix
    if region_whitelist == "None":
        reads_list = None
    if filt_bam_file_suffix == "None":
        filt_bam_file_suffix = None

    # print some starting info related to version, used program and to the input arguments
    print(f"[{datetime.now()}] filter_bam_file_by_perc_map.py program.", flush=True)
    print(f"[{datetime.now()}] Input arguments:", flush=True)
    for argument in args.__dict__.keys():
        print(f"\t- {argument} --> {args.__dict__[argument]}", flush=True)

    # launch main program
    # filter starting bam file
    filter_bam_file_by_perc_map(fasta_filepath=fasta_filepath, 
                                bam_filepath=bam_filepath, 
                                perc_threshold=perc_threshold, 
                                region_list=region_whitelist, 
                                filt_bam_file_suffix=filt_bam_file_suffix)
