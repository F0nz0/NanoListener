from __init__ import __version__, __author__
import pysam,sys
import random, os
from collections import defaultdict
from datetime import datetime

###### Description: A script to downsample overepresented regions in a BAM file
###### Generated with ChatGPT and modified where required #####

# I/O
input_bam = sys.argv[1]   #-------------------------> Input BAM file (only forward reads)
output_bam = sys.argv[2]  #-------------------------> Output downsampled BAM file
max_reads_per_contig = int(sys.argv[3]) #-----------> Maximum number of reads per-region to perform downsampling (e.g. if a regions has more than "max_reads_per_contig" so a random sample of "max_reads_per_contig" will be retained.
min_reads_per_contig = int(sys.argv[4]) #-----------> Minumum number of reads per-region to retain the region. If a region has less than "min_reads_per_contig" thus the region will be discarded.

# print logs
print(f"[{datetime.now()}] NanoListener {__version__}.", flush=True)
print(f"[{datetime.now()}] sample_bam_by_contig.py\nInput:")
print("BAM input:", input_bam)
print("BAM output:", output_bam)
print("Max Reads per Contig/Region:", max_reads_per_contig)
print("Min Reads per Contig/Region:", min_reads_per_contig)

# read input bam file
bamfile = pysam.AlignmentFile(input_bam, "rb")

# group read per contig/region
contig_reads = defaultdict(list)

print(f"[{datetime.now()}] Collecting contigs/regions and corresponding reads...")
for read in bamfile.fetch(until_eof=True):
    if not read.is_unmapped:
        contig = bamfile.get_reference_name(read.reference_id)
        contig_reads[contig].append(read)

bamfile.close()

print(f"[{datetime.now()}] Perform downsampling on overrepresented regions...")
# if required, take a random sample of reads for contigs/regions with more than
# max_reads_per_contig
sampled_reads = []

for contig, reads in contig_reads.items():
    if len(reads) >= min_reads_per_contig:
        if len(reads) <= max_reads_per_contig:
            sampled_reads.extend(reads)
        else:
            sampled_reads.extend(random.sample(reads, max_reads_per_contig))

# write on disk the downsampled BAM file...
bamfile = pysam.AlignmentFile(input_bam, "rb")  # open input bam file to retrieve header...
out_bam = pysam.AlignmentFile(output_bam, "wb", header=bamfile.header)

for read in sampled_reads:
    out_bam.write(read)

bamfile.close()
out_bam.close()

# make index on output downsampled BAM file 
print(f"[{datetime.now()}] Indexing output BAM file...")
os.system(f"samtools sort -O BAM {output_bam} > {output_bam}.sorted.bam")
os.system(f"samtools index {output_bam}.sorted.bam")
os.remove(output_bam)

print(f"[{datetime.now()}] Filtering and Dowsampling completed.")
print(f"[{datetime.now()}] Final downsampled, sorted and indexed BAM file: {output_bam}.sorted.bam")
print(f"[{datetime.now()}] Computation Finished.")