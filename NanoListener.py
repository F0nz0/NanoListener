# import basic modules
import os, sys, pysam, argparse, shutil
import pandas as pd
from datetime import datetime
from multiprocessing import Queue, Process, Value
import numpy as np
from ont_fast5_api.fast5_interface import get_fast5_file

# define functions
def decode_seq(poss, refs):
    '''
    Function to decode reference bases in function of the positions in search for movements into events.
    '''
    seq = []
    last_pos = None
    for pos,kmer in zip(poss, refs):
        if last_pos == None:
            last_pos = pos
            seq.append(kmer[0])
        else:
            if pos != last_pos:
                last_pos = pos
                seq.append(kmer[0])
    seq = "".join(seq)
    return seq


# defining useful functions to handle with fast5 files
def raw_to_pA(f5):
    '''
    Function to transform back from raw signal to pA scale.
    '''
    try:
        raw_unit = f5.get_channel_info()["range"] / f5.get_channel_info()["digitisation"]
        pA_signal = (f5.get_raw_data() + f5.get_channel_info()["offset"]) * raw_unit
        return pA_signal
    except Exception as e:
        print("AN EXCEPTION HAS OCCURRED!\n", e, flush=True)


def retrieve_read_pA_fastq_from_fast5(fast5_fullpath, read_name_id):
    '''
    Retrieve fastq and pA converted data related to a given readname_id from a fast5 file 
    '''
    with get_fast5_file(fast5_fullpath) as f5:
        r = f5.get_read(read_name_id)
        if r.read_id == read_name_id:
            read_name = r.read_id
            #print(f"One putative read example found with id {read_name} in fast5 file: {fast5_fullpath}", flush=True)
            pA_data = raw_to_pA(r)
            fastq = r.get_analysis_dataset("Basecall_1D_000/BaseCalled_template", "Fastq")
    return pA_data, fastq


# defining useful functions to handle with pod5 files
def retrieve_read_pA_from_pod5(pod5_fullpath, read_name_id):
    with pod5.Reader(pod5_fullpath) as reader:
        for read in reader.reads([read_name_id]):
            if str(read.read_id) == read_name_id:
                signal = read.signal
                raw_unit = read.calibration.scale
                offset = read.calibration.offset
                pA_signal = (offset + signal) * raw_unit
                return pA_signal


def producer(eventalign_path, q, threads_n, n_reads_to_process_per_region, regions):
    '''
    Function that process the input eventalign file splitting event per-read and then put these into
    a queue which will be processed by consumers.
    '''
    print(f"[{datetime.now()}] [Producer Message] Starting Iteration on the input eventalign file...\n", flush=True)
    start = datetime.now()
    
    # start computation for each region
    total_processed_reads = 0
    for region in regions:
        print(f"[{datetime.now()}] [Producer Message] Producer is starting to process region {region}.", flush=True)
        with open(eventalign_path, "rt") as f:
            # create a list of events for each read and append to the queue
            read = []
            read_name = None
            processed_reads = 0    # processed reads counter for current region
            for l in f:
                if l.startswith("contig"):
                    continue
                else:
                    line = l.rstrip().split("\t")
                    if line[0] == region:
                        if read_name == None:
                            read_name = line[3]
                            processed_reads += 1
                            total_processed_reads += 1
                            read.append(line)
                        elif line[3] == read_name:
                            read.append(line)
                        else:
                            # append to queue and restart for another read
                            q.put([read_name, read])
                            read_name = line[3]
                            processed_reads += 1
                            total_processed_reads += 1
                            read = []
                            # stop if number of reads to process has been reached
                            if processed_reads > n_reads_to_process_per_region:
                                break

    # append end signal for every consumer
    for t in range(threads_n):
        q.put(["None", None])

    # print producer's statistics.
    stop = datetime.now()
    print(f"[{datetime.now()}] [Producer Message] Producer finished to load reads into queue.", flush=True)
    print(f"[{datetime.now()}] [Producer Message] Total processed reads: {total_processed_reads-1}", flush=True)
    print(f"[{datetime.now()}] [Producer Message] Producer Elapsed time: {stop-start}", flush=True)
    
    # close queue and join all queue threads
    q.close()
    q.join_thread()
    
    print(f"[{datetime.now()}] [Producer Message] Finished.", flush=True)

def consumer_worker(q, id_consumer, summary_table, out_dir, min_l, max_l, counter_r, n_reads_to_process, mod_pos_dir=None, clip_outliers=[40,165], ref_filepath=None):
    # open global output file
    out_filepath = os.path.join(out_dir, f"cons_{id_consumer}.tsv")
    out_file = open(out_filepath, "w")
    min_l = min_l
    len_chunks = max_l
    average_len = np.mean([min_l, len_chunks])
    print(f"[{datetime.now()}] [Consumer {id_consumer} Message] Start processing reads extracted by producers.", flush=True)
    if not ref_filepath != None:
        sys.exit("ERROR!!! Reference fasta file have to be provided! Exiting.")
    ref = pysam.FastaFile(ref_filepath)
    while True:
        # get eventalign read
        input_evread = q.get()
        read_name = input_evread[0]
        eventalign_read = input_evread[1]
        if eventalign_read != None:
            # try load summary for current read
            if not summary_table.empty:
                try:
                    summary_table_curr_read = summary_table.query(f"read_name == '{read_name}'")
                except:
                    print(f"[{datetime.now()}] [Consumer {id_consumer} Message] Error on read {read_name} while loading fast5_filepath, scale and shift parameters.", flush=True)
                    print(summary_table_curr_read, flush=True)
            # try to load corresponding mod_pos_table from mod_pos_directory
            if mod_pos_dir != None:
                try:
                    mod_pos_table = pd.read_table(os.path.join(mod_pos_dir, f"{read_name}.tsv"), names=["region", "pos_0_based", "alt_base"])
                except Exception as error:
                    print(f"[{datetime.now()}] Impossible to load mod_pos_table for read {read_name}. Working in unmodified mode.\nError: {error}")
                    mod_pos_table = pd.DataFrame()
            else:
                # if mod_pos_dir is not provided working in normal mode (no modification correction for alternative base)
                mod_pos_table = pd.DataFrame()
            # concatenate read events and measurements into 3 lists: 1) currents, 2) refs (kmers) and 3) poss (positions)
            currents = []
            refs = []
            poss = []
            start_idxs = []
            stop_idxs = []
            region = ""
            for _,e in enumerate(eventalign_read):
                if _ == 0:
                    region = e[0] # region
                e_currents = [float(c) for c in e[-1].split(",")]
                currents += e_currents
                if not summary_table.empty:
                    try:
                        summary_table_curr_read = summary_table.query(f"read_name == '{read_name}'")
                        scale = summary_table_curr_read["scale"].values[0]
                        shift = summary_table_curr_read["shift"].values[0]
                        fast5_filpath = summary_table_curr_read["fast5_path"].values[0]
                    except:
                        print(f"[{datetime.now()}] [Consumer {id_consumer} Message] Error on read {read_name} while loading scale and shift parameters.", flush=True)
                        print(summary_table_curr_read, flush=True)
                if mod_pos_table.empty:
                    refs += [e[2]] * len(e_currents) # kmers
                else:
                    # detect if, for the current read, the current site is in the mod_pos_table list
                    detect = mod_pos_table[(mod_pos_table["region"]==e[0])&(mod_pos_table["pos_0_based"]==int(e[1]))]
                    if detect.empty:
                        refs += [e[2]] * len(e_currents) # kmers
                    else:
                        alt_base = detect["alt_base"].iloc[0]
                        refs += [alt_base + e[2][1:]] * len(e_currents)
                poss += [int(e[1])] * len(e_currents) # positions
                start_idxs += [int(e[-3])] * len(e_currents) # start_idx positions
                stop_idxs += [int(e[-2])] * len(e_currents) # stop_idx positions
            # scale to pA scale if summary table is available
            currents = [round(sa,3) for sa in map(lambda x: (scale*x)+shift, currents)]
            # assess if the length of currents is at least two times the required length of chunks
            if len(currents) > len_chunks*2:
                if len(poss) == len(refs) == len(currents) == len(start_idxs) == len(stop_idxs):
                    # start processing read and extracting chunks and printing progress message if needed
                    with counter_r.get_lock():
                        counter_r.value += 1
                    if counter_r.value % (n_reads_to_process / 10) == 0:
                        print(f"[{datetime.now()}] [Consumer {id_consumer} Message] Progress processed_reads/total_reads: {counter_r.value} / {n_reads_to_process}.\n", flush=True)
                    # take n number of chunks with the following formula
                    # OLD VERSION pos_to_extract = set([0] + np.random.randint(0, len(currents), int(len(currents)/(len(currents)/150))).tolist())
                    pos_to_extract = set([0] + np.random.randint(0, len(currents)-int(len_chunks*0.8), int(len(currents)/average_len)*50).tolist())
                    ### produce second version directly from fast5 using fast5 measuements indexes
                    pA_data, fastq = retrieve_read_pA_fastq_from_fast5(fast5_filpath, read_name)
                    for start in pos_to_extract:
                        # take chunks with random length between min_l and max_l
                        random_l = np.random.randint(min_l, int(len_chunks*0.8), 1)[0]
                        currents_chunk = currents[start:start+random_l]
                        if len(currents_chunk) > min_l: # position close to the end of the read signals could be very short...discarding these
                            p = poss[start:start+random_l]
                            s = start_idxs[start:start+random_l]
                            S = stop_idxs[start:start+random_l]
                            r = refs[start:start+random_l]
                            # retrieve reference
                            seq = decode_seq(p, r) # from eventalign
                            ref_chunk = ref.fetch(region, min(p), max(p)+1)
                            # create alignment positions / ref_chunk in case of mod_pos_table
                            if not mod_pos_table.empty:
                                df_r_p = pd.DataFrame([[region for e in range(min(p), max(p)+1)], [i for i in range(min(p), max(p)+1, 1)], list(ref_chunk)]).T
                                df_r_p.columns = ["region", "pos_0_based", "ref_base"]
                                ref_chunk_alt = ""
                                for i in df_r_p.itertuples():
                                    # detect if, for the current read, the current site is in the mod_pos_table list
                                    detect = mod_pos_table[(mod_pos_table["region"]==i.region)&(mod_pos_table["pos_0_based"]==i.pos_0_based)]
                                    if not detect.empty:
                                        alt_base = detect["alt_base"].iloc[0]
                                        ref_chunk_alt += alt_base
                                    else:
                                        ref_chunk_alt += i.ref_base
                                ref_chunk = ref_chunk_alt
                            pA_data_chunk = pA_data[min(s):max(S)][::-1]
                            # clipping of outliers values (
                            # mean +/- 4*stdev by default) if requested
                            if clip_outliers:
                                currents_chunk_df = pd.Series(pA_data_chunk) # convert to pandas series
                                currents_chunk_df[(currents_chunk_df<clip_outliers[0])|(currents_chunk_df>clip_outliers[1])] = round(currents_chunk_df.mean(), 3) # clip outliers to the average values of the current chunk
                                pA_data_chunk = currents_chunk_df.tolist() # come back to list for subsequent operation
                            if len(pA_data_chunk) <= max_l:
                                # padding to fixed max chunks length
                                padding_to_max_l = len_chunks - len(pA_data_chunk)
                                out_line = "\t".join([str(i) for i in pA_data_chunk])+"\t"*padding_to_max_l+"\t"+"\t".join(str(x) for x in [region, min(p), max(p), min(s), max(S), read_name, fast5_filpath, ref_chunk])+"\n"
                                out_file.write(out_line)
        else:
            # Stopping the loop of the consumer if found a end-signal (None) in the Queue.
            print(f"[{datetime.now()}] [Consumer {id_consumer} Message] Found end of Queue.\n", flush=True)
            out_file.close()
            ref.close()
            break


def eventalign_to_chunks(eventalign_path, summary_table_path=None, 
                         threads_n=1, min_l=300, max_l=300, mod_pos_dir=None, 
                         n_reads_to_process=100, clip_outliers=[40,165], outdir_prefix=None, ref_filepath=None):
    start_global = datetime.now()
    eventalign_path = eventalign_path
    summary_table_path = summary_table_path # to scale to pA scale eventalign currents measurements
    threads_n = threads_n
    min_l = min_l
    max_l = max_l
    mod_pos_dir = mod_pos_dir
    n_reads_to_process = n_reads_to_process
    clip_outliers = clip_outliers
    # generate output directory
    if outdir_prefix == None:
        out_dir = f"{eventalign_path}.{min_l}to{max_l}.chunks"
    else:
        out_dir = f"{outdir_prefix}.{min_l}to{max_l}.chunks"
    # print some starting messages
    print(f"[{datetime.now()}] [Main Process Message] NanoListener program.", flush=True)
    print(f"[{datetime.now()}] [Main Process Message] Input eventalign file: {eventalign_path}", flush=True)
    print(f"[{datetime.now()}] [Main Process Message] Input summary eventalign file: {summary_table_path}", flush=True)
    print(f"[{datetime.now()}] [Main Process Message] Output folder for eventalign collapsed reads: {out_dir}", flush=True)
    print(f"[{datetime.now()}] [Main Process Message] Thread used: {threads_n}", flush=True)
    print(f"[{datetime.now()}] [Main Process Message] Min chunks length used: {min_l}", flush=True)
    print(f"[{datetime.now()}] [Main Process Message] Max chunks length used: {max_l}", flush=True)
    print(f"[{datetime.now()}] [Main Process Message] Provided modified positions tables directory: {mod_pos_dir}", flush=True)
    print(f"[{datetime.now()}] [Main Process Message] Requested number of reads to be processed: {n_reads_to_process}", flush=True)
    print(f"[{datetime.now()}] [Main Process Message] Reference fasta file: {ref_filepath}", flush=True)
    print(f"[{datetime.now()}] [Main Process Message] Outliers clipping limits (lower, higher) to change with mean: {clip_outliers}", flush=True)
    
    # creating a new folder if it doesn't exist and if it exist eliminate that and create a new directory
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)
    
    # open summary file as pandas dataframe if provided
    if summary_table_path:
        print(f"[{datetime.now()}] [Main Process Message] Loading summaryfile table.", flush=True)   
        summary_table = pd.read_table(summary_table_path)[["read_name", "shift", "scale", "fast5_path"]]
    else:
        # if not provided set to None
        summary_table = pd.DataFrame()
        
    # detect regions in order to extract a balanced number of reads from each region
    print(f"[{datetime.now()}] [Main Process Message] Detecting the number of regions...", flush=True)
    regions = set()
    with open(eventalign_path, "r") as ev:
        for l in ev:
            if not l.startswith("contig"):
                r = l.split("\t")[0]
                if not r in regions:
                    regions.add(r)
    print(f"[{datetime.now()}] [Main Process Message] Number of regions detected into Eventalign file: {len(regions)}", flush=True)
    print(f"[{datetime.now()}] [Main Process Message] Regions detected into Eventalign file:\n {regions}", flush=True)

    # calulate number of reads to process per region (at least one read per region)
    if n_reads_to_process > len(regions):
        n_reads_to_process_per_region = int(n_reads_to_process / len(regions))
    else:
        n_reads_to_process_per_region = 1
        print(f"[{datetime.now()}] [Main Process Message] Warning: the number of reads to process is less than the number of regions...", flush=True)        
    print(f"[{datetime.now()}] [Main Process Message] Number of reads to exctract for each region: {n_reads_to_process_per_region}", flush=True)
    n_reads_to_process = n_reads_to_process_per_region * len(regions)
    print(f"[{datetime.now()}] [Main Process Message] Number of total reads to be processed based on n. of regions detected: {n_reads_to_process}", flush=True)
    
    # create a counter for the progress of number of reads elaborated
    counter_r = Value("i", 0)
    
    q = Queue(maxsize=threads_n*3)
    
    # create consumers processes
    consumers = []
    for t in range(threads_n):
        # q, id_consumer, summary_table, out_dir, min_l, max_l, counter_r, n_reads_to_process, mod_pos_dir=None, clip_outliers, ref_filepath
        consumers.append( Process(target=consumer_worker, args=(q,t+1,summary_table,out_dir,min_l,max_l,counter_r,n_reads_to_process,mod_pos_dir,clip_outliers,ref_filepath)))
    print(f"[{datetime.now()}] [Main Process Message] Generating requested consumers. N° of Consumers: {len(consumers)}", flush=True)
    # start consumers
    for c in consumers:
        c.start()
    
    # create a producer process
    # eventalign_path, q, threads_n, n_reads_to_process
    p = Process(target=producer, args=(eventalign_path,q,threads_n,n_reads_to_process_per_region,regions))
    p.start()
    
    # join consumers
    for c in consumers:
        c.join()
    # join producer
    p.join()
    
    stop_global = datetime.now()
    print(f"[{datetime.now()}] [Main Process Message] Computation finished. Global Elapsed time: {stop_global - start_global}", flush=True)
    print(f"[{datetime.now()}] [Main Process Message] EXITING...Queue final size is:", q.qsize(), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="NanoListener. A features extractor from eventalign and fast5 files: chunks of a given length will be extracted along with related sequence and other useful informations. These chunks can be used to train a model from ionic currents measurements.")
    
    parser.add_argument("-e",
                        "--eventalign_filepath",
                        required=True,
                        type=str,
                        help="--eventalign_filepath: \t a <str> indicating the fullpath for the eventalign table to be used to extract ionic currents features/chunks.")

    
    parser.add_argument("-s",
                        "--summary_table_filepath",
                        required=False,
                        type=str,
                        default=None,
                        help="--summary_table_filepath: \t a <str> with the full_path for the summary file produced by the eventalign sofware. It is needed if you want to go back to pA values. [None]")

    
    parser.add_argument("-t",
                        "--threads_n", 
                        required=False,
                        type=int,
                        default=1,
                        help="--threads_n: \t a <int> indicating the number of threads/consumers to generate in parallel (do not use more threads than available CPUs). [1]")


    parser.add_argument("-ml",
                        "--min_length_chunks",
                        required=False,
                        default=300,
                        type=int,
                        help=f"--min_length_chunks: \t a <int> indicating the minimum length of the chunks to extract. [300]")
    
    
    parser.add_argument("-l",
                        "--max_length_chunks",
                        required=False,
                        default=300,
                        type=int,
                        help=f"--max_length_chunks: \t a <int> indicating the maximum length of the chunks to extract. [300]")


    parser.add_argument("-m",
                        "--mod_pos_tables_dir",
                        required=False,
                        default=None,
                        type=str,
                        help=f"-mod_pos: \t a <str> indicating the fullpath for the directory containing tsv tables with the expected modified positions for each read (region pos0based alt_base). [None]")


    parser.add_argument("-n",
                        "--n_reads_to_process",
                        required=False,
                        default=100,
                        type=int,
                        help=f"-n_reads_to_process: \t a <int> indicating the total number of reads to be processed. [100]")

    
    parser.add_argument("-c",
                        "--clip_outliers",
                        required=False,
                        default="40-165",
                        type=str,
                        help=f"-clip_outliers: \t a <str> indicating the min and max pA values to clip with mean. [40-165]")

    parser.add_argument("-r",
                        "--ref_filepath",
                        required=True,
                        type=str,
                        help=f"--ref_filepath: \t a <str> indicating the reference fasta file used to align reads and to perform the f5c eventalign procedure.")

    
    parser.add_argument("-o",
                        "--outdir_prefix",
                        required=False,
                        default=None,
                        type=str,
                        help=f"-outdir_prefix: \t a <str> indicating the ouput directory prefix to use. By default the eventalign basename will be used as prefix. [None]")

    args = parser.parse_args()
    eventalign_filepath = args.eventalign_filepath
    summary_table_filepath = args.summary_table_filepath
    if summary_table_filepath == "None":
        summary_table_filepath = None
    threads_n = args.threads_n
    min_length_chunks = args.min_length_chunks
    max_length_chunks = args.max_length_chunks
    mod_pos_tables_dir = args.mod_pos_tables_dir
    if mod_pos_tables_dir == "None":
        mod_pos_tables_dir = None
    n_reads_to_process = args.n_reads_to_process
    clip_outliers = args.clip_outliers
    if not clip_outliers == "None":
        clip_outliers = [int(i) for i in clip_outliers.split("-")]
    else:
        clip_outliers = None
    ref_filepath = args.ref_filepath
    outdir_prefix = args.outdir_prefix

    # launch main command
    eventalign_to_chunks(eventalign_path=eventalign_filepath, 
                        summary_table_path=summary_table_filepath, 
                        threads_n=threads_n,
                        min_l=min_length_chunks,
                        max_l=max_length_chunks,
                        mod_pos_dir=mod_pos_tables_dir,
                        n_reads_to_process=n_reads_to_process,
                        clip_outliers=clip_outliers,
                        ref_filepath=ref_filepath,
                        outdir_prefix=outdir_prefix)
