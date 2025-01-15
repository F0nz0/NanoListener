# import basic modules
import os, argparse
from glob import glob
import pandas as pd
import numpy as np
from datetime import datetime
from sklearn.model_selection import train_test_split

# define a useful function to make datasets of diffent chunks uniform in length
def lengthen_chunks(df, targent_len=5000, input_chunk_len=2500):
    '''
    A function to make NanoListener dataframes longer to a target length (for instance, from 2.5k to 5k).
    '''
    if not targent_len > input_chunk_len:
        raise Exception("Error! Please target len have to be longer than starting chunk length.")
    print(f"[{datetime.now()}] Making the input dataset from chunks with length {input_chunk_len} to {targent_len}.", flush=True)
    print(f"[{datetime.now()}] Starting dataframe shape (n_rows, n_cols): {df.shape}.", flush=True)
    # reset index
    df.reset_index(inplace=True, drop=True)
    # split currents chunks from metadata on the last columns
    chunks = df.iloc[:,0:input_chunk_len]
    metadata = df.iloc[:,input_chunk_len:]
    # create a df with NaNs
    NaNs = pd.DataFrame(np.full((df.shape[0],int(targent_len-input_chunk_len)), np.nan), columns=[i for i in range(input_chunk_len, targent_len)])
    chunks = pd.concat([chunks, NaNs], axis=1)
    # delete df and NaNs to spare memory
    del(df)
    del(NaNs)
    # concat and return final longer df
    df = pd.concat([chunks, metadata], axis=1)
    df.columns = [i for i in range(df.shape[1])]
    print(f"[{datetime.now()}] Dataframe shape after the 'elongation' of chunks (n_rows, n_cols): {df.shape}.", flush=True)
    return df

# define main function
def nanolistener_cons_out_filt(nanolistener_chunks_dir, target_max_len=5000, mod_base = None, 
                             kmer_len_interval = [5,100], balancing_strategy = 0, 
                             train_test_ratio = 0.03, sample_n = None, cons_out_n = None, outsuffix=None):
    '''
    A function to filter and balance consumers outputs from NanoListener before 
    NanoSpeech model training. It will merge selected cons. outputs into X and y_meta
    dataframes and it will split them into X_train, X_test, y_train_meta and y_test_meta
    datasets on disk at the same folderpath of the nanolistener chunks directory.
    '''
    start_time = datetime.now()
    if not outsuffix:
        # define outputs filepaths basing on NanoListener chunks directory
        X_train_output_path = os.path.join(nanolistener_chunks_dir, "X_train.tsv")
        y_train_meta_output_path = os.path.join(nanolistener_chunks_dir, "y_train_meta.tsv")
        X_test_output_path = os.path.join(nanolistener_chunks_dir, "X_test.tsv")
        y_test_meta_output_path = os.path.join(nanolistener_chunks_dir, "y_test_meta.tsv")
    else:
        # define outputs filepaths basing on NanoListener chunks directory
        X_train_output_path = os.path.join(nanolistener_chunks_dir, f"X_train_{outsuffix}.tsv")
        y_train_meta_output_path = os.path.join(nanolistener_chunks_dir, f"y_train_meta_{outsuffix}.tsv")
        X_test_output_path = os.path.join(nanolistener_chunks_dir, f"X_test_{outsuffix}.tsv")
        y_test_meta_output_path = os.path.join(nanolistener_chunks_dir, f"y_test_meta_{outsuffix}.tsv")

    # create the four new empty output files (delete if exist)
    with open(X_train_output_path, 'w') as out:
        pass
    with open(y_train_meta_output_path, 'w') as out:
        pass
    with open(X_test_output_path, 'w') as out:
        pass
    with open(y_test_meta_output_path, 'w') as out:
        pass

    # detect consumer outputs into the provided folder
    # all consumer outputs
    cons_list = glob(os.path.join(nanolistener_chunks_dir, "cons_*.tsv"))
    # if required work on a subset of the consumer tsv output tables
    if cons_out_n:
        print(f"[{datetime.now()}] Working on the first {cons_out_n} consumer outputs.")
        cons_list = [os.path.join(nanolistener_chunks_dir, f"cons_{id}.tsv") for id in range(1,cons_out_n+1)]
    print(f"[{datetime.now()}] Consumer output tsv used:\n{[os.path.basename(_) for _ in cons_list]}")
    for c,cons_out in enumerate(cons_list):
        try:
            print(f"\n##########################################",flush=True)
            print(f"##########################################",flush=True)
            print(f"[{datetime.now()}] Consumer Output tsv {c+1} out of {len(cons_list)}.\n\t- Elapsed time: {datetime.now()-start_time}.")
            print(f"[{datetime.now()}] Processing consumers output at: {cons_out}.", flush=True)
            # detect chunks length
            lens = [int(i) for i in os.path.basename(nanolistener_chunks_dir).split(".")[-2].split("to")]
            print(f"[{datetime.now()}] Chunks min and Max lengths detected:", lens, flush=True)
            # load data chuncks
            print(f"[{datetime.now()}] Loading consumer output in memory.", flush=True)
            df = pd.read_table(cons_out, header=None)
            # scale to 5k-long chunks if needed
            if lens[1] < target_max_len:
                print(f"[{datetime.now()}] Lengthening of chunks needed from {lens[1]} to {target_max_len}.", flush=True)
                df = lengthen_chunks(df, targent_len=target_max_len, input_chunk_len=lens[1])
            # select only chunks containing modified base if requested
            if mod_base:
                print(f"[{datetime.now()}] Retaining only chunks with at least one {mod_base}. Starting shape: {df.shape}.")
                df = df[df.iloc[:,-1].str.contains(mod_base)].copy()
                print(f"[{datetime.now()}] df shape after mod_base containing chunks selection: {df.shape}.")
            # remove duplicates
            df.drop_duplicates(inplace=True)
            df.reset_index(inplace=True, drop=True)
            print(f"[{datetime.now()}] df shape after duplicates removal:", df.shape, flush=True)
            # retain only lengths within an interval
            if kmer_len_interval:
                print(f"[{datetime.now()}] Filtering out chunks outsite the given kmer length [min,max]: {kmer_len_interval}.", flush=True)
                df = df[(df.iloc[:,-1].apply(len)>=kmer_len_interval[0])&(df.iloc[:,-1].apply(len)<=kmer_len_interval[1])].copy()
                print(f"[{datetime.now()}] df shape after seqs length filtering:", df.shape, flush=True)
            # balancing
            # 0: no balancing, 1: balancing on regions
            if balancing_strategy != 0:
                print(f"[{datetime.now()}] Balancing of chunks requested with strategy: {balancing_strategy}", flush=True)
                if balancing_strategy == 1:
                    print(f"[{datetime.now()}] Balancing chunks on regions...", flush=True)
                    df = df.groupby(df.shape[1]-8)
                    df = df.apply(lambda x: x.sample(df.size().min())).reset_index(drop=True)
            print(f"[{datetime.now()}] df stats after filtering (start coords stats):\n", df.groupby(df.shape[1]-8)[df.shape[1]-7].describe(), flush=True)
            print(f"[{datetime.now()}] df shape after filtering:", df.shape, flush=True)
            # if needed select a ranodom sample with n chunks
            if sample_n:
                #print(type(sample_n)) #### DEVELOPMENT/TEST
                if sample_n < df.shape[0]:
                    print(f"[{datetime.now()}] Sampling a random sample of chunks: {sample_n}")
                    df = df.sample(n=sample_n, ignore_index=True)
                    print(f"[{datetime.now()}] df shape after random sampling: {df.shape}")
                else:
                    print(f"[{datetime.now()}] Provided 'sample_n' {sample_n} less than chunks for this consumer. All the chunks will be retained.")
            # split X and y_meta
            print(f"[{datetime.now()}] Shuffling and splitting into data and metadata dataframes.", flush=True)
            X = df.iloc[:,:target_max_len].values
            y_meta = df.iloc[:,target_max_len:]
            y_meta = y_meta[y_meta.columns[[0,1,2,5,6,7]]].copy()
            y_meta.columns = [c for c in range(target_max_len,target_max_len+6)]
            # delete df to spare memory
            del(df)
            # split training and test sets
            print(f"[{datetime.now()}] Splitting into training and test datasets.", flush=True)
            X_train, X_test, y_train_meta, y_test_meta = train_test_split(X, y_meta, test_size=train_test_ratio, shuffle=True, random_state=24)
            # delete to spare memory
            del(X, y_meta)
            print("X_train, X_test, y_train_meta, y_test_meta shapes", X_train.shape, X_test.shape, y_train_meta.shape, y_test_meta.shape)
            # save to disk training a test datasets
            print(f"[{datetime.now()}] Appending to disk at paths:", flush=True)
            print(f"\tX_train --> {X_train_output_path}", flush=True)
            print(f"\tX_test --> {X_test_output_path}", flush=True)
            print(f"\ty_train_meta --> {y_train_meta_output_path}", flush=True)
            print(f"\ty_test_meta --> {y_test_meta_output_path}", flush=True)
            # save Xtrain
            with open(X_train_output_path, 'a') as out:
                np.savetxt(out, X_train, delimiter='\t')
            # save y_train_meta (header only first consumer)
            if c == 0:
                y_train_meta.to_csv(y_train_meta_output_path, sep='\t', mode="a")
            else:
                y_train_meta.to_csv(y_train_meta_output_path, sep='\t', mode="a", header=False)
            # save Xtest
            with open(X_test_output_path, 'a') as out:
                np.savetxt(out, X_test, delimiter='\t')
            # save y_test_meta (header only first consumer)
            if c == 0:    
                y_test_meta.to_csv(y_test_meta_output_path, sep='\t', mode="a")
            else:
                y_test_meta.to_csv(y_test_meta_output_path, sep='\t', mode="a", header=False)
        except Exception as e:
            print(f"[{datetime.now()}] ERROR!!! PROBLEM WITH CURRENT NanoListener CONSUMER {c+1} THIS CONSUMER OUTPUT WILL BE SKIPPED.", flush=True) 
            print(f"[{datetime.now()}] Expetion Raised:", e, flush=True)
            continue
    print(f"[{datetime.now()}] Computation finished. Elapsed time: {datetime.now()-start_time}.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="nanolistener_cons_out_filt.py. A program to filter and balance consumers outputs produced by NanoListener before NanoSpeech model training. It merges all the selected consumer outputs into X and y_meta dataframes and it then splits them into X_train, X_test, y_train_meta and y_test_meta datasets saving these on disk at the same folder-path of the NanoListener chunks directory.")

    parser.add_argument("-nld",
                        "--nanolistener_chunks_dir",
                        required=True,
                        type=str,
                        help="--nanolistener_chunks_dir: \t a <str> indicating the fullpath for the NanoListener output folder.")

    
    parser.add_argument("-tml",
                        "--target_max_len",
                        required=True,
                        type=int,
                        help="--target_max_len: \t a <int> with target max length of the chunks. It will try to fit to this value the input chunks.")

    
    parser.add_argument("-mb",
                        "--mod_base", 
                        required=False,
                        default=None,
                        help="--mod_base: \t a <str> indicating the modified base to be searched into output kmer (at least one mod_base needed to mantain the chunks). [None]")


    parser.add_argument("-kli",
                        "--kmer_len_interval",
                        required=False,
                        type=str,
                        default="5-100",
                        help=f"--kmer_len_interval: \t a <str> indicating the min and maximum allowed length of output kmers. [5-100]")
    
    
    parser.add_argument("-bal",
                        "--balancing_strategy",
                        required=False,
                        default=1,
                        type=int,
                        help=f"--balancing_strategy: \t a <int> indicating (0: no balancing, 1: balancing on regions). [1]")


    parser.add_argument("-ttr",
                        "--train_test_ratio",
                        required=False,
                        default=0.03,
                        type=float,
                        help=f"-train_test_ratio: \t a <float> indicating the ratio of chunks to be splitted into test set. [0.03]")


    parser.add_argument("-n",
                        "--sample_n",
                        required=False,
                        default=None,
                        type=int,
                        help=f"-sample_n: \t a <int> indicating the number of chunks to sample randomly. [None]")

    
    parser.add_argument("-c",
                        "--cons_out_n",
                        required=False,
                        default=None,
                        type=int,
                        help=f"-cons_out_n: \t a <int> indicating how many consumer outputs to considers. [None]")


    parser.add_argument("-o",
                        "--outsuffix",
                        required=False,
                        default=None,
                        type=str,
                        help=f"-outsuffix: \t a <str> with a suffix to be added to the output X and y merged and filtered datasets. [None]")

    args = parser.parse_args()
    nanolistener_chunks_dir = args.nanolistener_chunks_dir
    target_max_len = args.target_max_len
    mod_base = args.mod_base
    if mod_base == "None":
        mod_base = None
    kmer_len_interval = args.kmer_len_interval
    if not kmer_len_interval == "None":
        kmer_len_interval = [int(i) for i in kmer_len_interval.split("-")]
    else:
        kmer_len_interval = None
    balancing_strategy = args.balancing_strategy
    train_test_ratio = args.train_test_ratio
    sample_n = args.sample_n
    if sample_n == "None":
        sample_n = None
    cons_out_n = args.cons_out_n
    if cons_out_n == "None":
        cons_out_n = None
    outsuffix = args.outsuffix
    if outsuffix == "None":
        outsuffix = None

    # print some starting info related to version, used program and to the input arguments
    print(f"[{datetime.now()}] nanolistener_cons_out_filt.py program.", flush=True)
    print(f"[{datetime.now()}] Input arguments:", flush=True)
    for argument in args.__dict__.keys():
        print(f"\t- {argument} --> {args.__dict__[argument]}", flush=True)

    # launch main command
    nanolistener_cons_out_filt(nanolistener_chunks_dir, target_max_len=target_max_len, mod_base = mod_base, 
                               kmer_len_interval = kmer_len_interval, balancing_strategy = balancing_strategy, 
                               train_test_ratio = train_test_ratio, sample_n = sample_n, cons_out_n = cons_out_n,
                               outsuffix=outsuffix)