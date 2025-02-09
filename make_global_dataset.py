# %%
# import needed modules
import os,glob,sys, shutil
from glob import glob
from datetime import datetime
import pandas as pd

# %% [markdown]
# **STEPS:**
# 1) merge all X and y_meta datasets into a single file for each datasets group (train or test)
# 2) concat all the chunks from the merged datasets into a single file balancing every single dataset
# 3) shuffling concatenated merged datasets
# %%
# INPUT
single_datasets_list_fp = sys.argv[1]
outdir = sys.argv[2]
dataset_group = sys.argv[3]
dataset_suffix = sys.argv[4]
overwrite = sys.argv[5]

print(f"Dataset group: {dataset_group}")

df_out_folder = os.path.join(outdir, f"df_{dataset_group}")


if os.path.exists(df_out_folder):
    if str(overwrite) == "True":
        print(f"Output folder already exists, removing...")
        shutil.rmtree(df_out_folder)
    else:
        print("Provided output folder already exists...please force ovrerwriting! Exiting.")
        sys.exit()
print(f"Creating output folder for the dataset: {df_out_folder}.")
os.mkdir(df_out_folder)

# load single datasets list
single_datasets_list = pd.read_csv(single_datasets_list_fp, names=["sample_name", "dataset_path"])

print(f"Loaded single datasets list:")
print(single_datasets_list)
print()


# %%

start_time = datetime.now()
for sample in single_datasets_list.itertuples():
    sample_name = sample.sample_name
    dataset_fp = sample.dataset_path
    print(f"\nProcessing sample: {sample_name}")
    print(f"Dataset filepath: {dataset_fp}")
    # deducing c_min and c_max
    dataset_fp_basename = os.path.basename(dataset_fp)
    c_min, c_max = dataset_fp_basename.split(".")[-2].split("to")

    sample_names = [sample_name]
    run_names = [f"{sample_name}_{c_min}to{c_max}"]

    # deduce X and y NanoListener datasets
    Xs = [os.path.join(dataset_fp, f"X_{dataset_group}{dataset_suffix}.tsv")]

    ys = [os.path.join(dataset_fp, f"y_{dataset_group}_meta{dataset_suffix}.tsv")]

    output_merged_filepath = os.path.join(df_out_folder, list(set(sample_names))[0])+f".x_y_meta.{dataset_group}.tsv"
    print("Writing on the output file:", output_merged_filepath)
    
    written_lines = 0
    with open(output_merged_filepath, "w") as output_merged:
        for s,r,x_filepath,y_filepath in zip(sample_names,run_names, Xs, ys):
            print("\n#######################")
            print("Sample name:", s)
            print("Run name:", r)
            print("X:", x_filepath)
            print("y_meta:", y_filepath)
            print("Some Stats:")
            print("X rows:")
            os.system(f"ls -lh {x_filepath}")
            os.system(f"cat {x_filepath} | wc -l")
            print("y_meta rows:")
            os.system(f"ls -lh {y_filepath}")
            os.system(f"cat {y_filepath} | wc -l")
            print("Writing on output file...")
            with open(x_filepath) as X:
                with open(y_filepath) as y:
                    next(y) # skip header of metadata
                    for audio,t in zip(X,y):
                        # generate new line with mereged X and y data
                        line = audio.rstrip()+"\t"+"\t".join(t.split("\t")[1:]) # exclude for the moment index column
                        output_merged.write(line)
                        written_lines += 1

# evaluate merged dataset
print("\nWriting of intermediate output files finished...performing some checks:")
os.system(f"ls -lh {output_merged_filepath}")
os.system(f"cat {output_merged_filepath} | wc -l")
print("\nElapsed time:", datetime.now()-start_time)

# %% [markdown]
# #### Merging into a uniq dataset and shuffling...then splitting into X and y removing not useful files

# %%
concat_merged_x_y_meta_filepath = os.path.join(df_out_folder, f"concat.x_y_meta.{dataset_group}.tsv")
print("Concatenating merged X-y_meta datasets into a uniq dataframe at:", concat_merged_x_y_meta_filepath)

# %%
os.system(f"ls {os.path.join(df_out_folder, '*tsv')}")

# %%
print("Concatenating x_y datasets from all the runs...")
os.system(f"cat {os.path.join(df_out_folder, '*tsv')} > {concat_merged_x_y_meta_filepath}")
print("Finished.")

# %%
# random shuffling of concatenated and merged X-y_meta dataset
start_time = datetime.now()
concat_merged_x_y_meta_shuffled_filepath = os.path.splitext(concat_merged_x_y_meta_filepath)[0] + ".shuffled.tsv"
print("Shuffling concatenated and merged dataset...saving to:", concat_merged_x_y_meta_shuffled_filepath)
print("chunks before shuffling:")
os.system(f"cat {concat_merged_x_y_meta_filepath} | wc -l")
os.system(f"shuf {concat_merged_x_y_meta_filepath} > {concat_merged_x_y_meta_shuffled_filepath}")
print("chunks after shuffling:")
os.system(f"cat {concat_merged_x_y_meta_shuffled_filepath} | wc -l")
print("Elapsed time:", datetime.now()-start_time)

# %%
# print metadata part head before and after shuffling
print("First 3 rows before shuffling:")
os.system(f"cat {concat_merged_x_y_meta_filepath} | cut -f {int(c_max)+1}-{int(c_max)+7} | head -3")
print()
print("First 3 rows after shuffling:")
os.system(f"cat {concat_merged_x_y_meta_shuffled_filepath} | cut -f {int(c_max)+1}-{int(c_max)+7} | head -3")

# %%
print("Region distribution before shuffling:")
os.system(f"cat {concat_merged_x_y_meta_filepath} | cut -f {int(c_max)+1} | sort | uniq -c")

# %%
print("Region distribution after shuffling:")
os.system(f"cat {concat_merged_x_y_meta_shuffled_filepath} | cut -f {int(c_max)+1} | sort | uniq -c")

# %%
print("Removing concat not shuffled x_y merged dataset...")
os.system(f"rm {concat_merged_x_y_meta_filepath}")
print("Computation finished.")