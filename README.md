# NanoListener

<p align="center">
<img src="https://github.com/F0nz0/NanoListener/blob/master/NatMet_NanoListener_NanoSpeech_Fig1_HD_GitHub.png" align="center">
</p>

<p align="justify">
NanoListener is a small suite of Python scripts to create training datasets for training direct-RNA modification-aware basecaller models. 
Raw signals stored in fast5, or converted pod5 files, are basecalled via modification-unaware basecalling models and mapped against reference sequences. 
Alignments are filtered collecting reads falling on regions of interest and respecting strict criteria, using NanoListener accessory scripts. 
Ionic currents needs to be re-squiggled onto the reference via f5c eventalign. 
  NanoListener takes advantage of these re-squiggled signals to retrieve context information and extract random chunks of electric measurements trying to avoid problematic 
  regions (in red) due to modified nucleotides. 
  It randomly exploits as anchors, low-noise “unmodified” flanking regions to extract whole raw signals from fast5, filling in turn, re-squiggling gaps. 
  For every extracted chunk, the deduced k-mer is marked for modified nucleotides in accordance with tables of per-read modification positions. 
Finally, a training dataset is returned, consisting in pairs of chunks/annotated k-mers where a padding is added to currents measurements making these, uniform in length. 
NanoListener needs positional information, at a per-read level, to annotate all the output k-mers. 
  These can be obtained either using synthetic in-vitro transcribed molecules or using other orthogonal methods. 
</p>

## **Required Softwares**:
NanoSpeech uses internally (and not) some software that should be installed preferably into a new conda enviroment. \
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

		# create virtual environment inside Conda NanoSpeech env
		python3 -m venv NanoListener_venv

    4) Activate the venv:
	
	    source NanoListener_venv/bin/activate

    5) Upgrade pip version:
    	
    	    python3 -m pip install --upgrade pip
    
    6) Install wheel package via pip:
    	
    	    pip install wheel
    
    7) Install required Python packages using the requirements_latest_tf_2_7_0.txt file:
    
            python -m pip install -r requirements.txt

## **Basic Usage**:

