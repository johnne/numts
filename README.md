# Numts removal

This repository contains a workflow and utility scripts for removing [NUMTs](https://en.wikipedia.org/wiki/Nuclear_mitochondrial_DNA_segment) from ASV data.

The snakemake workflow found in `Snakefile` can run on any fasta file matching `data/{dataset}.fasta` and will perform the following:

1. Downloads HMM database files from https://github.com/Hajibabaei-Lab/SCVUC_COI_metabarcode_pipeline/releases/tag/v4.3.0
2. Runs ORFfinder on the input fasta and selects the longest ORF for each sequence
3. Runs hmmscan on the longest ORFs against the HMM database
4. Filters the input sequences based on length and hmmscan bitscore according to method described in [Porter & Hajibabei, 2021](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04180-x)


## Installation

1. Clone this repository:

```bash
git clone git@github.com:johnne/numts.git
```

2. Create and activate the conda environment
```bash
conda env create -f environment.yml
conda activate numts
```

## Run the workflow
Here we assume a fasta file named `data/asvs.fasta` is present:

```bash
snakemake -c 1 results/filter/asvs/{length,bitscore}_filtered.txt
```

This will run the orffinder and hmmscan steps and output two files:
```
-results/filter/asvs/
  - length_filtered.txt
  - bitscore_filtered.txt
```

Alternatively, you can enter the name of the dataset in a yaml-formatted config file, _e.g._ `config.yml`:

```yaml
dataset: asvs
```

then run the workflow as:

```bash
snakemake -c 1 
```