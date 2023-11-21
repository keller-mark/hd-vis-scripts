# hd-vis-scripts

## Setup

```sh
conda env create -f environment.yml
```

### Run

```sh
conda activate hd-vis

export S2_API_KEY="MY_KEY_HERE"
snakemake -j 4 --rerun-triggers mtime --keep-incomplete
```

