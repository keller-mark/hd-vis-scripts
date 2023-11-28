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

### Troubleshooting

If presigned URLs in `data/raw/citations_meta.json` and `data/raw/papers_meta.json` have expired, then:

```sh
snakemake -j 1 --until download_s2_papers_meta
snakemake -j 1 --until download_s2_citations_meta
snakemake -j 1 --touch
```

Execute download rules on different nodes:

```sh
snakemake --until download_s2_papers_bulk_part --slurm \
    --default-resources slurm_account=mk596 slurm_partition=short

```

Resources:
- https://stackoverflow.com/questions/80801/how-can-i-merge-many-sqlite-databases
