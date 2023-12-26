# hd-vis-scripts

This repository contains scripts for bibliometric analysis of literature which cites dimensionality reduction method papers.
These scripts download the entire Semantic Scholar Academic Graph (S2AG) and insert all paper and citation records into a PostgreSQL database.
Jupyter notebooks plot the results of queries against this database.
The scripts to download Semantic Scholar bulk files and insert records into the SQL database (everything upstream of querying the database) are general and could be reused for other bibliometric analyses.


## Setup

```sh
conda env create -f environment.yml
```

Also see [DB setup](./db-root/README.md)

Request a semantic scholar API key from https://www.semanticscholar.org/product/api#api-key

### Troubleshooting

To install `altair_saver`, you will need to have NodeJS (tested with v20) and NPM (v8 or earlier) installed.

If the installation fails, run the NPM install command manually:

```sh
npm install vega-lite vega-cli canvas
```

On macOS, the [canvas](https://github.com/Automattic/node-canvas#installation) package may require some Homebrew dependencies.

## Run

### 1. Start postgres server in long-running interactive job

Adjust the time (`-t`) that will be requested for the `srun` command based on the intended client jobs that will be run against the server.

```sh
srun -p medium --pty -t 2-12:00 -n 4 --mem 16G bash
cd ~/lab
module load postgresql
pg_ctl -D $(pwd)/hd-vis-db -l logfile start
hostname
```

Copy the hostname

### 2. Start client jobs from login node

This took approximately one week to run (with 50-100 jobs running simultaneously as specified with `-j`).

This command will download the bulk database files, unzip them, and then begin reading through the JSONL files and inserting records into the Postgres database. 

```sh
cd ~/research/hd-vis-scripts

conda activate hd-vis

export SLURM_ACCOUNT=$(sshare -u mk596 -U | cut -d ' ' -f 1 | tail -n 1)
export HD_VIS_DB_NAME=mk596
export HD_VIS_DB_USER=mk596
export HD_VIS_DB_HOST="HOST_HERE"
export HD_VIS_DB_PASSWORD="MY_PASSWORD_HERE"
export S2_API_KEY="MY_KEY_HERE"

snakemake -j 50 --rerun-triggers mtime --keep-incomplete --keep-going --latency-wait 30 --slurm \
    --default-resources slurm_account=$SLURM_ACCOUNT slurm_partition=short runtime=30
```

The `.jsonl.gz` files are 220 GB total. After unzipping the `.jsonl` files are 840 GB total.



#### 2.1. Repair errors

The insertion 

This took approximately one week to run (with 50-100 jobs running simultaneously as specified with `-j`).

```sh
cd ~/research/hd-vis-scripts

conda activate hd-vis

export SLURM_ACCOUNT=$(sshare -u mk596 -U | cut -d ' ' -f 1 | tail -n 1)
export HD_VIS_DB_NAME=mk596
export HD_VIS_DB_USER=mk596
export HD_VIS_DB_HOST="HOST_HERE"
export HD_VIS_DB_PASSWORD="MY_PASSWORD_HERE"
export S2_API_KEY="MY_KEY_HERE"

snakemake repair_all -j 50 --rerun-triggers mtime --keep-incomplete --keep-going --latency-wait 30 --slurm \
    --default-resources slurm_account=$SLURM_ACCOUNT slurm_partition=short runtime=180
```

After inserting and repairing errors the database directory should be approximately 142 GB:

```sh
du -sh ~/lab/hd-vis-db
# 142G	/home/mk596/lab/hd-vis-db
```

### 3. Run queries against running postgres server

These queries take approximately 3 hours to run.

```sh
srun -p interactive --pty -t 8:00:00 -n 4 --mem 32G bash

cd ~/research/hd-vis-scripts
conda activate hd-vis

export HD_VIS_DB_NAME=mk596
export HD_VIS_DB_USER=mk596
export HD_VIS_DB_HOST="HOST_HERE"
export HD_VIS_DB_PASSWORD="MY_PASSWORD_HERE"

cd src

python q01_update_preprints.py
python q02_update_papers.py
python q03_get_field_year_percentiles.py
python q04_get_citing_methods_l1.py
python q05_get_citing_methods_l2.py
python q06_get_l1_papers_fields.py
python q07_get_l1_papers_fields.py
python q08_summary_stats.py
```

The resulting extracted CSV files are approximately 2 GB total.

### More troubleshooting

If presigned URLs in `data/raw/citations_meta.json` and `data/raw/papers_meta.json` have expired, then:

```sh
snakemake -j 1 --until download_s2_papers_meta
snakemake -j 1 --until download_s2_citations_meta
snakemake -j 1 --touch
```

Execute download rules on different nodes:

```sh
snakemake -j 10 --until download_s2_papers_bulk_part --slurm \
    --default-resources slurm_account=$SLURM_ACCOUNT slurm_partition=short runtime=30
snakemake -j 10 --until download_s2_citations_bulk_part --slurm \
    --default-resources slurm_account=$SLURM_ACCOUNT slurm_partition=short runtime=30
```


Resources:
- https://stackoverflow.com/questions/80801/how-can-i-merge-many-sqlite-databases
- https://hcc.unl.edu/docs/applications/app_specific/running_postgres/
- https://stackoverflow.com/questions/1237725/copying-postgresql-database-to-another-server
