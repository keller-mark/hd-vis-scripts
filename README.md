# hd-vis-scripts

## Setup

```sh
conda env create -f environment.yml
```

Also see [DB setup](./db-root/README.md)

### Run

#### 1. Start postgres server in long-running interactive job

```
srun -p medium --pty -t 2-12:00 -n 4 --mem 16G bash
cd ~/lab
module load postgresql
pg_ctl -D $(pwd)/hd-vis-db -l logfile start
hostname
```

Copy the hostname

#### 2. Start client jobs from login node

```sh
cd ~/research/hd-vis-scripts

conda activate hd-vis

export SLURM_ACCOUNT=$(sshare -u mk596 -U | cut -d ' ' -f 1 | tail -n 1)
export HD_VIS_DB_NAME=mk596
export HD_VIS_DB_USER=mk596
export HD_VIS_DB_HOST="HOST_HERE"
export HD_VIS_DB_PASSWORD="MY_PASSWORD_HERE"
export S2_API_KEY="MY_KEY_HERE"
snakemake -j 4 --rerun-triggers mtime --keep-incomplete
```

### Run on cluster

```sh
snakemake -j 100 --rerun-triggers mtime --keep-incomplete --keep-going --latency-wait 30 --slurm \
    --default-resources slurm_account=$SLURM_ACCOUNT slurm_partition=short runtime=30
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
snakemake -j 10 --until download_s2_papers_bulk_part --slurm \
    --default-resources slurm_account=$SLURM_ACCOUNT slurm_partition=short runtime=30
snakemake -j 10 --until download_s2_citations_bulk_part --slurm \
    --default-resources slurm_account=$SLURM_ACCOUNT slurm_partition=short runtime=30

```

Navigate to the data directory:

```sh
cd /n/data1/hms/dbmi/gehlenborg/lab/hd-vis-star
```

List things in the data directory:

```sh
ls /n/data1/hms/dbmi/gehlenborg/lab/hd-vis-star/
ls -l /n/data1/hms/dbmi/gehlenborg/lab/hd-vis-star/processed/papers
ls -l /n/data1/hms/dbmi/gehlenborg/lab/hd-vis-star/processed/citations
ls -lh /n/data1/hms/dbmi/gehlenborg/lab/hd-vis-star/raw/papers
ls -lh /n/data1/hms/dbmi/gehlenborg/lab/hd-vis-star/raw/citations
```

Check file sizes in the data directory:

```sh
du -sh /n/data1/hms/dbmi/gehlenborg/lab/hd-vis-star/raw/papers/*
du -sh /n/data1/hms/dbmi/gehlenborg/lab/hd-vis-star/raw/citations/*
du -sh /n/data1/hms/dbmi/gehlenborg/lab/hd-vis-star/processed/s2.db
```

Resources:
- https://stackoverflow.com/questions/80801/how-can-i-merge-many-sqlite-databases
- https://hcc.unl.edu/docs/applications/app_specific/running_postgres/