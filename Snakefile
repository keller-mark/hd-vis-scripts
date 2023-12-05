include: "./common.smk"

from os.path import join
import os

NUM_PAPERS_FILES = 30
NUM_CITATIONS_FILES = 30

RELEASE_ID = "2023-11-07"

# 2.4B citations in 30 1.5GB files
# 200M papers in 30 8.5 GB files
PAPER_LIMIT = 1_000_000
PAPER_PARTS_PER_FILE = 10
PAPER_OFFSETS = [i * PAPER_LIMIT for i in range(PAPER_PARTS_PER_FILE)]
CITATION_LIMIT = 2_500_000
CITATION_PARTS_PER_FILE = 40
CITATION_OFFSETS = [i * CITATION_LIMIT for i in range(CITATION_PARTS_PER_FILE)]

envvars:
  "S2_API_KEY",
  "HD_VIS_DB_NAME",
  "HD_VIS_DB_HOST",
  "HD_VIS_DB_USER",
  "HD_VIS_DB_PASSWORD"

rule all:
  input:
    expand(
      join(PROCESSED_DIR, "papers", "part{file_i}_offset{offset}_complete.json"),
      file_i=range(NUM_PAPERS_FILES),
      offset=PAPER_OFFSETS,
    ),
    expand(
      join(PROCESSED_DIR, "citations", "part{file_i}_offset{offset}_complete.json"),
      file_i=range(NUM_CITATIONS_FILES),
      offset=CITATION_OFFSETS,
    )

# Create the SQLite DB
rule create_db:
  output:
    join(PROCESSED_DIR, "db.json")
  params:
    db_name=os.environ["HD_VIS_DB_NAME"],
    db_host=os.environ["HD_VIS_DB_HOST"],
    db_user=os.environ["HD_VIS_DB_USER"],
    db_password=os.environ["HD_VIS_DB_PASSWORD"]
  resources:
    slurm_partition="short",
    runtime=60*1, # 1 hour
    mem_mb=4_000, # 4 GB
    cpus_per_task=4
  script:
    join(SRC_DIR, "create_db.py")

# Insert papers into the SQLite DB
rule insert_papers:
  input:
    db=join(PROCESSED_DIR, "db.json"),
    papers_part=join(RAW_DIR, "papers", "part{file_i}.jsonl"),
  params:
    limit=PAPER_LIMIT,
    db_name=os.environ["HD_VIS_DB_NAME"],
    db_host=os.environ["HD_VIS_DB_HOST"],
    db_user=os.environ["HD_VIS_DB_USER"],
    db_password=os.environ["HD_VIS_DB_PASSWORD"]
  output:
    papers_part=join(PROCESSED_DIR, "papers", "part{file_i}_offset{offset}_complete.json"),
  resources:
    slurm_partition="medium",
    runtime=60*24, # 1 day
    mem_mb=2_000, # 2 GB
    cpus_per_task=1
  script:
    join(SRC_DIR, "insert_papers.py")

rule insert_citations:
  input:
    db=join(PROCESSED_DIR, "db.json"),
    citations_part=join(RAW_DIR, "citations", "part{file_i}.jsonl"),
  params:
    limit=CITATION_LIMIT,
    db_name=os.environ["HD_VIS_DB_NAME"],
    db_host=os.environ["HD_VIS_DB_HOST"],
    db_user=os.environ["HD_VIS_DB_USER"],
    db_password=os.environ["HD_VIS_DB_PASSWORD"]
  output:
    citations_part=join(PROCESSED_DIR, "citations", "part{file_i}_offset{offset}_complete.json"),
  resources:
    slurm_partition="medium",
    runtime=60*24, # 1 day
    mem_mb=2_000, # 2 GB
    cpus_per_task=1
  script:
    join(SRC_DIR, "insert_citations.py")

# Un-gz the bulk data
rule unzip_s2_papers_bulk_part:
  input:
    join(RAW_DIR, "papers", "part{file_i}.jsonl.gz")
  output:
    join(RAW_DIR, "papers", "part{file_i}.jsonl")
  resources:
    slurm_partition="short",
    runtime=60*1, # 1 hour
    mem_mb=16_000, # 16 GB
    cpus_per_task=4
  shell:
    """
    gunzip -c {input} > {output}
    """

rule unzip_s2_citations_bulk_part:
  input:
    join(RAW_DIR, "citations", "part{file_i}.jsonl.gz")
  output:
    join(RAW_DIR, "citations", "part{file_i}.jsonl")
  resources:
    slurm_partition="short",
    runtime=60*1, # 1 hour
    mem_mb=16_000, # 16 GB
    cpus_per_task=4
  shell:
    """
    gunzip -c {input} > {output}
    """

# Download the bulk data
rule download_s2_papers_bulk_part:
  input:
    join(RAW_DIR, "papers_meta.json")
  output:
    join(RAW_DIR, "papers", "part{file_i}.jsonl.gz")
  params:
    s2_api_key=os.environ["S2_API_KEY"]
  resources:
    slurm_partition="short",
    runtime=60*8, # 8 hours
    mem_mb=16_000, # 16 GB
    cpus_per_task=4
  shell:
    """
    curl -o {output} -H "x-api-key: {params.s2_api_key}" -L "$(cat {input} | jq -r '.files[{wildcards.file_i}]')"
    """

rule download_s2_citations_bulk_part:
  input:
    join(RAW_DIR, "citations_meta.json")
  output:
    join(RAW_DIR, "citations", "part{file_i}.jsonl.gz")
  params:
    s2_api_key=os.environ["S2_API_KEY"]
  resources:
    slurm_partition="short",
    runtime=60*8, # 8 hours
    mem_mb=16_000, # 16 GB
    cpus_per_task=4
  shell:
    """
    curl -o {output} -H "x-api-key: {params.s2_api_key}" -L "$(cat {input} | jq -r '.files[{wildcards.file_i}]')"
    """

# Download metadata containing URLs for bulk download
rule download_s2_papers_meta:
  output:
    join(RAW_DIR, "papers_meta.json")
  params:
    file_url=f"http://api.semanticscholar.org/datasets/v1/release/{RELEASE_ID}/dataset/papers",
    s2_api_key=os.environ["S2_API_KEY"]
  shell:
    '''
    curl -o {output} -H "x-api-key: {params.s2_api_key}" -L "{params.file_url}"
    '''

rule download_s2_citations_meta:
  output:
    join(RAW_DIR, "citations_meta.json")
  params:
    file_url=f"http://api.semanticscholar.org/datasets/v1/release/{RELEASE_ID}/dataset/citations",
    s2_api_key=os.environ["S2_API_KEY"]
  shell:
    '''
    curl -o {output} -H "x-api-key: {params.s2_api_key}" -L "{params.file_url}"
    '''