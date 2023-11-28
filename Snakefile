include: "./common.smk"

from os.path import join
import os

NUM_PAPERS_FILES = 2 # TODO: revert to 30
NUM_CITATIONS_FILES = 2 # TODO: revert to 30

RELEASE_ID = "2023-11-07"

# 2.4B citations in 30 1.5GB files
# 200M papers in 30 8.5 GB files
PAPER_LIMIT = 1_000_000
PAPER_PARTS_PER_FILE = 10
PAPER_OFFSETS = [i * PAPER_LIMIT for i in range(PAPER_PARTS_PER_FILE)]
CITATION_LIMIT = 10_000_000
CITATION_PARTS_PER_FILE = 10
CITATION_OFFSETS = [i * CITATION_LIMIT for i in range(CITATION_PARTS_PER_FILE)]


rule all:
  input:
    join(PROCESSED_DIR, "s2.db"),
    join(RAW_DIR, "papers_meta.json"),
    join(RAW_DIR, "citations_meta.json"),
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
    join(PROCESSED_DIR, "s2.db")
  script:
    join(SRC_DIR, "create_db.py")

# Insert papers into the SQLite DB
rule insert_papers:
  input:
    db=join(PROCESSED_DIR, "s2.db"),
    papers_part=join(RAW_DIR, "papers", "part{file_i}.jsonl"),
  params:
    limit=PAPER_LIMIT
  output:
    papers_part=join(PROCESSED_DIR, "papers", "part{file_i}_offset{offset}_complete.json"),
  script:
    join(SRC_DIR, "insert_papers.py")

rule insert_citations:
  input:
    db=join(PROCESSED_DIR, "s2.db"),
    citations_part=join(RAW_DIR, "citations", "part{file_i}.jsonl"),
  params:
    limit=CITATION_LIMIT
  output:
    citations_part=join(PROCESSED_DIR, "citations", "part{file_i}_offset{offset}_complete.json"),
  script:
    join(SRC_DIR, "insert_citations.py")

# Un-gz the bulk data
rule unzip_s2_papers_bulk_part:
  input:
    join(RAW_DIR, "papers", "part{file_i}.jsonl.gz")
  output:
    join(RAW_DIR, "papers", "part{file_i}.jsonl")
  shell:
    """
    gunzip -c {input} > {output}
    """

rule unzip_s2_citations_bulk_part:
  input:
    join(RAW_DIR, "citations", "part{file_i}.jsonl.gz")
  output:
    join(RAW_DIR, "citations", "part{file_i}.jsonl")
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
    partition="short",
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
    partition="short",
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