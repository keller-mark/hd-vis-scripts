import pandas as pd
import numpy as np
from os.path import join, basename, dirname
import requests
import sqlite3
import os
import json
from peewee import SqliteDatabase, chunked
from playhouse.reflection import generate_models, print_model, print_table_sql

LINE_BATCH_SIZE = 100
BULK_BATCH_SIZE = 100

def get_doi(p_info):
  return (None if ('externalids' not in p_info or p_info['externalids'] is None) else p_info['externalids'].get('DOI', None))

if __name__ == "__main__":
  db = SqliteDatabase(snakemake.input['db'], timeout = 120)
  models = generate_models(db)

  globals().update({
    "Paper": models['papers'],
    "PaperToField": models['paper_to_field'],
    "Citation": models['citations'],
  })

  with open(snakemake.input['citations_part'], 'r') as f:
    line_i = 0
    citations = []
    
    for line in f:
      citation_dict = json.loads(line)
      citation_obj = Citation(
        citation_id=citation_dict['citationid'],
        citing_corpus_id=citation_dict['citingcorpusid'],
        cited_corpus_id=citation_dict['citedcorpusid'],
        source_file=basename(snakemake.input['citations_part'])
      )

      citations.append(citation_obj)
      if line_i % LINE_BATCH_SIZE == 0:
        with db.atomic():
          Citation.bulk_create(citations, batch_size=BULK_BATCH_SIZE)
        citations = []
      line_i += 1

  with open(snakemake.output['citations_part'], 'w') as out_f:
    json.dump({ "line_i": line_i }, out_f)