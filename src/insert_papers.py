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
  db = SqliteDatabase(snakemake.input['db'])
  models = generate_models(db)

  globals().update({
    "Paper": models['papers'],
    "PaperToField": models['paper_to_field'],
    "Citation": models['citations'],
  })

  with open(snakemake.input['papers_part'], 'r') as f:
    line_i = 0
    papers = []
    fields = []

    for line in f:
      paper_dict = json.loads(line)
      paper_obj = Paper(
        corpus_id=paper_dict['corpusid'],
        doi=get_doi(paper_dict),
        title=paper_dict['title'],
        year=(int(paper_dict['year']) if paper_dict['year'] is not None else None),
        citation_count=paper_dict['citationcount'],
        venue=paper_dict['venue'],
        venue_id=paper_dict['publicationvenueid'],
        is_preprint=0,
        is_domain_full=0,
        is_domain_partial=0,
        is_method_primary=0,
        is_method_secondary=0,
        has_dr_vis=0,
        source_file=basename(snakemake.input['papers_part']),
      )
  
      if paper_dict['s2fieldsofstudy'] is not None:
        for field in paper_dict['s2fieldsofstudy']:
          ptf_obj = PaperToField(
            corpus_id=paper_dict['corpusid'],
            field=field['category'],
            source=field['source']
          )
          fields.append(ptf_obj)
          
      papers.append(paper_obj)
      if line_i % LINE_BATCH_SIZE == 0:
        with db.atomic():
          Paper.bulk_create(papers, batch_size=BULK_BATCH_SIZE)
        with db.atomic():
          PaperToField.bulk_create(fields, batch_size=BULK_BATCH_SIZE)
        papers = []
        fields = []
      line_i += 1
  
  with open(snakemake.output['papers_part'], 'w') as out_f:
    json.dump({ "line_i": line_i }), out_f