import pandas as pd
import numpy as np
from os.path import join, basename, dirname
import requests
import os
import json
from peewee import PostgresqlDatabase, chunked

from mod.utils import (
  get_doi,
  consume,
  LINE_BATCH_SIZE,
  BULK_BATCH_SIZE
)
from mod.models import get_models

if __name__ == "__main__":
  db = PostgresqlDatabase(
    snakemake.params['db_name'],
    host=snakemake.params['db_host'],
    user=snakemake.params['db_user'],
    password=snakemake.params['db_password']
  )
  models = get_models(db)

  offset = int(snakemake.wildcards['offset'])
  limit = int(snakemake.params['limit'])

  print(offset, limit)

  globals().update({
    "Paper": models['papers'],
    "PaperToField": models['paper_to_field'],
    "Citation": models['citations'],
  })

  with open(snakemake.input['papers_part'], 'r') as f:
    line_i = 0
    papers = []
    fields = []

    consume(f, offset)

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

      if line_i >= limit:
        break
  
  with open(snakemake.output['papers_part'], 'w') as out_f:
    json.dump({ "line_i": line_i }, out_f)