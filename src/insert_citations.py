import pandas as pd
import numpy as np
from os.path import join, basename, dirname
import requests
import os
import json
from peewee import PostgresqlDatabase, chunked, IntegrityError

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

  with open(snakemake.input['citations_part'], 'r') as f:
    line_i = 0
    error_lines = []
    citations = []

    consume(f, offset)
    
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
        try:
          with db.atomic():
            Citation.bulk_create(citations, batch_size=BULK_BATCH_SIZE)
        except IntegrityError as e:
          print(e)
          error_lines.append(line_i)
        citations = []
      line_i += 1
    
      if line_i >= limit:
        break

  with open(snakemake.output['citations_part'], 'w') as out_f:
    json.dump({ "line_i": line_i, "errors": error_lines }, out_f)