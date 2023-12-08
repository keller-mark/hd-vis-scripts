import pandas as pd
import numpy as np
from os.path import join, basename, dirname
import requests
import os
import json
from peewee import *

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

  globals().update({
    "Paper": models['papers'],
    "PaperToField": models['paper_to_field'],
    "Citation": models['citations'],
  })

  with open(snakemake.input['citations_errors'], 'r') as error_f:
    error_lines = json.load(error_f)["errors"]


  with open(snakemake.input['citations_part'], 'r') as f:
    more_error_lines = []

    consume(f, offset)
    n_consumed = offset
    
    for error_batch_line_i in error_lines:

      to_consume = error_batch_line_i - n_consumed
      consume(f, to_consume)
      n_consumed += to_consume

      i = 0
      for line in f:
        citation_dict = json.loads(line)
        citation_obj = Citation(
          citation_id=citation_dict['citationid'],
          citing_corpus_id=citation_dict['citingcorpusid'],
          cited_corpus_id=citation_dict['citedcorpusid'],
          source_file=basename(snakemake.input['citations_part'])
        )
        try:
          with db.atomic():
            citation_obj.save()
        except (DataError, IntegrityError) as e:
          print(e)
          more_error_lines.append(error_batch_line_i + i)
      
        i += 1
        if i >= LINE_BATCH_SIZE:
          break

  with open(snakemake.output['citations_part'], 'w') as out_f:
    json.dump({ "more_errors": more_error_lines }, out_f)