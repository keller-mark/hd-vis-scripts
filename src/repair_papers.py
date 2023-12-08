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

  with open(snakemake.input['papers_errors'], 'r') as error_f:
    error_lines = json.load(error_f)["errors"]

  with open(snakemake.input['papers_part'], 'r') as f:
    more_error_lines = []

    consume(f, offset)
    n_consumed = offset

    for error_batch_line_i in error_lines:

      to_consume = error_batch_line_i - n_consumed
      consume(f, to_consume)
      n_consumed += to_consume

      i = 0
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
        try:
          with db.atomic():
            paper_obj.save()
        except (DataError, IntegrityError) as e:
          # peewee.IntegrityError: null values where non-null expected
          # peewee.DataError: value too long for type character varying(255)
          print(e)
          more_error_lines.append(error_batch_line_i + i)
    
        if paper_dict['s2fieldsofstudy'] is not None:
          for field in paper_dict['s2fieldsofstudy']:
            try:
              with db.atomic():
                ptf_obj, created = PaperToField.get_or_create(
                  corpus_id=paper_dict['corpusid'],
                  field=field['category'],
                  source=field['source']
                )
            except (DataError, IntegrityError) as e:
              # peewee.IntegrityError: null values where non-null expected
              # peewee.DataError: value too long for type character varying(255)
              if error_batch_line_i + i not in more_error_lines:
                more_error_lines.append(error_batch_line_i + i)
        i += 1
        if i >= LINE_BATCH_SIZE:
          break
  
  with open(snakemake.output['papers_part'], 'w') as out_f:
    json.dump({ "more_errors": more_error_lines }, out_f)