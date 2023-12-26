import pandas as pd
import numpy as np
from os.path import join
import requests
import os
import json
from peewee import *

from mod.models import get_models

if __name__ == "__main__":
  db = PostgresqlDatabase(
    os.environ['HD_VIS_DB_NAME'],
    host=os.environ['HD_VIS_DB_HOST'],
    user=os.environ['HD_VIS_DB_USER'],
    password=os.environ['HD_VIS_DB_PASSWORD']
  )
  models = get_models(db)

  globals().update({
    "Paper": models['papers'],
    "PaperToField": models['paper_to_field'],
    "Citation": models['citations'],
  })

  # All papers
  all_papers_count = (
    Paper
      .select()
      .count()
  )

  # All citations count
  all_citations_count = (
    Citation
      .select()
      .count()
  )

  # All paper-to-field rows
  all_paper_to_field_count = (
    PaperToField
      .select()
      .count()
  )

  # All paper-to-field rows with source = 'external'
  all_paper_to_field_external_count = (
    PaperToField
      .select()
      .where(PaperToField.source == 'external')
      .count()
  )

  # All paper-to-field rows with source = 's2-fos-model'
  all_paper_to_field_s2_fos_model_count = (
    PaperToField
      .select()
      .where(PaperToField.source == 's2-fos-model')
      .count()
  )

  # Non-preprints
  non_preprints_count = (
    Paper
      .select()
      .where(Paper.is_preprint == False)
      .count()
  )

  # Since 2013
  all_since_2013_count = (
    Paper
      .select()
      .where(Paper.year >= 2013)
      .count()
  )

  # Non-preprints since 2013
  non_preprints_since_2013_count = (
    Paper
      .select()
      .where(Paper.is_preprint == False, Paper.year >= 2013)
      .count()
  )
    
  df = pd.DataFrame(data=[
    {
      "query": "All paper records",
      "count": all_papers_count
    },
    {
      "query": "All citation records",
      "count": all_citations_count
    },
    {
      "query": "Non-preprint paper records",
      "count": non_preprints_count
    },
    {
      "query": "Since 2013 paper records",
      "count": all_since_2013_count
    },
    {
      "query": "Non-preprint since 2013 paper records",
      "count": non_preprints_since_2013_count
    },
    {
      "query": "All paper-to-field records",
      "count": all_paper_to_field_count
    },
    {
      "query": "Paper-to-field records from external source",
      "count": all_paper_to_field_external_count
    },
    {
      "query": "Paper-to-field records from s2-fos-model",
      "count": all_paper_to_field_s2_fos_model_count
    }
  ])

  df.to_csv('q08_summary_stats.csv')
