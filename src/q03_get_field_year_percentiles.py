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

  # Create a table with the columns:
  # Citation count | Number of papers with citation count | Year | Field
  rows = (
    PaperToField
      .select(PaperToField.field, PaperToField.source, Paper.citation_count, Paper.year, fn.COUNT(PaperToField.corpus_id).alias('paper_count'))
      .where(Paper.is_preprint == False, Paper.year >= 2013)
      .join(Paper, on=(PaperToField.corpus_id == Paper.corpus_id))
      .group_by(PaperToField.field, PaperToField.source, Paper.year, Paper.citation_count)
  )
    
  df = pd.DataFrame(data=[
    {
      "field": row.field,
      "source": row.source,
      "year": row.paper.year,
      "citation_count": row.paper.citation_count,
      "paper_count": row.paper_count
    } for row in rows
  ])

  # TODO: sort by citation count increasing in each year and field
  # Calculate cumulative frequencies in percentages (Table 3 of Bornmann and Williams Scientometrics 2020)
  # - CP-IN is the percentage of papers with citation count at or below X.
  # - CP-EX is the percentage of papers with citation count below X.
  df.to_csv('q03_get_field_year_percentiles.csv')
