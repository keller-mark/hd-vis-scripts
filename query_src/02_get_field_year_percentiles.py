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
    snakemake.params['db_name'],
    host=snakemake.params['db_host'],
    user=snakemake.params['db_user'],
    password=snakemake.params['db_password']
  )
  models = get_models(db)

  # Create a table with the columns:
  # Citation count | Number of papers with citation count | Year | Field
  rows = PaperToField
    .select(PaperToField.field, PaperToField.corpus_id, Paper.is_preprint, Paper.citation_count, Paper.year, fn.COUNT(PaperToField.corpus_id).alias('paper_count'))
    .where(Paper.is_preprint == False, Paper.year >= 2013)
    .join(Paper, on=(PaperToField.corpus_id == Paper.corpus_id))
    .group_by(PaperToField.field, Paper.year, Paper.citation_count)
  
  df = pd.DataFrame(data=[
    {
      "field": row.field,
      "year": row.paper.year,
      "citation_count": row.paper.citation_count,
      "paper_count": row.paper_count
    } for row in rows
  ])

  # TODO: sort by citation count increasing in each year and field
  # Calculate cumulative frequencies in percentages (Table 3 of Bornmann and Williams Scientometrics 2020)
  # - CP-IN is the percentage of papers with citation count at or below X.
  # - CP-EX is the percentage of papers with citation count below X.
  df.to_csv(snakemake.output[0])
