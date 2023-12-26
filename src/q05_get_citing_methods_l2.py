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

  dr_df = pd.read_csv(join("..", "dr_algorithms_with_ids.csv"))
  dr_corpus_ids = dr_df['ss_corpus_id'].astype(str).tolist()

  # Do a second query which only gets the L1 "citing_corpus_id" column
  # to pass to the L2 query's .in_() call.
  l1_corpus_ids = (
    Citation
      .select(Citation.citing_corpus_id)
      .where(Citation.cited_corpus_id.in_(dr_corpus_ids))
  )

  l2_rows = (
    Citation
      .select(Citation.citing_corpus_id, Citation.cited_corpus_id)
      .where(Citation.cited_corpus_id.in_(l1_corpus_ids))
  )
  l2_corpus_df = pd.DataFrame(data=[
    {
      "citing_corpus_id": row.citing_corpus_id,
      "cited_corpus_id": row.cited_corpus_id
    } for row in l2_rows
  ])
  l2_corpus_df.to_csv('q05_get_citing_methods_l2.csv')