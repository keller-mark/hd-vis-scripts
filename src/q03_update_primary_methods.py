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

  # Set is_preprint field
  query = (
    Paper
      .update(is_method_primary=True)
      .where(Paper.corpus_id.in_(dr_corpus_ids))
  )
  query.execute()

  print("Successfully updated method status for primary methods.")