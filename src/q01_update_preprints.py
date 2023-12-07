import pandas as pd
import numpy as np
from os.path import join
import requests
import os
import json
from peewee import *

from mod.models import get_models

PREPRINT_VENUES = [
  "bioRxiv",
  "medRxiv",
  "arXiv.org",
  "PeerJ Preprints"
]

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

  # Set is_preprint field
  query = Paper
    .update(is_preprint=True)
    .where(Paper.venue.in_(PREPRINT_VENUES))
  query.execute()

  print("Successfully updated preprint status.")

  
