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
    snakemake.params['db_name'],
    host=snakemake.params['db_host'],
    user=snakemake.params['db_user'],
    password=snakemake.params['db_password']
  )
  models = get_models(db, force=True)

  # Set is_preprint field
  query = Paper
    .update(is_preprint=True)
    .where(Paper.venue.in_(PREPRINT_VENUES))
  query.execute()

  