import pandas as pd
import numpy as np
from os.path import join
import requests
import os
import json
from peewee import PostgresqlDatabase, chunked

from mod.models import get_models


if __name__ == "__main__":
  db = PostgresqlDatabase(
    snakemake.params['db_name'],
    host=snakemake.params['db_host'],
    user=snakemake.params['db_user'],
    password=snakemake.params['db_password']
  )
  models = get_models(db)

  with open(snakemake.output[0], 'w') as out_f:
    json.dump({ "models": list(models.keys()) }, out_f)
