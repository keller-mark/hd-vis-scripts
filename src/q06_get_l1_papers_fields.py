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

  rows = (
    PaperToField
      .select(PaperToField.field, PaperToField.source, Paper.citation_count, Paper.year, Paper.title, Paper.venue, Paper.corpus_id, Paper.doi, Citation.cited_corpus_id)
      .where(Paper.is_preprint == False, Paper.year >= 2013, Paper.corpus_id.in_(l1_corpus_ids))
      .join(Paper, on=(PaperToField.corpus_id == Paper.corpus_id))
      .join(Citation, on=(Paper.corpus_id == Citation.citing_corpus_id))
  )

  df = pd.DataFrame(data=[
    {
      "field": row.field,
      "source": row.source,
      "year": row.paper.year,
      "citation_count": row.paper.citation_count,
      "title": row.paper.title,
      "venue": row.paper.venue,
      "corpus_id": row.paper.corpus_id,
      "doi": row.paper.doi,
      "cited_corpus_id": row.paper.citation.cited_corpus_id
    } for row in rows
  ])
  df.to_csv('q06_get_l1_papers_fields.csv')
