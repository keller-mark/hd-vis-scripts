import pandas as pd
import numpy as np
from os.path import join
import requests
import sqlite3
import os
import json
from peewee import SqliteDatabase, chunked
from playhouse.reflection import generate_models, print_model, print_table_sql


if __name__ == "__main__":
  connection = sqlite3.connect(snakemake.output[0])

  cursor = connection.cursor()
  try:
    cursor.execute("DROP TABLE papers")
    cursor.execute("DROP TABLE paper_to_field")
    cursor.execute("DROP TABLE citations")
  except:
    pass

  cursor.execute("CREATE TABLE papers (corpus_id TEXT, doi TEXT, title TEXT, year INTEGER, citation_count INTEGER, venue TEXT, venue_id TEXT, is_domain_full INTEGER, is_domain_partial INTEGER, is_preprint INTEGER, is_method_primary INTEGER, is_method_secondary INTEGER, has_dr_vis INTEGER, source_file TEXT)")
  cursor.execute("CREATE TABLE paper_to_field (corpus_id TEXT, field TEXT, source TEXT)")
  cursor.execute("CREATE TABLE citations (citation_id TEXT, citing_corpus_id TEXT, cited_corpus_id TEXT, source_file TEXT)")
