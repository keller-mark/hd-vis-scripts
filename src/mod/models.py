import os
from peewee import *


def get_models(db):  
  class Paper(Model):
    corpus_id = CharField(unique=True, primary_key=True)
    doi = CharField(null=True)
    title = CharField()
    year = IntegerField(null=True)
    citation_count = IntegerField()
    venue = CharField()
    venue_id = CharField()
    is_domain_full = BooleanField()
    is_domain_partial = BooleanField()
    is_preprint = BooleanField()
    is_method_primary = BooleanField()
    is_method_secondary = BooleanField()
    has_dr_vis = BooleanField()
    source_file = CharField()
  
    class Meta:
      database = db
  
  class PaperToField(Model):
    corpus_id = CharField()
    field = CharField()
    source = CharField()
  
    class Meta:
      database = db
  
  class Citation(Model):
    citation_id = CharField(unique=True, primary_key=True)
    citing_corpus_id = CharField()
    cited_corpus_id = CharField()
    source_file = CharField()
  
    class Meta:
      database = db

  db.connect()

  existing_tables = db.get_tables()
  if len(existing_tables) == 0:
    db.create_tables([Paper, PaperToField, Citation])

  return ({
    "papers": Paper,
    "paper_to_field": PaperToField,
    "citations": Citation,
  })