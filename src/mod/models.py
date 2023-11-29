import os
from peewee import *


def get_models(db, force=False):  
  class Paper(Model):
    corpus_id = CharField(unique=True, primary_key=True)
    doi = CharField(null=True)
    title = CharField(null=True)
    year = IntegerField(null=True)
    citation_count = IntegerField(null=True)
    venue = CharField(null=True)
    venue_id = CharField(null=True)
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
    citing_corpus_id = CharField(null=True)
    cited_corpus_id = CharField(null=True)
    source_file = CharField()
  
    class Meta:
      database = db

  db.connect()

  existing_tables = db.get_tables()
  if len(existing_tables) == 0 or force:
    if len(existing_tables) > 0 and force:
      db.drop_tables([Paper, PaperToField, Citation])
    db.create_tables([Paper, PaperToField, Citation])

  return ({
    "papers": Paper,
    "paper_to_field": PaperToField,
    "citations": Citation,
  })