import itertools


LINE_BATCH_SIZE = 100
BULK_BATCH_SIZE = 100

def get_doi(p_info):
  return (None if ('externalids' not in p_info or p_info['externalids'] is None) else p_info['externalids'].get('DOI', None))

def consume(it, n):
  if n > 0:
    return next(itertools.islice(it, n-1, n), None)