#!/usr/bin/env python
from __future__ import division
import re
import os
import sys
import argparse
import pipeline

OPT_DEFAULTS = {}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""
FASTQ_NAME_REGEX = r'(.+)_[12]\.(fastq|fq)(\.gz)?$'

def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('ref', metavar='reference.fa',
    help='The reference genome.')
  parser.add_argument('-f', '--fastq-dir', required=True,
    help='Directory containing the fastq files. The filenames should be in the format '
         '{sampleid}_1.fastq.gz.')
  parser.add_argument('-t', '--families-table', required=True,
    help='A file containing the list of families. Format: A tab-delimited file with one line per '
         'each family. The first column is the id of the family, and the following columns are the '
         'ids of the samples in that family.')
  parser.add_argument('-F', '--family-id')
  parser.add_argument('pipeargs', metavar='...', nargs=argparse.REMAINDER,
    help='The remaining arguments will be passed directly to pipeline.py. Make sure to put them '
         'at the end, after the arguments for this script.')

  args = parser.parse_args(argv[1:])

  fastqs = get_fastqs(args.fastq_dir)

  families = read_families(args.families_table)
  for family_id, sample_ids in families.items():
    print family_id
    for sample_id in sample_ids:
      print '\t'+sample_id
      if sample_id in fastqs:
        fastq_pair = fastqs[sample_id]
        if len(fastq_pair) != 2:
          raise PipefamError('Wrong number of fastq\'s associated with sample "{}": {}'
                             .format(sample_id, ', '.join(fastq_pair)))
      else:
        raise PipefamError('No fastq\'s found for sample ""'.format(sample_id))
      fastq1 = os.path.join(args.fastq_dir, fastq_pair[0])
      fastq2 = os.path.join(args.fastq_dir, fastq_pair[1])
      print '\t\t'+fastq1
      print '\t\t'+fastq2



def read_families(families_filepath):
  """Read the families tsv file and return a dict mapping family ID's to lists of sample ID's in
  that family."""
  families = {}
  with open(families_filepath) as families_file:
    for line in families_file:
      fields = line.rstrip('\r\n').split('\t')
      if len(fields) < 2:
        continue
      family_id = None
      sample_ids = []
      for field in fields:
        if family_id is None:
          family_id = field
        else:
          sample_ids.append(field)
      families[family_id] = sample_ids
  return families


def get_fastqs(fastq_dir):
  """Read the files in a directory and return a mapping of sample ID's to fastq filenames.
  The fastqs are found by matching FASTQ_NAME_REGEX and sample_id's are extracted that way.
  Returns a dict mapping each sample_id to a list of the filenames of fastq's that match that
  sample_id."""
  fastqs = {}
  # Read the files, build dict of fastqs.
  for filename in os.listdir(fastq_dir):
    match = re.search(FASTQ_NAME_REGEX, filename)
    if not match:
      return
    sample_id = match.group(1)
    fastq_pair = fastqs.get(sample_id, [])
    fastq_pair.append(filename)
    fastqs[sample_id] = fastq_pair
  # Check and sort the fastq filename lists.
  for sample_id, fastq_pair in fastqs.items():
    if len(fastq_pair) != 2:
      raise PipefamError('More or less than 2 files match the same sample ID "{}": {}'
                         .format(sample_id, ', '.join(fastq_pair)))
    fastq_pair.sort()
  return fastqs


class PipefamError(Exception):
  def __init__(self, message=None):
    if message:
      Exception.__init__(self, message)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
