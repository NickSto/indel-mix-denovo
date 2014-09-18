#!/usr/bin/env python
from __future__ import division
import re
import os
import sys
import random
import shutil
import string
import logging
import argparse
import subprocess
import lavintervals
import fastareader
import lavreader

OPT_DEFAULTS = {'output_name':'Cleaned', 'score_by':'length', 'fasta_width':70}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""
EPILOG = """"""
TMP_DIR_BASE = 'cleaning'
SPADES_NAME_PATTERN = r'^(NODE_\d+)_length_\d+_cov_\d+\.?\d+_ID_\d+$'


def main():

  parser = argparse.ArgumentParser(description=DESCRIPTION, epilog=EPILOG)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('ref', metavar='reference.fa',
    help='The FASTA file of the reference. N.B.: Must contain only one '
      'sequence. If the full reference contains multiple chromosomes, you will '
      'have to break it into multiple files and run this once per chromsome.')
  parser.add_argument('asm', metavar='assembly.fa',
    help='The FASTA file of the raw assembly.')
  parser.add_argument('-s', '--score-by', choices=['length', 'id', 'support'],
    help='The method used to determine which of two overlapping contigs to use '
      'in constructing the final sequence. Default: %(default)s.')
  parser.add_argument('-O', '--orient', action='store_true',
    help='Only orient the contigs, but do not alter them otherwise.')
  parser.add_argument('-o', '--output', metavar='assembly-unified.fa',
    help='Write the processed assembly FASTA to this file instead of stdout.')
  parser.add_argument('-n', '--output-name', metavar='SeqName',
    help='Use this name for the output sequence in the FASTA header line. '
      'Default: "%(default)s".')
  parser.add_argument('-l', '--log', metavar='logfile.txt',
    help='A log file to use for writing details of the process, if one is '
      'desired.')
  parser.add_argument('-v', '--verbosity', type=int, metavar='level',
    help='How verbose the log file printing should be. If -v is given but -l '
      'is not, log printing will be turned on, to stderr. Give a number from 0 '
      '(silent) to 3 (most verbose). Default: 2 when printing to a file, '
      '1 when printing to stderr (-l -).')
  parser.add_argument('-W', '--fasta-width', metavar='characters', type=int,
    help='Line width of the output FASTA file. Default: %(default)s.')

  args = parser.parse_args()

  if args.score_by not in ('length'):
    raise NotImplementedError('--score-by option not implemented yet.')

  if args.output:
    outfile = open(args.output, 'w')
  else:
    outfile = sys.stdout

  # Set up logger
  # logging level
  if args.verbosity == 0:
    loglevel = logging.WARNING
  elif args.verbosity == 1:
    loglevel = logging.INFO
  elif args.verbosity == 2:
    loglevel = logging.DEBUG
  elif args.verbosity >= 3:
    loglevel = logging.NOTSET
  elif args.verbosity is None:
    # default loglevel when printing to screen
    if args.log == '-':
      loglevel = logging.INFO
    # default loglevel when printing to file
    else:
      loglevel = logging.DEBUG
  else:
    fail('Error: Invalid verbosity "'+str(args.verbosity)+'"')
  # open logger, set to correct output destination
  if args.log == '-' or args.verbosity is not None:
    logging.basicConfig(stream=sys.stderr, level=loglevel, format='%(message)s')
  elif args.log:
    logging.basicConfig(filename=args.log, filemode='w', level=loglevel,
      format='%(message)s')
  else:
    logging.disable(logging.CRITICAL)

  # open LAV and FASTA files, make temporary directory
  tmpdir = get_tmp_path(TMP_DIR_BASE)
  try:
    os.makedirs(tmpdir)
  except OSError:
    fail('Error: temporary directory "'+tmpdir+'" exists.')

  # on any exception, first close the outfile and remove the temp directory
  def cleanup_excepthook(exceptype, value, traceback):
    cleanup(outfile, tmpdir)
    sys.__excepthook__(exceptype, value, traceback)
  sys.excepthook = cleanup_excepthook

  # Orient all contigs in forward direction
  pre_lav = align(args.ref, args.asm, tmpdir) # LASTZ align asm to ref
  asm_fasta_path = orient(args.asm, pre_lav, tmpdir, args.fasta_width)
  if args.orient:
    with open(asm_fasta_path) as asm_fasta:
      for line in asm_fasta:
        sys.stdout.write(line)
    cleanup(outfile, tmpdir)
    sys.exit(0)
  asm_fasta = fastareader.FastaLineGenerator(asm_fasta_path)

  # Now align the oriented assembly to the reference
  lav = align(args.ref, asm_fasta_path, tmpdir)

  # convert alignments to a set of intervals and map each to its alignment
  #TODO 3: System to avoid collisions of intervals with identical start/ends
  #        Maybe include an identifier of the origin alignment in the tuple?
  #        (A simple third field of "query name" helps but won't be sufficient.)
  interval_to_aln = lavintervals.alignments_to_intervals(lav)
  # extract the intervals into a simple list
  intervals = interval_to_aln.keys()

  # Construct a final, non-redundant sequence out of the original, by walking
  # along the reference, adding sequence from the assembly.
  # Algorithm:
  # While intervals left:
  #   Sort list of intervals by starting coordinate
  #   Pop the first interval off the list
  #   If the interval is unique (does not overlap any other interval):
  #     Take the corresponding sequence from the assembly & add it to the output
  #   Else (it overlaps):
  #     Break the interval in two: a unique region, followed by the rest
  #     Add the unique region (if any) to the list of intervals
  #     Decide between the two overlapping intervals:
  #       For each, align the overlapping region to the original assembly
  #       Use the one with the highest LASTZ score
  #     If the larger interval wins:
  #       Pop the smaller interval from the list
  #     Else (smaller interval wins):
  #       Split the larger interval and pop it from the list
  #       Add the non-overlapping segment to the list
  final_sequence = ''
  while len(intervals) > 0:
    # Pop the next interval, peek the one after it (if it exists)
    #TODO 3: add new intervals intelligently, to preserve order & avoid sorting
    intervals.sort(key=lambda interval: interval[0]) # sort by start coord
    interval = intervals[0]
    if len(intervals) > 1:
      next_interval = intervals[1]
    else:
      next_interval = None
    del(intervals[0])
    assert next_interval is None or interval[0] <= next_interval[0]

    # Is this interval unique (no overlap between this interval and the next)?
    if next_interval is None or interval[1] < next_interval[0]:
      # Add query sequence of the interval to the output sequence
      seq_name = interval_to_aln[interval].parent.query['name']
      interval_on_query = convert_with_alignment(interval,
                                                 interval_to_aln[interval],
                                                 fail='tryharder')
      logging.debug('Using {} for {:>5} to {:>5}'.format(nickname(seq_name),
                                                         interval[0],
                                                         interval[1]))
      final_sequence += asm_fasta.extract(*interval_on_query, chrom=seq_name)
      # Is there a gap between this interval and the next?
      if next_interval is not None and interval[1] - next_interval[0] > 1:
        #TODO 2: add N's (goal is to keep reference coordinates)
        #        not necessary for R33S10 (those gaps are in redundant seq)
        raise NotImplementedError('There must be no gaps between contigs.')

    # ..or do the intervals overlap?
    else:
      # The overlapping region
      overlap = (next_interval[0], min(interval[1], next_interval[1]))
      # Portion of interval before the overlapping region
      unique = (interval[0], next_interval[0]-1)
      # All parts of interval after the unique region
      rest = (next_interval[0], interval[1])
      # Add the unique sequence before the overlap (if any) to the list
      if length(unique) > 0:
        intervals.append(unique)
        interval_to_aln[unique] = interval_to_aln[interval]
      # Determine which sequence to use in the overlap section
      alignment1 = interval_to_aln[interval]
      alignment2 = interval_to_aln[next_interval]
      if choose_sequence(alignment1, alignment2, overlap, tmpdir, lav, args):
        (winner, loser) = (rest, next_interval)
      else:
        (winner, loser) = (next_interval, rest)
      assert winner[0] == loser[0]
      # Make sure the list contains the winning interval and not the losing one.
      if loser in intervals:
        intervals.remove(loser)
      if winner not in intervals:
        intervals.append(winner)
        if winner == next_interval:
          interval_to_aln[winner] = interval_to_aln[next_interval]
        else:
          interval_to_aln[winner] = interval_to_aln[interval]
      # If the losing interval is the larger one, split it and keep the part
      # after the overlap.
      if length(loser) > length(winner):
        loser_rest = (winner[1]+1, loser[1])
        intervals.append(loser_rest)
        if loser == next_interval:
          interval_to_aln[loser_rest] = interval_to_aln[next_interval]
        else:
          interval_to_aln[loser_rest] = interval_to_aln[interval]

  outfile.write(fasta_format(final_sequence, args.output_name, args.fasta_width))
  cleanup(outfile, tmpdir)


def length(interval):
  """Get length of interval.
  1-based: length((10, 10)) == 1"""
  return interval[1] - interval[0] + 1


def align(ref_path, asm_path, tmpdir):
  """LASTZ align two FASTA files and return an LavReader of the result.
  Executes "$ lastz ref_path asm_path" and writes the output to an LAV file in
  "tmpdir"."""
  basename = os.path.splitext(os.path.split(asm_path)[1])[0]
  lav_path = os.path.join(tmpdir, basename+'.lav')
  with open(lav_path, 'w') as lavfile:
    logging.info("$ "+" ".join(['lastz', ref_path, asm_path]))
    subprocess.call(['lastz', ref_path, asm_path], stdout=lavfile)
    #TODO 3: Check exit code for success or failure
  return lavreader.LavReader(lav_path)


def orient(in_fasta_path, lav, tmpdir, fasta_width):
  basename = os.path.splitext(os.path.split(in_fasta_path)[1])[0]
  out_fasta_path = os.path.join(tmpdir, basename+'.oriented.fa')

  # Read the LAV alignment to determine the orientation of each sequence.
  # If there are hits to both orientations, decide based on the sum of the
  # alignment scores for each.
  orientations = {}
  scores = {}
  for hit in lav:
    name = hit.query['name']
    #TODO 2: Check that summing alignment scores is a reasonable way of
    #        scoring a hit.
    score = 0
    for alignment in hit:
      score += alignment.score
    if score > scores.get(name, -1):
      scores[name] = score
      orientations[name] = hit.query['revcomp']

  # Read through the input FASTA, printing to each sequence in the correct
  # orientation to the output FASTA.
  in_fasta = fastareader.FastaLineGenerator(in_fasta_path)
  name = None
  revcomp = None
  seqbuffer = ''
  with open(out_fasta_path, 'w') as out_fasta:
    for line in in_fasta:
      # Started a new sequence; finish up the last one and print header.
      if in_fasta.name != name:
        if revcomp:
          revcomp_seq = get_revcomp(seqbuffer)
          out_fasta.write(fasta_format(revcomp_seq, name, width=fasta_width,
                                       header=False))
        out_fasta.write(">"+in_fasta.name+"\n")
        name = in_fasta.name
        revcomp = orientations.get(name, False)
        seqbuffer = ''
      # If it's reversed, save up the whole sequence so it can be revcomp'd as
      # a whole at the end.
      if revcomp:
        seqbuffer += line
      else:
        out_fasta.write(line+"\n")
    if revcomp:
      revcomp_seq = get_revcomp(seqbuffer)
      out_fasta.write(fasta_format(revcomp_seq, name, width=fasta_width,
                                   header=False))
  return out_fasta_path


def choose_sequence(alignment1, alignment2, overlap, tmpdir, lav, args):
  """Returns True if the first sequence is best, False otherwise.
  Decides based on the criteria specified in args.score_by:
    "length":  Chooses the longer contig.
    "id":      Chooses the contig with the highest % identity compared to the
               reference.
    "support": Chooses the contig with the highest coverage in supporting reads.
    """
  if args.score_by == 'length':
    return choose_sequence_length(alignment1, alignment2)
  elif args.score_by == 'id':
    return choose_sequence_id(alignment1, alignment2, lav)
  elif args.score_by == 'support':
    return choose_sequence_support(alignment1, alignment2, tmpdir, args)


def choose_sequence_length(alignment1, alignment2):
  id1 = nickname(alignment1.parent.query['id'])
  id2 = nickname(alignment2.parent.query['id'])
  logging.debug("choosing between {} and {}:".format(id1, id2))
  length1 = alignment1.parent.query['length']
  length2 = alignment2.parent.query['length']
  if length1 >= length2:
    logging.debug('  winner: {} ({} >= {})'.format(id1, length1, length2))
    return True
  else:
    logging.debug('  winner: {} ({} < {})'.format(id2, length1, length2))
    return False

def choose_sequence_id(alignment1, alignment2, lav):
  raise NotImplementedError


def choose_sequence_support(alignment1, alignment2, tmpdir, args):
  raise NotImplementedError


def choose_sequence_old(alignment1, alignment2, overlap, tmpdir, asm_raw_file):
  best_hits = []
  logging.debug("processing {}:".format(overlap))
  for (alignment, name) in zip((alignment1, alignment2), ('seq1', 'seq2')):
    fastapath = os.path.join(tmpdir, name+'.fa')
    overlap_on_query = convert_with_alignment(overlap, alignment,
                                              fail='tryharder')
    sequence = asm_merge.extract(*overlap_on_query)
    with open(fastapath, 'w') as fastafile:
      fastafile.write(fasta_format(sequence, name))
    # Perform LASTZ alignment
    lavpath = os.path.join(tmpdir, name+'.lav')
    with open(lavpath, 'w') as lavfile:
      # sys.stderr.write("$ "+" ".join(['lastz', fastapath, asm_raw_file])+"\n")
      subprocess.call(['lastz', fastapath, asm_raw_file], stdout=lavfile)
    lav = lavreader.LavReader(lavpath)
    # Get the top-scoring alignment
    logging.debug("  scores for {}:".format(overlap_on_query))
    (top_length, top_score) = lav_top_length(lav)
    best_hits.append({'length':top_length, 'score':top_score})
  # If both match the same contig the best, choose the highest alignment score
  # (if length is equal, use alignment score as a tiebreaker).
  if best_hits[0]['length'] == best_hits[1]['length']:
    return best_hits[0]['score'] > best_hits[1]['score']
  else:
    return best_hits[0]['length'] > best_hits[1]['length']


def lav_top_score(lav):
  """Score an LAV by the top score of all alignments."""
  top_score = 0
  for hit in lav:
    for alignment in hit:
      logging.debug("    {}".format(alignment.score))
      if alignment.score > top_score:
        top_score = alignment.score
  return top_score


def lav_top_id(lav):
  """Score an LAV by the top % identity of all alignments.
  % identity of an alignment is computed by an average of all the % identities
  of its blocks, weighted by block length."""
  top_id = 0
  for hit in lav:
    for alignment in hit:
      total_id = 0
      total_length = 0
      for block in alignment:
        length = block.query['length']
        total_length += length
        total_id += length * block.identity
      id_pct = total_id / total_length
      logging.debug("    {}".format(round(id_pct, 2)))
      if id_pct > top_id:
        top_id = id_pct
  return top_id


def lav_top_length(lav):
  """Score an LAV by the length of query sequence of the best hit.
  The best hit is determined by its top alignment score."""
  length_of_best = 0
  top_score = 0
  for hit in lav:
    log_msg = "    len: {}\tscores: ".format(hit.query['length'])
    for alignment in hit:
      log_msg += "{} ".format(alignment.score)
      if alignment.score > top_score:
        top_score = alignment.score
        length_of_best = hit.query['length']
    logging.debug(log_msg)
  logging.debug("    best: {}".format(length_of_best))
  return (length_of_best, top_score)


def fasta_format(sequence, name, width=70, header=True):
  """Turn a sequence and name into a FASTA-formatted string."""
  if header:
    output = '>'+name+'\n'
  else:
    output = ''
  for start in range(0, len(sequence), width):
    end = start + width
    if end > len(sequence):
      end = len(sequence)
    output += sequence[start:end]+'\n'
  return output


def get_revcomp(sequence):
  delete_chars = '\r\n '
  table = string.maketrans('acgtrymkbdhvACGTRYMKBDHV',
                           'tgcayrkmvhdbTGCAYRKMVHDB')
  return sequence.translate(table, delete_chars)[::-1]


def convert_with_alignment(interval, alignment, fail='throw'):
  """One-off interval conversion, using a pre-chosen alignment.
  Converts the interval's start/end coordinates from subject to query
  coordinates using the given alignment.
  Assumes the interval is actually contained in the alignment."""
  #TODO 1: What happens when a coordinate is in a small gap between blocks?
  #        Make sure the resulting sequence doesn't include 1bp duplications.
  table = lavintervals.alignments_to_conv_table([alignment],
                                                query_to_subject=False)
  try:
    begin = lavintervals.convert(table, interval[0], fail=fail)[0]
    end = lavintervals.convert(table, interval[1], fail=fail)[0]
  except Exception as e:
    if len(e.args) > 1 and e.args[1] == 'fail':
      raise AssertionError('Interval must be contained in alignment.')
    else:
      raise
  if begin > end:
    (begin, end) = (end, begin)
  return (begin, end)


def get_tmp_path(base, max_tries=20):
  """Return an unoccupied path based on the supplied one.
  The returned path will be the argument plus ".tmp", or if that's taken, with a
  number up to max_tries, like ".3.tmp". Once max_tries has been reached, it
  will throw an exception.
  N.B.: It will be a relative path if the input is."""
  attempts = 0
  candidate = base + '.tmp'
  while os.path.exists(candidate):
    candidate = "{}.{}.{}".format(base, attempts, 'tmp')
    attempts+=1
    if attempts > 20:
      raise Exception('Cannot find an unoccupied temp directory name.')
  return candidate


# Turn some types of verbose contig names into more concise ones
def nickname(raw_name):
  # SPAdes
  new_name = raw_name
  match = re.search(SPADES_NAME_PATTERN, raw_name)
  if match:
    new_name = match.group(1)
  return new_name


def cleanup(outfile, tmpdir):
  logging.shutdown()
  if outfile is not sys.stdout:
    outfile.close()
  shutil.rmtree(tmpdir)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  main()
