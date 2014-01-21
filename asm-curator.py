#!/usr/bin/env python
import os
import sys
import random
import quicksect
import bioreaders
from optparse import OptionParser

OPT_DEFAULTS = {'lav':'', 'asm':'', 'ref':'', 'contig_limit':1000, 'slop':20,
  'test_output':False}
USAGE = ("USAGE: %prog [opts] "
  +"(-a asm.fa -r ref.fa|-l align.lav) [-o asm-new.fa]")
DESCRIPTION = """Analyze an assembly via its LASTZ alignment to the reference,
and curate the assembly sequence. Provide either a LASTZ alignment (in LAV
format) or the assembly and reference sequence so that an alignment can be
performed. It will print to stdout an analysis of the alignment and any assembly
issues it indicates. If a new assembly filename is provided with -o, it will
create a curated assembly, with some redundant contigs removed.
Uses curation algorithm version 2: Just remove contigs which are entirely
contained within another contig (in terms of their footprints on the reference).
"""
EPILOG = """"""

#TODO: Distinguish contigs with exactly the same start/end points.
#      Maybe use IntervalNode.linenum as unique identifier.
#TODO: When a contig is made of multiple alignments, check whether they look
#      incorrect (basically anything except a break at the reference edge).
def main():

  parser = OptionParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)

  parser.add_option('-l', '--lav', dest='lav',
    default=OPT_DEFAULTS.get('lav'),
    help='The LASTZ alignment of the assembly to the reference, in LAV format.')
  parser.add_option('-a', '--asm', dest='asm',
    default=OPT_DEFAULTS.get('asm'),
    help='The assembly itself, in a FASTA file.')
  parser.add_option('-r', '--ref', dest='ref',
    default=OPT_DEFAULTS.get('ref'),
    help='The reference genome, in a FASTA file.')
  parser.add_option('-o', '--out', dest='out',
    default=OPT_DEFAULTS.get('out'),
    help='Write a curated version of the assembly to this filename.')
  parser.add_option('-C', '--contig-limit', dest='contig_limit', type='int',
    default=OPT_DEFAULTS.get('contig_limit'),
    help="""Maximum allowed contigs in the assembly. If there are more than this
many contigs in the assembly, abort to avoid exceeding resources. Currently the
number of contigs is approximated by the number of LASTZ hits. Set to 0 for no
limit (at your own risk). Default: %default""")
  parser.add_option('-s', '--slop', dest='slop', type='int',
    default=OPT_DEFAULTS.get('slop'),
    help="""Leeway in the contig overlap comparison. Essentially, when deciding
whether contig A wholly contains contig B, contig A's flanks will be extended
by this much. A larger slop means more contigs are removed. Note: in a mutual
overlap situation, only the smaller contig is removed. Default: %default""")
  parser.add_option('-T', '--test-output', dest='test_output',
    action='store_true', default=OPT_DEFAULTS.get('test_output'),
    help='Print legacy test output.')

  (options, arguments) = parser.parse_args()

  files = []
  if options.lav:
    files.append(options.lav)
    if options.asm:
      files.append(options.asm)
  else:
    if options.asm and options.ref:
      files.append(options.asm)
      files.append(options.ref)
    else:
      parser.print_help()
      fail("\nError: Please provide either an LAV file or an assembly & "
        +"reference FASTA file.")
  if options.out:
    if not options.asm:
      parser.print_help()
      fail("\nError: Please provide the assembly FASTA file in order to create "
        +"a curated assembly.")

  for path in files:
    if not os.path.isfile(path):
      fail("\nError: cannot find input file "+path)

  if options.lav:
    lavpath = options.lav
  else:
    lavpath = align(options.asm, options.ref)

  lav = bioreaders.LavReader(lavpath)
  # Abort if assembly is too fragmented
  #TODO: count actual contigs, not hits, when assembly FASTA is provided
  if options.contig_limit and len(lav) > options.contig_limit:
    sys.stderr.write('Warning: Too many contigs in assembly ('+str(len(lav))
      +' LASTZ hits). Aborting.\n')
    sys.exit(0)
  lav.convert()

  intervals = lav_to_intervals(lav)
  all_overlaps = get_all_overlaps(intervals)
  all_overlaps = discard_redundant(all_overlaps, slop=options.slop)
  if options.test_output:
    test_output(all_overlaps, intervals)
    sys.exit(0)
  # remove discarded intervals from main dict
  for interval in intervals.keys():
    if interval not in all_overlaps:
      del(intervals[interval])

  #TODO: print alignment analysis

  if not (options.asm and options.out):
    sys.exit(0)

  curate(intervals, options.asm, options.out)



def length(interval):
  return interval[1] - interval[0] + 1


def test_output(all_overlaps, intervals):
  for interval in sorted(intervals, key=lambda x: x[0]):
    if interval not in all_overlaps:
      continue
    print "overlapping",format_interval(interval, intervals)
    # for overlap in sorted(all_overlaps[interval], key=lambda x: x[0]):
    #   print format_interval(overlap, intervals)
    # print "  unique:"
    # merged = merge(all_overlaps[interval])
    # merged.sort(key=lambda x: x[0])
    # for unique in get_uniques(merged, interval):
    #   print "    "+format_interval(unique)


def align(asmpath, refpath):
  fail("LASTZ ALIGNMENT NOT YET IMPLEMENTED")


def format_interval(interval, intervals=None):
  output = str(interval)
  output += ' '+str(interval[1]-interval[0]+1)+' bp'
  if intervals:
    name = intervals[interval].parent.query['id']
    if name.startswith('NODE_'):
      fields = name.split('_')
      if len(fields) > 2:
        name = '_'.join(fields[:2])
    output += ":\t"+name
  return output


def lav_to_intervals(lav):
  """Convert a LASTZ alignment to a series of intervals along the reference
  sequence.
  Give an LavReader and it will return a dict with the intervals as keys and
  the corresponding LavAlignments as values. Each interval is a tuple of the
  "begin" and "end" values of an LavAlignment."""
  intervals = {}
  for hit in lav:
    for alignment in hit:
      interval = (alignment.subject['begin'], alignment.subject['end'])
      intervals[interval] = alignment
  return intervals


def get_all_overlaps(intervals):
  """Compare each interval to the rest, finding ones with any overlap.
  Returns a dict mapping each interval to a list of intervals which overlap.
  The original interval is always excluded from the list."""
  # build tree
  tree = None
  random.seed(1)
  for interval in intervals:
    if tree is None:
      tree = quicksect.IntervalNode(interval[0], interval[1])
    else:
      tree = tree.insert(interval[0], interval[1])
  # find ones that intersect each interval
  all_overlaps = {}
  for interval in intervals:
    overlaps = []
    # reporter function: is given the matching interval when one is found
    add_overlap = lambda node: overlaps.append((node.start, node.end))
    tree.intersect(interval[0], interval[1], add_overlap)
    # remove the query interval from the results
    overlaps = [overlap for overlap in overlaps if overlap != interval]
    all_overlaps[interval] = overlaps
  return all_overlaps


def discard_redundant(all_overlaps, slop=0):
  """Remove redundant intervals from all_overlaps.
  The discarded intervals will be removed from both the
  set of keys and the lists of overlapping intervals."""
  redundant = set()
  # gather redundant intervals
  for interval in all_overlaps:
    if is_redundant(interval, all_overlaps[interval], slop=slop):
      redundant.add(interval)
  # go through all_overlaps again, discarding redundant intervals
  for interval in all_overlaps.keys():
    if interval in redundant:
      del(all_overlaps[interval])
    else:
      all_overlaps[interval] = [overlap for overlap in all_overlaps[interval]
        if overlap not in redundant]
  return all_overlaps


def is_redundant(interval, overlaps, slop=0):
  """Return True if the interval is redundant."""
  # print "is_redundant on "+format_interval(interval)
  for overlap in overlaps:
    # if interval is wholly contained within overlap (+ slop)
    if (overlap[0] - slop <= interval[0] and interval[1] <= overlap[1] + slop
        and length(interval) < length(overlap)):
      # print "  wholly contained within "+format_interval(overlap)
      return True
  return False


def merge(intervals):
  """Merge a list of intervals into a new list of non-overlapping intervals.
  Merging is 1-based (will merge if end == start).
  Implementation is totally naive: O(n^2)."""
  #TODO: replace with an actually efficient merge, e.g. from bx-python
  merged = []
  # go through input list, adding 1 interval at a time to a growing merged list
  for interval in intervals:
    i = 0
    while i < len(merged):
      target = merged[i]
      # if any overlap
      if interval[0] <= target[1] and interval[1] >= target[0]:
        # replace the target with a combination of it and the input interval
        del(merged[i])
        # merge the two intervals by taking the widest possible region
        interval = (min(interval[0], target[0]), max(interval[1], target[1]))
      else:
        i+=1
    merged.append(interval)
  return merged


def get_uniques(merged, int_of_interest):
  """Find all regions in overlaps that do not overlap with any other interval.
  merged is a merge of all the intervals overlapping the one of interest. It
  must be a set of non-overlapping intervals, in order from lowest to highest
  coordinate.
  interval is the interval of interest.
  Returns a list of intervals representing the unique spans."""
  begin = int_of_interest[0]
  end = int_of_interest[1]
  uniques = []
  for i in range(len(merged)):
    # add region before first overlap (if unique)
    if i == 0 and begin < merged[i][0]:
      uniques.append((begin, merged[i][0]-1))
    if i > 0:
      unique = (merged[i-1][1]+1, merged[i][0]-1)
      if unique[0] <= unique[1]:
        uniques.append(unique)
    # add region after last overlap (if unique)
    if i == len(merged) - 1 and end > merged[i][1]:
      uniques.append((merged[i][1]+1, end))
  if not merged:
    uniques.append((begin, end))
  return uniques


def curate(intervals, asmpath, outpath):
  """Remove all contigs not present in intervals and write FASTA to outpath."""
  contigs = get_contig_names(intervals)
  fasta = bioreaders.FastaLineGenerator(asmpath)
  last_name = None
  skipping = False
  with open(outpath, 'w') as outfile:
    for line in fasta:
      if fasta.name != last_name:
        if fasta.name in contigs:
          skipping = False
          outfile.write('>'+fasta.name+'\n')
        else:
          skipping = True
        last_name = fasta.name
      if not skipping:
        outfile.write(line+'\n')


def get_contig_names(intervals):
  """Get the FASTA sequence names for each contig represented in the intervals
  dict. Returns the names in a set."""
  names = set()
  for alignment in intervals.values():
    hit = alignment.parent
    names.add(hit.query['name'])
  return names


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
