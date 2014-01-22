#!/usr/bin/env python
# requires Python 2.7
from __future__ import division
import os
import sys
import random
import quicksect
import bioreaders
import collections
from optparse import OptionParser

OPT_DEFAULTS = {'lav':'', 'asm':'', 'ref':'', 'report':'', 'fragmented':500,
  'minor_len':600, 'min_flank':75, 'min_gap':10, 'contig_limit':10000,
  'slop':20, 'test_output':False}
USAGE = ("USAGE: %prog [opts] "
  +"(-a asm.fa -r ref.fa|-l align.lav) [-o asm-new.fa]")
DESCRIPTION = """Analyze an assembly via its LASTZ alignment to the reference,
and curate the assembly sequence. Provide either a LASTZ alignment (in LAV
format) or the assembly and reference sequence so that an alignment can be
performed. It will print to stdout an analysis of the alignment and any assembly
issues it indicates. If a new assembly filename is provided with -o, it will
write a curated assembly there, with some redundant contigs removed. If LASTZ
finds no hits, the curated assembly file will be created, but empty. On LAV
input: this is only written to work with output of LASTZ 1.02.00 run with
default options (only arguments are the reference and assembly filenames, then
"--format=lav").
Uses curation algorithm version 2: Just remove contigs which are entirely
contained within another contig (in terms of their footprints on the reference).
"""
EPILOG = """"""

DUPLICATE_THRES = 1.8
REPORT_TEXT = {
  'hits':['h','contigs with hits to the reference'],
  'curated':['c','contigs kept after curation'],
  'failure':['F','assembly failure: no contigs with similarity to the reference'],
  'fragmented':['f','fragmented assembly: a very high number of contigs'],
  'duplicated':['d','whole-chromosome duplications'],
  'nonref_flanks':['n','non-reference contig flanks'],
}

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
    help='Write a curated version of the assembly to this file.')
  parser.add_option('-R', '--report', dest='report',
    default=OPT_DEFAULTS.get('report'),
    help='Write an analysis of the assembly to this file instead of stdout.')
  parser.add_option('-f', '--fragmented', dest='fragmented', type='int',
    default=OPT_DEFAULTS.get('fragmented'),
    help="""Assemblies with more than this number of contigs will be marked as
"fragmented" in the report. Currently the contig number is approximated by the
number of LASTZ hits. Set to 0 to never mark as fragmented. Default: %default""")
  parser.add_option('-m', '--minor-len', dest='minor_len', type='int',
    default=OPT_DEFAULTS.get('minor_len'),
    help="""Minor contig threshold. Contigs longer than this will be considered
"major" contigs, and the rest are "minor" contigs. Default: %default""")
  parser.add_option('-F', '--min-flank', dest='min_flank', type='int',
    default=OPT_DEFAULTS.get('min_flank'),
    help="""Minimum length a non-reference flank must be to count it in the
report. Only major contigs are considered. Default: %default""")
  parser.add_option('-g', '--min-gap', dest='min_gap', type='int',
    default=OPT_DEFAULTS.get('min_gap'),
    help="""Minimum length for an assembly gap to be counted in the report.
Default: %default""")
  parser.add_option('-C', '--contig-limit', dest='contig_limit', type='int',
    default=OPT_DEFAULTS.get('contig_limit'),
    help="""Maximum allowed contigs in the assembly. If there are more than this
many contigs in the assembly, abort to avoid exceeding resources. Currently the
number of contigs is approximated by the number of LASTZ hits. Set to 0 for no
limit. Default: %default""")
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

  # set up I/O
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
  if options.report:
    report_handle = open(options.report, 'w')
  else:
    report_handle = sys.stdout
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

  report = get_report(lav, intervals, options)
  report_handle.write(format_report(report, REPORT_TEXT))

  if isinstance(report_handle, file):
    report_handle.close()

  if not (options.asm and options.out):
    sys.exit(0)
  curate(intervals, options.asm, options.out)



def length(interval):
  return interval[1] - interval[0] + 1


def align(asmpath, refpath):
  fail("LASTZ ALIGNMENT NOT YET IMPLEMENTED")


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
  for overlap in overlaps:
    # if interval is wholly contained within overlap (+ slop)
    if (overlap[0] - slop <= interval[0] and interval[1] <= overlap[1] + slop
        and length(interval) < length(overlap)):
      return True
  return False


def get_report(lav, intervals, options):
  """Generate report on the assembly quality.
  Call after redundants have been removed from intervals.
  Assumption: The entire query and subject sequences have been used in alignment
  (the end coordinate in the LAV h stanza == the sequence length)."""
  report = collections.OrderedDict()
  # number of "contigs"
  report['hits'] = len(lav)
  # number of kept contigs
  report['curated'] = len(get_contig_names(intervals))
  # failed assembly? (zero good contigs)
  if len(lav) == 0:
    report['failure'] = True
    return report
  else:
    report['failure'] = False
  # failed assembly? (many contigs)
  if options.fragmented and len(lav) > options.fragmented:
    report['fragmented'] = True
  else:
    report['fragmented'] = False
  # is there a contig that looks like a double of the reference?
  report['duplicated'] = 0
  for hit in lav:
    # length of contig over length of reference
    if hit.query['end'] / hit.subject['end'] > DUPLICATE_THRES:
      #TODO: check that there are actually two alignments ~ as big as reference
      report['duplicated'] += 1
  # are there contigs with long flanks of non-reference sequence?
  report['nonref_flanks'] = len(get_nonref_flanks(intervals, options))
  # how many gaps are in the assembly?
  #TODO: take chromosome into account
  merged = merge(intervals.keys(), sort=True)
  gaps = subtract(merged, (1, lav.hits[0].subject['end']))
  report['gaps'] = len(gaps)
  return report


def format_report(report, report_text):
  always_print = ['hits', 'curated']
  output = ''
  for issue in report:
    if report[issue] or issue in always_print:
      messages = report_text[issue]
      output += "\t".join([messages[0], str(report[issue]), messages[1]]) + "\n"
  return output


def merge(intervals, sort=False):
  """Merge a list of intervals into a new list of non-overlapping intervals.
  Optionally sort the output list by starting coordinate. Merging is 1-based
  (will merge if end == start). Totally naive implementation: O(n^2)."""
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
  if sort:
    merged.sort(key=lambda x: x[0])
  return merged


def subtract(merged, interval):
  """Subtract intervals in merged from interval. Merged must be a list of
  "merged" intervals: non-overlapping, in order from lowest to highest
  coordinate. Returns a list of intervals representing the regions of interval
  not overlapping any intervals in merged."""
  begin = interval[0]
  end = interval[1]
  results = []
  for i in range(len(merged)):
    # add region before first overlap (if any)
    if i == 0 and begin < merged[i][0]:
      results.append((begin, merged[i][0]-1))
    if i > 0:
      result = (merged[i-1][1]+1, merged[i][0]-1)
      if result[0] <= result[1]:
        results.append(result)
    # add region after last overlap (if any)
    if i == len(merged) - 1 and end > merged[i][1]:
      results.append((merged[i][1]+1, end))
  if not merged:
    results.append((begin, end))
  return results


def get_nonref_flanks(intervals, options):
  """Return a list of the names of contigs which have non-reference flanking
  sequences. There will be one entry for every flank, so a contig with two will
  appear twice (and len() of the list == number of flanks)."""
  # build a map of contigs to their regions that align to the reference
  asm_ref_regions = {}
  asm_lengths = {}
  for interval in intervals:
    hit = intervals[interval].parent
    queryname = hit.query['name']
    if queryname in asm_ref_regions:
      continue
    asm_lengths[queryname] = hit.query['end']
    ref_regions = []
    for alignment in hit:
      interval = (alignment.query['begin'], alignment.query['end'])
      ref_regions.append(interval)
    asm_ref_regions[queryname] = merge(ref_regions, sort=True)
  # check contigs for non-reference flanks
  nonref_flanks = []
  for queryname in asm_ref_regions:
    length = asm_lengths[queryname]
    ref_regions = asm_ref_regions[queryname]
    if length <= options.minor_len:
      continue
    # a contig with no alignment shouldn't appear, but deal with it anyway
    # by counting it as having two nonref flanks
    if not ref_regions:
      if length/2 >= options.min_flank:
        nonref_flanks.append(queryname)
        nonref_flanks.append(queryname)
      continue
    lflank = ref_regions[0][0] - 1
    rflank = length - ref_regions[-1][1]
    if lflank > options.min_flank:
      nonref_flanks.append(queryname)
    if rflank > options.min_flank:
      nonref_flanks.append(queryname)
  return nonref_flanks


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


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
