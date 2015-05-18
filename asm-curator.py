#!/usr/bin/env python
# requires Python 2.7
from __future__ import division
import os
import sys
import argparse
import subprocess
import collections
import lavreader
import fastareader
import lavintervals
import distutils.spawn

EXPECTED_VERSIONS = {
  'lavreader':'0.71', 'fastareader':'0.7', 'lavintervals':'0.81'
}

OPT_DEFAULTS = {'lav':'', 'asm':'', 'ref':'', 'report':'', 'fragmented':500,
  'minor_len':600, 'min_flank':75, 'min_gap':10, 'contig_limit':10000,
  'slop':20}
USAGE = """%(prog)s -l align.lav [-r report.tsv] [-a asm.fa -o new-asm.fa]
       %(prog)s -r ref.fa -a asm.fa [-l align.lav] [-R report.tsv] [-o new-asm.fa]"""
DESCRIPTION = """Analyze an assembly via its LASTZ alignment to the reference.
In the most basic mode, you can provide a LASTZ alignment (in LAV format) and it
will print a report on the quality of the assembly. If the assembly FASTA is
provided, it will produce a curated version, with redundant contigs removed.
If the reference genome is also provided, it will produce the LASTZ alignment
itself. The curation algorithm just removes contigs which are entirely contained
within another contig (in terms of their footprints on the reference)."""
EPILOG = """"""

REPORT_TEXT = {
  'hits':('h','contigs with hits to the reference'),
  'curated':('c','contigs kept after curation'),
  'failure':('F','assembly failure: no contigs with similarity to the reference'),
  'fragmented':('f','fragmented assembly: a very high number of contigs'),
  'duplicated':('d','whole-chromosome duplications'),
  'nonref_flanks':('n','non-reference contig flanks'),
  'gaps':('g','assembly gaps (between contigs)')
}
DUPLICATE_THRES = 1.8
# lastz ref.fa asm.fa --format=lav > out.lav
LASTZ_OPTS = ['--format=lav']


#TODO: Distinguish contigs with exactly the same start/end points.
#      Maybe use IntervalNode.linenum as unique identifier.
#TODO: When a contig is made of multiple alignments, check whether they look
#      incorrect (basically anything except a break at the reference edge).
def main():
  version_check(EXPECTED_VERSIONS)

  parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION,
    epilog=EPILOG)
  parser.set_defaults(**OPT_DEFAULTS)

  io_group = parser.add_argument_group('Input and output files')
  io_group.add_argument('-l', '--lav',
    help='The LAV input, if no reference is provided. If a reference (and '
      'assembly) is provided, the alignment will be performed and the LAV '
      'output will be written to this filename. It will not overwrite an '
      'existing file. If an output LAV filename is not given, it will create '
      'one using the assembly filename base. About the LASTZ alignment: it is '
      'done with the default options, and the "lastz" command must be on the '
      'PATH. It is only guaranteed to work with LASTZ version 1.02.00.')
  io_group.add_argument('-a', '--asm',
    help='The assembly, in a FASTA file. Only required to produce a curated '
      'assembly or to perform an alignment with the reference.')
  io_group.add_argument('-r', '--ref',
    help='The reference genome, in a FASTA file. Only required to perform an '
      'alignment with the reference.')
  io_group.add_argument('-o', '--out',
    help='Write a curated version of the assembly to this file. If no hits are '
      'found in the alignment, the file will be empty.')
  io_group.add_argument('-R', '--report',
    help='Write an analysis of the assembly to this file. If omitted, it will '
      'be printed to stdout.')
  params_group = parser.add_argument_group('Parameters')
  params_group.add_argument('-f', '--fragmented', type=int,
    help='Assemblies with more than this number of contigs will be marked as '
      '"fragmented" in the report. Currently the contig number is approximated '
      'by the number of LASTZ hits. Set to 0 to never mark as fragmented. '
      'Default: %(default)s')
  params_group.add_argument('-m', '--minor-len', type=int,
    help='Minor contig threshold. Contigs longer than this will be considered '
      '"major" contigs, and the rest are "minor" contigs. Default: %(default)s')
  params_group.add_argument('-F', '--min-flank', type=int,
    help='Minimum length a non-reference flank must be to count it in the '
      'report. Only major contigs are considered. Default: %(default)s')
  params_group.add_argument('-g', '--min-gap', type=int,
    help='Minimum length for an assembly gap to be counted in the report. '
      'Default: %(default)s')
  params_group.add_argument('-C', '--contig-limit', type=int,
    help='Maximum allowed contigs in the assembly. If there are more than this '
      'many contigs in the assembly, abort to avoid exceeding resources. '
      'Currently the number of contigs is approximated by the number of LASTZ '
      'hits. Set to 0 for no limit. Default: %(default)s')
  params_group.add_argument('-s', '--slop', type=int,
    help='Leeway in the contig overlap comparison. Essentially, when deciding '
      'whether contig A wholly contains contig B, contig A\'s flanks will be '
      'extended by this much. A larger slop means more contigs are removed. '
      'Note: in a mutual overlap situation, only the smaller contig is '
      'removed. Default: %(default)s')
  params_group.add_argument('-T', '--test-output', action='store_true',
    help='Print legacy test output.')

  args = parser.parse_args()

  # An LAV file is required unless we're producing the alignment ourselves,
  # in which case we need to create the output LAV filename.
  if not args.lav:
    if args.ref and args.asm:
      args.lav = os.path.splitext(args.asm)[0]+'.lav'
    else:
      parser.print_help()
      fail("\nError: Either an input LAV file or the assembly and reference is "
        "required.")
  # If the reference is provided, but not the assembly, assume a failed attempt
  # to produce an alignment.
  if args.ref and not args.asm:
    parser.print_help()
    fail("\nError: To run the LASTZ alignment, please provide both the "
      "reference and assembly sequences.")
  # If an assembly output name is given, but no assembly
  if args.out and not args.asm:
    parser.print_help()
    fail("\nError: To produce a curated assembly, an input assembly is required.")
  # Don't overwrite existing output assembly files
  if args.out and args.out != os.devnull and os.path.exists(args.out):
    fail('\nError: Curated assembly output file "'+args.out+'" already exists.')

  if args.report:
    report_handle = open(args.report, 'w')
  else:
    report_handle = sys.stdout

  # Do LASTZ alignment if both a reference and assembly were given
  if args.ref and args.asm:
    if os.path.exists(args.lav):
      fail('\nError: Attempted to perform alignment and write LAV output to "'+
        args.lav+'", but the file already exists.')
    if not distutils.spawn.find_executable('lastz'):
      fail('\nError: "lastz" command not found on PATH.')
    exit_code = align(args.ref, args.asm, args.lav)
    if exit_code != 0:
      fail('\nError: lastz alignment failure.')

  lav = lavreader.LavReader(args.lav)
  # Abort if assembly is too fragmented
  #TODO: count actual contigs, not hits, when assembly FASTA is provided
  if args.contig_limit and len(lav) > args.contig_limit:
    sys.stderr.write('Warning: Too many contigs in assembly ('+str(len(lav))
      +' LASTZ hits). Aborting.\n')
    sys.exit(0)

  intervals = lavintervals.alignments_to_intervals(lav)
  all_overlaps = lavintervals.get_all_overlaps(intervals)
  all_overlaps = lavintervals.discard_redundant(all_overlaps, slop=args.slop)
  if args.test_output:
    test_output(all_overlaps, intervals)
    sys.exit(0)
  # remove discarded intervals from main dict
  for interval in intervals.keys():
    if interval not in all_overlaps:
      del(intervals[interval])

  report = get_report(lav, intervals, args)
  report_handle.write(format_report(report, REPORT_TEXT))

  if report_handle is not sys.stdout:
    report_handle.close()

  if args.asm and args.out:
    curate(intervals, args.asm, args.out)



def align(refpath, asmpath, lavpath):
  """Perform LASTZ alignment, output LAV to lavpath.
  Returns the exit code from LASTZ.
  Can raise OSError, subprocess.CalledProcessError (at least).
  """
  lavhandle = open(lavpath, 'w')
  command = ['lastz', refpath, asmpath]
  command.extend(LASTZ_OPTS)
  exit_code = subprocess.call(command, stdout=lavhandle)
  lavhandle.close()
  return exit_code


def get_report(lav, intervals, args):
  """Generate report on the assembly quality.
  Call after redundants have been removed from intervals.
  Assumption: The entire query and subject sequences have been used in alignment
  (the end coordinate in the LAV h stanza == the sequence length)."""
  report = collections.OrderedDict()
  # number of contigs (with hits to reference)
  contigs = set()
  for hit in lav:
    contigs.add(hit.query['name'])
  report['hits'] = len(contigs)
  # number of kept contigs
  report['curated'] = len(get_contig_names(intervals))
  # failed assembly? (zero good contigs)
  if report['hits'] == 0:
    report['failure'] = True
    return report
  else:
    report['failure'] = False
  # failed assembly? (many contigs)
  if args.fragmented and report['hits'] > args.fragmented:
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
  report['nonref_flanks'] = len(get_nonref_flanks(intervals, args))
  # how many gaps are in the assembly?
  #TODO: take chromosome into account
  merged = lavintervals.merge(intervals.keys(), sort=True)
  gaps = lavintervals.subtract(merged, (1, lav.hits[0].subject['end']))
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


def get_nonref_flanks(intervals, args):
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
    asm_ref_regions[queryname] = lavintervals.merge(ref_regions, sort=True)
  # check contigs for non-reference flanks
  nonref_flanks = []
  for queryname in asm_ref_regions:
    length = asm_lengths[queryname]
    ref_regions = asm_ref_regions[queryname]
    if length <= args.minor_len:
      continue
    # a contig with no alignment shouldn't appear, but deal with it anyway
    # by counting it as having two nonref flanks
    if not ref_regions:
      if length/2 >= args.min_flank:
        nonref_flanks.append(queryname)
        nonref_flanks.append(queryname)
      continue
    lflank = ref_regions[0][0] - 1
    rflank = length - ref_regions[-1][1]
    if lflank > args.min_flank:
      nonref_flanks.append(queryname)
    if rflank > args.min_flank:
      nonref_flanks.append(queryname)
  return nonref_flanks


def curate(intervals, asmpath, outpath):
  """Remove all contigs not present in intervals and write FASTA to outpath."""
  contigs = get_contig_names(intervals)
  fasta = fastareader.FastaLineGenerator(asmpath)
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


def version_check(expected):
  actual = {}
  for module_name in expected:
    module = sys.modules[module_name]
    for version_name in ['version', 'VERSION', '__version__']:
      if version_name in dir(module):
        actual[module_name] = getattr(module, version_name)
  for module_name in actual:
    assert actual[module_name] == expected[module_name], (
      "Wrong version of "+module_name+". Expected: "+expected[module_name]
      +", actual: "+actual[module_name]
    )


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)


if __name__ == "__main__":
  main()
