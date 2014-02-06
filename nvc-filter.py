#!/usr/bin/env python
from __future__ import division
import os
import re
import sys
import errno
import vcfreader
from collections import deque
from optparse import OptionParser

# Command line options

OPT_DEFAULTS = {'output':'-', 'freq_thres':0.0, 'covg_thres':0, 'minor':False,
  'strand':False, 'keep_types':'', 'rm_types':'', 'compliant':False,
  'split_file':False, 'drop_end_len':0, 'ends':''}
USAGE = """USAGE: %prog [options] nvc-results.vcf -o filtered-variants.vcf
       cat nvc-results.vcf | %prog [options] - > filtered-variants.vcf"""
DESCRIPTION = """Filters the output of the Naive Variant Caller, removing
variant calls according to user-set criteria. It will use the raw variant
counts in the sample columns to determine which variants satisfy the criteria.
Ones that do not will be removed from the ALT column. If no variant at a
position passes the filters, the position will be omitted in the output."""
EPILOG = """"""

def main():

  parser = OptionParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)

  parser.add_option('-o', '--output', dest='output',
    default=OPT_DEFAULTS.get('output'), help="""Write output to the given file
instead of stdout.""")
  parser.add_option('-f', '--freq-thres', dest='freq_thres', type='float',
    default=OPT_DEFAULTS.get('freq_thres'),
    help="""Frequency threshold. Variants less abundant than this will be
filtered out. Given in %. Default: %default%""")
  parser.add_option('-c', '--covg-thres', dest='covg_thres', type='int',
    default=OPT_DEFAULTS.get('covg_thres'),
    help="""Coverage threshold. Sites with fewer total reads than this will be
filtered out. Default: %default%""")
  parser.add_option('-k', '--keep-types', dest='keep_types',
    default=OPT_DEFAULTS.get('keep_types'), help="""Preserve these types of
variants always, regardless of filters. Give a string of letters, e.g. "SID" to
keep SNVs (S), insertions (I), and deletions (D). Default: "%default" """)
  parser.add_option('-r', '--rm-types', dest='rm_types',
    default=OPT_DEFAULTS.get('rm_types'), help="""Remove these types of
variants always, regardless of filters. Use same notation as for keep-types.
Default: "%default" """)
  parser.add_option('-m', '--minor', dest='minor', action='store_const',
    const=not OPT_DEFAULTS.get('minor'), default=OPT_DEFAULTS.get('minor'),
    help="""***NOT YET IMPLEMENTED*** Use this to find minority variants only.
When selected, variants that are the majority will be filtered (whether the
same as the reference or not""")
  parser.add_option('-s', '--strand', dest='strand', action='store_const',
    const=not OPT_DEFAULTS.get('strand'), default=OPT_DEFAULTS.get('strand'),
    help="""***NOT YET IMPLEMENTED*** Attempt to filter out indels which show
strand bias.""")
  parser.add_option('-S', '--split-file', dest='split_file',
    action='store_const', const=not OPT_DEFAULTS.get('split_file'),
    default=OPT_DEFAULTS.get('split_file'),
    help="""Create an output VCF file for every sample in the input. Currently,
this leaves some fields in an inconsistent state. The ones that will be correct
are CHROM, POS, REF, ALT, and the variant counts. The names for the new files
will be based on the either the output filename, the input filename, or the
name of the first sample, preferred in that order.""")
  parser.add_option('-C', '--compliant', dest='compliant',
    action='store_const', const=not OPT_DEFAULTS.get('compliant'),
    default=OPT_DEFAULTS.get('compliant'),
    help="""***NOT YET IMPLEMENTED*** Produce standard-compliant VCF output.
Currently this means the presence of the REF allele is indicated with an INFO
flag instead of putting it in the ALT column.""")
  parser.add_option('-e', '--ends', dest='ends',
    default=OPT_DEFAULTS.get('ends'),
    help="""Specify the start and end coordinates of the reference. Give a
comma-separated pair of numbers like "1,23000". Leave one coordinate
unspecified by omitting it, but leaving the comma, like "100," """)
  parser.add_option('-E', '--drop-end-len', dest='drop_end_len', type='int',
    default=OPT_DEFAULTS.get('drop_end_len'),
    help="""Exclude (filter out) sites this many bases from either end. If no
--ends are given, """)

  (options, arguments) = parser.parse_args()

  ends = options.ends
  output = options.output
  strand = options.strand
  rm_types = options.rm_types
  keep_types = options.keep_types
  covg_thres = options.covg_thres
  freq_thres = options.freq_thres
  split_file = options.split_file
  drop_end_len = options.drop_end_len

  start = None
  end = None
  if ',' in ends:
    (start_str, end_str) = ends.split(',')
    if start_str:
      try:
        start = int(start_str)
      except ValueError:
        fail("Error: Bad --ends value \""+start_str+"\"")
    if end_str:
      try:
        end = int(end_str)
      except ValueError:
        fail("Error: Bad --ends value \""+end_str+"\"")

  if not arguments:
    parser.print_help()
    fail("Please provide an input VCF file or \"-\" for stdin.")
  else:
    vcfpath = arguments[0]

  # Setup I/O

  # global to allow cleanup() to close open filehandle on premature exit
  global infile
  if vcfpath == '-':
    infile = sys.stdin
  else:
    if os.path.exists(vcfpath):
      infile = open(vcfpath, 'rU')
    else:
      fail("Error: Input VCF file "+vcfpath+" not found.")

  if output == '-':
    outfile = sys.stdout
  elif not split_file:
    outfile = open(output, 'w')

  vcffile = vcfreader.VCFReader(infile)
  meta_header = vcffile.get_meta_header()
  column_header = vcffile.get_column_header()
  sample_names = vcffile.get_sample_names()

  # global to allow cleanup() to close open filehandles on premature exit
  global outfiles
  if split_file:
    (outfile_base, outfile_ext) = get_outfile_base(output, vcfpath)
    outfiles = get_outfiles(outfile_base, sample_names, outfile_ext)
    headers = []
    column_header_pre = '\t'.join(column_header.split('\t')[:9])
    for sample_name in sample_names:
      headers.append(meta_header+column_header_pre+'\t'+sample_name+'\n')
  else:
    outfiles = [outfile]
    headers = [meta_header+column_header]
  
  for (outfile, header) in zip(outfiles, headers):
    try:
      outfile.write(header)
    except IOError as ioerror:
      handle_broken_pipe(ioerror)

  buffers = []
  for file_ in outfiles:
    buffers.append(deque())

  # Main loop: parse VCF line-by-line

  first_pos = last_pos = 0
  if start is not None:
    first_pos = start
  for site in vcffile:
    pos = site.get_pos()
    if not first_pos:
      first_pos = pos
    if pos - first_pos < drop_end_len:
      continue
    last_pos = pos

    if split_file:
      sites = site.split()
    else:
      sites = [site]

    # sys.stderr.write(str(len(site.get_varcounts()))+"\t")
    for (site, buffer_, outfile) in zip(sites, buffers, outfiles):
      varcounts = site.get_varcounts(stranded=False)
      coverages = site.get_coverages()
      passed_vars = set()

      # sys.stderr.write(str(site.get_pos())+":"+str(len(varcounts))+"\t")
      for sample_name in varcounts:
        varcount = varcounts[sample_name]
        coverage = coverages[sample_name]
        for variant in varcount:
          # frequency threshold
          passed = varcount[variant]/coverage * 100 > freq_thres
          # coverage threshold
          passed = passed and coverage >= covg_thres
          # always keep types
          vartype = get_vartype(variant)
          if vartype in keep_types:
            passed = True
          # always remove types
          if vartype in rm_types:
            passed = False
          # final verdict
          if passed:
            passed_vars.add(variant)

      passed_alts = set()
      for variant in passed_vars:
        passed_alts.add(site.variant_to_alt(variant))
      alts_old = site.get_alt()
      alts_new = []
      #TODO: see if VCF allows the REF allele in the ALT column
      # want alts_new to be in same order as original, with additions on the end
      for alt in alts_old:
        if alt in passed_alts:
          alts_new.append(alt)
          passed_alts.discard(alt)
      for alt in passed_alts:
        if options.compliant and alt == site.get_ref():
          #TODO: add INFO flag
          continue
        alts_new.append(alt)    

      if alts_new:
        site.set_alt(alts_new)
        #TODO: recompute INFO and genotype values to be consistent with new ALTs
        #TODO: revise the REF to the shortest sequence necessary to represent
        #      the new set of ALTs
        line = str(site)+'\n'
        # sys.stderr.write(sample_name+" ")
        buffer_.append((line, pos))
        while len(buffer_) > drop_end_len:
          (line, pos) = buffer_.popleft()
          try:
            # sys.stderr.write('printing')
            outfile.write(line)
          except IOError as ioerror:
            handle_broken_pipe(ioerror)
        # sys.stderr.write("\t")
    # sys.stderr.write("\n")

  # flush buffers
  if end is not None:
    last_pos = end
  for (buffer_, outfile) in zip(buffers, outfiles):
    for (line, pos) in buffer_:
      if last_pos - pos > drop_end_len:
        try:
          outfile.write(line)
        except IOError as ioerror:
          handle_broken_pipe(ioerror)

  cleanup()


def get_outfile_base(outfile, infile):
  pre = post = ''
  if outfile and outfile != '-':
    pre_base = outfile
  elif infile and infile != '-':
    pre_base = infile
  else:
    return (pre, post)
  match = re.search(r'^(.+)\.vcf([.-][^/]+)?$', pre_base)
  if match:
    pre = match.group(1)
    post = match.group(2)
    if post is None:
      post = ''
  else:
    parts = pre_base.split('.')
    if len(parts) == 1:
      pre = pre_base
    else:
      pre = '.'.join(parts[0:-1])
      post = parts[-1]
  return (pre, post)


def get_outfiles(pre, sample_names, post):
  outfiles = []
  if not sample_names:
    fail("Error: No samples")
  if not pre:
    pre = sample_names[0]
  for sample in sample_names:
    base = pre+'-'+sample
    ext = '.vcf'+post
    final_name = get_unique_filename(base, ext)
    outfiles.append(open(final_name, 'w'))
  return outfiles


def get_unique_filename(base, ext):
  tries = 1
  proposed = base+ext
  while os.path.exists(proposed):
    tries+=1
    proposed = base+'-'+str(tries)+ext
    if tries >= 1000:
      fail("Error finding unused filename like "+proposed)
  return proposed


def get_vartype(variant):
  """Return 'S', 'I', or 'D' based on whether it's an SNV, insertion or
  deletion, respectively. If it is unrecognized, return None"""
  if len(variant) == 1:
    return 'S'
  elif variant[0] == 'd':
    try:
      int(variant[1:])
      return 'D'
    except ValueError:
      return None
  elif len(variant) > 1:
    return 'I'


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)


def handle_broken_pipe(ioerror):
  if ioerror.errno == errno.EPIPE:
    cleanup()
    sys.exit(0)


def cleanup():
  if isinstance(infile, file):
    infile.close()
  for outfile in outfiles:
    if isinstance(outfile, file):
      outfile.close()


def test():
  vcfhandle = open('../testdata/run19.vcf', 'rU')
  reader = vcfreader.VCFReader(vcfhandle)
  position = reader.next()
  position = reader.next()
  varcounts = position.get_varcounts()['S1']
  print varcounts

if __name__ == "__main__":
  main()
  # test()