#!/usr/bin/env python
from __future__ import division
import os
import re
import sys
import errno
from optparse import OptionParser
from bioreaders import VCFReader
from bioreaders import VCFSite

OPT_DEFAULTS = {'output':'-', 'freq_thres':0.0, 'covg_thres':0, 'minor':False,
  'strand':False, 'keep_types':'', 'rm_types':'', 'compliant':False,
  'split_file':False}
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
  parser.add_option('-C', '--compliant', dest='compliant',
    action='store_const', const=not OPT_DEFAULTS.get('compliant'),
    default=OPT_DEFAULTS.get('compliant'),
    help="""***NOT YET IMPLEMENTED*** Produce standard-compliant VCF output.
Currently this means the presence of the REF allele is indicated with an INFO
flag instead of putting it in the ALT column.""")
  parser.add_option('-S', '--split-file', dest='split_file',
    action='store_const', const=not OPT_DEFAULTS.get('split_file'),
    default=OPT_DEFAULTS.get('split_file'),
    help="""Create an output VCF file for every sample in the input. Currently,
this leaves some fields in an inconsistent state. The ones that will be correct
are CHROM, POS, REF, ALT, and the sample genotype columns.""")

  (options, arguments) = parser.parse_args()

  output = options.output
  stranded = options.strand
  covg_thres = options.covg_thres
  freq_thres = options.freq_thres
  split_file = options.split_file
  keep_types = options.keep_types
  rm_types = options.rm_types

  if not arguments:
    parser.print_help()
    fail("Please provide an input VCF file or \"-\" for stdin.")
  else:
    vcfpath = arguments[0]

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
  else:
    outfile = open(output, 'w')

  vcfreader = VCFReader(infile)
  meta_header = vcfreader.get_meta_header()
  column_header = vcfreader.get_column_header()
  sample_names = vcfreader.get_sample_names()

  global outfiles
  if split_file:
    (outfile_base, outfile_ext) = get_outfile_base(output, vcfpath)
    outfiles = get_outfiles(outfile_base, sample_names, outfile_ext)
    headers = []
    column_header_pre = '\t'.join(column_header.split('\t')[:8])
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

  for site in vcfreader:
    if split_file:
      sites = site.split()
    else:
      sites = [site]

    for (site, outfile) in zip(sites, outfiles):
      varcounts = site.get_varcounts(stranded=False)
      coverages = site.get_coverages()
      passed_vars = set()

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
        try:
          outfile.write(str(site)+'\n')
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
  deletion, respectively."""
  if len(variant) == 1:
    return 'S'
  elif re.match(r'^d\d+$', variant):
    return 'D'
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
  reader = VCFReader(vcfhandle)
  position = reader.next()
  position = reader.next()
  varcounts = position.get_varcounts()['S1']
  print varcounts

if __name__ == "__main__":
  main()
  # test()