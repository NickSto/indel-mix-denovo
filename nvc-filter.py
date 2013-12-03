#!/usr/bin/env python
from __future__ import division
import os
import sys
from optparse import OptionParser
from bioreaders import VCFReader
from bioreaders import VCFSite

OPT_DEFAULTS = {'output':'-', 'freq_thres':0.0, 'minor':False, 'strand':False,
  'keep_types':'', 'rm_types':'', 'compliant':False, 'split_file':False}
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
  parser.add_option('-f', '--freq_thres', dest='freq_thres', type='float',
    default=OPT_DEFAULTS.get('freq_thres'),
    help="""Frequency threshold. Variants less abundant than this will be
filtered out. Given in %. Default: %default%""")
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
    help="""Use this to find minority variants only. When selected, variants
that are the majority will be filtered (whether the same as the reference or
not""")
  parser.add_option('-s', '--strand', dest='strand', action='store_const',
    const=not OPT_DEFAULTS.get('strand'), default=OPT_DEFAULTS.get('strand'),
    help="""***NOT YET IMPLEMENTED*** Attempt to filter out indels which show
strand bias.""")
  parser.add_option('-c', '--compliant', dest='compliant',
    action='store_const', const=not OPT_DEFAULTS.get('compliant'),
    default=OPT_DEFAULTS.get('compliant'),
    help="""***NOT YET IMPLEMENTED*** Produce standard-compliant VCF output.
Currently this means the presence of the REF allele is indicated with an INFO
flag instead of putting it in the ALT column.""")
  parser.add_option('-S', '--split-file', dest='split_file',
    action='store_const', const=not OPT_DEFAULTS.get('split_file'),
    default=OPT_DEFAULTS.get('split_file'),
    help="""***NOT YET IMPLEMENTED*** Create an output VCF file for every
sample in the input.""")

  (options, arguments) = parser.parse_args()

  stranded = options.strand
  freq_thres = options.freq_thres

  if not arguments:
    parser.print_help()
    fail("Please provide an input VCF file or \"-\" for stdin.")
  else:
    vcfpath = arguments[0]

  if vcfpath == '-':
    infile = sys.stdin
  else:
    if os.path.exists(vcfpath):
      infile = open(vcfpath, 'rU')
    else:
      fail("Error: Input VCF file "+vcfpath+" not found.")

  if options.output == '-':
    outfile = sys.stdout
  else:
    outfile = open(options.output, 'w')

  #TODO: first implement it all with only one filter: freq_thres
  vcfreader = VCFReader(infile)
  header = vcfreader.get_header()
  try:
    outfile.write(header)
  except IOError:
    sys.exit()
  for site in vcfreader:

    varcounts = site.get_varcounts(stranded=False)
    coverages = site.get_coverages()
    passed = set()
    for sample_name in varcounts:

      varcount = varcounts[sample_name]
      coverage = coverages[sample_name]
      for variant in varcount:
        if varcount[variant]/coverage * 100 > freq_thres:
          passed.add(variant)

    passed_alts = []
    for variant in passed:
      passed_alts.append(site.variant_to_alt(variant))
    alts_old = site.get_alt()
    alts_new = []
    #TODO: add reference allele if it passes thresholds
    # - might require making alts_new by adding passed_alts instead of
    #   subtracting from alts_old
    #TODO: see if VCF allows the REF allele in the ALT column
    for alt in alts_old:
      if alt in passed_alts:
        alts_new.append(alt)

    #TODO: recompute INFO and genotype values to be consistent with new ALTs

    if alts_new:
      site.set_alt(alts_new)
      try:
        outfile.write(str(site)+'\n')
      except IOError:
        break

  if isinstance(infile, file):
    infile.close()
  if isinstance(outfile, file):
    outfile.close()



def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)


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