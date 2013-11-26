#!/usr/bin/env python
import os
import sys
from optparse import OptionParser
from FormatReaders import VCFReader
from FormatReaders import VCFPosition

OPT_DEFAULTS = {'types':'ID', 'freq_thres':0.0, 'minor':False, 'strand':False}
USAGE = "USAGE: %prog [options] input.vcf"
DESCRIPTION = """Filters the output of the Naive Variant Caller, removing 
spurious indels that fail to satisfy the given criteria."""
EPILOG = """"""

def main():

  parser = OptionParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)

  parser.add_option('-t', '--types', dest='types',
    default=OPT_DEFAULTS.get('types'), help='The types of indels to find. '
    +'Give a string of letters, "I" for insertions, "D" for deletions. Only '
    +'indels of the selected type(s) will be preserved in the output. '
    +'Default: %default')
  parser.add_option('-f', '--freq_thres', dest='freq_thres', type='float',
    default=OPT_DEFAULTS.get('freq_thres'),
    help='Frequency threshold. Indels less abundant than this will be '
    +'filtered out.')
  parser.add_option('-m', '--minor', dest='minor', action='store_const',
    const=not OPT_DEFAULTS.get('minor'), default=OPT_DEFAULTS.get('minor'),
    help='When selected, indels that are the majority variant will be filtered '
    +'out. This allows looking for minority indels only.')
  parser.add_option('-s', '--strand', dest='strand', action='store_const',
    const=not OPT_DEFAULTS.get('strand'), default=OPT_DEFAULTS.get('strand'),
    help='Attempt to filter out indels which show strand bias. ***NOT YET '
    +'IMPLEMENTED***')

  (options, arguments) = parser.parse_args()

  if not arguments:
    parser.print_help()
    fail("Please provide an input VCF file or \"-\" for stdin.")
  else:
    vcfpath = arguments[0]

  if vcfpath == '-':
    vcfhandle = sys.stdin
  else:
    if os.path.exists(vcfpath):
      vcfhandle = open(vcfpath, 'rU')
    else:
      fail("Error: Input VCF file "+vcfpath+" not found.")

  vcfreader = VCFReader(vcfhandle)
  header = vcfreader.get_header()
  for position in vcfreader:

    varcounts = position.get_varcounts()
    for sample_name in varcounts:

      if options.strand:
        varcount = varcounts[sample_name]
        #TODO
      else:
        varcount = strandless(varcounts[sample_name])


def strandless(varcount):
  varcount_new = {}
  for variant in varcount:
    variant_new = variant.lstrip('+-')
    varcount_new[variant_new] = (varcount[variant]
      + varcount_new.get(variant_new, 0))
  return varcount_new


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)


def test():
  from FormatReaders import VCFReader
  from FormatReaders import VCFPosition
  vcfhandle = open('../testdata/real-mit.vcf.in', 'rU')
  reader = VCFReader(vcfhandle)
  position = reader.next()
  position = reader.next()
  varcounts = position.get_varcounts()['S1']
  print varcounts
  print strandless(varcounts)

if __name__ == "__main__":
  #main()
  test()