#!/usr/bin/env python
import bamslicer
import optparse
import sys
import os

def main():

  FUNCTIONS = {
    'bamslicer.get_reads_and_stats':bamslicer_get_reads_and_stats,
    'inspect-reads.variants_from_vcf':inspect_reads_variants_from_vcf,
  }

  OPT_DEFAULTS = {'bam':'', 'vcf':'', 'variants':'', 'no_reads':False}
  USAGE = "USAGE: %prog [options] function.to.test"
  DESCRIPTION = """Run test on a given function and input BAM and print results.
    Give one of the following function names: """+', '.join(FUNCTIONS)
  EPILOG = """"""

  parser = optparse.OptionParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)

  parser.add_option('-b', '--bam', dest='bam', default=OPT_DEFAULTS.get('bam'),
    help="""Input BAM file. Necessary for some functions.""")
  parser.add_option('-V', '--vcf', dest='vcf', default=OPT_DEFAULTS.get('vcf'),
    help="""Input VCF file. Necessary for some functions.""")
  parser.add_option('-v', '--variants', dest='variants',
    default=OPT_DEFAULTS.get('variants'),
    help="""Use these variants. Give a comma-separated list, in the format
"chrom:pos-type[:alt]" e.g. "chr1:2345-D:2". "pos" is the 1-based coordinate
of the SNV, or the base before the insertion or deletion. "type" is "S", "I", or
"D", for SNV, insertion, or deletion. "alt" is optional, and specifies the
alternate allele: the SNV alternate base, the length of the deletion, or the
inserted sequence (including the preceding base, as in VCF). If "alt" is not
provided, it will select any read with that type of variant at that location.""")
  parser.add_option('-R', '--no-reads', dest='no_reads', action='store_const',
    const=not OPT_DEFAULTS.get('no_reads'), default=OPT_DEFAULTS.get('no_reads'),
    help="Don't print info on individual reads.")

  (options, arguments) = parser.parse_args()

  if len(arguments) == 1:
    function = arguments[0]
  else:
    parser.print_help()
    fail("Error: Please only give a function name as a single positional "
      +"argument.")

  if function not in FUNCTIONS:
    fail('Error: function "'+function+'" not supported. Please pick one from '
      +'the list: '+', '.join(FUNCTIONS))

  FUNCTIONS[function](options)


##### UNIT TEST FUNCTIONS #####

def bamslicer_get_reads_and_stats(options):
  # kludge to get around dash in name
  inspect_reads = __import__("inspect-reads")
  if options.bam:
    if not os.path.exists(options.bam):
      fail('Error: Cannot find BAM file "'+options.bam+'"')
  else:
    fail("Error: A BAM file is required for this function.")
  if options.variants:
    variants = inspect_reads.variants_from_str(options.variants)
  else:
    fail("Error: A list of variants is required for this function.")

  #TODO: make sure BAM is sorted
  #TODO: take chrom in to account when sorting variants: use order in BAM header
  variants.sort(key=lambda variant: variant['coord'])

  (read_sets, stat_sets) = bamslicer.get_reads_and_stats(options.bam, variants)

  for (reads, stats, variant) in zip(read_sets, stat_sets, variants):
    print variant['chrom'], variant['coord'], variant['type']
    print "  flags:", stats['flags']
    sys.stdout.write("  mapqs:")
    for (quality, total) in enumerate(stats['mapqs']):
      if total > 0:
        sys.stdout.write(' %s=%d' % (quality, total))
    print
    if not options.no_reads:
      print "  reads:"
      for read in reads:
        print "   ", read.get_read_name().split(':')[-1]


def inspect_reads_variants_from_vcf(options):
  # kludge to get around dash in name
  inspect_reads = __import__("inspect-reads")
  if options.vcf:
    if not os.path.exists(options.vcf):
      fail('Error: Cannot find VCF file "'+options.vcf+'"')
  else:
    fail('Error: A VCF file is required for this function.')
  variants = inspect_reads.variants_from_vcf(options.vcf)
  for variant in variants:
    sys.stdout.write(str(variant.get('chrom')))
    sys.stdout.write(':'+str(variant.get('coord')))
    sys.stdout.write('-'+str(variant.get('type')))
    if variant.get('alt') is not None:
      sys.stdout.write(':'+str(variant.get('alt')))
    print


##### UTILITY FUNCTIONS #####


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
