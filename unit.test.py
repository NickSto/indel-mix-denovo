#!/usr/bin/env python
"""The guts of the read selection and statistics calculation will be done by the
bamslicer module, so that it can also be used by nvc-filter.py. Then
nvc-filter.py can use the statistics as a filter."""
import bamslicer
import optparse
import sys
import os

def main():

  FUNCTIONS = {
    'bamslicer.get_reads_and_stats':bamslicer_get_reads_and_stats,
  }

  OPT_DEFAULTS = {'variants':'', 'no_reads':False}
  USAGE = "USAGE: %prog [options] function.to.test reads.bam"
  DESCRIPTION = """Run test on a given function and input BAM and print results.
    Give one of the following function names: """+', '.join(FUNCTIONS)
  EPILOG = """"""

  parser = optparse.OptionParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)

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

  if len(arguments) == 2:
    (function, bamfilepath) = arguments
  else:
    parser.print_help()
    fail("Error: Please provide a BAM file and a list of variants.")

  if function not in FUNCTIONS:
    fail('Error: function "'+function+'" not supported. Please pick one from '
      +'the list: '+', '.join(FUNCTIONS))
  if not os.path.exists(bamfilepath):
    fail('Error: cannot find BAM file "'+bamfilepath+'"')

  FUNCTIONS[function](bamfilepath, options)


def bamslicer_get_reads_and_stats(bamfilepath, options):
  if options.variants:
    variants = variants_from_str(options.variants)
  else:
    fail("Error: A list of variants is required for this function.")

  #TODO: make sure BAM is sorted
  #TODO: take chrom in to account when sorting variants: use order in BAM header
  variants.sort(key=lambda variant: variant['coord'])

  (read_sets, stat_sets) = bamslicer.get_reads_and_stats(bamfilepath, variants)

  for (reads, stats, variant) in zip(read_sets, stat_sets, variants):
    print variant['chrom'], variant['coord'], variant['type']
    print "  flags:", stats['flags']
    sys.stdout.write("  mapqs:")
    for (quality, total) in enumerate(stats['mapqs']):
      if total > 0:
        sys.stdout.write(" "+str(quality)+"="+str(total))
    print
    if not options.no_reads:
      print "  reads:"
      for read in reads:
        print "   ", read.get_read_name().split(':')[-1]


def variants_from_str(variants_str):
  """Parse the list of variants passed in from the command line.
  Example list: "chrM:310-S,1:2345-D:2,pUC18:4210-I:GAT"
  Will return: [
    {'chrom':'chrM', 'coord':310, 'type':'S', 'alt':None},
    {'chrom':'1', 'coord':2345, 'type':'D', 'alt':'2'},
    {'chrom':'pUC18', 'coord':4210, 'type':'I', 'alt':'GAT'},
  ]
  """
  variants = []
  for variant_str in variants_str.split(','):
    if '-' not in variant_str:
      fail('Error: Incorrect format in variants list: "'+variant_str+'"')
    fields = variant_str.split('-')
    location = '-'.join(fields[:-1])
    var_details = fields[-1]
    try:
      (chrom, coord) = location.split(':')
      coord = int(coord)
    except ValueError:
      fail('Error: Incorrect format in variants list: "'+variant_str+'"')
    if ':' in var_details:
      (vartype, alt) = var_details.split(':')
      if not (vartype in 'SID' and valid_variant(vartype, alt)):
        fail('Error: Incorrect format in variants list: "'+variant_str+'"')
    else:
      vartype = var_details
      alt = None
    variants.append({'chrom':chrom, 'coord':coord, 'type':vartype, 'alt':alt})
  return variants


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
