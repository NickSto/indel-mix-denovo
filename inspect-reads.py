#!/usr/bin/env python
"""The guts of the read selection and statistics calculation will be done by the
bamslicer module, so that it can also be used by nvc-filter.py. Then
nvc-filter.py can use the statistics as a filter."""
import bamslicer
import sys
import os
import optparse

OPT_DEFAULTS = {'variants':'', 'vcf':'', 'output_bam':'', 'all':False,
  'opposing':False, 'containing':False}
USAGE = "USAGE: %prog [options] reads.bam"
DESCRIPTION = """Retrieve the reads supporting a given set of variants, and
report statistics on them. Provide the variants in a VCF or in a list on the
command line."""
EPILOG = """WARNING: Work in progress. Not yet implemented features:
--containing option, --vcf option, considering the "alt" information of variants
(the only thing considered right now is the location and type). And most
important of all, ***SNVS ARE NOT SUPPORTED!*** Right now only indels are
recognized."""

def main():

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
  parser.add_option('-V', '--vcf', dest='vcf',
    default=OPT_DEFAULTS.get('vcf'),
    help='Use the variants in the ALT column of this VCF file.')
  parser.add_option('-o', '--output-bam', dest='output_bam',
    default=OPT_DEFAULTS.get('output_bam'),
    help='Output the selected reads to this BAM file.')
  parser.add_option('-a', '--all', dest='all', action='store_const',
    const=not OPT_DEFAULTS.get('all'), default=OPT_DEFAULTS.get('all'),
    help="""Select all the reads covering the variants, not just those
supporting the variant. The statistics reported will be on these reads, not the
supporting reads.""")
  parser.add_option('-O', '--opposing', dest='opposing', action='store_const',
    const=not OPT_DEFAULTS.get('opposing'),default=OPT_DEFAULTS.get('opposing'),
    help="""Select the reads OPPOSING the variants. The statistics reported will
be on these reads, not the supporting reads.""")
  parser.add_option('-C', '--containing', dest='containing',
    action='store_const', const=not OPT_DEFAULTS.get('containing'),
    default=OPT_DEFAULTS.get('containing'),
    help="""For deletions, select any read with a deletion which *contains* the
specified one. So if chr1:2345-D:2 is a variant of interest, and a read has the
variant chr1:2344-D:3, it will select that read even though the variants don't
have the same start coordinates.""")

  (options, arguments) = parser.parse_args()

  if arguments:
    #TODO: support multiple BAMs
    bamfilepath = arguments[0]
  else:
    parser.print_help()
    fail("Error: Please provide a BAM file and a list of variants.")

  if options.vcf:
    variants = variants_from_vcf(options.vcf)
  elif options.variants:
    variants = variants_from_str(options.variants)
  else:
    parser.print_help()
    fail("Error: Please provide a list of variants in either a VCF file or a "
      +"command line option.")

  #TODO: make sure BAM is sorted
  #TODO: take chrom in to account when sorting variants: use order in BAM header
  variants.sort(key=lambda variant: variant['coord'])

  (read_sets, stat_sets) = bamslicer.get_reads_and_stats(bamfilepath, variants)




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
    try:
      (location, var_details) = variant_str.split('-')
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

#TODO
def valid_variant(vartype, alt):
  """Make sure the alt is the correct format, based on the vartype"""
  return True

#TODO
def variants_from_vcf(vcffilename):
  variants = []
  return variants


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
