#!/usr/bin/env python
"""The guts of the read selection and statistics calculation will be done by the
bamslicer module, so that it can also be used by nvc-filter.py. Then
nvc-filter.py can use the statistics as a filter."""
from __future__ import division
import bioreaders
import bamslicer
import optparse
import sys
import os
# ./inspect-reads.py ~/backuphide/bfx/R19S5-new-nm.bam -H -v chrM-R19S5-new-nm-dedup:15873-I,chrM-R19S5-new-nm-dedup:15873-D,chrM-R19S5-new-nm-dedup:13718-D,chrM-R19S5-new-nm-dedup:11571-D,chrM-R19S5-new-nm-dedup:3757-D
# ./inspect-reads.py tests/cigar-tests.bam -H -v chrM:5-I,chrM:199-D,chrM:199-I,chrM:3106-D,chrM:6110-D,chrM:16568-I

OPT_DEFAULTS = {'variants':'', 'vcf':'', 'output_bam':'', 'all':False,
  'opposing':False, 'containing':False, 'strand_bias':0, 'mate_bias':0}
USAGE = "USAGE: %prog [options] (-v variantslist|-V variants.vcf) reads.bam"
DESCRIPTION = """Retrieve the reads supporting a given set of variants, and
report statistics on them. Provide the variants in a VCF or in a list on the
command line."""
EPILOG = """WARNING: Work in progress. Not yet implemented features:
--containing, --output-bam, --all, and --opposing options, and
considering the "alt" information of variants (the only thing considered right
now is the location and type). And most important of all, ***SNVS ARE NOT
SUPPORTED!*** Right now only indels are recognized."""

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
  parser.add_option('-H', '--human', dest='human', action='store_const',
    const=not OPT_DEFAULTS.get('human'),default=OPT_DEFAULTS.get('human'),
    help="""Print statistics in a human-readable format.""")
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
  parser.add_option('-s', '--strand-bias', dest='strand_bias', type='float',
    default=OPT_DEFAULTS.get('strand_bias'),
    help="""Filter out variants with a strand bias value greater than or equal
to this. The strand bias statistic is method 1 (SB) of Guo et al., 2012. A
typical cutoff is 1.0.""")
  parser.add_option('-m', '--mate-bias', dest='mate_bias', type='float',
    default=OPT_DEFAULTS.get('mate_bias'),
    help="""Filter out variants with a higher mate bias value than this. The
statistic is calculated identically to the strand bias one, but replacing
forward/reverse strand with first/second mate in the pair.""")

  (options, arguments) = parser.parse_args()

  if arguments:
    #TODO: support multiple BAMs
    bamfilepath = arguments[0]
  else:
    parser.print_help()
    fail("\nError: Please provide a BAM file and a list of variants.")

  if options.vcf:
    variants = variants_from_vcf(options.vcf)
  elif options.variants:
    variants = variants_from_str(options.variants)
  else:
    parser.print_help()
    fail("\nError: Please provide a list of variants in either a VCF file or a "
      +"command line option.")

  #TODO: make sure BAM is sorted
  #TODO: take chrom in to account when sorting variants: use order in BAM header
  # sort variants by coordinate
  variants.sort(key=lambda variant: variant['coord'])

  (read_sets, stat_sets) = bamslicer.get_reads_and_stats(bamfilepath, variants)

  for (variant, reads, stats) in zip(variants, read_sets, stat_sets):
    if filter_out(variant, reads, stats, options):
      continue
    if options.human:
      sys.stdout.write(human_stats(variant, reads, stats))


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

#TODO
def valid_variant(vartype, alt):
  """Make sure the alt is the correct format, based on the vartype"""
  return True


def variants_from_vcf(vcffilename):
  #TODO: support multiple samples
  variants = []
  vcfreader = bioreaders.VCFReader(open(vcffilename, 'rU'))
  for site in vcfreader:
    chrom = site.get_chrom()
    pos = site.get_pos()
    ref = site.get_ref()
    for alt in site.get_alt():
      #TODO: make sure this is the desired handling of reference allele
      if alt == ref:
        continue
      varstr = site.alt_to_variant(alt)
      variants.append(parse_varstr(varstr, chrom, pos))
  return variants


def parse_varstr(varstr, chrom, coord):
  """Take in a variant string like in the NVC sample columns and return a data
  structure of the format specified in variants_from_str().
  Input string examples:
    SNV:        'A', 'C'
    insertion:  'AG', 'GAATC'
    deletion:   'd1', 'd4'
  """
  #TODO: handle reference allele (currently classified as an SNV)
  if len(varstr) == 1:
    return {'chrom':chrom, 'coord':coord, 'type':'S', 'alt':varstr}
  elif varstr[0] == 'd':
    try:
      int(varstr[1:])
      return {'chrom':chrom, 'coord':coord, 'type':'D', 'alt':varstr[1:]}
    except ValueError:
      return None
  elif len(varstr) > 1:
    return {'chrom':chrom, 'coord':coord, 'type':'I', 'alt':varstr}


def filter_out(variant, reads, stats, options):
  filter_out = False
  if options.strand_bias and stats['strand_bias'] is not None:
    if stats['strand_bias'] > options.strand_bias:
      filter_out = True
      return filter_out
  if options.mate_bias and stats['mate_bias'] is not None:
    if stats['mate_bias'] > options.mate_bias:
      filter_out = True
      return filter_out
  return filter_out


def human_stats(variant, reads, stats):
  output = ""
  output += variant['chrom']+':'+str(variant['coord'])+'-'+variant['type']+"\n"
  output += "  coverage:         "+str(stats['coverage'])+"\n"
  total = stats['supporting']
  freq = str(round(100*stats['freq'],2))+"%"
  output += "  supporting reads: "+str(total)+" ("+freq+")\n"
  if total == 0:
    return output
  flags = stats['flags']
  output += ("  unmapped:            "+pct(flags[2], total)+"\n")
  output += ("  not proper pair:     "+pct(total-flags[1], total)+"\n")
  output += ("  forward/reverse:     "+pct(total-flags[4], total)+"/"
    +pct(flags[4], total, const_width=False)+"\n")
  output += ("  1st/2nd mate:        "+pct(flags[6], total)+"/"
    +pct(flags[7], total, const_width=False)+"\n")
  output += ("  secondary alignment: "+pct(flags[8], total)+"\n")
  output += ("  marked duplicate:    "+pct(flags[11], total)+"\n")
  mapqs = stats['mapqs']
  output += "  MAPQ == 0:  "+pct(mapqs[0], total)+"\n"
  output += "  MAPQ >= 20: "+pct(mapq_ge_thres(mapqs, 20), total)+"\n"
  output += "  MAPQ >= 30: "+pct(mapq_ge_thres(mapqs, 30), total)+"\n"
  output += "  MAPQ == "+str(len(mapqs)-1)+": "+pct(mapqs[-1], total)+"\n"
  output += "  strand bias: "+str(stats['strand_bias'])+"\n"
  output += "  mate bias:   "+str(stats['mate_bias'])+"\n"
  return output


def pct(count, total, decimals=1, const_width=True):
  """Calculate a percentage and return a formatted, constant-width and constant-
  precision string."""
  if const_width:
    width = str(decimals+3)
    if count <= 0 or total <= 0:
      return ("%"+width+"d%%") % 0
    if count == total:
      return ("%"+width+"d%%") % 100
    fraction = 100*count/total
    return ("%"+width+"."+str(decimals)+"f%%") % fraction
  else:
    if count <= 0 or total <= 0:
      return "0%"
    if count == total:
      return "100%"
    fraction = 100*count/total
    return str(round(fraction, decimals))+"%"


def mapq_ge_thres(mapqs, thres):
  """Return the number of reads greater than or equal to a given MAPQ."""
  count = 0
  i = thres
  while i < len(mapqs):
    count += mapqs[i]
    i+=1
  return count


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
