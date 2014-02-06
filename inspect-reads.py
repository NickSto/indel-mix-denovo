#!/usr/bin/env python
# requires Python 2.7
"""The guts of the read selection and statistics calculation will be done by the
bamslicer module, so that it can also be used by nvc-filter.py. Then
nvc-filter.py can use the statistics as a filter."""
from __future__ import division
import collections
import simplewrap
import vcfreader
import bamslicer
import argparse
import sys
import os
# ./inspect-reads.py ~/backuphide/bfx/R19S5-new-nm.bam -H -v chrM-R19S5-new-nm-dedup:15873-I,chrM-R19S5-new-nm-dedup:15873-D,chrM-R19S5-new-nm-dedup:13718-D,chrM-R19S5-new-nm-dedup:11571-D,chrM-R19S5-new-nm-dedup:3757-D
# ./inspect-reads.py tests/cigar-tests.bam -H -v chrM:5-I,chrM:199-D,chrM:199-I,chrM:3106-D,chrM:6110-D,chrM:16568-I

"""Need to use argparse.RawTextHelpFormatter to preserve formatting in the
description of columns in the tsv output. But to still accommodate different
terminal widths, dynamic wrapping with simplewrap will be necessary."""
WIDTH = simplewrap.termwidth(70)
wrapper = simplewrap.Wrapper(width=WIDTH)
wrap = wrapper.wrap

OPT_DEFAULTS = {'human':True, 'tsv':False, 'vartypes':'ID',
  'strand_bias':0, 'mate_bias':0}
USAGE = "%(prog)s [options] -v variantslist -V variants.vcf reads.bam"
DESCRIPTION = wrap("Retrieve the reads supporting a given set of "
  "variants, and report statistics on them. Provide the variants in a VCF "
  "and/or in a list on the command line.")
EPILOG = wrap("WARNING: Work in progress. Not yet implemented "
  "features:\n"
  "--containing --output-bam --all --opposing options\n"
  "Considering the \"alt\" information of variants (the only thing considered "
  "right now is the location and type).\n"
  "And most important of all, ***SNVS ARE NOT SUPPORTED!*** Right now only "
  "indels are recognized.")

VALID_BASES = 'GATCNgatcn'

def main():

  wrapper.width = WIDTH-24

  parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION,
    epilog=EPILOG, formatter_class=argparse.RawTextHelpFormatter)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('bamfilepath', metavar='reads.bam',
    help='the input reads')
  parser.add_argument('-v', '--variants',
    help=wrap('Use these variants. Give a comma-separated list, in the format '
      '"chrom:pos-type[:alt]" e.g. "chr1:2345-D:2". "pos" is the 1-based '
      'coordinate of the SNV, or the base before the insertion or deletion. '
      '"type" is "S", "I", or "D", for SNV, insertion, or deletion. "alt" is '
      'optional, and specifies the alternate allele: the SNV alternate base, '
      'the length of the deletion, or the inserted sequence (including the '
      'preceding base, as in VCF). If "alt" is not provided, it will select '
      'any read with that type of variant at that location.'))
  parser.add_argument('-V', '--vcf',
    help='Use the variants in the ALT column of this VCF file.')
  parser.add_argument('-o', '--output-bam',
    help='Output the selected reads to this BAM file.')
  parser.add_argument('-H', '--human', action='store_true',
    help='Print statistics in a human-readable format (default).')
  parser.add_argument('-t', '--tsv', action='store_true',
    help=wrap('Print statistics in a tab-delimited columnar format. The first '
      'four columns are the same data as described in the --variants format. '
      'The strand and mate bias statistics are described in the -s and -m '
      'options. The null value is "None".\n'
      'Columns:')
      +'\n'+wrap(
      '1:  chrom\n'
      '2:  coord\n'
      '3:  variant type\n'
      '4:  alt allele\n'
      '5:  coverage, in # of reads\n'
      '6:  # of variant-supporting reads\n'
      '7:  %% frequency\n'
      '8:  %% of supporting reads that are unmapped\n'
      '9:  %% not mapped in proper pair\n'
      '10: %% on the forward strand\n'
      '11: %% that are the 1st mate in the pair\n'
      '12: %% that are a secondary alignment\n'
      '13: %% that are marked duplicates\n'
      '14: %% with a MAPQ == 0\n'
      '15: %% with a MAPQ >= 20\n'
      '16: %% with a MAPQ >= 30\n'
      '17: %% with the highest MAPQ\n'
      '18: the highest MAPQ in the supporting reads\n'
      '19: strand bias\n'
      '20: mate bias\n'
      '21: The total number of supporting reads with each SAM flag. This is a '
          'comma-separated list of the total for each flag, from lowest to '
          'highest bit value.',
      lspace=4, indent=-4))
  parser.add_argument('-L', '--no-labels', action='store_true',
    help=wrap('If tsv output is selected, do not print column labels (normally '
      'the first line, begins with #).'))
  parser.add_argument('-T', '--vartypes', metavar='TYPES',
    help=wrap('Only consider these variant types. Give a string of letters, '
      'e.g. "SID" to keep SNVs (S), insertions (I), and deletions (D). '
      'Default: "%(default)s"'))
  parser.add_argument('-a', '--all', action='store_true',
    help=wrap('Select all the reads covering the variants, not just those '
      'supporting the variant. The statistics reported will be on these reads, '
      'not the supporting reads.'))
  parser.add_argument('-O', '--opposing', action='store_true',
    default=OPT_DEFAULTS.get('opposing'),
    help=wrap('Select the reads OPPOSING the variants. The statistics reported '
      'will be on these reads, not the supporting reads.'))
  parser.add_argument('-C', '--containing', action='store_true',
    help=wrap('For deletions, select any read with a deletion which *contains* '
      'the specified one. So if chr1:2345-D:2 is a variant of interest, and a '
      'read has the variant chr1:2344-D:3, it will select that read even '
      'though the variants don\'t have the same start coordinates.'))
  parser.add_argument('-s', '--strand-bias', type=float, metavar='SBIAS',
    help=wrap('Filter out variants with a strand bias value greater than or '
      'equal to this. The strand bias statistic is method 1 (SB) of Guo et '
      'al., 2012. A typical cutoff is 1.0.'))
  parser.add_argument('-m', '--mate-bias', type=float, metavar='MBIAS',
    help=wrap('Filter out variants with a higher mate bias value than this. '
      'The statistic is calculated identically to the strand bias one, but '
      'replacing forward/reverse strand with first/second mate in the pair.'))

  args = parser.parse_args()

  variants = []
  if args.vcf:
    variants.extend(variants_from_vcf(args.vcf))
  if args.variants:
    variants.extend(variants_from_str(args.variants))
  if not variants:
    parser.print_help()
    fail("\nError: Please provide a list of variants in either a VCF file or a "
      +"command line option.")
  # eliminate variants that aren't one of the specified types
  variants[:] = filter(lambda var: var['type'] in args.vartypes, variants)

  if args.tsv:
    outformat = 'tsv'
  else:
    outformat = 'human'

  #TODO: make sure BAM is sorted
  #TODO: take chrom in to account when sorting variants: use order in BAM header
  # sort variants by coordinate
  variants.sort(key=lambda variant: variant['coord'])

  #TODO: make readgroup-aware
  #      probably have to return a dict of reads and stats for every variant
  (read_sets, stat_sets) = bamslicer.get_reads_and_stats(
    args.bamfilepath, variants)

  for (variant, reads, stats) in zip(variants, read_sets, stat_sets):
    if filter_out(variant, reads, stats, args):
      continue
    output_stats = get_output_stats(variant, stats)
    if outformat == 'human':
      sys.stdout.write(humanize(output_stats))
    elif outformat == 'tsv':
      sys.stdout.write("\t".join(map(str, output_stats.values()))+"\n")



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
      if not valid_variant(vartype, alt):
        fail('Error: Incorrect format in variants list: "'+variant_str+'"')
    else:
      vartype = var_details
      alt = None
    variants.append({'chrom':chrom, 'coord':coord, 'type':vartype, 'alt':alt})
  return variants


def valid_variant(vartype, alt):
  """Make sure the vartype is valid, and that the alt is the correct format for
  the vartype."""
  if vartype == 'S':
    if len(alt) == 1 and alt in VALID_BASES:
      return True
  elif vartype == 'I':
    if len(alt) < 2:
      return False
    for base in alt:
      if base not in VALID_BASES:
        return False
    return True
  elif vartype == 'D':
    if len(alt) < 2:
      return False
    if alt[0] != 'd':
      return False
    try:
      int(alt[1:])
      return True
    except ValueError:
      return False
  return False


def variants_from_vcf(vcffilename):
  #TODO: support multiple samples
  variants = []
  vcf = vcfreader.VCFReader(open(vcffilename, 'rU'))
  for site in vcf:
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


def filter_out(variant, reads, stats, args):
  filter_out = False
  if args.strand_bias and stats['strand_bias'] is not None:
    if stats['strand_bias'] > args.strand_bias:
      filter_out = True
      return filter_out
  if args.mate_bias and stats['mate_bias'] is not None:
    if stats['mate_bias'] > args.mate_bias:
      filter_out = True
      return filter_out
  return filter_out


def get_output_stats(variant, stats):
  output = collections.OrderedDict()
  output['chrom']       = variant.get('chrom')
  output['coord']       = variant.get('coord')
  output['type']        = variant.get('type')
  output['alt']         = variant.get('alt')
  output['coverage']    = stats.get('coverage')
  total = stats.get('supporting')
  output['supporting']  = total
  output['freq']        = round(100*stats['freq'], 2)
  if total == 0:
    # the rest of the fields aren't meaningful if no supporting reads
    for i in range(14):
      output[i] = None
    return output
  flags = stats['flags']
  output['unmapped']    = pct(flags[2], total)
  output['improper']    = pct(total-flags[1], total)
  output['forward']     = pct(total-flags[4], total)
  output['1st-mate']    = pct(flags[6], total)
  output['2ndary']      = pct(flags[8], total)
  output['duplicate']   = pct(flags[11], total)
  mapqs = stats['mapqs']
  output['mapq0']       = pct(mapqs[0], total)
  output['mapq20']      = pct(mapq_ge_thres(mapqs, 20), total)
  output['mapq30']      = pct(mapq_ge_thres(mapqs, 30), total)
  output['mapq-top']    = pct(mapqs[-1], total)
  output['mapq-best']   = len(mapqs)-1
  output['strand_bias'] = stats.get('strand_bias')
  output['mate_bias']   =  stats.get('mate_bias')
  output['flags'] = ','.join([str(flag) for flag in flags])
  return output


def humanize(output_stats):
  """Format variant statistics into a human-readable printout.
  Input: output of get_output_stats()
  Output: a string ready to print (including ending newline)
  """
  # print 'N/A' as null value
  for key in output_stats:
    if output_stats[key] is None and key != 'alt':
      output_stats[key] = 'N/A'
  output = ""
  output += (output_stats['chrom']+':'+str(output_stats['coord'])+'-'
    +output_stats['type'])
  if output_stats['alt']:
    output += ':'+output_stats['alt']
  output += "\n"
  output += "  coverage:         %s\n" % output_stats['coverage']
  output += "  supporting reads: %s (%s%%)\n" % (output_stats['supporting'],
    output_stats['freq'])
  if output_stats['supporting'] == 0:
    return output
  output += "  unmapped:            %s\n" % pct_str(output_stats['unmapped'])
  output += "  not proper pair:     %s\n" % pct_str(output_stats['improper'])
  output += "  forward strand:      %s\n" % pct_str(output_stats['forward'])
  output += "  1st mate:            %s\n" % pct_str(output_stats['1st-mate'])
  output += "  secondary alignment: %s\n" % pct_str(output_stats['2ndary'])
  output += "  marked duplicate:    %s\n" % pct_str(output_stats['duplicate'])
  output += "  MAPQ == 0:  %s\n" % pct_str(output_stats['mapq0'])
  output += "  MAPQ >= 20: %s\n" % pct_str(output_stats['mapq20'])
  output += "  MAPQ >= 30: %s\n" % pct_str(output_stats['mapq30'])
  output += "  MAPQ == %2s: %s\n" % (output_stats['mapq-best'],
    pct_str(output_stats['mapq-top']))
  output += "  strand bias: %s\n" % output_stats['strand_bias']
  output += "  mate bias:   %s\n" % output_stats['mate_bias']
  return output


def mapq_ge_thres(mapqs, thres):
  """Return the number of reads greater than or equal to a given MAPQ."""
  count = 0
  i = thres
  while i < len(mapqs):
    count += mapqs[i]
    i+=1
  return count


def pct(count, total, decimals=1):
  try:
    return round(100*count/total, decimals)
  except ZeroDivisionError:
    return None


def pct_str(percent, decimals=1, const_width=True):
  """Format a percentage into a constant-width, constant precision string."""
  if const_width:
    width = str(decimals+3)
    if percent <= 0:
      return ("%"+width+"d%%") % 0
    elif percent == 100:
      return ("%"+width+"d%%") % 100
    elif percent is None:
      return ("%"+width+"s%%") % "N/A"
    else:
      return ("%"+width+"."+str(decimals)+"f%%") % percent
  else:
    if percent <= 0:
      return "0%"
    elif percent == 100:
      return "100%"
    elif percent is None:
      return "N/A"
    else:
      return str(round(percent, decimals))+"%"


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
