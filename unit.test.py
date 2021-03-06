#!/usr/bin/env python
#TODO: incorporate improvements from version in library repo
import bamslicer
import optparse
import sys
import os

def main():

  FUNCTIONS = {
    'bamslicer.get_reads_and_stats':bamslicer_get_reads_and_stats,
    'inspect-reads.variants_from_vcf':inspect_reads_variants_from_vcf,
    'pipeline-family.choose_asm':pipeline_family_choose_asm,
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
  parser.add_option('-r', '--reports-dir', dest='reports_dir')

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
    for sample in reads:
      print variant['chrom'], variant['coord'], variant['type']
      print "   flags:", stats[sample]['flags']
      sys.stdout.write("   mapqs:")
      for (quality, total) in enumerate(stats[sample]['mapqs']):
        if total > 0:
          sys.stdout.write(' %s=%d' % (quality, total))
      print
      print 'pos-dist: '+','.join(map(str,stats[sample]['var_pos_dist']))

      if not options.no_reads:
        sys.stdout.write("   reads:")
        for read in reads[sample]:
          sys.stdout.write(" "+str(read.get_read_name().split(':')[-1]))
        print


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


def pipeline_family_choose_asm(options):
  pipeline_family = __import__("pipeline-family")
  # Find, read, and parse reports.
  # Then print their parsed representation and store in a dict.
  reports = {}
  for filename in os.listdir(options.reports_dir):
    if not filename.endswith('.tsv'):
      continue
    report_path = os.path.join(options.reports_dir, filename)
    if not os.path.isfile(report_path):
      continue
    with open(report_path) as report_file:
      report = pipeline_family.parse_report(report_file.read())
    print filename
    for key, value in report.items():
      print "\t{}\t{}".format(key, value)
    reports[filename] = report
  max_contigs = pipeline_family.OPT_DEFAULTS['max_contigs']
  include_dup = False
  best_asm = pipeline_family.choose_asm(reports, max_contigs, include_dup)
  print 'chose: '+best_asm


##### UTILITY FUNCTIONS #####


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
