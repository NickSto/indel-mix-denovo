#!/usr/bin/env perl
# Warning: This grew organically from a quick script and is a mess.
# It's kept because it's useful, but needs a total rewrite.
use warnings;
use strict;
use File::Basename;
my $LIMIT_DEFAULT = 1000000;
my $LIMIT_MAJMIN = 2;

my $script_name = basename($0);
my $USAGE = <<USAGE;
USAGE:
  \$ $script_name [options] naive-variant-caller.vcf
  \$ cat naive-variant-caller.vcf | $script_name [options]
For summarizing the output of the Naive Variant Caller (or nvc-filter.py). Does
a simple parsing of an input VCF and prints a summary of the variants in the
sample columns of each line (Note: it combines all samples on each line).
  -r: Print the reference allele from the REF column
  -c: Give read counts for each variant
  -f: Give frequency percentages for each variant
  -a: Order variants alphabetically instead of by frequency, for consistency
  -p: Give more precise frequencies (two decimal points instead of one)
  -i, -d: Print only certain variants: -i for insertions, -d for deletions
  -l: Limit the number of variants printed to this number
  -F: Print the filename as the last column.
  -m: Print a different, tab-delimited format: the chromosome, the location,
      the top variant (major allele), its frequency, the #2 variant (minor
      allele), and its frequency (or more variants, if requested with -l):
        chrM    12417   C       94.94   d1      4.53
USAGE

# homegrown option system because getopt has problems
my $ins = 0;
my $del = 0;
my $refs = 0;
my $freqs = 0;
my $alpha = 0;
my $counts = 0;
my $majmin = 0;
my $precise = 0;
my $filename = 0;
my $limit = $LIMIT_DEFAULT;
my $last_was_limit = 0;
my @new_args = ();
for my $arg (@ARGV) {
  my $opt = 0;
  if ($arg =~ m/^-\w*h\w*$/ || $arg eq '--help') {
    print $USAGE;
    exit(0);
  }
  if ($last_was_limit) {
    $limit = $arg;
    $last_was_limit = 0;
    $opt = 1;
  }
  if ($arg =~ m/^-\w*r\w*$/) {
    $refs = 1; $opt = 1;
  }
  if ($arg =~ m/^-\w*f\w*$/) {
    $freqs = 1; $opt = 1;
  }
  if ($arg =~ m/^-\w*a\w*$/) {
    $alpha = 1; $opt = 1;
  }
  if ($arg =~ m/^-\w*c\w*$/) {
    $counts = 1; $opt = 1;
  }
  if ($arg =~ m/^-\w*p\w*$/) {
    $precise = 1; $opt = 1;
  }
  if ($arg =~ m/^-\w*m\w*$/) {
    $majmin = 1; $opt = 1;
  }
  if ($arg =~ m/^-\w*i\w*$/) {
    $ins = 1; $opt = 1;
  }
  if ($arg =~ m/^-\w*d\w*$/) {
    $del = 1; $opt = 1;
  }
  if ($arg =~ m/^-\w*F\w*$/) {
    $filename = 1; $opt = 1;
  }
  if ($arg =~ m/^-\w*l\w*$/) {
    $last_was_limit = 1;
    $opt = 1;
  }
  unless ($opt) {
    push(@new_args, $arg);
  }
}
if ($majmin) {
  $freqs = 1;
  $precise = 1;
  if ($limit == $LIMIT_DEFAULT) {
    $limit = $LIMIT_MAJMIN;
  }
}
@ARGV = @new_args;

while (<>) {
  chomp();
  my $line = $_;
  next if (substr($line, 0, 1) eq '#');
  my @columns = split("\t", $line);
  next if (@columns < 10);
  my $ref = $columns[3];
  my $coverage = 0;
  my %variants;
  for my $column (@columns[9..$#columns]) {
    my @fields = split(':', $column);
    next if (@fields < 4);
    my @variants = split(',', $fields[3]);
    for my $variant (@variants) {
      if ($variant =~ m/^[+-]?([^=]+)=(\d+)/) {
        $variants{$1} += $2;
        $coverage += $2;
      }
    }
  }
  my $buffer = "";
  # print (to buffer) chrom, coordinate, and maybe ref allele
  if ($majmin) {
    $buffer .= "$columns[0]\t$columns[1]\t";
  } else {
    print "$columns[0] $columns[1]";
    if ($refs) {
      $buffer .= " $ref:\t";
    } else {
      $buffer .= ": \t";
    }
  }
  # print variants
  my @vartypes;
  if ($alpha && ! $majmin) {
    @vartypes = sort(keys %variants);
  } else {
    @vartypes = sort_by_value(%variants);
  }
  my $varcount = 0;
  for my $variant (@vartypes) {
    if ($ins || $del) {
      my $print = 0;
      $print = 1 if ($ins && $variant =~ m/^[GATCN][GATCN]+$/i);
      $print = 1 if ($del && $variant =~ m/^d\d+$/i);
      next unless ($print);
    }
    $buffer .= $variant;
    if ($counts && ! $majmin) {
      $buffer .= ":$variants{$variant}";
    }
    if ($freqs || $majmin) {
      my $percent = 100*$variants{$variant}/$coverage;
      my $format = "%.1f";
      $format = "%.2f" if ($precise);
      if ($majmin) {
        $buffer .= sprintf("\t".$format, $percent);
      } else {
        $buffer .= ":".sprintf($format, "$percent")."%";
      }
    }
    $buffer .= "\t";
    $varcount++;
    if ($varcount >= $limit) {
      last;
    }
  }
  if ($varcount) {
    print $buffer;
    # always print a constant number of columns, even empty ones
    while ($varcount < $limit) {
      print "\t\t";
      $varcount++;
    }
    if ($filename) {
      print "$ARGV\t";
    }
    print "\n";
  }
}

# sorts keys by their associated values, in descending order
sub sort_by_value {
  my %hash = @_;
  my @keys = keys(%hash);
  return
    map  { $_->[0] }
    sort { $b->[1] <=> $a->[1] }
    map  { [ $_, $hash{$_} ] }
    @keys;
}