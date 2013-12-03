#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;

my $script_name = basename($0);
my $USAGE = <<USAGE;
USAGE:
  \$ $script_name [options] naive-variant-caller.vcf
  \$ cat naive-variant-caller.vcf | $script_name [options]
For testing nvc-filter.py. Does a simple parsing of an input VCF and prints a
summary of the variants in the sample columns of each line.
  -r: Print the reference allele from the REF column
  -c: Give read counts for each variant
  -f: Give frequency percentages for each variant
  -p: Give more precise frequencies (two decimal points instead of one)
USAGE

my $refs = 0;
my $freqs = 0;
my $counts = 0;
my $precise = 0;
my @new_args = ();
for my $arg (@ARGV) {
  my $opt = 0;
  if ($arg =~ m/^-\w*h\w*$/ || $arg eq '--help') {
    print $USAGE;
    exit(0);
  }
  if ($arg =~ m/^-\w*r\w*$/) {
    $refs = 1;
    $opt = 1;
  }
  if ($arg =~ m/^-\w*f\w*$/) {
    $freqs = 1;
    $opt = 1;
  }
  if ($arg =~ m/^-\w*c\w*$/) {
    $counts = 1;
    $opt = 1;
  }
  if ($arg =~ m/^-\w*p\w*$/) {
    $precise = 1;
    $opt = 1;
  }
  unless ($opt) {
    push(@new_args, $arg);
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
  print "$columns[0] $columns[1]";
  if ($refs) {
    print " $ref:\t";
  } else {
    print ": \t";
  }
  my @vartypes = sort(keys %variants);
  for my $variant (@vartypes) {
    print $variant;
    if ($counts) {
      print ":$variants{$variant}";
    }
    if ($freqs) {
      my $percent = 100*$variants{$variant}/$coverage;
      my $format = "%.1f";
      $format = "%.2f" if ($precise);
      print ":".sprintf($format, "$percent")."%";
    }
    print "\t";
  }
  print "\n";
}