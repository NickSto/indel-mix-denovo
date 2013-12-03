#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;

if (@ARGV && $ARGV[0] =~ m/^-?-h/) {
  my $script_name = basename($0);
  print <<USAGE;
USAGE:
  \$ $script_name naive-variant-caller.vcf
  \$ cat naive-variant-caller.vcf | $script_name
For testing nvc-filter.py. Does a simple parsing of an input VCF and prints a
summary of the variants in the sample columns of each line.
USAGE
  exit(0);
}

while (<>) {
  chomp();
  my $line = $_;
  next if (substr($line, 0, 1) eq '#');
  my @columns = split("\t", $line);
  next if (@columns < 10);
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
  print "$columns[0] $columns[1]:\t";
  my @vartypes = sort(keys %variants);
  print "@vartypes\n";
}