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
Finds all the types of variants output in the sample columns of the Naive
Variant Caller and prints them, sorted by frequency of occurrence (not read
count).
USAGE
  exit(0);
}

my %types;
while (<>) {
  chomp();
  my $line = $_;
  next if (substr($line, 0, 1) eq '#');
  my @fields = split("\t", $line);
  next if (@fields < 10);
  my @subfields = split(':', $fields[9]);
  next if (@subfields < 4);
  my @variants = split(',', $subfields[3]);
  for my $variant (@variants) {
    if ($variant =~ m/^[+-]?([^=]+)=(\d+)/) {
      $types{$1} += $2;
    }
  }
}

# sort patterns by frequency of occurrence
my @type_list = keys(%types);
my @types_sorted =
  map  { $_->[0] }
  sort { $b->[1] <=> $a->[1] }
  map  { [ $_, $types{$_} ] }
  @type_list;

for my $type (@types_sorted) {
  print "$types{$type}\t$type\n";
}