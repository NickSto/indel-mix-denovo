#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;
use Getopt::Std;

my $ARG_TAKING_OPTS = 'nr';
my $script_name = basename($0);
my $USAGE = <<USAGE;
USAGE:
  \$ $script_name [options] reads.sam
Prints each unique CIGAR "pattern", from most to least common, with totals and
example full CIGAR strings. A "pattern" is a CIGAR string with all digits
removed, leaving the sequence of operations without regard to their lengths.
  e.g. "4M2I245M" has the pattern "MIM"
Output format is tab delimited:
  Columns:  OCCURRENCES  PATTERN  EXAMPLE
  e.g.:     10946        MIM      4M2I245M
  -s: Print a full SAM file as output instead. Instead of printing the normal 3
      columns for each pattern, instead it will print a read with that pattern
      (still one per unique pattern, and in order of abundance).
  -H: Don't print SAM headers if using SAM output (ignored otherwise).
  -n: The maximum number of patterns (or reads) to print
  -r: Which read to pick as the example for the pattern. The default is 1,
      meaning it will pick the first read it encounters with that pattern. 2
      would make it pick the 2nd (or the 1st if there is only one), and so on.
USAGE

my $opt_arg = 0;
my $infilehandle;
for my $arg (@ARGV) {
  if ($infilehandle) {
    print $USAGE;
    print STDERR "Error: options must come before the filename.\n";
    exit(1);
  }
  if ($arg eq '-' || substr($arg, 0, 1) ne '-') {
    if (! $opt_arg) {
      open($infilehandle, "<$arg") or
        die "Error opening file $arg: $!";
    }
  }
  if ($arg =~ m/^-[$ARG_TAKING_OPTS]$/) {
    $opt_arg = 1;
  } else {
    $opt_arg = 0;
  }
}

my $skip_header = 0;
my $which_read  = 1;
my $maximum     = 0;
my $help        = 0;
my $sam         = 0;

my %opt;
getopts('n:r:shH', \%opt);
$skip_header = $opt{H} if (defined($opt{H}));
$which_read  = $opt{r} if (defined($opt{r}));
$maximum     = $opt{n} if (defined($opt{n}));
$help        = $opt{h} if (defined($opt{h}));
$sam         = $opt{s} if (defined($opt{s}));

# print "-s: $sam\n";
# print "-h: $help\n";
# print "-n: $maximum\n";
# print "-r: $which_read\n";
# print "-H: $skip_header\n";
# exit();

if ($help) {
  print $USAGE;
  exit(0);
}

if (! $infilehandle) {
  print $USAGE;
  print STDERR "Error: must provide a filename or \"-\" for stdin.\n";
  exit(1);
}


my %reads = ();
my %cigars = ();
my %patterns = ();
my $header = "";
while (<$infilehandle>) {
  my $line = $_;
  if (substr($line, 0, 1) eq '@' && $line =~ m/^@\w\w\t/) {
    $header .= $line;
    next;
  }
  my @fields = split("\t", $line);
  my $cigar = $fields[5];
  my $pattern = $cigar;
  $pattern =~ s/\d//g;
  if (! $patterns{$pattern} || $patterns{$pattern} < $which_read) {
    if ($sam) {
      $reads{$pattern} = $line;
    } else {
      $cigars{$pattern} = $cigar;
    }
    # print "$cigar:\t$pattern\n";
  }
  $patterns{$pattern}++;
}

# sort patterns by frequency of occurrence
my @pattern_list = keys(%patterns);
my @patterns_sorted =
  map  { $_->[0] }
  sort { $b->[1] <=> $a->[1] }
  map  { [ $_, $patterns{$_} ] }
  @pattern_list;



if ($sam && ! $skip_header) {
  print $header;
}
my $printed = 0;
for my $pattern (@patterns_sorted) {
  $printed++;
  if ($maximum && $printed > $maximum) {
    last;
  }
  if ($sam) {
    print $reads{$pattern};
  } else {
    print $patterns{$pattern}."\t$pattern\t".$cigars{$pattern}."\n";
  }
}
