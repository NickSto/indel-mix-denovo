#TODO: Replace with Python script that can detect the columns from the headers, and also have a
#      better interface for excluded regions.
BEGIN {
  # Command line options (set with -v)
  if (! covg) {
    covg = 1000;
  }
  if (! freq) {
    freq = 0.01;
  }
  if (! strand) {
    strand = 1.0;
  }
  if (! mate) {
    mate = 1.0;
  }
  FS="\t";
  OFS="\t";
}

$6 > covg && $8 > freq && $20 < strand && $21 < mate { print $0; }
