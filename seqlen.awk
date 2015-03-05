#!/usr/bin/awk -f
# Usage: $ awk -f seqlen.awk [-v target=chr_id] seq.fa
# Print the lengths of each sequence in the FASTA file.
# Output format: one sequence per line, two tab-delimited columns:
# (1) sequence id, (2) sequence length.
# If target is specified with -v, only chromosomes with an id matching the
# given one will be output.
# Originally from https://stackoverflow.com/a/23992773/726773

BEGIN {
  OFS = "\t";
}

substr($0, 1, 1) == ">" {
  if (seqlen) {
    if (! target || id == target) {
      print id, seqlen;
    }
  }
  split($0, fields);
  id = substr(fields[1], 2);
  seqlen = 0;
  next;
}

{
  seqlen = seqlen + length($0);
}

END {
  if (! target || id == target) {
    print id, seqlen;
  }
}