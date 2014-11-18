#!/usr/bin/env python
import os
import sys
import pysam

infile = sys.argv[1]
(dirname, filename) = os.path.split(infile)
outfile = os.path.join(dirname, 'dechim.rlen.'+filename)

def check_chim(read):
  tags = dict(read.tags)
  if 'SA' in tags.keys():
    sa = tags['SA'].split(';')[:-1]
    if len(sa) == 1:
      (chrom, pos, strand, cigar, mapq, nm) = sa[0].split(',')
      if chrom == 'chrM' and (int(pos) >= 16000 or int(pos) <= 600):
        return read
      else:
        pass
    else:
      pass
  else:
    return read

sam = pysam.Samfile(infile, 'rb')
out = pysam.Samfile(outfile, 'wb', template=sam)
for read in sam:
  try:
    read1 = read
    read2 = sam.next()
    if (check_chim(read1) and check_chim(read2) and
        read1.rlen >= 100 and read2.rlen >= 100):
      out.write(read1)
      out.write(read2)
    else:
      pass
  except StopIteration:
    sam.close()
    out.close()

