#!/usr/bin/python
import os
import sys
import pysam


def get_nm(read):
    tags=dict(read.tags)
    return float(tags['NM'])


infile = sys.argv[1]
(dirname, filename) = os.path.split(infile)
outfile = os.path.join(dirname, 'nm-ratio.'+filename)

sam=pysam.Samfile(infile,'rb')
out=pysam.Samfile(outfile,'wb',template=sam)


for read in sam:
    try:
        read1=read
        read2=sam.next()
        if get_nm(read1)<=(0.02*read1.rlen) and get_nm(read2)<=(0.02*read2.rlen):
            out.write(read1)
            out.write(read2)
        else:
            pass

    except StopIteration:
            sam.close()
            out.close()
