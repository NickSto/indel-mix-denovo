#!/usr/bin/env bash
dirname=$(dirname $0)
set -ue
echo -e "\tbamslicer.get_reads_and_stats ::: cigar-tests.bam:"
$dirname/../unit.test.py bamslicer.get_reads_and_stats $dirname/cigar-tests.bam -v chrM:5-I,chrM:199-D,chrM:199-I,chrM:3106-D,chrM:6110-D,chrM:16568-I | diff -s - $dirname/cigar-tests.bamslicer.out
