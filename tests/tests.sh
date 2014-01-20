#!/usr/bin/env bash
dirname=$(dirname $0)
set -ue

# unit tests
echo -e "\tbamslicer.get_reads_and_stats ::: cigar-tests.bam:"
$dirname/../unit.test.py bamslicer.get_reads_and_stats -b $dirname/cigar-tests.bam -v chrM:5-I,chrM:199-D,chrM:199-I,chrM:3106-D,chrM:6110-D,chrM:16568-I | diff -s - $dirname/cigar-tests.bamslicer.out
echo -e "\tinspect-reads.variants_from_vcf ::: R19S5-head.vcf.in:"
$dirname/../unit.test.py inspect-reads.variants_from_vcf -V $dirname/R19S5-head.vcf.in | diff -s - $dirname/R19S5-head.var_from_vcf.out 

# functional tests
echo -e "\tasm-curator.py ::: R19S11.lav/R19S11.fa:"
$dirname/../asm-curator.py -l $dirname/R19S11.lav -a $dirname/R19S11.fa -o $dirname/tmp-R19S11.fa.test
diff -s $dirname/tmp-R19S11.fa.test $dirname/R19S11.fa.out
rm $dirname/tmp-R19S11.fa.test
echo -e "\tasm-curator.py ::: R20S11.lav/R20S11.fa:"
$dirname/../asm-curator.py -l $dirname/R20S11.lav -a $dirname/R20S11.fa -o $dirname/tmp-R20S11.fa.test
diff -s $dirname/tmp-R20S11.fa.test $dirname/R20S11.fa.out
rm $dirname/tmp-R20S11.fa.test