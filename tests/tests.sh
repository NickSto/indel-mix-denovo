#!/usr/bin/env bash
dirname=$(dirname $0)


# unit tests
BAMSLICER_VARS='chrM:5-I:AG,chrM:199-D:1,chrM:199-I,chrM:310-I:A,chrM:1900-D:2,chrM:3106-D,chrM:6110-D,chrM:16568-I'
echo -e "\tbamslicer.get_reads_and_stats ::: cigar-tests.bam:"
$dirname/../unit.test.py bamslicer.get_reads_and_stats -b $dirname/cigar-tests.bam -v $BAMSLICER_VARS | diff -s - $dirname/cigar-tests.bamslicer.out

echo -e "\tinspect-reads.variants_from_vcf ::: R19S5-head.vcf.in:"
$dirname/../unit.test.py inspect-reads.variants_from_vcf -V $dirname/R19S5-head.vcf.in | diff -s - $dirname/R19S5-head.var_from_vcf.out 


# functional tests
echo -e "\tasm-curator.py ::: R19S11.lav/R19S11.fa:"
$dirname/../asm-curator.py -l $dirname/R19S11.lav -a $dirname/R19S11.fa -o $dirname/tmp-R19S11.fa.test >/dev/null
diff -s $dirname/tmp-R19S11.fa.test $dirname/R19S11.fa.out
rm $dirname/tmp-R19S11.fa.test

echo -e "\tasm-curator.py ::: R20S11.lav/R20S11.fa:"
$dirname/../asm-curator.py -l $dirname/R20S11.lav -a $dirname/R20S11.fa -o $dirname/tmp-R20S11.fa.test >/dev/null
diff -s $dirname/tmp-R20S11.fa.test $dirname/R20S11.fa.out
rm $dirname/tmp-R20S11.fa.test

echo -e "\tclean-merge.py ::: R19S11/R20S1/R31S4/R33S10:"
for sample in R19S11 R20S1 R31S4 R33S10; do
  $dirname/../clean-merge.py $dirname/clean-merge/$sample.lav $dirname/clean-merge/$sample-merge.fa $dirname/clean-merge/$sample-asm.fa > $dirname/clean-merge/$sample.tmp.fa
  if [[ $(crc32 $dirname/clean-merge/$sample.tmp.fa) != $(crc32 $dirname/clean-merge/$sample.fa.out) ]]; then
    echo "FAILURE! Files $dirname/clean-merge/$sample.fa.out and $dirname/clean-merge/$sample.tmp.fa are different"
  else
    echo "Files $dirname/clean-merge/$sample.fa.out and $dirname/clean-merge/$sample.tmp.fa are identical"
  fi
  rm $dirname/clean-merge/$sample.tmp.fa
done

echo -e "\tasm-unifier.py ::: R33S10.fa ::: --orient:"
$dirname/../asm-unifier.py -O $dirname/chrM-rCRS.fa $dirname/asm-unifier/R33S10.fa | diff -s - $dirname/asm-unifier/R33S10-oriented.fa.out
echo -e "\tasm-unifier.py ::: R19S6.fa ::: --orient:"
$dirname/../asm-unifier.py -O $dirname/chrM-rCRS.fa $dirname/asm-unifier/R19S6.fa | diff -s - $dirname/asm-unifier/R19S6-oriented.fa.out

echo -e "\tgroup-filter.py ::: R20S9.tsv/R24S8.tsv/R35S11.tsv/R35S2.tsv:"
$dirname/../group-filter.py -e test -s 1 -m 1 $dirname/tsv-vars/R20S9.tsv $dirname/tsv-vars/R24S8.tsv $dirname/tsv-vars/R35S11.tsv $dirname/tsv-vars/R35S2.tsv
diff -s $dirname/tsv-vars/R20S9-test.tsv $dirname/tsv-vars/R20S9.tsv.out
diff -s $dirname/tsv-vars/R24S8-test.tsv $dirname/tsv-vars/R24S8.tsv.out
diff -s $dirname/tsv-vars/R35S11-test.tsv $dirname/tsv-vars/R35S11.tsv.out
diff -s $dirname/tsv-vars/R35S2-test.tsv $dirname/tsv-vars/R35S2.tsv.out
[[ -f $dirname/tsv-vars/R20S9-test.tsv ]] && rm $dirname/tsv-vars/R20S9-test.tsv
[[ -f $dirname/tsv-vars/R24S8-test.tsv ]] && rm $dirname/tsv-vars/R24S8-test.tsv
[[ -f $dirname/tsv-vars/R35S11-test.tsv ]] && rm $dirname/tsv-vars/R35S11-test.tsv
[[ -f $dirname/tsv-vars/R35S2-test.tsv ]] && rm $dirname/tsv-vars/R35S2-test.tsv
