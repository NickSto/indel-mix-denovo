version: 0569ad9

Description
===========

This might not qualify as a bug, but this is more to document the effects of a known deficiency in the inspect-reads.py.

The deficiency is that it is unable to differentiate between insertions of the same length. The
effect is that it's showing the same stats for different insertions of the same length (duh).
The effect seems to be to count every read that contains either insertion as supporting both.
So the final stats are identical for both insertions.

An example is site 8271 in /nfs/brubeck.bx.psu.edu/scratch5/nick/family/all4/misc/8272align/pipeline/SC8C1-bl.bt2.--rdg5,2--rfg5,2.
As you can see in the NVC output, there are 3 insertions, 2 of which have the same length:

    $ awkt '$2 == 8271 {print $10}' nvc.vcf | cut -d : -f 4 | tr ,= '\n\t' | awkt 'length($1) > 2 {print $2, substr($1,1,1) substr($1,3)}' | sort -k 2
    142	-ACCCCCTCT
    70	+ACCCCCTCT
    1	-ACTCCCTCT
    2	-CCCTCT
    2	+CCCTCT

The output from inspect-reads.py is:

    $ head -n 1 vars_asm.tsv && bioawk -c header '$Coord == 8271' vars_asm.tsv | sort
    #Sample	Chrom	Coord	Type	Alt	Covg	Reads	Freq	Unmap	Improp	Fwd	First	Sndary	Dup	Mapq0	Mapq20	Mapq30	MapqMax	MapqMaxReads	SBias	MBias	Flags	PosDist	Seq
    __NONE__	SC8-ch	8271	I	ACCCCCTCT	3347	244	7.29	0.0	0.0	37.754.9	0.0	0.0	0.0	100.0	92.2	50.8	42	0.0979	0.1317	244,244,0,0,152,92,134,110,0,0,0,0	127,37,3,0,0,0,0,4,32,41	TATTTACCCTATAGCaCCCCCTCTACCCCC
    __NONE__	SC8-ch	8271	I	ACTCCCTCT	3347	244	7.29	0.0	0.0	37.754.9	0.0	0.0	0.0	100.0	92.2	50.8	42	0.0979	0.1317	244,244,0,0,152,92,134,110,0,0,0,0	127,37,3,0,0,0,0,4,32,41	TATTTACCCTATAGCaCCCCCTCTACCCCC
    __NONE__	SC8-ch	8271	I	CCCTCT	3347	4	0.12	0.0	0.0	50.0	50.00.0	0.0	0.0	100.0	50.0	50.0	40	0.6351	0.0652	4,4,0,0,2,2,2,2,0,0,0,0	2,0,0,0,0,0,0,0,0,2	TATTTACCCTATAGCaCCCCCTCTACCCCC

The first two indels are given identical stats.

This issue was acceptable when I was only giving it a select list of indels to inspect. It would
be very unlikely that some would be insertions of the same length at the same location. But now,
I'm giving it a list of all variants in the entire alignment, so collisions are expected.

Solution
========

Unfortunately, the solution can't come from altering inspect-reads.py or even bamslicer.py.
The limitation is in pyBamParser.bam.BAMRead.get_indels(), which only gives the lengths of indels.
And even if I were to alter pyBamParser, it wouldn't be easy. get_indels() just parses the CIGAR
string, which only includes alignment positioning. I'd have to get the read sequence and walk along
it as I parse the CIGAR string at the same time.
