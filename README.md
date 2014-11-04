indel-mix-denovo
================

Purpose
-------

This set of tools forms a pipeline for discovering short indels at non-diploid allele frequencies in short genomes. It was developed for mitochondrial DNA, but this situation also applies to viral sequencing data.

Usage
-----

The input to this pipeline is a set of filtered reads and a *de novo* assembly of those reads. These tools were developed with the output of the [SPAdes](http://bioinf.spbau.ru/spades) assembler (version 2.5.1).

Once you have an assembly, you can use `asm-unifier.py` or `asm-curator.py` to clean it up. Both are reference-based methods, and `asm-curator.py` is a more conservative approach, only discarding seemingly redundant contigs, while `asm-unifier.py` is more aggressive, rearranging and merging contigs. Both require [LASTZ](http://www.bx.psu.edu/~rsharris/lastz/) to align to the reference.

The next step is to align the original reads to your *de novo* assembly. Again, the method used in development was [BWA-MEM](http://bio-bwa.sourceforge.net/). Then, quality filtering of the alignment is recommended. The filtering pipeline this was developed with can be found in [this study](http://dx.doi.org/10.1073/pnas.1409328111).

Once you have a filtered BAM, you must count the allele frequencies with the [Naive Variant Caller](https://toolshed.g2.bx.psu.edu/view/blankenberg/naive_variant_caller) tool, either on [Galaxy](https://usegalaxy.org/) or on the command line (the command is found in the tool source at `tools/naive_variant_caller.py`). The end result will be a VCF file with counts for every allele at every position.

Then filter this VCF with `nvc-filter.py`. Add the arguments `-r S` to remove SNPs and extract indels only. You can also add your coverage and allele frequency thresholds as filters.

Next, use `inspect-reads.py` to further filter the indels, based on read and alignment data from the assembly-aligned BAM. Use the `-t` option to produce simple tab-delimited output, one per variant (per sample, if multiple read groups).

These are your final variants, but they are still in the coordinate system of your *de novo* assembly, not the reference. You can use `quick-liftover.py` to convert them to reference coordinates, though beware the inherent caveat that some assembly coordinates (like those in an insertion) will not exist in the reference and so cannot be converted.