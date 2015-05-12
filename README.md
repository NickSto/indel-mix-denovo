![icon](http://v.nstoler.com/img/hg19-32.png) indel-mix-denovo
================

Purpose
-------

This set of tools forms a pipeline for discovering short indels at non-diploid allele frequencies in short genomes. It was developed for mitochondrial DNA, but it also applies to other systems like viral populations.

Automated Pipeline
------------------

### Requirements

The pipeline itself requires no installation. Simply clone the repository and run the scripts from the top-level directory. The requirements must be installed, though, and must be placed on the `PATH`, except the Python modules and Picard tools. The path to the Picard jars should be set as the `PICARD_DIR` environment variable.

Version numbers in parentheses are what the development environment uses. Most are not hard requirements, but results may differ with different versions.

Some versions which **are** hard requirements are Python 2.7+, bash 4.0+, and PySAM, which must be 0.7.7 exactly, as later versions (specifically 0.8.1) broke APIs that at least one script uses. If you already have 0.8.1 or later installed, you can use [virtualenv](https://virtualenv.pypa.io/en/latest/) to run the pipeline with a different version without disturbing the global installation.

* [bash](https://www.gnu.org/software/bash/bash.html) (**4.0**)
* [Python](https://www.python.org/) (**2.7**)
* [Java](https://www.java.com/) (1.6.0_26)
* [LASTZ](http://www.bx.psu.edu/~rsharris/lastz/) (1.02.00)
* [SPAdes](http://bioinf.spbau.ru/spades) (2.5.1)
* [SAMtools](http://sourceforge.net/projects/samtools/) (0.1.18)
* [BWA-MEM](http://bio-bwa.sourceforge.net/) (0.7.12-r1039)
* [Naive Variant Caller](https://toolshed.g2.bx.psu.edu/view/blankenberg/naive_variant_caller) (7:3c3ef33b9d56)
* [pyBamParser](https://bitbucket.org/dan/pybamparser) (10:e0d3642d2648)
* [pyBamTools](https://bitbucket.org/dan/pybamtools) (21:4a1bcf3d85fa)
* [PySAM](https://pysam.readthedocs.org/en/latest/) (**0.7.7**)
* [Picard](https://broadinstitute.github.io/picard/) (1.100(1571))
* [bamtools](https://github.com/pezmaster31/bamtools) (2.3.0)
* `bamleftalign` from [freebayes](https://github.com/ekg/freebayes) (v9.9.2-27-g5d5b8ac-dirty)

### Usage

Cheatsheet:  
`$ pipeline.py -s sampleid -r refid -l 100 -m 0 reference.fa reads_1.fq reads_2.fq output/dir`  
(100bp paired-end reads, linear (non-circular) genome.)

The automated pipeline script will filter the input reads, create an assembly, align to that assembly, filter that alignment, and call high-quality indels.

See details of the individual arguments by running `./pipeline.py -h`. One option worth mentioning explicitly is `-m`, which should be set to `0` when the genome is linear. The default is to assume it is circular, which affects filtering of reads that appear chimeric because they wrap around the end.

Manual Pipeline
---------------

### Requirements

(See above [explanation](#requirements) about version numbers.)

* [Python](https://www.python.org/) (**2.7**)
* [SPAdes](http://bioinf.spbau.ru/spades) (2.5.1)
* [LASTZ](http://www.bx.psu.edu/~rsharris/lastz/) (1.02.00)
* [BWA-MEM](http://bio-bwa.sourceforge.net/) (0.7.12-r1039)
* [Naive Variant Caller](https://toolshed.g2.bx.psu.edu/view/blankenberg/naive_variant_caller) (7:3c3ef33b9d56)
* [pyBamParser](https://bitbucket.org/dan/pybamparser) (10:e0d3642d2648)
* [pyBamTools](https://bitbucket.org/dan/pybamtools) (21:4a1bcf3d85fa)

### Usage

If you'd like to alter the pipeline steps, you can use the individual tools directly. Following is an outline of the steps; for details of the commands, you can inspect them in the source of `pipeline.py`.

The input to this pipeline is a set of filtered reads and a *de novo* assembly of those reads. These tools were developed with the output of the SPAdes assembler.

Once you have an assembly, you can use `asm-unifier.py` or `asm-curator.py` to clean it up. Both are reference-based methods, and `asm-curator.py` is a more conservative approach, only discarding seemingly redundant contigs, while `asm-unifier.py` is more aggressive, rearranging and merging contigs. Both require LASTZ to align to the reference.

The next step is to align the original reads to your *de novo* assembly. Again, the method used in development was BWA-MEM. Then, quality filtering of the alignment is recommended. The filtering pipeline this was developed with can be found in [this study](http://dx.doi.org/10.1073/pnas.1409328111). You can also see the steps in the automated pipeline scripts, `pipeline.py`, and especially `pre-process-mt.sh`.

Once you have a filtered BAM, you must count the allele frequencies with the Naive Variant Caller tool, either on [Galaxy](https://usegalaxy.org/) or on the command line (the script is found in the tool source at `tools/naive_variant_caller.py`). The end result will be a VCF file with counts for every allele at every position.

Then filter this VCF with `nvc-filter.py`. Add the arguments `-r S` to remove SNPs and extract indels only. You can also add your coverage and allele frequency thresholds as filters.

Next, use `inspect-reads.py` to further filter the indels, based on read and alignment data from the assembly-aligned BAM. Use the `-t` option to produce simple tab-delimited output, one per variant (per sample, if multiple read groups).

These are your final variants, but they are still in the coordinate system of your *de novo* assembly, not the reference. You can use `quick-liftover.py` to convert them to reference coordinates, though beware the inherent caveat that some assembly coordinates (like those in an insertion) will not exist in the reference and so cannot be converted. Use LASTZ to align your assembly to the reference, and then give `quick-liftover.py` the alignment and your indels file.