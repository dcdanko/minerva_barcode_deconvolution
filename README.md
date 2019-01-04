# Minerva Barcoded Read Deconvolution

[![CircleCI](https://circleci.com/gh/dcdanko/minerva_barcode_deconvolution.svg?style=svg)](https://circleci.com/gh/dcdanko/minerva_barcode_deconvolution)

[![CodeFactor](https://www.codefactor.io/repository/github/dcdanko/minerva_barcode_deconvolution/badge)](https://www.codefactor.io/repository/github/dcdanko/minerva_barcode_deconvolution)

Emerging linked-read technologies (aka Read-Cloud or barcoded short-reads) have revived interest in short-read technology as a viable way to understand large-scale structure in genomes and metagenomes. Linked-read technologies, such as the 10x Chromium system, use a microfluidic system and a specialized set of barcodes to tag short DNA reads sourced from the same long fragment of DNA. Subsequently, the tagged reads are sequenced on standard short read platforms.

This approach results in interesting compromises. Each long fragment of DNA is only sparsely covered by reads, no information about the ordering of reads from the same fragment is preserved, and barcodes match reads from roughly 2-20 long fragments of DNA. However, compared to long read technologies the cost per base to sequence is far lower, far less input DNA is required, and the per base error rate is that of Illumina short-reads.

In the accompanying paper, we formally describe a particular algorithmic issue for linked-read technology: the deconvolution of reads with a single barcode into clusters that represent single long fragments of DNA. We also present Minerva, an algorithm which approximately solves the barcode deconvolution problem for metagenomic data. This codebase implements Minerva.

[Minerva: An Alignment and Reference Free Approach to Deconvole Linked-Reads for Metagenomics](https://genome.cshlp.org/content/early/2018/12/06/gr.235499.118.full.pdf+html)

## Installation

From PyPi
```
pip install minerva_deconvolve
```

From source
```
git clone <url>   
cd minerva_barcode_deconvolution
python setup.py install
```

## Deconvolving Reads

Use the following command to run barcode deconvolution. `<fastq>` should be an interleaved fastq file where reads have a `BX` tag designating barcode (this is the default output of [longranger basic](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/other-pipelines))
```
cat <fastq> | minerva_deconvolve -k 20 -w 40 -d 8 -a 20 --remove-stopwords > ebc_assignments.tsv
```

For more options
```
minerva_deconvolve --help
```

### Output

Minerva assigns barcoded reads to clusters within each barcode called deconvolved barcodes. The `minerva_deconvolve` command outputs a tsv file with three columns: read id, barcode, and cluster id. The deconvolved barcode for a read is a tuple of (barcode, cluster id).

```
$ head <minerva_output_file>
BX:Z:GTGCCTTAGTCCGTAT-1 D00547:847:HYHNTBCXX:1:1207:20627:25951 0
BX:Z:GTGCCTTAGTCCGTAT-1 D00547:847:HYHNTBCXX:1:1113:11082:83578 0
BX:Z:GTGCCTTAGTCCGTAT-1 D00547:847:HYHNTBCXX:1:2206:4393:100450 0
BX:Z:GTGCCTTAGTCCGTAT-1 D00547:847:HYHNTBCXX:2:1111:2014:28730  1
BX:Z:GTGCCTTAGTCCGTAT-1 D00547:847:HYHNTBCXX:1:2216:17277:16384 2
BX:Z:GTGCCTTAGTCCGTAT-1 D00547:847:HYHNTBCXX:1:1201:19163:82220 2
BX:Z:GTGCCTTAGTCCGTAT-1 D00547:847:HYHNTBCXX:1:1202:16780:78102 0
BX:Z:GTGCCTTAGTCCGTAT-1 D00547:847:HYHNTBCXX:1:1210:7460:13722  2
```

### Performance

This is a demonstration program and is not intended to be performant. Runtimes over 10 hours are common even on small datasets.
RAM usage is typically 50-100Gb.

## Datasets

The datasets used in the paper may be downloaded from AWS.
 - [Dataset 1](https://s3.us-east-2.amazonaws.com/minerva-datasets/10M.data1_atgctgaaq.fq.gz)
 - [Dataset 2](https://s3.us-east-2.amazonaws.com/minerva-datasets/10M.data2_accctcct.fq.gz)


## Credits

This algorithm was devloped and tested with help from Dmitrii Meleshko, Daniela Bezdan, Chris Mason, and Iman Hajirasouliha.

This package is written and maintained by [David C. Danko](mailto:dcdanko@gmail.com)