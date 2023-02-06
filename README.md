# De novo discovery of viral families in virome data
Current (as of 2023) workflows for analysis of virome data involve mapping reads to a public virus database, thereby ignoring the vast amounts of viral dark matter that exists within such data sets. Here we provide the code that we used for de novo discovery of new viral families in a virome data set, as outlined in https://doi.org/10.1101/2021.07.02.450849

## Assembly of reads into contigs
### Read QC
Read QC was performed as shown below. The vsearch step can be skipped if the virome was unamplified (MDA).
```
zcat input.fqz | fastq_quality_trimmer -t 13 -l 32 -Q 33 | fastq_quality_filter -p 90 -q 13 -Q 33 | cutadapt -a CTGTCTCTTATACACATCT -m 32 - | vsearch --derep_prefix /dev/stdin --output /dev/stdout)  > filtered.fqc
```
We used cutadapt to remove residual illumina adapters. Whether this step is necessary for other data sets depends on the quality of your sequences as well as the sequence of the adapters used.

### assembly
```
spades.py -1 1.fq.gz -2 2.fq.gz -s 3.fq.gz --meta -t 48 -m 200 --only-assembler -o out/$1
```

## clustering of contigs into species-level vOTUs
Assemled contigs from all samples must be pooled into a single FASTA file, making sure that sequence headers are able to distinguish contigs from different samples. Then use BLAT (https://github.com/djhshih/blat) to do an all-against-all alignment:
```
blat contigs.all.fna contigs.all.fna contigs.all.blat -out=blast8
```
The output from BLAT can be used to build ~95% sequence clusters as follows:
```
cut -f1,2,12 contigs.all.blat | hashsums | joincol contigs.all.lengths 2 | sort -k4,4nr -k1,1 | awk '{if ($3/$NF >= .90) print $1 "\t" $2}'
```
The information in the resulting output can be used to boil down `contigs.all.fna` into `vOTUs.fna` like so:
cat contigs.all.fna | f2s | 

## vOTU gene calling, and protein clustering
Calling genes on the vOTUs and submitting the resulting proteins to a sensitive all-against-all sequence search allows for two things:
 - uncovering of deeper evolutionary relationships so taxonomic groups like genera and families are revealed
 - grouping of viral proteins into de novo viral ortholog groups (VOGs) that can be used alongside the vOTUs for further downstream analysis

This is done lile 
