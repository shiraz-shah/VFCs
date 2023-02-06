# VFCs
Code for de novo delieation of viral families in virome data

## Assembly og reads into contigs
### Read QC
 zcat input.fqz | fastq_quality_trimmer -t 13 -l 32 -Q 33 | fastq_quality_filter -p 90 -q 13 -Q 33

