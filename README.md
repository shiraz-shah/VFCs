# De novo discovery of viral families in virome data
Current (as of 2023) workflows for analysis of virome data involve mapping reads to a public virus database, thereby ignoring the vast amounts of viral dark matter that exists within such data sets. Here we provide the code that we used for de novo discovery of new viral families in a virome data set, as outlined in [Shah et al. 2021](https://doi.org/10.1101/2021.07.02.450849).

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
cat contigs.all.fna | f2s | seqlengths > contigs.all.lengths
cut -f1,2,12 contigs.all.blat | hashsums | joincol contigs.all.lengths 2 | sort -k4,4nr -k1,1 | awk '{if ($3/$NF >= .90) print $1 "\t" $2}' > vOTUs.tsv
```
The information in the resulting output can be used to boil down `contigs.all.fna` into `vOTUs.fna` like so:
```
cat contigs.all.fna | f2s | joincol <(cut -f2 vOTUs.tsv) | awk '$NF == 1' | cut -f1,2 > s2f > vOTUs.fna
```

## vOTU gene calling, and protein comparison
Calling genes on the vOTUs and submitting the resulting proteins to a sensitive all-against-all sequence search allows for two things:
 - uncovering of deeper evolutionary relationships so taxonomic groups like genera and families are revealed
 - grouping of viral proteins into de novo viral ortholog groups (VOGs) that can be used alongside the vOTUs for further downstream analysis

This is done with prodigal (https://github.com/hyattpd/Prodigal) and fasta36 (https://github.com/wrpearson/fasta36) like follows:
```
cat vOTUs.fna | prodigal -a vOTUs.faa -p meta > vOTUs.gbk
fasta36 vOTUs.faa vOTUs.faa -m 8 > vOTUs.fasta36
```

## Defining VOGs using protein comparison result:
Applying an orthology-filter (based on length, coverage and alignment positioning cutoffs as first described in [Shah et al. 2018](https://doi.org/10.1101/262675)) to the above comparison and submitting that to Markov Clustering (https://github.com/micans/mcl) enables clustering of viral proteins into VOGs:
```
cat vOTUs.faa | f2s | seqlengths > vOTUs.faa.lengths
cat vOTUs.fasta36 | joincol vOTUs.faa.lengths | joincol vOTUs.faa.lengths 2 | awk '{print $1 "\t" $2 "\t" $11 "\t" $13/$14 "\t" ($8-$7)/(2*$13)+($10-$9)/(2*$14) "\t" ($7+$8-$9-$10)/($13+$14)}' | awk '{if ($3 <= 0.05) print}' | awk '{if ($5 >= 0.4) print}' | awk '{if (sqrt(($4-1)^2) - (sqrt(sqrt($5))-.8) + sqrt($6^2) <= 0.1) print $1 "\t" $2}' | mcl - -o - --abc | awk '{j++; for (i = 1; i <= NF; i++) {print $i "\t" j}}' > vOTUs.VOGs.tsv
```

## building an Aggregate Protein Similarity (APS) tree for taxonomic delineation:
The all-against-all protein comparison can also be used to construct a distance-matrix between the vOTUs, which can itself be used to construct a tree:
```
awk '{if ($11 <= 0.05) print $1 "\t" $2 "\t" $12}' vOTUs.fasta36 | rev | sed 's/\t[[:digit:]]\+_/\t/' | rev | sed 's/_[[:digit:]]\+\t/\t/' | sort | hashsums | tree_bray > vOTUs.mat
```
The distance-matrix can be input into phylip (https://bio.tools/PHYLIP) or rapidnj (https://github.com/somme89/rapidNJ) to construct a neighbour-joining tree, e.g. like follows:
```
rapidnj -i pd vOTUs.mat > vOTUs.nwk
```

## using phylotreelib for cutting the tree to obtain viral genera, subfamilies and families
After manually rooting the tree (using FigTree), we used the following cutoffs with  to obtain viral genera, subfamilies and family-level clusters (VFCs) and order-level clusters (VOCs):
```
treetool -I newick --clustcut=0.025 vOTUs.rooted.nwk > vOTUs.VOCs.tsv
treetool -I newick --clustcut=0.04 vOTUs.rooted.nwk > vOTUs.VFCs.tsv
treetool -I newick --clustcut=0.125 vOTUs.rooted.nwk > vOTUs.subfamilies.tsv
treetool -I newick --clustcut=0.250 vOTUs.rooted.nwk > vOTUs.genera.tsv
```
The above commands require installation of https://github.com/agormp/treetool and https://github.com/agormp/phylotreelib
