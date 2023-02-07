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
Left and right reads (1, and 2), as well as unpaired reads left over from read QC (3) were used as input for assembly with [spades](https://github.com/ablab/spades) as follows:
```
spades.py -1 1.fq.gz -2 2.fq.gz -s 3.fq.gz --meta -t 48 -m 200 --only-assembler -o out/$1
```
We disabled "read hamming" as reads were already QC'd, and this substantially accellerated assembly speed without compromising its quality.

## clustering of contigs into species-level vOTUs
For clustering similar viruses accross samples into species-level clusters, the assemled contigs from all samples were first pooled into a single FASTA file. Then we used BLAT (https://github.com/djhshih/blat) to do an all-against-all alignment:
```
blat contigs.all.fna contigs.all.fna contigs.all.blat -out=blast8
```
The output from BLAT was used to build ~95% sequence clusters as follows:
```
cat contigs.all.fna | f2s | seqlengths > contigs.all.lengths
cut -f1,2,12 contigs.all.blat | hashsums | joincol contigs.all.lengths 2 | sort -k4,4nr -k1,1 | awk '{if ($3/$NF >= .90) print $1 "\t" $2}' > vOTUs.tsv
```
The information in the resulting output can be used to boil down `contigs.all.fna` into `vOTUs.fna` like so:
```
cat contigs.all.fna | f2s | joincol <(cut -f2 vOTUs.tsv) | awk '$NF == 1' | cut -f1,2 > s2f > vOTUs.fna
```

## decontamination of viral species
Since most virome extractions contain some amount of bacterial contaminating DNA, some of the vOTUs from above may represent contaminant species. In our study we decontaminated the vOTUs manually by clustering them by encoded protein similarity (as shown below) and examining each cluster for viral signatures.

We do not recommend the manual approach now as tools have since been developed for this task. We recommend [CheckV](https://bio.tools/checkv). We do not recommend older tools such as DeepVirFinder as their performance is sub-par. In the end, a subset the vOTUs will be deemed viral, with the rest being likely contaminants. All subsequent steps should be limited to this decontaminated vOTU subset. We have not provided code for subsetting as our approach was different.

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

## Estimation of viral relative abundances and generation of OTU table
For calculation of relative abundances we mapped QC'd reads from each sample to assembled contigs from that sample using [BWA](https://github.com/lh3/bwa) followingly:
```
bwa mem -t 6 -a -x intractg sample.contigs.fna <(cat sample_?.fqc) | samtools view -F4 -b | samtools sort -m 20G -@ 7 -T sample.bam | samtools view -C -T sample.contigs.fna > sample.cram
```
BWA apparently doesn't support mixing paired and upaired reads (left over from read QC), so all three fq files (left, right and unpaired) were concatenated. Bwa mem's `intractg` option was used to ensure mapping did not cross the species boundary. We're not sure whether this option is needed or even makes sense to use here. Mappings were saved as CRAM using [samtools](https://github.com/samtools/) to save space.

Relative abundances per vOTU were calculated using [msamtools](https://github.com/arumugamlab/msamtools), which normalises for length and interatively redistributes ambiguous mappings:
```
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample.cram | msamtools profile --label={} -o sample.profile.txt -
```

Abundance profiles from each sample was joined with the master vOTU list using some convoulted bash code.
```
(cat samples.list | tr '\n' '\t' | sed 's/\t$/\n/'; paste <(cut -f2 vOTUs.tsv | uniq) <((echo -n "paste"; cat samples.list | while read sample; do echo -n " <(cut -f2 vOTUs.tsv | uniq | joincol <(tail -n +2 $sample.profile.txt | joincol vOTUs.tsv | awk '{print \$3 \"\t\" \$2}') | cut -f2)"; done) | bash)) > samples.OTU.mat
```
There may be an easier way to do this in R.

## merging everything in R inside a phyloseq object
A sample data table was prepared containing sample metadata, including sequencing batches, sequencing depths, no. of reads passing QC, no. of reads mapping to vOTUs, number of reads deemed bacterial contaminants as per [ViromeQC](https://github.com/SegataLab/viromeqc) etc.

A "taxoonmy table" was also prepared outlining which genus, subfamily, VFC and VOC each vOTU belonged to. This table had additional per-vOTU info about vOTU length, the predicted bacterial host for phages (we used [CrisprOpenDB](https://github.com/edzuf/CrisprOpenDB) with a custom database, along with [WiSH](https://github.com/soedinglab/wish)).

The above two tables along with the tree and the OTU table were merged in R using phyloseq, enabling further downstream analyses:
```
library(tidyr)
library(phyloseq)
library(ape)
phyloseq(sample_data(read.table("sample_data.tsv", sep = "\t"), \
 tax_table(read.table("tax_table.tsv") %>% as.matrix, \
 read
```

