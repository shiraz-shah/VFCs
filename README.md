# De novo discovery of viral families in virome data
Current (as of 2023) workflows for analysis of virome data involve mapping reads to a public virus database, thereby ignoring the vast amounts of viral dark matter that exists within such data sets. Here we provide the code that we used for de novo discovery of new viral families in a virome data set, as outlined in [Shah et al. 2023](https://www.nature.com/articles/s41564-023-01345-7). The code has been modified since so it works better with newer software dependencies and alternative data configurations.

## Assembly of reads into contigs
### Read QC
Read QC for our study was performed as shown below. Your viromes sequences may benefit from a more updated read qc pipeline based on [trimmomatic](https://github.com/usadellab/Trimmomatic) or [FastP](https://github.com/OpenGene/fastp). The vsearch read deduplication step was used here because our viromes were MDA amplified.
```
zcat input.fqz | fastq_quality_trimmer -t 13 -l 32 -Q 33 | fastq_quality_filter -p 90 -q 13 -Q 33 | cutadapt -a CTGTCTCTTATACACATCT -m 32 - | vsearch --derep_prefix /dev/stdin --output /dev/stdout  > filtered.fqc
```
We used cutadapt to remove residual illumina adapters. You may want to use [trimmomatic](https://github.com/usadellab/Trimmomatic) or [FastP](https://github.com/OpenGene/fastp) instead. Residual illumina adapters are a pervasive problem across many different virome extraction protocols due to abnormally small insert sizes. Whether this step is necessary for your data set depends on the quality of your sequences as well as the sequence of the adapters used.

### assembly
Left and right reads (1, and 2), as well as unpaired reads left over from read QC (3) were used as input for assembly with [spades](https://github.com/ablab/spades) as follows:
```
spades.py -1 sample_1.fq.gz -2 sample_2.fq.gz -s sample_3.fq.gz --meta -t 48 -m 200 --only-assembler -o sample.assembly
```
We disabled "read hamming" as reads were already QC'd, and this substantially accellerated assembly speed without compromising its quality.

We have nice switched to [megahit](https://github.com/voutcn/megahit) for assembly:
```
megahit -1 sample_1.fq.gz -2 sample_2.fq.gz -r sample_3.fq.gz -o sample.assembly
```
Using this assembler we obtain very similar results to spades, but with more efficent use of CPU and memory making assemblies faster.

## clustering of contigs into species-level vOTUs
For clustering similar viruses accross samples into species-level clusters, the assemled contigs from all samples were first pooled into a single FASTA file. Then we used [BLAT](https://github.com/djhshih/blat) to do an all-against-all alignment:
```
blat contigs.all.fna contigs.all.fna contigs.all.blat -out=blast8
```
The output from BLAT was used to build ~95% sequence clusters as follows, while avoiding the selection of chimeric assemblies as OTU representatives:
```
cat contigs.all.fna | f2s | seqlengths | joincol <(cat contigs.all.blat | awk '{if ($1 == $2) print $1 "\t" $12}' | hashsums | tail -n +2) > contigs.all.lengths
cat contigs.all.lengths | awk '$3/$2 > 2.15' | cut -f1 > contigs.all.chimeras.list
cut -f1,2,12 contigs.all.blat | hashsums | tail -n +2 | joincol contigs.all.chimeras.list | awk '{if ($NF == 0) print $1 "\t" $2 "\t" $3}' | joincol contigs.all.lengths | joincol contigs.all.lengths 2 | sort -k4,4nr -k1,1 | awk '{if ($3/$NF >= .90) print $1 "\t" $2}' | perl -lane 'unless (exists($clusters{$F[1]})) {$clusters{$F[1]} = $F[0]; print "$F[1]\t$F[0]"}' > OTUs.tsv
```
The information in the resulting output can be used to boil down `contigs.all.fna` into `OTUs.fna` like so:
```
cat contigs.all.fna | f2s | joincol <(cut -f2 OTUs.tsv) | awk '$NF == 1' | cut -f1,2 | s2f > OTUs.fna
```

## decontamination of viral species
Since most virome extractions contain some amount of bacterial contaminating DNA, some of the OTUs from above may represent contaminant species. In our study we decontaminated the OTUs manually by clustering them by encoded protein similarity (as shown below) and examining each cluster for viral signatures.

We do not recommend the manual approach now as tools have since been developed for this task. Instead we recommend using [geNomad](https://portal.nersc.gov/genomad/) followed by [CheckV](https://bio.tools/checkv), because we find that geNomad is sensitive and CheckV adds specificity. We have had mixed results with [VirSorter2](https://github.com/jiarong/VirSorter2). We definitely do not recommend older tools such as DeepVirFinder, or so-called deep-learning-based tools such as PPR-Meta as their performance is sub-par. In the end, a subset the OTUs will be deemed viral, with the rest being likely contaminants. All subsequent steps should be limited to the decontaminated viral subset of `OTUs.fna`, henceforth refered to as `vOTUs.fna`. We have not provided code for subsetting as our method was manual. In our case we had more than 300k OTUs of which only 15k were vOTUs.

## vOTU gene calling, and protein comparison
Calling genes on the vOTUs and submitting the resulting proteins to a sensitive all-against-all sequence search allows for two things:
 - uncovering of deeper evolutionary relationships so taxonomic groups like genera and families are revealed
 - grouping of viral proteins into de novo viral ortholog groups (VOGs) that can be used alongside the vOTUs for statistical analyses against sample metadata.

This is done with [Prodigal](https://github.com/hyattpd/Prodigal) and [fasta36](https://github.com/wrpearson/fasta36) like follows:
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
The distance-matrix can be input into [PHYLIP](https://bio.tools/PHYLIP) or [rapidNJ](https://github.com/somme89/rapidNJ) to construct a neighbour-joining tree, e.g. like follows:
```
rapidnj -i pd vOTUs.mat > vOTUs.nwk
```

## using PhyloTreeLib for cutting the tree to obtain viral genera, subfamilies and families
The above tree can be cut at certain cutoffs to generate viral genera, subfamilies and family-level clusters (VFCs).

After installing [treetool](https://github.com/agormp/treetool) and [PhyloTreeLib](https://github.com/agormp/phylotreelib) you can find appropriate cutoffs for various taxonomic levels, using `treetool`'s `cladeinfo` option. To invoke this fuction one first saves a list of viruses from the tree that one knows belong to the same genus, subfamily or family into a file called e.g. `cladefile.txt`. `treetool --cladeinfo=cladefile.txt` will then return the cutoff to cut the tree and reproduce that taxon.

The `cladeinfo` option can be used to find appropriate cutoffs for all three taxonomic levels above as long as there is a group of viruses within the tree for which the taxonomy is fully resolved. This is why it is useful to spike in your vOTUs.faa file with proteins from additional viruses from [ICTV](ictv.global). In fact, we used the viral family _Herelleviridae_ in order to determine the very cutoffs shown below. At the time, _Herelleviridae_ was the only viral family with a fully resovled viral taxonomy according to ICTVs new criteria. Now there are multiple viral families that are fully resolved on [ICTV](ictv.global) and we recommend running `cladeinfo` on multiple families, subfamilies and genera, and then determining the correct cutoff for each taxonomic level by taking the median.
```
python3 /path/to/treetool.py -I newick --clustcut=0.04 vOTUs.rooted.nwk; mv clusterdir/clusterinfo.txt vOTUs.VFCs.tsv; rmdir clusterdir
python3 /path/to/treetool.py -I newick --clustcut=0.125 vOTUs.rooted.nwk; mv clusterdir/clusterinfo.txt vOTUs.subfamilies.tsv; rmdir clusterdir
python3 /path/to/treetool.py -I newick --clustcut=0.250 vOTUs.rooted.nwk; mv clusterdir/clusterinfo.txt vOTUs.genera.tsv; rmdir clusterdir
```
The cutoffs we used above for our study were based on a manually rooted tree, and can thus be directly interpreted in AAI terms. In practice however, rooting a tree manually is time consuming and clustcut should work fine even if the tree is not properly rooted, as long as you cut the tree using the `cladeinfo` cutoffs generated from your specific tree using valid cladefiles. In your case, if the tree is not properly rooted, the cutoffs you find will not be as close to 0 as ours. They could be something like 0.54, 0.63 and 0.78 for the three taxonomic levels respectively. Their numerical values will not be meaningful outside the scope of your own tree. But they should work correctly as a means to define families, subfamilies and genera based on your tree none the less.

## Estimation of viral relative abundances and generation of OTU table
For calculation of relative abundances we mapped QC'd reads from each sample to assembled contigs from that sample using [BWA](https://github.com/lh3/bwa) followingly:
```
bwa mem -a sample.contigs.fna <(cat sample_?.fqc) | samtools view -F 4 -C -T sample.contigs.fna > sample.cram
```
BWA apparently doesn't support mixing paired and upaired reads (left over from read QC), so all three fq files (left, right and unpaired) were concatenated. Mappings were saved as CRAM using [samtools](https://github.com/samtools/) to save space.

Relative abundances per vOTU were calculated using [msamtools](https://github.com/arumugamlab/msamtools), which normalises for length and iteratively redistributes ambiguous mappings. Providing the total number of reads before mapping to msamtools provides an additional "unknown" relative abundance row that is a good estimation of virome contamination:
```
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample.cram | msamtools profile --label=sample --total=$(cat sample_?.fqc | f2s | tail -n +2 | wc -l) -o sample.profile.gz -
```

Abundance profiles from each sample were joined with vOTUs.tsv and pasted together into a matrix using some convoulted bash code. There may be an easier way to do this in R.
```
(cat samples.list | tr '\n' '\t' | sed 's/\t$/\n/'; paste <(cut -f2 vOTUs.tsv | uniq) <((echo -n "paste"; cat samples.list | while read sample; do echo -n " <(cut -f2 vOTUs.tsv | uniq | joincol <(tail -n +2 $sample.profile.txt | joincol vOTUs.tsv | awk '{print \$3 \"\t\" \$2}') | cut -f2)"; done) | bash)) > samples.OTU.mat
```
Mapping to the sample's own contigs (i.e. local mapping) instead of all vOTUs in vOTUs.fna (global mapping) provides better specificity at the cost of sensitivity, and local mapping is what we did for this study. However global mapping is equally valid and seems to be the norm, but will require some denoising of the resulting OTU-table. Also with global mapping you will not need to translate contig names into vOTU names, so the above code for generating the OTU table will be simpler. Global mapping will however result in 10 times the alpha-diversity per sample than local mapping. Much of that is noise, and different studies have employed different measures for denoising the OTU table resulting from global mapping. Minimum prevalence and abundance cutoffs per OTU are widely used to filter the OTU table. Roux et al. recommend requiring a minimum vOTU read mapping coverage 75% for a vOTU to be considered present in a particular sample. We have been using a combination of coverage and average depth of at least 50% and 1x respectively to consider a vOTU present. Such stats per vOTU per sample can be obtained using msamtools too:
```
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample.cram | msamtools coverage --summary -o sample.coverage.gz -
```
Filtering a globally mapped `sample.profile.gz` with 1x and 50% depth and coverage cutoffs from `sample.coverage.gz` will result in an OTU table that is conservative and will dramatically reduce overall alpha diversity. Local mapping is equally conservative and will definitely get rid of noise, while also certainly ignoring any viruses that were not abundant enough to generate assembled contigs. A global mapping without any denoising has too much noise to be useful in most situations. Technical variables such as virome extraction batch and sequencing lane will completely overshadown any biological signals in the OTU table when no denoising is done on a global mapping. Which approach works for your data may require some benchmarking. We recommend either local mapping, or a global mapping denoised with the above cutoffs.

## merging everything in R inside a phyloseq object
Statistical analyses comparing viral species-counts against sample meta-data are made particularly practical using the [PhyloSeq](https://github.com/joey711/phyloseq) framework. PhyloSeq is designed for bacterial 16S data, where bacterial taxonomy can be used to meaningfully agglomorate the data for added statistical power. We can do the same now for viromics data because we have a rooted tree covering all vOTUs, as well as viral taxonomy at different levels from the PhyloTreeLib step.

A sample data table was prepared containing sample metadata, including sequencing batches, sequencing depths, no. of reads passing QC, no. of reads mapping to vOTUs, number of reads deemed bacterial contaminants as per [ViromeQC](https://github.com/SegataLab/viromeqc) etc. Our sample data table looked something like this:
```
sampleId batch lane     seqDepth   propOTU    viromeQC  
sample1  1     1-1      10876064   0.402      0.099
sample2  1     1-1      22717700   0.164      0.261
sample3  1     1-1      27739264   0.13       0.345
sample4  1     1-1      13398052   0.65       0.1
sample5  1     1-1      14690890   0.432      0.247
sample6  1     1-1      12316230   0.463      0.505
sample7  1     1-1      14045568   0.42       0.313
```

A "taxonomy table" was also prepared outlining which genus, subfamily, VFC and VOC each vOTU belonged to. This table had additional per-vOTU info about vOTU length, the predicted bacterial host for phages (we used [CrisprOpenDB](https://github.com/edzuf/CrisprOpenDB) with a custom database, along with [WiSH](https://github.com/soedinglab/wish)), viral lifestyle, genome completion etc. Our taxonomy table looked like this:
```
           category    class               order          family          subfamily  genus    species    contaminant  virulence  complete  length  nkids  gbkNspacers  mtgNspacers wishPvalue   hostTaxid    hostKingdom  hostPhylum             hostClass                hostOrder                 hostFamily                 hostGenus                hostSpecies 
OTU_9913   satellite   fragment            NA             NA              2739       3473     OTU_9913   1            NA         0         9243    33     0            0           0.0360257    997877       "Bacteria"  "Bacteroidetes"        "Bacteroidia"            "Bacteroidales"           "Bacteroidaceae"           "Bacteroides"             NA          
OTU_9904   satellite   fragment            NA             NA              2739       3473     OTU_9904   1            NA         0         9285    2      0            0           0.0287685    435590       "Bacteria"  "Bacteroidetes"        "Bacteroidia"            "Bacteroidales"           "Bacteroidaceae"           "Bacteroides"             NA          
OTU_4377   contaminant otherClass          NA             NA              5169       9759     OTU_4377   1            NA         NA        26924   26     0            2           0.00739772   816          "Bacteria"  "Bacteroidetes"        "Bacteroidia"            "Bacteroidales"           "Bacteroidaceae"           "Bacteroides"             NA          
OTU_901    virus       Caudoviricetes      Crassvirales   Frejaviridae    2738       3474     OTU_901    0            0          1         55095   1      30           15          0.00690197   816          "Bacteria"  "Bacteroidetes"        "Bacteroidia"            "Bacteroidales"           "Bacteroidaceae"           "Bacteroides"             NA          
OTU_830    virus       Caudoviricetes      Crassvirales   Frejaviridae    2738       3474     OTU_830    0            0          1         56814   4      24           15          0.00652455   816          "Bacteria"  "Bacteroidetes"        "Bacteroidia"            "Bacteroidales"           "Bacteroidaceae"           "Bacteroides"             NA          
OTU_1002   virus       Caudoviricetes      Crassvirales   Frejaviridae    2738       3474     OTU_1002   0            0          1         52037   18     31           14          NA           816          "Bacteria"  "Bacteroidetes"        "Bacteroidia"            "Bacteroidales"           "Bacteroidaceae"           "Bacteroides"             NA          
```

The above two tables along with the tree and the OTU table were merged in R using phyloseq, for statistical analyses against sample metadata:
```
library(phyloseq)
library(ape)
phyloseq(phy_tree(read.tree("all.rooted.nwk")), sample_data(read.table("samples.data.tab", sep = "\t")), otu_table(read.table("samples.OTU.mat", sep = "\t"), taxa_are_rows = T), tax_table(as.matrix(read.table("finalCuration.taxtable.v2names.tab", sep = "\t"))))
```

