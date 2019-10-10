# RSeQC: An RNA-seq Quality Control Package


## infer_experiment.py
This program is used to “guess” how RNA-seq sequencing were configured, particulary how reads were stranded for strand-specific RNA-seq data, through comparing the “strandness of reads” with the “standness of transcripts”.

The “strandness of reads” is determiend from alignment, and the “standness of transcripts” is determined from annotation.

For non strand-specific RNA-seq data, “strandness of reads” and “standness of transcripts” are independent.

For strand-specific RNA-seq data, “strandness of reads” is largely determined by “standness of transcripts”. See below 3 examples for details.

You don’t need to know the RNA sequencing protocol before mapping your reads to the reference genome. Mapping your RNA-seq reads as if they were non-strand specific, this script can “guess” how RNA-seq reads were stranded.

### For pair-end RNA-seq, there are two different ways to strand reads (such as Illumina ScriptSeq protocol):

* 1++,1–,2+-,2-+

read1 mapped to ‘+’ strand indicates parental gene on ‘+’ strand

read1 mapped to ‘-‘ strand indicates parental gene on ‘-‘ strand

read2 mapped to ‘+’ strand indicates parental gene on ‘-‘ strand

read2 mapped to ‘-‘ strand indicates parental gene on ‘+’ strand

* 1+-,1-+,2++,2–

read1 mapped to ‘+’ strand indicates parental gene on ‘-‘ strand

read1 mapped to ‘-‘ strand indicates parental gene on ‘+’ strand

read2 mapped to ‘+’ strand indicates parental gene on ‘+’ strand

read2 mapped to ‘-‘ strand indicates parental gene on ‘-‘ strand

### For single-end RNA-seq, there are also two different ways to strand reads:

* ++,–

read mapped to ‘+’ strand indicates parental gene on ‘+’ strand

read mapped to ‘-‘ strand indicates parental gene on ‘-‘ strand

* +-,-+

read mapped to ‘+’ strand indicates parental gene on ‘-‘ strand

read mapped to ‘-‘ strand indicates parental gene on ‘+’ strand

* Options:  

--version  
show program’s version number and exit

-h, --help  
show this help message and exit

-i INPUT_FILE, --input-file=INPUT_FILE  
Input alignment file in SAM or BAM format

-r REFGENE_BED, --refgene=REFGENE_BED  
Reference gene model in bed fomat.

-s SAMPLE_SIZE, --sample-size=SAMPLE_SIZE  
Number of reads sampled from SAM/BAM file. default=200000

-q MAP_QUAL, --mapq=MAP_QUAL  
Minimum mapping quality (phred scaled) for an alignment to be considered as “uniquely mapped”. default=30

### Example 1: Pair-end non strand specific:
```sh
infer_experiment.py -r hg19.refseq.bed12 -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam
```

* Output::
```sh
This is PairEnd Data
Fraction of reads failed to determine: 0.0172
Fraction of reads explained by "1++,1--,2+-,2-+": 0.4903
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4925
```

Interpretation: 1.72% of total reads were mapped to genome regions that we cannot determine the “standness of transcripts” (such as regions that having both strands transcribed). For the remaining 98.28% (1 - 0.0172 = 0.9828) of reads, half can be explained by “1++,1–,2+-,2-+”, while the other half can be explained by “1+-,1-+,2++,2–”. We conclude that this is NOT a strand specific dataset because “strandness of reads” was independent of “standness of transcripts”

### Example 2: Pair-end strand specific:
```sh
infer_experiment.py -r hg19.refseq.bed12 -i Pairend_StrandSpecific_51mer_Human_hg19.bam
```
* Output::
```sh
This is PairEnd Data
Fraction of reads failed to determine: 0.0072
Fraction of reads explained by "1++,1--,2+-,2-+": 0.9441
Fraction of reads explained by "1+-,1-+,2++,2--": 0.0487
```

Interpretation: 0.72% of total reads were mapped to genome regions that we cannot determine the “standness of transcripts” (such as regions that having both strands transcribed). For the remaining 99.28% (1 - 0.0072 = 0.9928) of reads, the vast majority was explained by “1++,1–,2+-,2-+”, suggesting a strand-specific dataset.

### Example 3: Single-end strand specific:
```sh
infer_experiment.py -r hg19.refseq.bed12 -i SingleEnd_StrandSpecific_36mer_Human_hg19.bam
```
* Output::

```sh
This is SingleEnd Data
Fraction of reads failed to determine: 0.0170
Fraction of reads explained by "++,--": 0.9669
Fraction of reads explained by "+-,-+": 0.0161
```
Interpretation: This is single-end, strand specific RNA-seq data. Strandness of reads are concordant with strandness of reference gene.


## read_distribution.py

Provided a BAM/SAM file and reference gene model, this module will calculate how mapped reads were distributed over genome feature (like CDS exon, 5’UTR exon, 3’ UTR exon, Intron, Intergenic regions). When genome features are overlapped (e.g. a region could be annotated as both exon and intron by two different transcripts) , they are prioritize as: CDS exons > UTR exons > Introns > Intergenic regions, for example, if a read was mapped to both CDS exon and intron, it will be assigned to CDS exons.

“Total Reads”: This does NOT include those QC fail,duplicate and non-primary hit reads

“Total Tags”: reads spliced once will be counted as 2 tags, reads spliced twice will be counted as 3 tags, etc. And because of this, “Total Tags” >= “Total Reads”

“Total Assigned Tags”: number of tags that can be unambiguously assigned the 10 groups (see below table).

Tags assigned to “TSS_up_1kb” were also assigned to “TSS_up_5kb” and “TSS_up_10kb”, tags assigned to “TSS_up_5kb” were also assigned to “TSS_up_10kb”. Therefore, “Total Assigned Tags” = CDS_Exons + 5’UTR_Exons + 3’UTR_Exons + Introns + TSS_up_10kb + TES_down_10kb.

When assign tags to genome features, each tag is represented by its middle point.

RSeQC cannot assign those reads that:

hit to intergenic regions that beyond region starting from TSS upstream 10Kb to TES downstream 10Kb.

hit to regions covered by both 5’UTR and 3’ UTR. This is possible when two head-to-tail transcripts are overlapped in UTR regions.

hit to regions covered by both TSS upstream 10Kb and TES downstream 10Kb.


* Options:
> --version  
show program’s version number and exit

> -h, --help  
show this help message and exit

> -i INPUT_FILE, --input-file=INPUT_FILE  
Alignment file in BAM or SAM format.

> -r REF_GENE_MODEL, --refgene=REF_GENE_MODEL  
Reference gene model in bed format.

Example:
```sh
read_distribution.py  -i Pairend_StrandSpecific_51mer_Human_hg19.bam -r hg19.refseq.bed12
```
* Output:

```sh
#!/bin/bash
export PATH=~/miniconda3/bin/:$PATH
name=(RC15  RC3  RW15  RW3 TC15  TC3  TW15  TW3)
gtf=/Share/home/tiangeng/Database/Reference_genome/Mus_musculus_Ensembl_GRCm38_star_genome-index/Mus_musculus.GRCm38.95.gtf
bed=/Share/home/tiangeng/reference_genome/mm10_bed/mm10_Gencode_VM18.bed
for i in ${name[@]}
do
        echo "read_distribution.py -i ../STAR/${i}_STAR/${i}Aligned.sortedByCoord.out.bam  -r ${bed}"
nohup read_distribution.py -i ../STAR/${i}_STAR/${i}Aligned.sortedByCoord.out.bam -r ${bed} > ${i}.log 2>&1 &
done

```