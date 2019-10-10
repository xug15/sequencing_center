# 6.3 A quick start for featureCounts in SourceForge Subread

You need to provide read mapping results (in either SAM or BAM format) and an annotation
file for the read summarization. The example commands below assume your annotation file
is in GTF format.

* Summarize SAM format single-end reads using 5 threads:  

```sh
featureCounts -T 5 -a annotation.gtf -t exon -g gene_id -o counts.txt mapping_results_SE.sam
```

* Summarize BAM format single-end read data:  

```sh
featureCounts -a annotation.gtf -t exon -g gene_id -o counts.txt mapping_results_SE.bam
```

* Summarize multiple libraries at the same time:

```sh
featureCounts -a annotation.gtf -t exon -g gene_id -o counts.txt mapping_results1.bam mapping_results2.bam
```

* Summarize paired-end reads and count fragments (instead of reads):

```sh
featureCounts -p -a annotation.gtf -t exon -g gene_id -o counts.txt mapping_results_PE.bam
```

* Count fragments satisfying the fragment length criteria, eg. [50bp, 600bp]:

```sh
featureCounts -p -P -d 50 -D 600 -a annotation.gtf -t exon -g gene_id -o counts.txt mapping_results_PE.bam
```

* Count fragments which have both ends successfully aligned without considering the fragment length constraint:

```sh
featureCounts -p -B -a annotation.gtf -t exon -g gene_id -o counts.txt mapping_results_PE.bam
```

* Exclude chimeric fragments from the fragment counting:

```sh
featureCounts -p -C -a annotation.gtf -t exon -g gene_id -o counts.txt mapping_results_PE.bam
```

