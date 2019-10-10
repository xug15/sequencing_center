# RSeQC: An RNA-seq Quality Control Package

* Options:
--version
show program’s version number and exit

-h, --help
show this help message and exit

-i INPUT_FILE, --input-file=INPUT_FILE
Alignment file in BAM or SAM format.

-r REF_GENE_MODEL, --refgene=REF_GENE_MODEL
Reference gene model in bed format.

Example:
```sh
read_distribution.py  -i Pairend_StrandSpecific_51mer_Human_hg19.bam -r hg19.refseq.bed12
```
Output:

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