

```sh
path=/home/xugang/zhangxiaolin20181015/2.cleandata/r6m310_NDSW43010/
out=/home/xugang/zhangxiaolin20181015/2.cleandata/r6m310_NDSW43010/

bowtie2 -p 6 -x /home/xugang/index/tair10.bowtie2/tair10.fa -1 $path/r6m310_NDSW43010_1.clean.fq.gz -2 $path/r6m310_NDSW43010_2.clean.fq.gz --seed 0 --very-sensitive -S $path/r6m310_NDSW43010.sam

```

```sh
path=/home/xugang/zhangxiaolin20181015/2.cleandata/r6m310_NDSW43010/
genome=/home/xugang/index/tair10.bowtie2/tair10.fa


samtools sort $path/r6m310_NDSW43010.bam > $path/r6m310_NDSW43010.sort.bam
samtools mpileup -E -gu -t DP -f $genome $path/r6m310_NDSW43010.sort.bam > $path/r6m310_NDSW43010.bcf
bcftools view $path/r6m310_NDSW43010.bcf > $path/r6m310_NDSW43010.vcf
bcftools call -m $path/r6m310_NDSW43010.vcf > $path/r6m310_NDSW43010.all.vcf
rm $path/r6m310_NDSW43010.vcf
rm $path/r6m310_NDSW43010.bcf
rm $path/r6m310_NDSW43010.sort.bam


#
docker run -it -v /home/xugang/zhangxiaolin20181015/2.cleandata/r6m310_NDSW43010:/data danielwells/shoremap
SHOREmap convert --marker r6m310_NDSW43010.all.vcf --folder r6m310_NDSW43010.shormap.all
SHOREmap backcross --chrsizes chrSize.txt --marker r6m310_NDSW43010.shormap.all/1_converted_variant.txt --consen r6m310_NDSW43010.shormap.all/1_converted_consen.txt --folder SHOREmap_backcross2/ -plot-bc --marker-score 20 --marker-freq 0.1 --min-coverage 5 --max-coverage 200 --bg 1_converted_variant.txt --bg-cov 1 --bg-freq 0.1 --bg-score 1 -non-EMS  --cluster 1 --marker-hit 1 -verbose

```

```sh
SHOREmap backcross --chrsizes chrSize.txt --marker r6m310_NDSW43010.shormap.all/1_converted_variant.txt --consen r6m310_NDSW43010.shormap.all/1_converted_consen.txt --folder SHOREmap_backcross2/ -plot-bc --marker-score 20 --marker-freq 0.1 --min-coverage 5 --max-coverage 200 --bg 1_converted_variant.txt --bg-cov 1 --bg-freq 0.1 --bg-score 1 -non-EMS  --cluster 1 --marker-hit 1 -verbose

SHOREmap annotate --chrsizes chrsizes.txt --snp markers.txt --chrom 4 --start 1600000 --end 1700000 --genome reference_seq.fa --gff gene_annotations.gff --peaks SHOREmap_boosted_peaks.txt --folder path/to/output_folder
```
