# Luboxun

## 1. Unzip file

> /Share/home/tiangeng/project_result/Riboseq/project_190814_luboxun/a1-data

**a1.ungz.sh**
```sh
for i in `ls|grep ten`;do
echo $i;
    for j in `ls ${i}|grep _1.fq.gz$`;do
    echo ${i}/${j};
    gunzip ${i}/${j};
    done;

done;
```

And move the file into 
/Share/home/tiangeng/project_result/Riboseq/project_190814_luboxun/a2-rowdata
```sh
7-111-R_1.fq  
7-111-T_1.fq  
7-7-R_1.fq  
7-7-T_1.fq
```
## Remove adapter.

> a3-cutadapt

**a1.cutadaptRPF.sh**
```sh
export PATH=$PATH:/Share/home/tiangeng/software/cap-miRNA/bin
rpf=(7-111-R 7-7-R)
adapt=CTGTAGGCACCATCAAT

for i in ${rpf[@]}
do
        nohup cutadapt -m 25 -M 35 --match-read-wildcards -a $adapt -o ${i}_trimmed.fastq ../a2-rawdata/${i}_1.fq > ${i}_trimmed.log 2>&1 &
        done
```
**a2.cutadaptTotal.sh**
```sh
#!/bin/zsh

export PATH=$PATH:/Share/home/tiangeng/software/cap-miRNA/bin
rpf=(7-111-T 7-7-T)
adapt=CTGTAGGCACCATCAAT

for i in ${rpf[@]}
do
        nohup cutadapt -m 17 --match-read-wildcards -a $adapt -o ${i}_trimmed.fastq ../a2-rawdata/${i}_1.fq > ${i}_trimmed.log 2>&1 &
done
```
## Before QC

> a4-beforeQC

**a1.beforqc.sh**

```sh
#!/bin/bash

rpf=(7-111-T 7-7-T 7-111-R 7-7-R)
for i in ${rpf[@]}
do
        echo -e " nohup /Share/home/tiangeng/software/FastQC/fastqc ../a3-cutadapt/${i}_trimmed.fastq -o ./ > ${i}.log 2>&1 &";
        nohup /Share/home/tiangeng/software/FastQC/fastqc ../a3-cutadapt/${i}_trimmed.fastq -o ./ > ${i}.log 2>&1 & 
done
```

## Step3 : filter
> a5-filter

**a1.filter.sh**

```sh
#!/bin/bash

rpf=(7-111-T 7-7-T 7-111-R 7-7-R)
for i in ${rpf[@]}
do
        nohup /Share/app/fastx_toolkit0.14/bin/fastq_quality_filter -Q33 -v -q 25 -p 75 -i ../a3-cutadapt/${i}_trimmed.fastq -o ${i}_trimmedQfilter.fastq > ${i}_Qfilter.log 2>&1 & 
done
```

## Step4 : contam (remove ribosome sequence)

> a6-contam
**a1.contam.sh**
```sh
#!/bin/bash

bowtieindex=/Share/home/tiangeng/Database/Reference_genome/Mus-musculus_rRNA_bowtie-index/musRibosomal
name=(7-111-T 7-7-T 7-111-R 7-7-R)

for i in ${name[@]}
do
        nohup bowtie -n 0 -norc --best -l 15 -p 8 --un=nocontam_${i}.fastq $bowtieindex -q ../a5-filter/${i}_trimmedQfilter.fastq ${i}.alin > ${i}.err 2>&1 & 
done
```

**a2.merge.sh**
```sh
name=(7-111-R 7-7-R 7-111-T 7-7-T)
head='Iterm'
for i in ${name[@]};
do
 head+="\t${i}";
done
echo -e $head >merge.counter;

for i in ${name[@]};
do
echo ${i}.err;
head -n 3 ${i}.err > ${i}.err.tmp;
sed -i 's/#//g' ${i}.err.tmp;
sed -i 's/:/\t/g' ${i}.err.tmp;
sed -i 's/^ //g' ${i}.err.tmp;
done
#
begin1=${name[0]};
begin2=${name[1]};
name2=("${name[@]:2}");
join -t $'\t' ${begin1}.err.tmp ${begin2}.err.tmp >merge.tmp

for i in ${name2[@]};
do 
echo ${i}.err.tmp;
join -t $'\t' merge.tmp ${i}.err.tmp >>merge.tmp2;
mv merge.tmp2 merge.tmp
done
#
cat merge.counter merge.tmp > merge2.tmp;
cut -f 2- merge2.tmp > summary.txt
rm merge.counter merge.tmp *.err.tmp
echo -e "Iterm\nTotal\nrRNA\nnon rRNA">name.txt;
paste -d $'\t' name.txt summary.txt > summary2.txt
mv summary2.txt summary.txt
rm name.txt merge2.tmp
```

## Step5: afterQC (overrepresent reads)

> a7-afterqc
**a1.afterqc.sh**
```sh
#!/bin/bash

rpf=(7-111-T 7-7-T 7-111-R 7-7-R)
for i in ${rpf[@]}
do
        nohup /Share/home/tiangeng/software/FastQC/fastqc ../a6-contam/nocontam_${i}.fastq -o ./ >${i}.log 2>&1 & 
done
```

## step6: Tophat + readsNumCal_intron_v3.py

> a8-tophat
**a1.tophat.sh**
```sh
ssh node-0-11
cd /Share/home/tiangeng/project_result/Riboseq/project_190814_luboxun/a8-tophat
gtf=/Share/home/tiangeng/Database/Reference_genome/Mus_musculus_Ensembl_GRCm38_star_genome-index/Mus_musculus.GRCm38.95.gtf
bowtieIndex=/Share/home/tiangeng/Database/Reference_genome/Mus_musculus_Ensembl_GRCm38_bowtie_genome-index/Mus_musculus.GRCm38.dna.primary_assembly
transIndex=/Share/home/tiangeng/Database/Reference_genome/Mus_musculus_Ensembl_GRCm38_tophat_trans-index/Mus_musculus.GRCm38.95
export PATH=/Share/home/tiangeng/software/bowtie-1.1.2:$PATH
export PATH=/Share/home/tiangeng/anaconda2/bin:$PATH
fileName=(7-111-T 7-7-T)
for i in ${fileName[@]}
do
    nohup tophat -p 8 -o ${i}_out -g 10 --bowtie1 --read-realign-edit-dist 0 --library-type fr-secondstrand -G $gtf --transcriptome-index=$transIndex --no-novel-juncs --segment-length=15 $bowtieIndex ../a6-contam/nocontam_${i}.fastq > ${i}.log 2>&1 & 
done
```

**a2.tophat.sh**
```sh
ssh node-0-12
cd /Share/home/tiangeng/project_result/Riboseq/project_190814_luboxun/a8-tophat
gtf=/Share/home/tiangeng/Database/Reference_genome/Mus_musculus_Ensembl_GRCm38_star_genome-index/Mus_musculus.GRCm38.95.gtf
bowtieIndex=/Share/home/tiangeng/Database/Reference_genome/Mus_musculus_Ensembl_GRCm38_bowtie_genome-index/Mus_musculus.GRCm38.dna.primary_assembly
transIndex=/Share/home/tiangeng/Database/Reference_genome/Mus_musculus_Ensembl_GRCm38_tophat_trans-index/Mus_musculus.GRCm38.95
export PATH=/Share/home/tiangeng/software/bowtie-1.1.2:$PATH
export PATH=/Share/home/tiangeng/anaconda2/bin:$PATH
fileName=(7-111-R 7-7-R)
for i in ${fileName[@]}
do
    nohup tophat -p 8 -o ${i}_out -g 10 --bowtie1 --read-realign-edit-dist 0 --library-type fr-secondstrand -G $gtf --transcriptome-index=$transIndex --no-novel-juncs --segment-length=15 $bowtieIndex ../a6-contam/nocontam_${i}.fastq > ${i}.log 2>&1 & 
done
```
**a3.readsnum.sh**
```sh
#!/bin/bash
export PATH=/Share/home/tiangeng/anaconda2/bin:$PATH
name=(7-111-R 7-7-R 7-111-T 7-7-T)
gtf=/Share/home/tiangeng/Database/Reference_genome/Mus_musculus_Ensembl_GRCm38_star_genome-index/Mus_musculus.GRCm38.95.gtf
for i in ${name[@]}
do
       nohup python ./readsNumCal_intron_v3.py ./${i}_out/accepted_hits.bam $gtf nocontam_${i}_mappedNum_intron.txt nocontam_${i} > read.${i}.log 2>&1 &
done
```

**a4.merge.sh**
```sh
name=(7-111-R 7-7-R 7-111-T 7-7-T)
head='Iterm'
for i in ${name[@]};
do
 head+="\t${i}";
done
echo -e $head >merge.counter;

for i in ${name[@]};
do
echo nocontam_${i}_mappedNum_intron.txt;
cp nocontam_${i}_mappedNum_intron.txt ${i}_num;
sed -i 's/:/\t/g' ${i}_num;
done
#
begin1=${name[0]};
begin2=${name[1]};
name2=("${name[@]:2}");
join -t $'\t' ${begin1}_num ${begin2}_num >merge.tmp

for i in ${name2[@]};
do 
echo ${i}_num;
join -t $'\t' merge.tmp ${i}_num >>merge.tmp2;
mv merge.tmp2 merge.tmp
done
#
cat merge.counter merge.tmp > merge2.tmp;
mv merge2.tmp summary.txt
rm merge.counter merge.tmp *_num
```

## Step7 : STAR + periodicity

> a9-STAR

run with screen.

**a1.STAR.sh**
```sh
#!/bin/bash
export PATH=$PATH:~/software/STAR-master/bin/Linux_x86_64_static
genomeFile=/Share/home/tiangeng/Database/Reference_genome/Mus_musculus_Ensembl_GRCm38_star_genome-index
fileName=(7-111-R 7-7-R 7-111-T 7-7-T)
for i in ${fileName[@]}
do
        mkdir -p ${i}_STAR
        cd ${i}_STAR
        STAR --runThreadN 8 --alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 --genomeDir $genomeFile --readFilesIn ../../a6-contam/nocontam_${i}.fastq --outFileNamePrefix $i --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts
        cd ../
done
```
## Step8: readlength

> a10-readlength

**a1.cutadapter.sh**
```sh
#!/bin/bash
export PATH=$PATH:/Share/home/tiangeng/software/cap-miRNA/bin
rpf=(7-111-R 7-7-R 7-111-T 7-7-T)
adapt=CTGTAGGCACCATCAAT

for i in ${rpf[@]}
do
        nohup cutadapt --match-read-wildcards -a $adapt -o ${i}_trimmed.fastq ../a2-rawdata/${i}_1.fq > ${i}_trimmed.log 2>&1 &
done
```

**a2.readlength.sh**
```sh
#!/bin/bash
export PATH=/Share/home/tiangeng/anaconda2/bin:$PATH
rpf=(7-111-R 7-7-R 7-111-T 7-7-T)

for i in ${rpf[@]}
do
        nohup python ./fq_len_stat.py ./${i}_trimmed.fastq ${i}_length.png > len${i}.out 2>&1 &
done
```

**a3.readlength_filter.sh**
```sh
#!/bin/bash
export PATH=/Share/home/tiangeng/anaconda2/bin:$PATH
rpf=(7-111-R 7-7-R 7-111-T 7-7-T)

for i in ${rpf[@]}
do
        nohup python ./fq_len_stat.py ../a6-contam/nocontam_${i}.fastq ${i}_nontam.png > ${i}.nontam.out 2>&1 &
done
```

## Step9 : report
> a11-report

**a1.metaplots.sh**
```sh
#!/bin/bash
export PATH="/Share/home/tiangeng/anaconda3/bin:$PATH"
fileName=(7-111-R 7-7-R 7-111-T 7-7-T)
for i in ${fileName[@]}
do
        nohup metaplots -a /Share/home/tiangeng/Database/Reference_genome/prepare_transcripts_Mus-musculus -r ../a9-STAR/${i}_STAR/${i}Aligned.toTranscriptome.out.bam -o $i > ${i}.err  2>&1 &
done
```

## Step 10: calculate count
> a12-counter

**a1.rna_counter.sh**
```sh
#!/bin/bash
export PATH=/Share/home/tiangeng/anaconda2/bin:$PATH
name=(7-111-T 7-7-T)
gtf=/Share/home/tiangeng/Database/Reference_genome/Mus_musculus_Ensembl_GRCm38_star_genome-index/Mus_musculus.GRCm38.95.gtf
data_p=../a9-STAR
for i in ${name[@]}
do 
echo -e "htseq-count -q -f bam -s reverse $data_p/${i}_STAR/${i}Aligned.sortedByCoord.out.bam ${gtf} > ${i}.count";
nohup htseq-count -q -f bam -s reverse $data_p/${i}_STAR/${i}Aligned.sortedByCoord.out.bam ${gtf} > ${i}.count 2>&1 &
done;
#htseq-count [options] <alignment_files> <gff_file>
```
**a2.ribo_counter.sh**
```sh
export PATH=/Share/home/tiangeng/anaconda2/bin:$PATH
gtf=/Share/home/tiangeng/Database/Reference_genome/Mus_musculus_Ensembl_GRCm38_star_genome-index/Mus_musculus.GRCm38.95.gtf
data_p=../a9-STAR
name=(7-111-R 7-7-R)
for i in ${name[@]}
do 
echo -e "python RPF_count_CDS.py  $data_p/${i}_STAR/${i}Aligned.sortedByCoord.out.bam ${gtf} > ${i}.count";
nohup python RPF_count_CDS.py  $data_p/${i}_STAR/${i}Aligned.sortedByCoord.out.bam ${gtf} > ${i}.count 2>${i}.log
done;
```
**a3.merge.sh**
```sh
name=(7-111-R 7-7-R 7-111-T 7-7-T)
head='gene'
for i in ${name[@]};
do
 head+=" ${i}";
done
echo -e $head >merge.counter;

begin1=${name[0]};
begin2=${name[1]};
name2=("${name[@]:2}");
join ${begin1}.count ${begin2}.count >merge.tmp
commander='join';
for i in ${name2[@]};
do 
echo ${i}.count;
join merge.tmp ${i}.count >>merge.tmp2;
mv merge.tmp2 merge.tmp
done

cat merge.counter merge.tmp > merge2.tmp;
rm merge.tmp
mv merge2.tmp merge.counter
sed -i 's/ \+/\t/g' merge.counter

grep -v '^__' merge.counter > merge.counter2
mv merge.counter2 merge.counter 
```
## Step 11. Xtail

```sh
mv merge.counter lbx.merge.counter
```
```R
lbxd=read.table('/home/lbx.merge.counter',header=T,row.name=1)
mrna=lbxd[,c(3,4)]
rpf=lbxd[,c(1,2)]
condition=c("control","treat")
test.results=xtail(mrna,rpf,condition,bins=1000,threads=2)
summary(test.results)

#
test.tab=resultsTable(test.results);
head(test.tab,5)

write.table(test.tab,"/home/lbx_results.txt",quote=F,sep="\t");

# Visualization
pdf('lbxFC.pdf',width=6,height=4,paper='special')
lbxfc=plotFCs(test.results)
dev.off()
write.table(lbxfc$resultsTable,"/home/lbxfc_results.txt",quote=F,sep="\t");

pdf('lbxRs.pdf',width=6,height=4,paper='special')
lbxrs=plotRs(test.results)
dev.off()
write.table(lbxrs$resultsTable,"/home/lbxrs_results.txt",quote=F,sep="\t");

pdf('lbxvolcano.pdf',width=6,height=4,paper='special')
volcanoPlot(test.results)
dev.off()
```
## Step 12 Ribo-code
> 
**a1.trans_ann.sh**
```sh
#!/bin/bash
export PATH=~/anaconda3/bin:$PATH
gtf=/Share/home/tiangeng/Database/Reference_genome/Mus_musculus_Ensembl_GRCm38_star_genome-index/Mus_musculus.GRCm38.95.gtf
fa=/Share/home/tiangeng/Database/Reference_genome/Mus_musculus_Ensembl_GRCm38_star_genome-index/Mus_musculus.GRCm38.dna.primary_assembly.fa
prepare_transcripts -g $gtf -f $fa -o /Share/home/tiangeng/Database/Reference_genome/prepare_transcripts_Mus-musculus
```
**a2.mergefile.sh**
```sh
#!/bin/bash
name=(7-111-R 7-111-T 7-7-R 7-7-T)
rm b1.merge_transcriptome.txt
for i in ${name[@]};
do 
echo /Share/home/tiangeng/project_result/Riboseq/project_190814_luboxun/a9-STAR/${i}_STAR/${i}Aligned.toTranscriptome.out.bam >>b1.merge_transcriptome.txt;
done

```

**a3.metaplot.sh**
```sh
#!/bin/bash
export PATH=~/anaconda3/bin:$PATH
mkdir b3-metaplot
rico_ann=/Share/home/tiangeng/Database/Reference_genome/prepare_transcripts_Mus-musculus
name=7-111-R
nohup metaplots -a $rico_ann -i b1.merge_transcriptome.txt -o b3-metaplot/meta > b3-metaplot/metaplot.log 2>&1 &
```

**a4.Ribocode.sh**
```sh
#!/bin/bash
export PATH=~/anaconda3/bin:$PATH
mkdir b4-RiboCode
rico_ann=/Share/home/tiangeng/Database/Reference_genome/prepare_transcripts_Mus-musculus
config=b3-metaplot/meta_pre_config.txt
nohup RiboCode -a ${rico_ann} -c ${config} -l no -g -o b4-RiboCode/RiboCode_ORFs_result  > b4-RiboCode/ribocode.log 2>&1 &
``` 

**a5.density.sh**
```sh
#!/bin/bash
export PATH=~/anaconda3/bin:$PATH

mkdir b5-plot_density_orf
rico_ann=/Share/home/tiangeng/Database/Reference_genome/prepare_transcripts_Mus-musculus
config=b3-metaplot/meta_pre_config.txt
tran_id=ENSMUST00000187405
orf_start=1
orf_end=900

plot_orf_density -a ${rico_ann} -c ${config} -t ${tran_id} -s ${orf_start} -e ${orf_end} -o b5-plot_density_orf/a1.transid

```

**a6.orfcount.sh**
```sh
#!/bin/bash
export PATH=~/anaconda3/bin:$PATH
#RiboCode_ORFs_result.gtf
mkdir b6-ORF_count
name=(7-111-R 7-7-R)

for i in ${name[@]};
do 
nohup ORFcount -g b4-RiboCode/RiboCode_ORFs_result.gtf -r /Share/home/tiangeng/project_result/Riboseq/project_190814_luboxun/a9-STAR/${i}_STAR/${i}Aligned.sortedByCoord.out.bam -f 15 -l 5 -e 100 -m 26 -M 34 -o b6-ORF_count/${i}.ORF.counts > b6-ORF_count/${i}.log 2>&1 &
done;
```

## RiboMiner.

```sh
sudo singularity exec -w -B /home/xugang/singularity_image:/home/share /home/xugang/singularity_image/ribocodeminer/ bash

pip uninstall RiboMiner

把包克隆下来，然后
cd RiboMiner-0.1

python setup.py install

 tAI --version



```

**Star docker.**
data

### Prepare sequences and annotaiton files on transcriptome level.

**a0-prepare.sh**

```sh
gtf=/home/share/riboseq/Mus-musculus_Ensembl_release-85/Mus_musculus.GRCm38.85.gtf
fa=/home/share/riboseq/Mus-musculus_Ensembl_release-85/Mus_musculus.GRCm38.dna.nonchromosomal-chromosomal.genome.fa
out=/home/share/riboseq/Ribocode

prepare_transcripts -g $gtf -f $fa -o $out
```

### Prepare the longest transcript annotaion files.

**a1-annotation.sh**

```sh
gtf=/home/share/riboseq/Mus-musculus_Ensembl_release-85/Mus_musculus.GRCm38.85.clean.gtf
fa=/home/share/riboseq/Mus-musculus_Ensembl_release-85/Mus_musculus.GRCm38.dna.nonchromosomal-chromosomal.genome.fa
cds=/home/share/riboseq/Ribocode/transcripts_cds.txt
trans=/home/share/riboseq/Ribocode/transcripts_sequence.fa
longest_info=/home/share/riboseq/RiboMiner/longest.transcripts.info.txt
all_info=/home/share/riboseq/RiboMiner/all.transcripts.info.txt
out=/home/share/riboseq/RiboMiner

mkdir $out

echo "OutputTranscriptInfo -c $cds -g $gtf -f $trans -o $longest_info -O $all_info"
OutputTranscriptInfo -c $cds -g $gtf -f $trans -o $longest_info -O $all_info
```


### Prepare the sequence file for the longest transcripts

a2-transcript.sh
```sh
transcripts_sequence=/home/share/riboseq/Ribocode/transcripts_sequence.fa
longest_info=/home/share/riboseq/RiboMiner/longest.transcripts.info.txt
transcript=/home/share/riboseq/RiboMiner/transcript

GetProteinCodingSequence -i $transcripts_sequence  -c $longest_info -o $transcript --mode whole --table 1 
```

### Prepare the UTR sequence

a3-utr.sh

```sh
transcripts_sequence=/home/share/riboseq/Ribocode/transcripts_sequence.fa
utr=/home/share/riboseq/RiboMiner/utr
transcript_cds=/home/share/riboseq/RiboMiner/transcript_cds_sequences.fa

GetUTRSequences -i $transcripts_sequence -o $utr -c $transcript_cds
```

### 

a4-metaplot.sh

```sh
metaplots -a /data/reference/RiboCode_annot -r /data/data/colAligned.toTranscriptome.out.bam -o /data/data/a4-col
metaplots -a /data/reference/RiboCode_annot -r /data/data/d14Aligned.toTranscriptome.out.bam -o /data/data/a4-d14
```

a5-periodicity.sh

```sh

Periodicity -i /data/data/colAligned.toTranscriptome.sort.bam -a /data/reference/RiboCode_annot -o /data/data/a5-col_periodicity -c /data/reference/tair_analy/longest.transcripts.info.txt -L 25 -R 35
Periodicity -i /data/data/d14Aligned.toTranscriptome.sort.bam -a /data/reference/RiboCode_annot -o /data/data/a5-d14_periodicity -c /data/reference/tair_analy/longest.transcripts.info.txt -L 25 -R 35
```

attributes.txt

```sh
/data/data/colAligned.toTranscriptome.sort.bam  30,31,32,33,34  12,12,13,13,13  col
/data/data/d14Aligned.toTranscriptome.sort.bam  30,31,32,33,34  11,12,12,13,13  d14
```

### Reads distribution among different reading frames.

a6-ribodensitydiffframe.sh
```sh
RiboDensityOfDiffFrames -f /data/data/attributes.txt -c /data/reference/tair_analy/longest.transcripts.info.txt -o /data/data/a6-ribo-density-diff-frame
```

### DNA contamination.

a7-dna-contamination.sh

```sh
StatisticReadsOnDNAsContam -i  /data/data/colAligned.sortedByCoord.out.bam  -g /data/reference/tair/Arabidopsis_thaliana.TAIR10.43.gtf -o /data/data/a7-dna-contamination.col 
StatisticReadsOnDNAsContam -i  /data/data/d14Aligned.sortedByCoord.out.bam  -g /data/reference/tair/Arabidopsis_thaliana.TAIR10.43.gtf -o /data/data/a7-dna-contamination.d14  
```

### Metagene Analysis (MA)

a8-metagene.sh

**Metagene analysis along the whole transcript region.**

```sh
MetageneAnalysisForTheWholeRegions -f /data/data/attributes.txt -c /data/reference/tair_analy/longest.transcripts.info.txt -o /data/data/a8-metagene -b 15,90,60 -l 100 -n 10 -m 1 -e 5 --plot yes
```

a9-PlotMetageneAnalysisForTheWholeRegions.sh

```sh
PlotMetageneAnalysisForTheWholeRegions -i /data/data/a8-metagene_scaled_density_dataframe.txt -o /data/data/a9-meta_gene_whole_regin -g group1,group2 -r group1,group2 -b 15,90,60 --mode all 

```

## Metagene analysis on CDS regions.

**b1-metagene_cds.sh**

```sh
## the first way
MetageneAnalysis -f /data/data/attributes.txt -c /data/reference/tair_analy/longest.transcripts.info.txt -o /data/data/b1-meat-cds -U codon -M RPKM -u 0 -d 500 -l 100 -n 10 -m 1 -e 5 --norm yes -y 100 --CI 0.95 --type CDS

```

b2-metagene_utr.sh

```sh
MetageneAnalysis -f /data/data/attributes.txt -c /data/reference/tair_analy/longest.transcripts.info.txt -o /data/data/b2-meat-utr -U nt -M RPKM -u 100 -d 100 -l 100 -n 10 -m 1 -e 5 --norm yes -y 50 --CI 0.95 --type UTR
```

b3-PolarityCalculation.sh

```sh
PolarityCalculation -f /data/data/attributes.txt -c /data/reference/tair_analy/longest.transcripts.info.txt -o /data/data/b3-polarity -n 64
```

b4-PlotPolarity.sh

```sh

PlotPolarity -i /data/data/b3-polarity_polarity_dataframe.txt -o /data/data/b4-plotpolarity -g col,d14 -r col__d14 -y 5 

```

### Feature Analysis (FA)

**Pick out transcripts enriched ribosomes on specific region.**

b5-transcript-enrich.sh
```sh
RiboDensityForSpecificRegion -f /data/data/attributes.txt -c /data/reference/tair_analy/longest.transcripts.info.txt -o /data/data/b5-transcript-enrich -U codon -M RPKM -L 25 -R 75
```

**Ribosome density at each kind of AA or codon**
b6-ribosome-aa.sh
```sh
RiboDensityAtEachKindAAOrCodon -f /data/data/attributes.txt -c /data/reference/tair_analy/longest.transcripts.info.txt -o /data/data/b6-ribosome-aa -M counts -S /data/select_trans.txt -l 100 -n 10 --table 1 -F /data/reference/tair_analy/transcript_cds_sequences.fa 
```

**Ribosome density on amino acids with positive or negative charge**

b7-PlotRiboDensityAtEachKindAAOrCodon.sh

```sh
PlotRiboDensityAtEachKindAAOrCodon -i /data/data/b6-ribosome-aa_all_codon_density.txt -o /data/data/b7-PlotRiboDensityAtEachKindAAOrCodon -g col,d14 -r col__d14 --level AA
```

**Pausing score of each triplete amino acid.**

b8-PausingScore.sh

```sh
PausingScore -f /data/data/attributes.txt -c /data/reference/tair_analy/longest.transcripts.info.txt -o  /data/data/b8-PausingScore -M counts -S /data/select_trans.txt  -l 100 -n 10 --table 1 -F  /data/reference/tair_analy/transcript_cds_sequences.fa
```

b9-ProcessPausingScore.sh
```sh
ProcessPausingScore -i /data/data/b8-PausingScore_col_pausing_score.txt,/data/data/b8-PausingScore_d14_pausing_score.txt -o /data/data/b9-ProcessPausingScore -g col,d14 -r col__d14 --mode raw --ratio_filter 2 --pausing_score_filter 0.5
```


**Ribosome density around the triplete amino acid (tri-AA) motifs.**

* As for a specific tri-AA motif, such as poly-proline (PPP)

c0-RiboDensityAroundTripleteAAMotifs.sh

```sh
RiboDensityAroundTripleteAAMotifs -f /data/data/attributes.txt -c /data/reference/tair_analy/longest.transcripts.info.txt -o /data/data/c0-RiboDensityAroundTripleteAAMotifs_PPP -M counts -S /data/select_trans.txt -l 100 -n 10 --table 1 -F /data/reference/tair_analy/transcript_cds_sequences.fa --type2 PPP --type1 PP
```

c1-PlotRiboDensityAroundTriAAMotifs.sh

```sh
PlotRiboDensityAroundTriAAMotifs -i /data/data/c0-RiboDensityAroundTripleteAAMotifs_PPP_motifDensity_dataframe.txt -o /data/data/c1-PPP_plot -g col,d14 -r col__d14 --mode mean
```

* Using following command to do such job:

c2-RiboDensityAroundTripleteAAMotifs.sh

```sh

RiboDensityAroundTripleteAAMotifs -f /data/data/attributes.txt -c /data/reference/tair_analy/longest.transcripts.info.txt -o  /data/data/c2-RiboDensityAroundTripleteAAMotifs -M counts -S /data/select_trans.txt -l 100 -n 10 --table 1 -F /data/reference/tair_analy/transcript_cds_sequences.fa --motifList1 /data/reference/tri_AA_motifs1.txt --motifList2 /data/reference/tri_AA_motifs2.txt

```
c2b-PlotRiboDensityAroundTriAAMotifs.sh
```sh
PlotRiboDensityAroundTriAAMotifs -i /data/data/c2-RiboDensityAroundTripleteAAMotifs_motifDensity_dataframe.txt -o /data/data/c2b-PPP_plot -g col,d14 -r col__d14 --mode mean


```


**RPFdist calculation.**

c3-RPFdist.sh

```sh
RPFdist -f /data/data/attributes.txt -c /data/reference/tair_analy/longest.transcripts.info.txt -o /data/data/c3-RPFdist -M counts -S /data/select_trans.txt -l 100 -n 10 -m 1 -e 5
```

**GC contents for sequences with a fasta format.**

c4-GCContent.sh

```sh
GCContent -i /data/reference/tair_analy/transcript_cds_sequences.fa -o /data/data/c4-GCContent-normal --mode normal
GCContent -i /data/reference/tair_analy/transcript_cds_sequences.fa -o /data/data/c4-GCContent-frames --mode frames
```

c5-PlotGCContent.sh

```sh
## normal mode
PlotGCContent -i /data/data/c4-GCContent-normal_GC_content.txt -o /data/data/c5-PlotGCContent-normal --mode normal
## frames mode
PlotGCContent -i /data/data/c4-GCContent-frames_GC_content_frames.txt -o /data/data/c5-PlotGCContent-frames --mode frames
```

**Local tRNA adaptation index and global tRNA adaptation index**

c6-tAI.sh

```sh
tAI -i /data/reference/tair_analy/transcript_cds_sequences_tAI.fa -t tair -o /data/data/c6-tAI -u 0 -d 500 --table 1 -N /data/aratha/araTha1-tRNAs-confidence-set.out


```
c7-tAIPlot.sh

```sh
tAIPlot -i /data/data/c6-tAI_tAI_dataframe.txt -o /data/data/c7-tAIPlot -u 0 -d 500 --mode all --start 5 --window 7 --step 1
```
**Local codon adaptation index and global codon adaptation index**

c8-cAI.sh

```sh
cAI -i /data/reference/tair_analy/transcript_cds_sequences_tAI.fa -o /data/data/c8-cAI -t tair -u 0 -d 500 --reference /data/reference/tair_analy/reference.fa

```

c9-cAIPlot.sh

```sh
cAIPlot -i /data/data/c8-cAI_local_cAI_dataframe.txt -o /data/data/c9-cAIPlot -u 0 -d 500 --mode all --start 5 --window 7 --step 1

```

```sh
GetProteinCodingSequence -i <transcrits_sequence.fa> -c <longest.trans.info.txt> -S /data/select_trans.txt -o <output_prefix> --mode whole --table 1 --id-type transcript-id
```

**Hydrophobicity calculation and Charge amino acids**

d1-hydropathyCharge.sh

```sh
## hydrophobicity calculation
hydropathyCharge  -i /data/reference/tair_analy/transcript_cds_sequences_tAI.fa -o /data/data/d1-hydropathyCharge -t select_gene --index /data/reference/hydropathy_index.txt -u 0 -d 500 --table 1
```

d2-charge.sh

```sh
##
hydropathyCharge  -i /data/reference/tair_analy/transcript_cds_sequences_tAI.fa -o /data/data/d2-charge -t select_gene --index /data/reference/AA_charge_index.txt -u 0 -d 500 --table 1
```

d3-PlotHydropathyCharge.sh

```sh
## hydrophobicity
PlotHydropathyCharge -i /data/data/d1-hydropathyCharge_values_dataframe.txt -o /data/data/d3-PlotHydropathyCharge  -u 0 -d 500 --mode all --ylab "Average Hydrophobicity"
```
d4-Plotcharges.sh

```sh
## charge
PlotHydropathyCharge -i /data/data/d2-charge_values_dataframe.txt -o /data/data/d4-Plotcharges -u 0 -d 500 --mode all --ylab "Average Charges"
```
### Enrichment Analysis (EA)

**Step 1: Calculate ribosome density at each position for each transcript.**

d5-RiboDensityAtEachPosition.sh

```sh
RiboDensityAtEachPosition -c /data/reference/tair_analy/longest.transcripts.info.txt -f /data/data/attributes.txt -o /data/data/d5-RiboDensityAtEachPosition -U codon
```
**Step 2: Calculate mean ribosome density for different replicates.**

d6-enrichmentMeanDensity.sh

```sh
enrichmentMeanDensity -i /data/data/d5-RiboDensityAtEachPosition_col_cds_codon_density.txt,/data/data/d5-RiboDensityAtEachPosition_d14_cds_codon_density.txt -o /data/data/d6-enrichmentMeanDensity
```
**Step 3: Enrichment analysis.**

d7-EnrichmentAnalysis.sh

```sh
## all transcripts
EnrichmentAnalysis --ctrl /data/data/d5-RiboDensityAtEachPosition_col_cds_codon_density.txt --treat /data/data/d5-RiboDensityAtEachPosition_d14_cds_codon_density.txt -c /data/reference/tair_analy/longest.transcripts.info.txt -o /data/data/d7-EnrichmentAnalysis -U codon -M RPKM -l 150 -n 10 -m 1 -e 30 --CI 0.95 -u 0 -d 500


```

```sh
## specific transcripts
EnrichmentAnalysis --ctrl <total-translatome.txt> --treat <IP-translatome.txt> -c /data/reference/tair_analy/longest.transcripts.info.txt -o <output_prefix> -U codon -M RPKM -l 150 -n 10 -m 1 -e 30 --CI 0.95 -u 0 -d 500 -S /data/select_trans.txt
```


**Step 4: Plot the enrichment ratio.**

d8-PlotEnrichmentRatio.sh

```sh
PlotEnrichmentRatio -i /data/data/d7-EnrichmentAnalysis_enrichment_dataframe.txt -o /data/data/d8-PlotEnrichmentRatio -u 0 -d 500 --unit codon --mode all

```
**Notes: if you want to see the enrichment ratio for a single transcript, the EnrichmentAnalysisForSingleTrans would be helpful.**

```sh
EnrichmentAnalysisForSingleTrans -i <output_prefix_codon_ratio.txt> -s <transcript_name> -o <output_prefix> -c <longest.trans.info.txt>  --id-type transcript_id --slide-window y --axhline 1
```



