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
