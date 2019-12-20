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

**Build star index**

a0.generate.sh

```sh
#!/bin/bash

fasta=/home/xugang/singularity_image/riboseq/Mus-musculus_Ensembl_release-85/genome/Mus_musculus.GRCm38.dna.primary_assembly.fa
gtf=/home/xugang/singularity_image/riboseq/Mus-musculus_Ensembl_release-85/genome/Mus_musculus.GRCm38.87.gtf
genomeout=/home/xugang/singularity_image/riboseq/Mus-musculus_Ensembl_release-85/genome/star

nohup STAR --runThreadN 8 --runMode genomeGenerate --limitGenomeGenerateRAM 162003700778 --genomeDir ${genomeout} --genomeFastaFiles ${fasta} --sjdbGTFfile ${gtf} >log.txt 2>&1 &

```

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

cd /home/test

pip uninstall RiboMiner

把包克隆下来，然后
cd RiboMiner-0.1

python setup.py install

 tAI --version



```

**Star docker.**
data

### Prepare sequences and annotaiton files on transcriptome level.


**a.datapreparation.sh**  

```sh
gtf=/home/share/riboseq/mus_ensemble/Mus_musculus.GRCm38.87.gtf
fa=/home/share/riboseq/mus_ensemble/Mus_musculus.GRCm38.dna.primary_assembly.fa
out=/home/share/riboseq/Ribocode

# Prepare sequences and annotaiton files on transcriptome level.
prepare_transcripts -g $gtf -f $fa -o $out

# Prepare the longest transcript annotaion files.

cds=/home/share/riboseq/Ribocode/transcripts_cds.txt
trans=/home/share/riboseq/Ribocode/transcripts_sequence.fa
longest_info=/home/share/riboseq/RiboMiner/longest.transcripts.info.txt
all_info=/home/share/riboseq/RiboMiner/all.transcripts.info.txt
out=/home/share/riboseq/RiboMiner
mkdir $out
echo "OutputTranscriptInfo -c $cds -g $gtf -f $trans -o $longest_info -O $all_info"
OutputTranscriptInfo -c $cds -g $gtf -f $trans -o $longest_info -O $all_info 

# Prepare the sequence file for the longest transcripts
transcripts_sequence=/home/share/riboseq/Ribocode/transcripts_sequence.fa
longest_info=/home/share/riboseq/RiboMiner/longest.transcripts.info.txt
transcript=/home/share/riboseq/RiboMiner/transcript

GetProteinCodingSequence -i $transcripts_sequence  -c $longest_info -o $transcript --mode whole --table 1

# Sometines, UTR sequences are needed. In this case, GetUTRSequences maybe helpful:
transcripts_sequence=/home/share/riboseq/Ribocode/transcripts_sequence.fa
utr=/home/share/riboseq/RiboMiner/utr
transcript_cds=/home/share/riboseq/Ribocode/transcripts_cds.txt

GetUTRSequences -i $transcripts_sequence -o $utr -c $transcript_cds
```

**b.qualitycontrol.sh**

```sh
#Periodicity checking
#Ribosome profiling data with a good quality tend to have a good 3-nt periodicity.
ribocode='/home/share/riboseq/Ribocode'
a111R='/home/share/riboseq/a9-STAR/7-111-R_STAR/7-111-RAligned.toTranscriptome.out.bam'
a111R_o='/home/share/riboseq/metaplot/7111r'
a7R='/home/share/riboseq/a9-STAR/7-7-R_STAR/7-7-RAligned.toTranscriptome.out.bam'
a7R_o='/home/share/riboseq/metaplot/77r'
[[ -d /home/share/riboseq/metaplot/ ]] || mkdir /home/share/riboseq/metaplot/
echo -e "metaplots -a $ribocode -r ${a111R} -o ${a111R_o}"
metaplots -a $ribocode -r ${a111R} -o ${a111R_o}
echo ""
echo -e "metaplots -a $ribocode -r ${a7R} -o ${a7R_o}"
metaplots -a $ribocode -r ${a7R} -o ${a7R_o}
echo ""
ribocode='/home/share/riboseq/Ribocode'
long='/home/share/riboseq/RiboMiner/longest.transcripts.info.txt'
a111R='/home/share/riboseq/a9-STAR/7-111-R_STAR/7-111-R.toTranscriptome.sort.bam'
a111R_o='/home/share/riboseq/a5-periodicity/7111r'
a7R='/home/share/riboseq/a9-STAR/7-7-R_STAR/7-7-R.toTranscriptome.sort.bam'
a7R_o='/home/share/riboseq/a5-periodicity/77r'
[ -d /home/share/riboseq/a5-periodicity ] || mkdir /home/share/riboseq/a5-periodicity
echo Periodicity -i $a111R -a $ribocode -o $a111R_o -c $long -L 25 -R 35
Periodicity -i $a111R -a $ribocode -o $a111R_o -c $long -L 25 -R 35
echo Periodicity -i $a7R -a $ribocode -o $a7R_o -c $long -L 25 -R 35
Periodicity -i $a7R -a $ribocode -o $a7R_o -c $long -L 25 -R 35

#Reads distribution among different reading frames.
RiboDensityOfDiffFrames -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/a6-ribo-density-diff-fram

#Length distribution.
LengthDistribution -i /home/share/riboseq/a9-STAR/7-111-R_STAR/7-111-RAligned.sortedByCoord.out.bam -o /home/share/riboseq/r111.length -f bam
LengthDistribution -i /home/share/riboseq/a9-STAR/7-7-R_STAR/7-7-RAligned.sortedByCoord.out.bam -o /home/share/riboseq/r7.length -f bam

#DNA contamination.
gtf=/home/share/riboseq/mus_ensemble/Mus_musculus.GRCm38.87.gtf
StatisticReadsOnDNAsContam -i  /home/share/riboseq/a9-STAR/7-111-R_STAR/7-111-RAligned.sortedByCoord.out.bam  -g $gtf -o /home/share/riboseq/a7-dna-contamination.111
StatisticReadsOnDNAsContam -i  /home/share/riboseq/a9-STAR/7-7-R_STAR/7-7-RAligned.sortedByCoord.out.bam  -g $gtf -o /home/share/riboseq/a7-dna-contamination.7
```

**c.metagene_analysis.sh**

```sh
#Metagene analysis along the whole transcript region.
MetageneAnalysisForTheWholeRegions -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/a8-metagene -b 15,90,60 -l 100 -n 10 -m 1 -e 5 --plot yes
PlotMetageneAnalysisForTheWholeRegions -i /home/share/riboseq/a8-metagene_scaled_density_dataframe.txt -o /home/share/riboseq/a9-meta_gene_whole_regin -g r7,r111 -r r7__r111 -b 15,90,60 --mode all

#Metagene analysis on CDS regions.
MetageneAnalysis -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/b1-meat-cds -U codon -M RPKM -u 0 -d 500 -l 100 -n 10 -m 1 -e 5 --norm yes -y 100 --CI 0.95 --type CDS

#Metagene analysis on UTR regions.
echo MetageneAnalysis -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/b2-meat-utr -U nt -M RPKM -u 100 -d 100 -l 100 -n 10 -m 1 -e 5 --norm yes -y 50 --CI 0.95 --type UTR
MetageneAnalysis -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/b2-meat-utr -U nt -M RPKM -u 100 -d 100 -l 100 -n 10 -m 1 -e 5 --norm yes -y 50 --CI 0.95 --type UTR

#Polarity calculation.
echo PolarityCalculation -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/b3-polarity -n 64
PolarityCalculation -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/b3-polarity -n 64
```

**d1.feature_analysis.sh**

```sh

# Pick out transcripts enriched ribosomes on specific region.
echo "Pick out transcripts enriched ribosomes on specific region.";
`date`
echo RiboDensityForSpecificRegion -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/b5-transcript-enrich -U codon -M RPKM -L 25 -R 75
RiboDensityForSpecificRegion -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/b5-transcript-enrich -U codon -M RPKM -L 25 -R 75
echo "Finished:Pick out transcripts enriched ribosomes on specific region."
`date`

echo "# Ribosome density at each kind of AA or codon."
date
# Ribosome density at each kind of AA or codon.
echo RiboDensityAtEachKindAAOrCodon -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/b6-ribosome-aa -M counts -S /home/share/riboseq/selec_trans_longest.txt -l 100 -n 10 --table 1 -F /home/share/riboseq/RiboMiner/transcript_cds_sequences.fa
RiboDensityAtEachKindAAOrCodon -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/b6-ribosome-aa -M counts -S /home/share/riboseq/selec_trans_longest.txt -l 100 -n 10 --table 1 -F /home/share/riboseq/RiboMiner/transcript_cds_sequences.fa
echo "Finished:Ribosome density at each kind of AA or codon"
date

echo "# Ribosome density on amino acids with positive or negative charge"
date
# Ribosome density on amino acids with positive or negative charge
echo PlotRiboDensityAtEachKindAAOrCodon -i /home/share/riboseq/b6-ribosome-aa_all_codon_density.txt -o /home/share/riboseq/b7-PlotRiboDensityAtEachKindAAOrCodon -g r7,r111 -r r7__r111 --level AA
PlotRiboDensityAtEachKindAAOrCodon -i /home/share/riboseq/b6-ribosome-aa_all_codon_density.txt -o /home/share/riboseq/b7-PlotRiboDensityAtEachKindAAOrCodon -g r7,r111 -r r7__r111 --level AA
echo "Finished: Ribosome density on amino acids with positive or negative charge"
date

echo "Pausing score of each triplete amino acid.";
date
# Pausing score of each triplete amino acid.
echo PausingScore -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o  /home/share/riboseq/b8-PausingScore -M counts -S /home/share/riboseq/selec_trans_longest.txt  -l 100 -n 10 --table 1 -F  /home/share/riboseq/RiboMiner/transcript_cds_sequences.fa
PausingScore -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o  /home/share/riboseq/b8-PausingScore -M counts -S /home/share/riboseq/selec_trans_longest.txt  -l 100 -n 10 --table 1 -F  /home/share/riboseq/RiboMiner/transcript_cds_sequences.fa

echo ProcessPausingScore -i /home/share/riboseq/b8-PausingScore_r7_pausing_score.txt,/home/share/riboseq/b8-PausingScore_r111_pausing_score.txt -o /home/share/riboseq/b9-ProcessPausingScore -g r7,r111 -r r7__r111 --mode raw --ratio_filter 2 --pausing_score_filter 0.5
ProcessPausingScore -i /home/share/riboseq/b8-PausingScore_r7_pausing_score.txt,/home/share/riboseq/b8-PausingScore_r111_pausing_score.txt -o /home/share/riboseq/b9-ProcessPausingScore -g r7,r111 -r r7__r111 --mode raw --ratio_filter 2 --pausing_score_filter 0.5
echo "Finshed: Pausing score of each triplete amino acid."
date

echo "# Ribosome density around the triplete amino acid (tri-AA) motifs. \
#As for a specific tri-AA motif, such as poly-proline (PPP)";
date

# Ribosome density around the triplete amino acid (tri-AA) motifs.
#As for a specific tri-AA motif, such as poly-proline (PPP)

echo "RPFdist calculation.";
date
# RPFdist calculation.
echo RPFdist -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/c3-RPFdist -M counts -S /home/share/riboseq/selec_trans_longest.txt -l 100 -n 10 -m 1 -e 5
RPFdist -f /home/share/riboseq/attributes.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/c3-RPFdist -M counts -S /home/share/riboseq/selec_trans_longest.txt -l 100 -n 10 -m 1 -e 5
echo "Finshed:RPFdist calculation."
date

echo "# GC contents for sequences with a fasta format.";
# GC contents for sequences with a fasta format.
date
echo GCContent -i /home/share/riboseq/RiboMiner/transcript_cds_sequences.fa -o /home/share/riboseq/c4-GCContent-normal --mode normal
echo GCContent -i /home/share/riboseq/RiboMiner/transcript_cds_sequences.fa -o /home/share/riboseq/c4-GCContent-frames --mode frames
GCContent -i /home/share/riboseq/RiboMiner/transcript_cds_sequences.fa -o /home/share/riboseq/c4-GCContent-normal --mode normal
GCContent -i /home/share/riboseq/RiboMiner/transcript_cds_sequences.fa -o /home/share/riboseq/c4-GCContent-frames --mode frames
echo PlotGCContent -i /home/share/riboseq/c4-GCContent-normal_GC_content.txt -o /home/share/riboseq/c5-PlotGCContent-normal --mode normal
PlotGCContent -i /home/share/riboseq/c4-GCContent-normal_GC_content.txt -o /home/share/riboseq/c5-PlotGCContent-normal --mode normal

echo "Finished:GC contents for sequences with a fasta format."
date

echo "## frames mode"
date
## frames mode
echo PlotGCContent -i /home/share/riboseq/c4-GCContent-frames_GC_content_frames.txt -o /home/share/riboseq/c5-PlotGCContent-frames --mode frames
PlotGCContent -i /home/share/riboseq/c4-GCContent-normal_GC_content.txt -o /home/share/riboseq/c5-PlotGCContent-normal --mode normal
PlotGCContent -i /home/share/riboseq/c4-GCContent-frames_GC_content_frames.txt -o /home/share/riboseq/c5-PlotGCContent-frames --mode frames
echo "Finished:## frames mode"
date
```

**d2.feature_analysis.sh**

We need tRNA file, that can be download from [GtRNAdb](http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/Mmusc10-gene-list.html)

[tair10 tRNA](./araTha1-tRNAs-confidence-set.txt)  
[mm10 tRNA](./mouse-tRNAs-confidence.txt)

|Chr|tRNA|Begin|End|Isotype|Anticodon|Upstream|Downstream|
| -| -| -| -|-| -|-|-|
|chr6	|95	|58141949	|58141877	|Ala	|AGC	|tttctccctc	|gtttcttgtc|
|chr6	|25	|26751918	|26751990	|Ala	|AGC	|agtgtagtgt	|gcttctttta|


```sh


echo "Local tRNA adaptation index and global tRNA adaptation index"
date
#Local tRNA adaptation index and global tRNA adaptation index
echo tAI -i /home/share/riboseq/RiboMiner/transcript_cds_sequences.fa -t mouse -o /home/share/riboseq/c6-tAI -u 0 -d 500 --table 1 -N /home/share/riboseq/mouse-tRNAs-confidence.txt
tAI -i /home/share/riboseq/RiboMiner/transcript_cds_sequences.fa -t mouse -o /home/share/riboseq/c6-tAI -u 0 -d 500 --table 1 -N /home/share/riboseq/mouse-tRNAs-confidence.txt

echo tAIPlot -i /home/share/riboseq/c6-tAI_tAI_dataframe.txt -o /home/share/riboseq/c7-tAIPlot -u 0 -d 500 --mode all --start 5 --window 7 --step 1
tAIPlot -i /home/share/riboseq/c6-tAI_tAI_dataframe.txt -o /home/share/riboseq/c7-tAIPlot -u 0 -d 500 --mode all --start 5 --window 7 --step 1

echo "Finished: Local tRNA adaptation index and global tRNA adaptation index"
date

echo "Local codon adaptation index and global codon adaptation index"
date
# Local codon adaptation index and global codon adaptation index
echo cAI -i /home/share/riboseq/RiboMiner/transcript_cds_sequences_tAI.fa -o /home/share/riboseq/c8-cAI -t mouse -u 0 -d 500 --reference /home/share/riboseq/RiboMiner/reference.fa
cAI -i /home/share/riboseq/RiboMiner/transcript_cds_sequences_tAI.fa -o /home/share/riboseq/c8-cAI -t mouse -u 0 -d 500 --reference /home/share/riboseq/RiboMiner/reference.fa

echo cAIPlot -i /home/share/riboseq/c8-cAI_local_cAI_dataframe.txt -o /home/share/riboseq/c9-cAIPlot -u 0 -d 500 --mode all --start 5 --window 7 --step 1
cAIPlot -i /home/share/riboseq/c8-cAI_local_cAI_dataframe.txt -o /home/share/riboseq/c9-cAIPlot -u 0 -d 500 --mode all --start 5 --window 7 --step 1
echo GetProteinCodingSequence -i <transcrits_sequence.fa> -c <longest.trans.info.txt> -S /home/share/riboseq/selec_trans_longest.txt -o <output_prefix> --mode whole --table 1 --id-type transcript-id
GetProteinCodingSequence -i <transcrits_sequence.fa> -c <longest.trans.info.txt> -S /home/share/riboseq/selec_trans_longest.txt -o <output_prefix> --mode whole --table 1 --id-type transcript-id
echo "Finished: Local codon adaptation index and global codon adaptation index"
date
# Hydrophobicity calculation and Charge amino acids
```

**e.Enrichment_Analysis.sh**

```sh
echo RiboDensityAtEachPosition -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -f /home/share/riboseq/attributes.txt -o /home/share/riboseq/d5-RiboDensityAtEachPosition -U codon
#RiboDensityAtEachPosition -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -f /home/share/riboseq/attributes.txt -o /home/share/riboseq/d5-RiboDensityAtEachPosition -U codon


echo enrichmentMeanDensity -i /home/share/riboseq/d5-RiboDensityAtEachPosition_r7_cds_codon_density.txt,/home/share/riboseq/d5-RiboDensityAtEachPosition_r111_cds_codon_density.txt -o /home/share/riboseq/d6-enrichmentMeanDensity
#enrichmentMeanDensity -i /home/share/riboseq/d5-RiboDensityAtEachPosition_r7_cds_codon_density.txt,/home/share/riboseq/d5-RiboDensityAtEachPosition_r111_cds_codon_density.txt -o /home/share/riboseq/d6-enrichmentMeanDensity



echo EnrichmentAnalysis --ctrl /home/share/riboseq/d5-RiboDensityAtEachPosition_r7_cds_codon_density.txt --treat /home/share/riboseq/d5-RiboDensityAtEachPosition_r111_cds_codon_density.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/d7-EnrichmentAnalysis -U codon -M RPKM -l 150 -n 10 -m 1 -e 30 --CI 0.95 -u 0 -d 500
#EnrichmentAnalysis --ctrl /home/share/riboseq/d5-RiboDensityAtEachPosition_r7_cds_codon_density.txt --treat /home/share/riboseq/d5-RiboDensityAtEachPosition_r111_cds_codon_density.txt -c /home/share/riboseq/RiboMiner/longest.transcripts.info.txt -o /home/share/riboseq/d7-EnrichmentAnalysis -U codon -M RPKM -l 150 -n 10 -m 1 -e 30 --CI 0.95 -u 0 -d 500

echo PlotEnrichmentRatio -i /home/share/riboseq/d7-EnrichmentAnalysis_enrichment_dataframe.txt -o /home/share/riboseq/d8-PlotEnrichmentRatio -u 0 -d 500 --unit codon --mode all
PlotEnrichmentRatio -i /home/share/riboseq/d7-EnrichmentAnalysis_enrichment_dataframe.txt -o /home/share/riboseq/d8-PlotEnrichmentRatio -u 0 -d 500 --unit codon --mode all
```


