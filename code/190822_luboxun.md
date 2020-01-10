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

echo "bam must be sort and index."


#Periodicity checking
#Ribosome profiling data with a good quality tend to have a good 3-nt periodicity.
ribocode='/home/share/riboseq/Ribocode'

out='/home/share/riboseq/bqualitycontrol'
quality_out=$out
meta_out=$quality_out'/metaplot'
a111R='/home/share/riboseq/a9-STAR/7-111-R_STAR/7-111-R.toTranscriptome.sort.bam'
a111Rgenome='/home/share/riboseq/a9-STAR/7-111-R_STAR/7-111-RAligned.sortedByCoord.out.bam'
a111R_o=$meta_out'/7111r'
a7R='/home/share/riboseq/a9-STAR/7-7-R_STAR/7-7-R.toTranscriptome.sort.bam'
a7Rgenome='/home/share/riboseq/a9-STAR/7-7-R_STAR/7-7-RAligned.sortedByCoord.out.bam'
a7R_o=$meta_out'/77r'
long='/home/share/riboseq/RiboMiner/longest.transcripts.info.txt'
periodicity_out=$quality_out'/a5-periodicity'
a111R_periodicity_o=$periodicity_out'/7111r'
a7R_periodicity_o=$periodicity_out'/77r'

echo $meta_out

gtf=/home/share/riboseq/mus_ensemble/Mus_musculus.GRCm38.87.gtf

[ -d $out ] || mkdir $out

[ -d $out/metaplot ] || mkdir $out/metaplot

metaplot()
{
echo "metaplot begin `date`"
echo -e "metaplots -a $ribocode -r ${a111R} -o ${a111R_o}"
metaplots -a $ribocode -r ${a111R} -o ${a111R_o}
echo "Begin `date`"
echo -e "metaplots -a $ribocode -r ${a7R} -o ${a7R_o}"
metaplots -a $ribocode -r ${a7R} -o ${a7R_o}
echo "Finished `date`"
echo "Begin to generate attribution`date`"
echo rm ${meta_out}/attributes.txt
rm ${meta_out}/attributes.txt
echo "SampleName\tAlignmentFile\tStranded\tP-siteReadLength\tP-siteLocations" > ${meta_out}/attributes2.txt;
echo -e " for i in `ls ${meta_out}|grep txt$`;do grep -v "#" ${meta_out}/${i} >> ${meta_out}/attributes2.txt;done;"
for i in `ls ${meta_out}|grep config.txt$`;do grep -v "#" ${meta_out}/${i} >> ${meta_out}/attributes2.txt;done;
perl -pi -e 's/^\n//g' ${meta_out}/attributes2.txt
awk '{print $2"\t"$4"\t"$5"\t"$1}' ${meta_out}/attributes2.txt > ${meta_out}/attributes.txt
echo "metaplot begin `end`"
}

periodicity()
{
echo "Periodicity start `date`"
[ -d ${periodicity_out} ] || mkdir ${periodicity_out}
echo -e " Periodicity -i $a111R -a $ribocode -o $a111R_periodicity_o -c $long -L 25 -R 35"
Periodicity -i $a111R -a $ribocode -o $a111R_periodicity_o -c $long -L 25 -R 35
echo -e " Periodicity -i $a7R -a $ribocode -o $a7R_periodicity_o -c $long -L 25 -R 35"
Periodicity -i $a7R -a $ribocode -o $a7R_periodicity_o -c $long -L 25 -R 35
echo "Periodicity end `date`"
}

ribodensity()
{
echo "RiboDensityOfDiffFrames begin `date`"
#Reads distribution among different reading frames.
echo -e "RiboDensityOfDiffFrames -f ${meta_out}/attributes.txt -c $long -o $out/a6-ribo-density-diff-fram"
RiboDensityOfDiffFrames -f ${meta_out}/attributes.txt -c $long -o $out/a6-ribo-density-diff-fram
echo "RiboDensityOfDiffFrames end `date`"
}

length_dis()
{
#Length distribution.
echo "LengthDistribution begin `date`"
echo -e "LengthDistribution -i $a111Rgenome -o $quality_out/r111.length -f bam"
LengthDistribution -i $a111Rgenome -o $quality_out/r111.length -f bam
echo -e "LengthDistribution -i $a7Rgenome -o $quality_out/r7.length -f bam"
LengthDistribution -i $a7Rgenome -o $quality_out/r7.length -f bam
echo "LengthDistribution end `date`"
}

DNA_contam()
{
#DNA contamination.
echo "StatisticReadsOnDNAsContam begin `date`"
echo "StatisticReadsOnDNAsContam -i $a111Rgenome  -g $gtf -o $out/a7-dna-contamination.111"
StatisticReadsOnDNAsContam -i $a111Rgenome  -g $gtf -o $out/a7-dna-contamination.111
echo -e "StatisticReadsOnDNAsContam -i  $a7Rgenome  -g $gtf -o $out/a7-dna-contamination.7"
StatisticReadsOnDNAsContam -i  $a7Rgenome  -g $gtf -o $out/a7-dna-contamination.7
echo "StatisticReadsOnDNAsContam end `date`"

}


#metaplot
#periodicity
#ribodensity
#length_dis
#DNA_contam

```

**c.metagene_analysis.sh**


```sh

meta_out='/home/share/riboseq/bqualitycontrol/metaplot'
ribocode='/home/share/riboseq/RiboMiner'
long='/home/share/riboseq/RiboMiner/longest.transcripts.info.txt'
out='/home/share/riboseq/c_metagene_analysis'
groupinfo='R111,R7'
replace='7-111-R.toTranscriptome.sort__7-7-R.toTranscriptome.sort'



[ -d $out ] || mkdir -p $out

wholeregion(){
echo "Name must from MetageneAnalysisForTheWholeRegions output"
echo "Start Metagene analysis along the whole transcript region. `date`"

echo "MetageneAnalysisForTheWholeRegions -f ${meta_out}/attributes.txt -c $long -o $out/a8-metagene -b 15,90,60 -l 100 -n 10 -m 1 -e 5 --plot yes"
MetageneAnalysisForTheWholeRegions -f ${meta_out}/attributes.txt -c $long -o $out/a8-metagene -b 15,90,60 -l 100 -n 10 -m 1 -e 5 --plot yes

echo "PlotMetageneAnalysisForTheWholeRegions -i $out/a8-metagene_scaled_density_dataframe.txt -o $out/a9-meta_gene_whole_regin -g $groupinfo -r $replace -b 15,90,60 --mode all
"
PlotMetageneAnalysisForTheWholeRegions -i $out/a8-metagene_scaled_density_dataframe.txt -o $out/a9-meta_gene_whole_regin -g $groupinfo -r $replace -b 15,90,60 --mode all
echo "End Metagene analysis along the whole transcript region. `date`"
}

cdsregion()
{
echo "Begin Metagene analysis on CDS regions. `date`";
#Metagene analysis on CDS regions.
echo "MetageneAnalysis -f ${meta_out}/attributes.txt -c $long -o $out/b1-meat-cds -U codon -M RPKM -u 0 -d 500 -l 100 -n 10 -m 1 -e 5 --norm yes -y 100 --CI 0.95 --type CDS"
MetageneAnalysis -f ${meta_out}/attributes.txt -c $long -o $out/b1-meat-cds -U codon -M RPKM -u 0 -d 500 -l 100 -n 10 -m 1 -e 5 --norm yes -y 100 --CI 0.95 --type CDS
echo "End Metagene analysis on CDS regions. `date`";

}

utrregion(){

echo "Begin Metagene analysis on UTR regions. `date`"
#Metagene analysis on UTR regions.
echo " MetageneAnalysis -f ${meta_out}/attributes.txt -c $long -o $out/b2-meat-utr -U nt -M RPKM -u 100 -d 100 -l 100 -n 10 -m 1 -e 5 --norm yes -y 50 --CI 0.95 --type UTR"
MetageneAnalysis -f ${meta_out}/attributes.txt -c $long -o $out/b2-meat-utr -U nt -M RPKM -u 100 -d 100 -l 100 -n 10 -m 1 -e 5 --norm yes -y 50 --CI 0.95 --type UTR

echo "End Metagene analysis on UTR regions. `date`"
}

polarity()
{
echo "Begin Polarity calculation. `date`"
#Polarity calculation.
echo "PolarityCalculation -f ${meta_out}/attributes.txt -c $long -o $out/b3-polarity -n 64"
PolarityCalculation -f ${meta_out}/attributes.txt -c $long -o $out/b3-polarity -n 64
echo "End Polarity calculation. `date`"
}

wholeregion
#cdsregion
#utrregion
#polarity
```

**d1.feature_analysis.sh**

```sh
out='/home/share/riboseq/d_feature_analysis'
meta_out='/home/share/riboseq/bqualitycontrol/metaplot'

RiboMiner='/home/share/riboseq/RiboMiner'
long=$RiboMiner'/longest.transcripts.info.txt'
select_gene='/home/share/riboseq/selec_trans_longest.txt'

groupinfo='R111,R7'
replace='7-111-R.toTranscriptome.sort__7-7-R.toTranscriptome.sort'

[  -d $out ] ||  mkdir d_feature_analysis/

densityspecific()
{
echo "Begin:Pick out transcripts enriched ribosomes on specific region `date`"
# Pick out transcripts enriched ribosomes on specific region.
echo "RiboDensityForSpecificRegion -f $meta_out/attributes.txt -c $long -o $out/b5-transcript-enrich -U codon -M RPKM -L 25 -R 75"
RiboDensityForSpecificRegion -f $meta_out/attributes.txt -c $long -o $out/b5-transcript-enrich -U codon -M RPKM -L 25 -R 75
echo "End: Pick out transcripts enriched ribosomes on specific region `date`"
}

densitycodon()
{
echo "Begin: Ribosome density at each kind of AA or codon. `date`"
# Ribosome density at each kind of AA or codon.
echo " RiboDensityAtEachKindAAOrCodon -f $meta_out/attributes.txt -c $long -o $out/b6-ribosome-aa -M counts -S $select_gene -l 100 -n 10 --table 1 -F $RiboMiner/transcript_cds_sequences.fa "
RiboDensityAtEachKindAAOrCodon -f $meta_out/attributes.txt -c $long -o $out/b6-ribosome-aa -M counts -S $select_gene -l 100 -n 10 --table 1 -F $RiboMiner/transcript_cds_sequences.fa 
echo "PlotRiboDensityAtEachKindAAOrCodon -i $out/b6-ribosome-aa_all_codon_density.txt -o $out/b7-PlotRiboDensityAtEachKindAAOrCodon -g $groupinfo -r $replace --level AA "
PlotRiboDensityAtEachKindAAOrCodon -i $out/b6-ribosome-aa_all_codon_density.txt -o $out/b7-PlotRiboDensityAtEachKindAAOrCodon -g $groupinfo -r $replace --level AA 
echo "End: Ribosome density at each kind of AA or codon. `date`"
}



triplete()
{
echo "Begin: Pausing score of each triplete amino acid`date`"
# Pausing score of each triplete amino acid.
echo "PausingScore -f $meta_out/attributes.txt -c $RiboMiner/longest.transcripts.info.txt -o  $out/b8-PausingScore -M counts -S $select_gene  -l 100 -n 10 --table 1 -F  $RiboMiner/transcript_cds_sequences.fa"
PausingScore -f $meta_out/attributes.txt -c $RiboMiner/longest.transcripts.info.txt -o  $out/b8-PausingScore -M counts -S $select_gene  -l 100 -n 10 --table 1 -F  $RiboMiner/transcript_cds_sequences.fa

echo "ProcessPausingScore -i $out/b8-PausingScore_7-7-R.toTranscriptome.sort_pausing_score.txt,$out/b8-PausingScore_7-111-R.toTranscriptome.sort_pausing_score.txt -o $out/b9-ProcessPausingScore -g $groupinfo -r $replace --mode raw --ratio_filter 1 --pausing_score_filter 0.01"
ProcessPausingScore -i $out/b8-PausingScore_7-7-R.toTranscriptome.sort_pausing_score.txt,$out/b8-PausingScore_7-111-R.toTranscriptome.sort_pausing_score.txt -o $out/b9-ProcessPausingScore -g $groupinfo -r $replace --mode raw --ratio_filter 1 --pausing_score_filter 0.01

echo "End: Pausing score of each triplete amino acid`date`"
}

RPFdistf()
{
echo "Begin: RPFdist calculation.`date`"

# RPFdist calculation.
echo "RPFdist -f $meta_out/attributes.txt -c $RiboMiner/longest.transcripts.info.txt -o $out/c3-RPFdist -M counts -S $select_gene -l 100 -n 10 -m 1 -e 5"
RPFdist -f $meta_out/attributes.txt -c $RiboMiner/longest.transcripts.info.txt -o $out/c3-RPFdist -M counts -S $select_gene -l 100 -n 10 -m 1 -e 5

echo "End: RPFdist calculation. `date`"
}

gccontents()
{
echo "Begin:GC contents for sequences with a fasta format `date`"

# GC contents for sequences with a fasta format.

echo "GCContent -i $RiboMiner/transcript_cds_sequences.fa -o $out/c4-GCContent-normal --mode normal"
echo "GCContent -i $RiboMiner/transcript_cds_sequences.fa -o $out/c4-GCContent-frames --mode frames"
GCContent -i $RiboMiner/transcript_cds_sequences.fa -o $out/c4-GCContent-normal --mode normal
GCContent -i $RiboMiner/transcript_cds_sequences.fa -o $out/c4-GCContent-frames --mode frames
echo "PlotGCContent -i $out/c4-GCContent-normal_GC_content.txt -o $out/c5-PlotGCContent-normal --mode normal"
PlotGCContent -i $out/c4-GCContent-normal_GC_content.txt -o $out/c5-PlotGCContent-normal --mode normal
echo "End: GC contents for sequences with a fasta format`date`"
}



#densityspecific
#densitycodon
#triplete
#RPFdistf
#gccontents
```
triplete 有问题。

**d2.feature_analysis.sh**

We need tRNA file, that can be download from [GtRNAdb](http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/Mmusc10-gene-list.html)

[tair10 tRNA](./araTha1-tRNAs-confidence-set.txt)  
[mm10 tRNA](./mouse-tRNAs-confidence.txt)

|Chr|tRNA|Begin|End|Isotype|Anticodon|Upstream|Downstream|
| -| -| -| -|-| -|-|-|
|chr6	|95	|58141949	|58141877	|Ala	|AGC	|tttctccctc	|gtttcttgtc|
|chr6	|25	|26751918	|26751990	|Ala	|AGC	|agtgtagtgt	|gcttctttta|

```sh
out='/home/share/riboseq/d_feature_analysis'
meta_out='/home/share/riboseq/bqualitycontrol/metaplot'

RiboMiner='/home/share/riboseq/RiboMiner'
long=$RiboMiner'/longest.transcripts.info.txt'
select_gene='/home/share/riboseq/selec_trans_longest.txt'

groupinfo='R111,R7'
replace='7-111-R.toTranscriptome.sort__7-7-R.toTranscriptome.sort'

tRNA_confidence='/home/share/riboseq/mouse-tRNAs-confidence.txt'

tAIf()
{
echo "Local tRNA adaptation index and global tRNA adaptation index `date`"
#Local tRNA adaptation index and global tRNA adaptation index
echo " tAI -i $RiboMiner/transcript_cds_sequences.fa -t mouse -o $out/c6-tAI -u 0 -d 500 --table 1 -N $tRNA_confidence"
tAI -i $RiboMiner/transcript_cds_sequences.fa -t mouse -o $out/c6-tAI -u 0 -d 500 --table 1 -N $tRNA_confidence

echo " tAIPlot -i $out/c6-tAI_tAI_dataframe.txt -o $out/c7-tAIPlot -u 0 -d 500 --mode all --start 5 --window 7 --step 1"
tAIPlot -i $out/c6-tAI_tAI_dataframe.txt -o $out/c7-tAIPlot -u 0 -d 500 --mode all --start 5 --window 7 --step 1

echo "Finished: Local tRNA adaptation index and global tRNA adaptation index `date`"

}

cAIf()
{
echo "Local codon adaptation index and global codon adaptation index `date`"
# Local codon adaptation index and global codon adaptation index
echo "cAI -i $RiboMiner/transcript_cds_sequences.fa -o $out/c8-cAI -t mouse -u 0 -d 500 --reference $RiboMiner/transcript_cds_sequences.fa"
cAI -i $RiboMiner/transcript_cds_sequences.fa -o $out/c8-cAI -t mouse -u 0 -d 500 --reference $RiboMiner/transcript_cds_sequences.fa

echo "cAIPlot -i $out/c8-cAI_local_cAI_dataframe.txt -o $out/c9-cAIPlot -u 0 -d 500 --mode all --start 5 --window 7 --step 1"
cAIPlot -i $out/c8-cAI_local_cAI_dataframe.txt -o $out/c9-cAIPlot -u 0 -d 500 --mode all --start 5 --window 7 --step 1 1
echo "Finished: Local codon adaptation index and global codon adaptation index `date`"

}

tAIf
cAIf

```

**d3.feature_analysis.sh**

[A_charge_index](./AA_charge_index.txt)

[hydropathy_index](./hydropathy_index.txt)

```sh
out='/home/share/riboseq/d_feature_analysis'
meta_out='/home/share/riboseq/bqualitycontrol/metaplot'

RiboMiner='/home/share/riboseq/RiboMiner'
long=$RiboMiner'/longest.transcripts.info.txt'
select_gene='/home/share/riboseq/selec_trans_longest.txt'


hydropathy_index='/home/share/riboseq/hydropathy_index.txt'
charge_index='/home/share/riboseq/A_charge_index.txt'

hydropathycharge()
{
echo "Start hydrophobicity calculation `date`"
## hydrophobicity calculation
echo "hydropathyCharge  -i $RiboMiner/transcript_cds_sequences.fa -o $out/d1_hydropathy -t total_cds --index $hydropathy_index -u 0 -d 500 --table 1 "
hydropathyCharge  -i $RiboMiner/transcript_cds_sequences.fa -o $out/d1_hydropathy -t total_cds --index $hydropathy_index -u 0 -d 500 --table 1

cho "hydropathyCharge  -i $RiboMiner/transcript_cds_sequences.fa -o $out/d2_charge -t total_cds --index $charge_index -u 0 -d 500 --table 1 "
hydropathyCharge  -i $RiboMiner/transcript_cds_sequences.fa -o $out/d2_charge -t total_cds --index $charge_index -u 0 -d 500 --table 1


echo "End hydrophobicity calculation `date`"
}

hydropathplot(){
echo "Start plot hydrophobicit `date`"
## hydrophobicity
echo "PlotHydropathyCharge -i $out/d1_hydropathy_values_dataframe.txt -o $out/d3_hydropathy -u 0 -d 500 --mode all --ylab "Average Hydrophobicity" "
PlotHydropathyCharge -i $out/d1_hydropathy_values_dataframe.txt -o $out/d3_hydropathy -u 0 -d 500 --mode all --ylab "Average Hydrophobicity"
## charge
echo "PlotHydropathyCharge -i $out/d2_charge_values_dataframe.txt -o $out/d4_charge -u 0 -d 500 --mode all --ylab "Average Charges" "
PlotHydropathyCharge -i $out/d2_charge_values_dataframe.txt -o $out/d4_charge -u 0 -d 500 --mode all --ylab "Average Charges" 

echo "End plot hydrophobicit `date`"
}

hydropathycharge
hydropathplot
```


## charge calculation
hydropathyCharge -i $workdir/2954_up_cds_sequences.fa,$workdir/1598_unblocked_cds_sequences.fa,$workdir/433_down_cds_sequences.fa -t 2954_up,1598_unblocked,433_down -o $workdir/featureAnalysis/test.yeast_charge -u 0 -d 500 --index $workdir/featureAnalysis/AA_charge.txt


* The input file seperate by ',', and the file name in the -t seperate by ','.

```sh
## hydrophobicity calculation
hydropathyCharge  -i <cds_sequence_1.fa,cds_sequence_2.fa...> -o d1_hydropathy_dataframe.txt -t <geneList1,geneList2...> --index <hydrophobicity_index.txt> -u 0 -d 500 --table 1
##
hydropathyCharge  -i <cds_sequence_1.fa,cds_sequence_2.fa...> -o d2_charge_values_dataframe.txt -t <geneList1,geneList2...> --index <charge_index.txt> -u 0 -d 500 --table 1
```



**e.Enrichment_Analysis.sh**

```sh

out='/home/share/riboseq/e_enrichmetn'
RiboMiner='/home/share/riboseq/RiboMiner'
meta_out='/home/share/riboseq/bqualitycontrol/metaplot'


[ -d $out  ] || mkdir $out

ribodensityeach()
{
echo "Start RiboDensityAtEachPosition `date`"
echo " RiboDensityAtEachPosition -c $RiboMiner/longest.transcripts.info.txt -f $meta_out/attributes.txt -o $out/d5-RiboDensityAtEachPosition -U codon "
 RiboDensityAtEachPosition -c $RiboMiner/longest.transcripts.info.txt -f $meta_out/attributes.txt -o $out/d5-RiboDensityAtEachPosition -U codon
echo "End RiboDensityAtEachPosition `date`"
}

enrichmeandensity()
{
echo "Start enrichmentMeanDensity `date`"
echo " enrichmentMeanDensity -i $out/d5-RiboDensityAtEachPosition_7-7-R.toTranscriptome.sort_cds_codon_density.txt,$out/d5-RiboDensityAtEachPosition_7-111-R.toTranscriptome.sort_cds_codon_density.txt -o $out/d6-enrichmentMeanDensity "
enrichmentMeanDensity -i $out/d5-RiboDensityAtEachPosition_7-7-R.toTranscriptome.sort_cds_codon_density.txt,$out/d5-RiboDensityAtEachPosition_7-111-R.toTranscriptome.sort_cds_codon_density.txt -o $out/d6-enrichmentMeanDensity
echo "End enrichmentMeanDensity `date`"
}

enrichanalysis()
{
    echo "Start EnrichmentAnalysis `date` "    
echo " EnrichmentAnalysis --ctrl $out/d5-RiboDensityAtEachPosition_7-7-R.toTranscriptome.sort_cds_codon_density.txt --treat $out/d5-RiboDensityAtEachPosition_7-111-R.toTranscriptome.sort_cds_codon_density.txt -c $RiboMiner/longest.transcripts.info.txt -o $out/d7-EnrichmentAnalysis -U codon -M RPKM -l 150 -n 10 -m 1 -e 30 --CI 0.95 -u 0 -d 500 "
EnrichmentAnalysis --ctrl $out/d5-RiboDensityAtEachPosition_7-7-R.toTranscriptome.sort_cds_codon_density.txt --treat $out/d5-RiboDensityAtEachPosition_7-111-R.toTranscriptome.sort_cds_codon_density.txt -c $RiboMiner/longest.transcripts.info.txt -o $out/d7-EnrichmentAnalysis -U codon -M RPKM -l 150 -n 10 -m 1 -e 30 --CI 0.95 -u 0 -d 500

echo "PlotEnrichmentRatio -i $out/d7-EnrichmentAnalysis_enrichment_dataframe.txt -o $out/d8-PlotEnrichmentRatio -u 0 -d 500 --unit codon --mode all"
PlotEnrichmentRatio -i $out/d7-EnrichmentAnalysis_enrichment_dataframe.txt -o $out/d8-PlotEnrichmentRatio -u 0 -d 500 --unit codon --mode all
echo "End EnrichmentAnalysis `date`"
}

ribodensityeach
enrichmeandensity
enrichanalysis
```


