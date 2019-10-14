# Luboxun

## step1 : unzip the data.

a1.unzip.sh  

```sh
gunzip Col-0_FKDL171663936-1A_1.clean.fq.gz
gunzip d14_FKDL171663937-1A_1.clean.fq.gz
```

a2.rename.sh
```sh
mv Col-0_FKDL171663936-1A_1.clean.fq col.fq
mv d14_FKDL171663937-1A_1.clean.fq d14.fq
```

## Step2 :remove 5' G*

a1.reG.sh

```sh
rpf=(col d14)

for i in ${rpf[@]}
do
        echo "perl a2-removeG.pl ../a1-raw/${i}.fq";
        nohup perl a2-removeG.pl ../a1-raw/${i}.fq > ${1}.rmg.log 2>&1 &
done
```

a2-removeG.pl

```pl
 open DATA, "<$ARGV[0]";
 open OUT, ">$ARGV[0].clean.fq";
 while(<DATA>){
     $seq=<DATA>;
     $qua=<DATA>;
     $quality=<DATA>;
     chomp($seq);
     chomp($quality);
             if($seq=~/^(G+)(.*)/){
             $len_g=length($1);
             $len_f=length($seq);
             $len_q=$len_f-$len_g;
             $start=$len_g;
             $seq=$2;
             $quality=substr $quality, $start, $len_q;
             }
     print OUT "$_";
     print OUT "$seq\n";
     print OUT "$qua";
     print OUT "$quality\n";
 }
```

a3.rename.sh

```sh
rpf=(col d14)

for i in ${rpf[@]}
do
        mv ../a1-raw/${i}.fq.clean.fq ${i}.reg.fq
done 
```

## step 3: remove ploy A

a1.rmploya.sh

```sh
export PATH=$PATH:/Share/home/tiangeng/software/cap-miRNA/bin
rpf=(col d14)
adapt=AAAAAA

for i in ${rpf[@]}
do
        echo "nohup cutadapt -m 25 -M 35 --match-read-wildcards -a $adapt -o ${i}_trimmed.fastq ../a2-removeG/${i}.reg.fq > ${i}_rmploya.log 2>&1 &";
        nohup cutadapt -m 25 -M 35 --match-read-wildcards -a $adapt -o ${i}_trimmed.fastq ../a2-removeG/${i}.reg.fq > ${i}_rmploya.log 2>&1 &
done
```

## step 4: remove low quality reads.

```sh
mkdir a4-filter-lowquality
```
a1.fliter.sh

```sh
#!/bin/bash
rpf=(col d14)
for i in ${rpf[@]}
do
        echo "nohup /Share/app/fastx_toolkit0.14/bin/fastq_quality_filter -Q33 -v -q 25 -p 75 -i ../a3-rmploya/${i}_trimmed.fastq -o ${i}_clean.fastq > ${i}_Qfilter.log 2>&1 &";
        nohup /Share/app/fastx_toolkit0.14/bin/fastq_quality_filter -Q33 -v -q 25 -p 75 -i ../a3-rmploya/${i}_trimmed.fastq -o ${i}_clean.fastq > ${i}_Qfilter.log 2>&1 & 
done
```

## Step 5 : QC

**a1.qc.sh**

```sh
#!/bin/bash
rpf=(col d14)
for i in ${rpf[@]}
do
        echo -e " nohup /Share/home/tiangeng/software/FastQC/fastqc ../a4-filter-lowquality/${i}_clean.fastq -o ./ > ${i}.log 2>&1 &";
        nohup /Share/home/tiangeng/software/FastQC/fastqc ../a4-filter-lowquality/${i}_clean.fastq -o ./ > ${i}.log 2>&1 &
done
```


## Step 6 : contam (remove ribosome sequence)

mkdir a6-contam

> a6-contam
**a1.contam.sh**
```sh
#!/bin/bash

bowtieindex=/Share/home/tiangeng/reference_genome/tair_rRNA_bowtie/tair.rRNA
name=(col d14)

for i in ${name[@]}
do
        echo "nohup bowtie -n 0 -norc --best -l 15 -p 7 --un=nocontam_${i}.fastq $bowtieindex -q ../a4-filter-lowquality/${i}_clean.fastq ${i}.alin > ${i}.err 2>&1 & ";
        nohup bowtie -n 0 -norc --best -l 15 -p 7 --un=nocontam_${i}.fastq $bowtieindex -q ../a4-filter-lowquality/${i}_clean.fastq ${i}.alin > ${i}.err 2>&1 & 
done
```

**a2.merge.sh**
```sh
name=(col d14)
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

## Step7: afterQC (overrepresent reads)

mkdir a7-afterQC

> a7-afterqc
**a1.afterqc.sh**
```sh
#!/bin/bash

rpf=(col d14)
for i in ${rpf[@]}
do
        nohup /Share/home/tiangeng/software/FastQC/fastqc ../a6-contam/nocontam_${i}.fastq -o ./ >${i}.log 2>&1 & 
done
```


## Step8 : STAR + periodicity

mkdir a8-star

> a8-STAR

run with screen.

**a1.STAR.sh**
```sh
#!/bin/bash
export PATH=$PATH:~/software/STAR-master/bin/Linux_x86_64_static
genomeFile=/Share/home/tiangeng/reference_genome/tair_star
fileName=(col d14)
for i in ${fileName[@]}
do
        mkdir -p ${i}_STAR
        cd ${i}_STAR
        STAR --runThreadN 8 --alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 --genomeDir $genomeFile --readFilesIn ../../a6-contam/nocontam_${i}.fastq --outFileNamePrefix $i --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts
        cd ../
done
```
a2.merge.sh
```sh
name=(col d14)
head='Iterm'
for i in ${name[@]};
do
 head+="\t${i}";
 done
 echo -e $head >merge.counter;

 for i in ${name[@]};
 do
     echo "${i}_STAR/${i}Log.final.out";

 cat ${i}_STAR/${i}Log.final.out > ${i}.err.tmp;
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

 cat merge.counter merge.tmp > merge2.tmp;
mv merge2.tmp summary.txt
rm *.err.tmp
rm merge.counter merge.tmp

sed -i 's/^[ ]*//g' summary.txt;
sed -i 's/|//g' summary.txt;

```

## step 9 : map statatistic.

mkdir a9-statistic

a1.readsnum.sh

```sh
#!/bin/bash
export PATH=/Share/home/tiangeng/anaconda2/bin:$PATH
name=(col d14)
gtf=/Share/home/tiangeng/reference_genome/tair/Arabidopsis_thaliana.TAIR10.43.gtf
for i in ${name[@]}
do
       nohup python /Share/home/tiangeng/apps/readsNumCal_intron_v3.py ../a8-star/${i}_STAR/${i}Aligned.sortedByCoord.out.bam $gtf nocontam_${i}_mappedNum_intron.txt nocontam_${i} > read.${i}.log 2>&1 &
done

```

**a2.merge.sh**

```sh
name=(col d14)
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


## Step10: readlength

mkdir b1-readlength

**a1.readlength_filter.sh**
```sh
#!/bin/bash
export PATH=/Share/home/tiangeng/anaconda2/bin:$PATH
rpf=(col d14)

for i in ${rpf[@]}
do
        echo "nohup python /Share/home/tiangeng/apps/fq_len_stat.py ../a6-contam/nocontam_${i}.fastq ${i}_nontam.png > ${i}.nontam.out 2>&1 &";
        nohup python /Share/home/tiangeng/apps/fq_len_stat.py ../a6-contam/nocontam_${i}.fastq ${i}_nontam.png > ${i}.nontam.out 2>&1 &
done
```


## Step 11: calculate count

mkdir b2-count

**a1.rna_counter.sh**
```sh
#!/bin/bash
export PATH=/Share/home/tiangeng/anaconda2/bin:$PATH
name=(col d14)
gtf=/Share/home/tiangeng/reference_genome/tair/Arabidopsis_thaliana.TAIR10.43.gtf
data_p=../a8-star
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
gtf=/Share/home/tiangeng/reference_genome/tair/Arabidopsis_thaliana.TAIR10.43.gtf
data_p=../a8-star
name=(col d14)
for i in ${name[@]}
do 
echo -e "python /Share/home/tiangeng/apps/RPF_count_CDS.py  $data_p/${i}_STAR/${i}Aligned.sortedByCoord.out.bam ${gtf} > ${i}.count";
nohup python /Share/home/tiangeng/apps/RPF_count_CDS.py  $data_p/${i}_STAR/${i}Aligned.sortedByCoord.out.bam ${gtf} > ${i}.count 2>${i}.log &
done;
```
**a3.merge.sh**
```sh
name=(col d14)
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
 
mkdir b3-ribocode-anno

**a1.trans_ann.sh**
```sh
#!/bin/bash
export PATH=~/anaconda3/bin:$PATH
gtf=/Share/home/tiangeng/reference_genome/tair/Arabidopsis_thaliana.TAIR10.43.gtf
fa=/Share/home/tiangeng/reference_genome/tair/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
prepare_transcripts -g $gtf -f $fa -o /Share/home/tiangeng/reference_genome/tair_ribocode/tair
```
**a2.mergefile.sh**
```sh
#!/bin/bash
name=(col d14)
rm b1.merge_transcriptome.txt
for i in ${name[@]};
do 
echo /Share/home/tiangeng/project_result/Riboseq/project_190920_guoruixin/a8-star/${i}_STAR/${i}Aligned.toTranscriptome.out.bam >>b1.merge_transcriptome.txt;
done

```

**a3.metaplot.sh**
```sh
#!/bin/bash
export PATH=~/anaconda3/bin:$PATH
mkdir b3-metaplot
rico_ann=/Share/home/tiangeng/reference_genome/tair_ribocode/tair
nohup metaplots -a $rico_ann -i b1.merge_transcriptome.txt -o b3-metaplot/meta > b3-metaplot/metaplot.log 2>&1 &
```

**a4.Ribocode.sh**
```sh
#!/bin/bash
export PATH=~/anaconda3/bin:$PATH
mkdir b4-RiboCode
rico_ann=/Share/home/tiangeng/reference_genome/tair_ribocode/tair
config=b3-metaplot/meta_pre_config.txt
nohup RiboCode -a ${rico_ann} -c ${config} -l no -g -o b4-RiboCode/RiboCode_ORFs_result  > b4-RiboCode/ribocode.log 2>&1 &
``` 

**a5.density.sh**
```sh
#!/bin/bash
export PATH=~/anaconda3/bin:$PATH

mkdir b5-plot_density_orf
rico_ann=/Share/home/tiangeng/reference_genome/tair_ribocode/tair
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
name=(col d14)

for i in ${name[@]};
do 
nohup ORFcount -g b4-RiboCode/RiboCode_ORFs_result.gtf -r /Share/home/tiangeng/project_result/Riboseq/project_190920_guoruixin/a8-star/${i}_STAR/${i}Aligned.sortedByCoord.out.bam -f 15 -l 5 -e 100 -m 26 -M 34 -o b6-ORF_count/${i}.ORF.counts > b6-ORF_count/${i}.log 2>&1 &
done;

```
