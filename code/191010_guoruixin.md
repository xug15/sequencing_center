# Guo Runxin

## 
```sh
pip install --upgrade RiboMiner
```

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
config=b3-metaplot/config.txt
nohup RiboCode -a ${rico_ann} -c ${config} -l no -g -o b4-RiboCode/RiboCode_ORFs_result  > b4-RiboCode/ribocode.log 2>&1 &
``` 

**a5.density.sh**
```sh
#!/bin/bash
export PATH=~/anaconda3/bin:$PATH

mkdir b5-plot_density_orf
rico_ann=/Share/home/tiangeng/reference_genome/tair_ribocode/tair
config=b3-metaplot/config.txt
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

## RiboAnalyzer

### Star docker.
data

**x1.rundocker.sh**

```sh
#/Users/xugang/Documents/data

docker run -dt --name ribo --restart unless-stopped -v /Users/xugang/Documents/data:/data gangxu/ribodocker:1.3
```

### Prepare sequences and annotaiton files on transcriptome level.
```sh
prepare_transcripts -g /data/reference/tair/Arabidopsis_thaliana.TAIR10.43.gtf -f /data/reference/tair/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -o /data/reference/RiboCode_annot
```

### Prepare the longest transcript annotaion files.
a1-annotation.sh
```sh
OutputTranscriptInfo -c /data/reference/RiboCode_annot/transcripts_cds.txt -g /data/reference/tair/Arabidopsis_thaliana.TAIR10.43.gtf -f /data/reference/RiboCode_annot/transcripts_sequence.fa -o /data/reference/tair_analy/longest.transcripts.info.txt -O /data/reference/tair_analy/all.transcripts.info.txt
```

### Prepare the sequence file for the longest transcripts

a2-transcript.sh
```sh
GetProteinCodingSequence -i /data/reference/RiboCode_annot/transcripts_sequence.fa  -c /data/reference/tair_analy/longest.transcripts.info.txt -o /data/reference/tair_analy/transcript --mode whole --table 1 
```

### Prepare the UTR sequence

a3-utr.sh

```sh
GetUTRSequences -i /data/reference/tair_analy/transcript_transcript_sequences.fa -o /data/reference/tair_analy/utr -c /data/reference/tair_ribocode/tair/transcripts_cds.txt
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
PlotMetageneAnalysisForTheWholeRegions -i /data/data/a8-metagene_scaled_density_dataframe.txt -o /data/data/a9-meta_gene_whole_regin -g group1,group2 -r group1__group2 -b 15,90,60 --mode all 

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










