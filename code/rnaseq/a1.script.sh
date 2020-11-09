
adaptor=CTGTCTCTTATACACATCT
data_path=/home/xugang/data/yangyiyi/b1-data
cutadaptf(){

name=$1
[[ -d $data_path/$name  ]] ||  mkdir -p $data_path/$name

echo cutadapt -m 20 -a ${adaptor} -A ${adaptor} -o ${data_path}/$name/$name.clean1.fastq -p ${data_path}/$name/$name.clean2.fastq ${data_path}/$name/$name_R1.fq.gz ${data_path}/${name}/$name_R2.fq.gz

nohup cutadapt -m 20 -a ${adaptor} -A ${adaptor} -o ${data_path}/$name/$name.clean1.fastq -p ${data_path}/$name/$name.clean2.fastq ${data_path}/$name/${name}_R1.fq.gz ${data_path}/${name}/${name}_R2.fq.gz > ${data_path}/$name/$name-R1.log  2>&1 &

}

#cutadaptf YXXKO2-1
#cutadaptf YXXKO2-2 
#cutadaptf YXXKO2-3
#cutadaptf YXXKO2-4
#cutadaptf YXXKO3-1
#cutadaptf YXXKO3-2
#cutadaptf YXXKO3-3
#cutadaptf YXXKO3-4
#cutadaptf YXXKO8-1
#cutadaptf YXXKO8-2
#cutadaptf YXXKO8-3
#cutadaptf YXXKO8-4
#cutadaptf YXXWT10-1
#cutadaptf YXXWT10-2
#cutadaptf YXXWT10-3
#cutadaptf YXXWT10-4
#cutadaptf YXXWT4-1
#cutadaptf YXXWT4-2
#cutadaptf YXXWT4-3
#cutadaptf YXXWT4-4
#cutadaptf YXXWT9-1
#cutadaptf YXXWT9-2
#cutadaptf YXXWT9-3
#cutadaptf YXXWT9-4

rRNAf(){
cd /home/xugang/data/reference/mouse
grep rRNA Mus_musculus.GRCm38.100.gff3 > rRNA.gff
bedtools getfasta -fi Mus_musculus.GRCm38.dna.primary_assembly.fa -bed rRNA.gff >rRNA.fa
mkdir rRNA
mkdir rRNA-bowtie
#bowtie2-build rRNA.fa rRNA/rRNA
bowtie-build rRNA.fa rRNA-bowtie/rRNA
cd -
}
#rRNAf

rmrRNA()
{
name=$1
[[ -d $data_path/b2-rmrRNA/${name} ]] ||  mkdir -p $data_path/b2-rmrRNA/${name}
index=/home/xugang/data/reference/mouse/rRNA-bowtie/rRNA
ls
echo " bowtie -n 0 -norc --best -l 15 -p 38 ${index} -1 ${data_path}/$name/$name.clean1.fastq -2 ${data_path}/$name/$name.clean2.fastq  $data_path/b2-rmrRNA/${name}/${name}.alin --un $data_path/b2-rmrRNA/${name}/${name}.ready  > $data_path/b2-rmrRNA/${name}/${name}.err  "

bowtie -n 0 -norc --best -l 15 -p 38 ${index} -1 ${data_path}/$name/$name.clean1.fastq -2 ${data_path}/$name/$name.clean2.fastq  $data_path/b2-rmrRNA/${name}/${name}.alin --un $data_path/b2-rmrRNA/${name}/${name}.fq  > $data_path/b2-rmrRNA/${name}/${name}.err 2> $data_path/b2-rmrRNA/${name}/${name}.err
rm $data_path/b2-rmrRNA/${name}/${name}.alin

}

#rmrRNA YXXKO2-1 
#rmrRNA YXXKO2-2
#rmrRNA YXXKO2-3
#rmrRNA YXXKO2-4
#rmrRNA YXXKO3-1
#rmrRNA YXXKO3-2
#rmrRNA YXXKO3-3
#rmrRNA YXXKO3-3
#rmrRNA YXXKO3-4
#rmrRNA YXXKO8-1
#rmrRNA YXXKO8-2
#rmrRNA YXXKO8-3
#rmrRNA YXXKO8-4
#rmrRNA YXXWT10-1
#rmrRNA YXXWT10-2
#rmrRNA YXXWT10-3
#rmrRNA YXXWT10-4
#rmrRNA YXXWT4-1
#rmrRNA YXXWT4-2
#rmrRNA YXXWT4-3
#rmrRNA YXXWT4-4
#rmrRNA YXXWT9-1
#rmrRNA YXXWT9-2
#rmrRNA YXXWT9-3
#rmrRNA YXXWT9-4

starp(){
name=$1
[[ -d $data_path/b3-STAR/${name} ]] ||  mkdir -p $data_path/b3-STAR/${name}
genomeFile=/home/xugang/data/reference/mouse/star

echo STAR --runThreadN 8 --alignEndsType EndToEnd --outFilterMismatchNmax 2 --outFilterMultimapNmax 8 --genomeDir $genomeFile --readFilesIn $data_path/b2-rmrRNA/${name}/${name}_1.fq $data_path/b2-rmrRNA/${name}/${name}_2.fq  --outFileNamePrefix $data_path/b3-STAR/${name}/${name} --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts

STAR --runThreadN 18 --alignEndsType EndToEnd --outFilterMismatchNmax 2 --outFilterMultimapNmax 8 --genomeDir $genomeFile --readFilesIn $data_path/b2-rmrRNA/${name}/${name}_1.fq $data_path/b2-rmrRNA/${name}/${name}_2.fq  --outFileNamePrefix $data_path/b3-STAR/${name}/${name} --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts

}

#starp YXXKO2-1
#starp YXXKO2-2
#starp YXXKO2-3
#starp YXXKO2-4
#starp YXXKO3-1
#starp YXXKO3-2
#starp YXXKO3-3
#starp YXXKO3-3
#starp YXXKO3-4
#starp YXXKO8-1
#starp YXXKO8-2
#starp YXXKO8-3
#starp YXXKO8-4
#starp YXXWT10-1
#starp YXXWT10-2
#starp YXXWT10-3
#starp YXXWT10-4
#starp YXXWT4-1
#starp YXXWT4-2
#starp YXXWT4-3
#starp YXXWT4-4
#starp YXXWT9-1
#starp YXXWT9-2
#starp YXXWT9-3
#starp YXXWT9-4

gfoldf(){

name=$1
[[ -d $data_path/b4-gfold/ ]] || mkdir -p $data_path/b4-gfold/  
gtf=/home/xugang/data/reference/mouse/gencode.vM25.annotation.gtf
echo -e " samtools view $data_path/b3-STAR/${name}/${name}Aligned.sortedByCoord.out.bam| gfold count -ann ${gtf} -tag stdin -o $data_path/b3-gfold/${name}.read_cnt"
nohup samtools view $data_path/b3-STAR/${name}/${name}Aligned.sortedByCoord.out.bam | gfold count -ann ${gtf} -tag stdin -o $data_path/b4-gfold/${name}.read_cnt  >$data_path/b4-gfold/${name}.err 2>&1 &

}

#gfoldf YXXKO2-1
#gfoldf YXXKO2-2
#gfoldf YXXKO2-3
#gfoldf YXXKO2-4
#gfoldf YXXKO3-1
#gfoldf YXXKO3-2
#gfoldf YXXKO3-3
#gfoldf YXXKO3-4
#gfoldf YXXKO8-1
#gfoldf YXXKO8-2
#gfoldf YXXKO8-3
#gfoldf YXXKO8-4
#gfoldf YXXWT10-1
#gfoldf YXXWT10-2
#gfoldf YXXWT10-3
#gfoldf YXXWT10-4
#gfoldf YXXWT4-1
#gfoldf YXXWT4-2
#gfoldf YXXWT4-3
#gfoldf YXXWT4-4
#gfoldf YXXWT9-1
#gfoldf YXXWT9-2
#gfoldf YXXWT9-3
#gfoldf YXXWT9-4

mergef(){
for i in `ls $data_path/b4-gfold/ |grep .read_cnt$`;
do
file=`echo $i| cut -f 1 -d '.'`;
echo -e "gene\t$file">$data_path/b4-gfold/$i.txt ;
cut -f 1,3 $data_path/b4-gfold/$i >> $data_path/b4-gfold/$i.txt;
done
cd $data_path/b4-gfold/
file=(`ls $data_path/b4-gfold/ |grep .read_cnt.txt$`);
begin1=${file[0]};
cp $begin1 tmp;
echo $begin1;
file2=("${file[@]:1}");
echo ${file2[0]};
for i in ${file2[@]};
do echo $i;
join tmp $i > tmp2
mv tmp2 tmp
done
mv tmp merge.tsv
sed -i 's/ \+/\t/g' merge.tsv
cd -
}
#mergef
mergegf(){
for i in `ls $data_path/b4-gfold/ |grep .read_cnt$`;
do
file=`echo $i| cut -f 1 -d '.'`;
echo -e "gene\t$file">$data_path/b4-gfold/$i.txt ;
cut -f 2,3 $data_path/b4-gfold/$i >> $data_path/b4-gfold/$i.txt;
done
cd $data_path/b4-gfold/
file=(`ls $data_path/b4-gfold/ |grep .read_cnt.txt$`);
begin1=${file[0]};
cp $begin1 tmp;
echo $begin1;
file2=("${file[@]:1}");
echo ${file2[0]};
for i in ${file2[@]};
do echo $i;
join tmp $i > tmp2
mv tmp2 tmp
done
mv tmp merge.gene.tsv
sed -i 's/ \+/\t/g' merge.gene.tsv
cd -
}
#mergegf

exname(){
perl a3.name.pl b4-gfold/df.up.tsv
perl a3.name.pl b4-gfold/df.down.tsv
perl a3.name.pl b4-gfold/cpm.tsv

}
exname





