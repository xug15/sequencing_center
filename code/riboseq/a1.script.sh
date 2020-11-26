adapter=CTGTAGGCACCATCAAT

outdir=output
gtf=/home/xugang/data/reference/mouse/gencode.vM25.annotation.gtf
genome=/home/xugang/data/reference/mouse/gencode.fa
rRNA_bowtie=/home/xugang/data/reference/mouse/rRNA-bowtie/rRNA
star_index=/home/xugang/data/reference/mouse/star2
[[ -d $outdir ]] || mkdir -p $outdir
cutadapterf(){
[[ -d $outdir/a2-cutadapter/ ]] || mkdir $outdir/a2-cutadapter/
cutadapt -m 18 --match-read-wildcards -a ${adapter} -o ${outdir}/a2-cutadapter/${name}_trimmed.fastq  ${fastq} > ${outdir}/a2-cutadapter/log.${name}.log

}
#cutadapterf
fastq_filterf(){
[[ -d $outdir/a3-filter/ ]] || mkdir $outdir/a3-filter/
fastq_quality_filter \
        -Q33 -v -q 25 -p 75 \
        -i ${outdir}/a2-cutadapter/${name}_trimmed.fastq \
        -o ${outdir}/a3-filter/${name}_trimmedQfilter.fastq > ${outdir}/a3-filter/log.${name}.log
}

#fastq_filterf

fastqcf(){
[[ -d $outdir/a4-qc/ ]] || mkdir -p $outdir/a4-qc/
	fastqc \
        ${outdir}/a3-filter/${name}_trimmedQfilter.fastq \
        -o ${outdir}/a4-qc/ > ${outdir}/a4-qc/log.${name}.log
}
#fastqcf
rmrRNA(){
[[ -d $outdir/a5-rmrRNA/nonrRNA ]] || mkdir -p $outdir/a5-rmrRNA/nonrRNA
homedir=`pwd`
bowtie \
        -n 0 -norc --best -l 15 -p 8 \
         --un=${homedir}/${outdir}/a5-rmrRNA/nonrRNA/nocontam_${name} ${rRNA_bowtie} \
         -q ${outdir}/a3-filter/${name}_trimmedQfilter.fastq \
         ${outdir}/a5-rmrRNA/${name}.alin > \
         ${outdir}/a5-rmrRNA/${name}.err 2>${outdir}/a5-rmrRNA/${name}.err && \
         rm -rf ${outdir}/a5-rmrRNA/${name}.alin
}
#rmrRNA
starf(){
[[ -d $outdir/a6-map ]] || mkdir -p $outdir/a6-map

STAR --runThreadN 8 --alignEndsType EndToEnd \
         --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 \
         --genomeDir ${star_index} \
         --readFilesIn ${outdir}/a5-rmrRNA/nonrRNA/nocontam_${name} \
         --outFileNamePrefix ${outdir}/a6-map/${name} \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts

}
#starf

ribocodeannf(){
[[ -d $outdir/a7-ribocode_annotation ]] || mkdir -p $outdir/a7-ribocode_annotation


prepare_transcripts -g ${gtf} -f ${genome} -o ${outdir}/a7-ribocode_annotation

}
#ribocodeannf
metaplotf(){

[[ -d ${outdir}/a8-ribocode  ]] || mkdir -p ${outdir}/a8-ribocode
metaplots -a ${outdir}/a7-ribocode_annotation \
        -r ${outdir}/a6-map/${name}Aligned.toTranscriptome.out.bam -o ${outdir}/a8-ribocode/${name}
       
}
#metaplotf
ribocodef(){

mkdir -p /home/sfs/${JobName}/a9-ribocode-result && /root/miniconda3/bin/RiboCode -a /home/sfs/${JobName}/a7-ribocode_annotation -c /home/sfs/${JobName}/a8-ribocode/a_pre_config.txt -l no -g -o
}


run_one_step(){
cutadapterf
fastq_filterf
fastqcf
rmrRNA
starf
metaplotf

}
# 1run one sample
name=racr
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/1.rawdata/ribosome-1-AC-R_FKDL202613506-1a-AK412-AK4544/ribosome-1-AC-R_FKDL202613506-1a-AK412-AK4544_1.fq.gz
run_one_step

# 2run one sample
name=ract
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/2.cleandata/ribosome-1-AC-T_FKDL202613506-1a-AK2126-AK4544/ribosome-1-AC-T_FKDL202613506-1a-AK2126-AK4544_1.clean.fq.gz
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/1.rawdata/ribosome-1-AC-T_FKDL202613506-1a-AK2126-AK4544/ribosome-1-AC-T_FKDL202613506-1a-AK2126-AK4544_1.fq.gz
run_one_step

# 3run one sample
name=rsbr
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/2.cleandata/ribosome-1-SB-R_FKDL202613506-1a-AK410-AK4544/ribosome-1-SB-R_FKDL202613506-1a-AK410-AK4544_1.clean.fq.gz
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/1.rawdata/ribosome-1-SB-R_FKDL202613506-1a-AK410-AK4544/ribosome-1-SB-R_FKDL202613506-1a-AK410-AK4544_1.fq.gz
run_one_step

# 4run one sample
name=rsbt
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/2.cleandata/ribosome-1-SB-T_FKDL202613506-1a-AK2124-AK4544/ribosome-1-SB-T_FKDL202613506-1a-AK2124-AK4544_1.clean.fq.gz
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/1.rawdata/ribosome-1-SB-T_FKDL202613506-1a-AK2124-AK4544/ribosome-1-SB-T_FKDL202613506-1a-AK2124-AK4544_1.fq.gz
run_one_step

# 5run one sample
name=rac2r
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/2.cleandata/ribosome-2-AC-R_FKDL202613506-1a-AK411-AK4544/ribosome-2-AC-R_FKDL202613506-1a-AK411-AK4544_1.clean.fq.gz
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/1.rawdata/ribosome-2-AC-R_FKDL202613506-1a-AK411-AK4544/ribosome-2-AC-R_FKDL202613506-1a-AK411-AK4544_1.fq.gz
run_one_step

# 6run one sample
name=rac2t
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/2.cleandata/ribosome-2-AC-T_FKDL202613506-1a-AK2127-AK4544/ribosome-2-AC-T_FKDL202613506-1a-AK2127-AK4544_1.clean.fq.gz
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/1.rawdata/ribosome-2-AC-T_FKDL202613506-1a-AK2127-AK4544/ribosome-2-AC-T_FKDL202613506-1a-AK2127-AK4544_1.fq.gz
run_one_step

# 7run one sample
name=rsb2r
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/2.cleandata/ribosome-2-SB-R_FKDL202613506-1a-AK401-AK4544/ribosome-2-SB-R_FKDL202613506-1a-AK401-AK4544_1.clean.fq.gz
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/1.rawdata/ribosome-2-SB-R_FKDL202613506-1a-AK401-AK4544/ribosome-2-SB-R_FKDL202613506-1a-AK401-AK4544_1.fq.gz
run_one_step

# 8run one sample
name=rsb2t
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/2.cleandata/ribosome-2-SB-T_FKDL202613506-1a-AK398-AK4544/ribosome-2-SB-T_FKDL202613506-1a-AK398-AK4544_1.clean.fq.gz
fastq=/home/xugang/data/xiqiaoran/X101SC20102772-Z01-J013/1.rawdata/ribosome-2-SB-T_FKDL202613506-1a-AK398-AK4544/ribosome-2-SB-T_FKDL202613506-1a-AK398-AK4544_1.fq.gz
run_one_step



