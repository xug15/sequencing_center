# GCS 流程设计

-----------------------------

## [Huawei Homepage](20200118_huawei.md)

## [Huawei Docker images](20200118_huawei_dockerimage.md)

## [Huawei 文件拷贝](20200118_huawei_copy.md)

## [Huawei GCS 流程设计](20200118_huawei_gcspipeline.md)

-----------------------------

## Content


* [进行数据拷贝](#进行数据拷贝) 
* [GCS 流程示意原理](#GCS流程示意原理)
* [两种流程执行的方式](#两种流程执行的方式)
* [文件挂载](#文件挂载)
* [流程描述](#流程描述)
* [RNA-seq toal](#RNA-seq-toal)
* [RNA-seq small](#RNA-seq-small)
* [RNA-seq ploy A](#RNA-seq-ploy-A)
* [Ribo-seq-Xtail](#Ribo-seq-Xtail) 
* [Ribo-seq Ribocode](#Ribo-seq-Ribocode) 
* [bash-with-parameters](#bash-with-parameters)
* [对字符串进行分割成数组](#对字符串进行分割成数组) 





## 进行数据拷贝
```sh
obsutil config -i=5ULAGR0CWKBAEDV57Y6P -k=gvroYZE9uUmp3igpEPAEQRfuQzUjcVQn9kBoHz02 -e=https://obs.cn-north-4.myhuaweicloud.com && obsutil cp -r -f -u obs://gene-container-xugang/gcs/huawei_file/a1-fastq/ /home/sfs && obsutil cp -r -f -u obs://gene-container-xugang/gcs/huawei_file/refrence/ /home/sfs && ls /home/sfs
```
https://support.huaweicloud.com/tr-gcs/gcs_tr_04_0002.html

## GCS流程示意原理
```sh
job-a:
  commands_iter:
    command: echo ${1} ${item}
    vars_iter:
      - [A, B, C]
```

```sh
如下示例中，commands有四行，则表示容器并发数量为4，每个容器分别执行不同的命令

commands:
  - sh /obs/gcscli/run-xxx/run.sh 1 a
  - sh /obs/gcscli/run-xxx/run.sh 2 a
  - sh /obs/gcscli/run-xxx/run.sh 1 b
  - sh /obs/gcscli/run-xxx/run.sh 2 b

  如果命令行是由多行组成，可以使用yaml语法中的“|”（保留换行符，整个字符串当做yaml中一个key的value）格式。这样就可以把大篇幅的命令行原封不动的拷贝过来，如：

commands:
  - |
    samtools merge -f -@ ${nthread} -b ${volume-path-sfs}/${sample}/mergelist.txt \
    ${volume-path-sfs}/${sample}/${sample}.sort.bam && \
    samtools flagstat ${volume-path-sfs}/${sample}/${sample}.sort.bam > ${volume-path-sfs}/${sample}/${sample}.sort.flagstat
```

## 两种流程执行的方式

看一下创建自定义流程的第5点
support.huaweicloud.com/bestpractice-gcs/gcs_bestpractice_001.html
看一下GCS_DATA_PVC这个变量的使用
support.huaweicloud.com/tr-gcs/gcs_tr_04_0004.html

```sh
   commands_iter:
      command: |
        bwa mem -t ${nthread} -M -R "@RG\tID:Sample\tPL:illumina\tSM:${sample}\tCN:GATK4" \
                ${volume-path-ref}/${reference-path}/${fastafile} \
                ${volume-path-sfs}/${sample}/${1}/R0.${1}.fastq  \
                ${volume-path-sfs}/${sample}/${1}/R1.${1}.fastq |\
        samtools view -F 4 -q 10 -bS /dev/stdin \
                >${volume-path-sfs}/${sample}/${sample}.${1}.bam && \
        samtools sort -@ ${nthread} \
                -o ${volume-path-sfs}/${sample}/${sample}.${1}.sort.bam \
                ${volume-path-sfs}/${sample}/${sample}.${1}.bam && \
        samtools index ${volume-path-sfs}/${sample}/${sample}.${1}.sort.bam
      vars_iter:
        - 'range(0, ${npart})'

    commands_iter:
      command: |
        if [ ! -f  ${volume-path-obs}/${1} ]; then echo "File ${volume-path-obs}/${1} not found" && exit 1; fi && \
        for i in `seq 0 ${npart}`; do mkdir -p ${volume-path-sfs}/${sample}/$i; done && \
        zlibfq -b 100000 -t ${npart} ${volume-path-obs}/${1} ${volume-path-sfs}/${sample} R${item}
      vars_iter:
        - '${fastq-files}'
```

## 流程描述
    description: '将数据拷贝到制定目录下/home/sfs'
    description: 将原始的下机数据去除接头
    description: 去除低质量的数据    
    description: 对数据进行质量评估
    description: 去除核糖体的reads
    description: 比对到基因组上
    description: 将结果拷贝回obsvolumn中

## RNA-seq toal
```sh

# cutadapter
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a2-cutadapter && \
        echo ${1} begin `date` /root/miniconda3/bin/cutadapt -m 18 \
        --match-read-wildcards -a ${adapter} \
        -o /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
        /home/obs/${obs_data_path}/${1} && \
        /root/miniconda3/bin/cutadapt -m 18 \
        --match-read-wildcards -a ${adapter} \
        -o /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
         /home/obs/${obs_data_path}/${1}
      vars_iter:
        - '${fastq_files}'

# Fastx
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a3-filter && \
        echo ${1} begin `date` && \
        /home/test/bin/fastq_quality_filter \
        -Q33 -v -q 25 -p 75 \
        -i /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
        -o /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq 
      vars_iter:
        - '${fastq_files}'

# fastQC
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a4-qc && \
        echo ${1} begin `date` && \
        /home/test/FastQC/fastqc \
        /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq \
        -o /home/sfs/${JobName}/a4-qc
      vars_iter:
        - '${fastq_files}'

# remove rRNA
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a5-rmrRNA && \
        mkdir -p /home/sfs/${JobName}/a5-rmrRNA/nonrRNA && \
        echo ${1} begin `date` && \
        bash /root/.bashrc && \
        /home/test/bowtie-1.2.3-linux-x86_64/bowtie \
        -n 0 -norc --best -l 15 -p 8 \
         --un=/home/sfs/${JobName}/a5-rmrRNA/nonrRNA/nocontam_${1} /home/obs/${obs_reference_rRNA_bowtie} \
         -q /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq \
         /home/sfs/${JobName}/a5-rmrRNA/${1}.alin > \
         /home/sfs/${JobName}/a5-rmrRNA/${1}.err && \
         rm -rf /home/sfs/${JobName}/a5-rmrRNA/${1}.alin 
         
      vars_iter:
        - '${fastq_files}'

# mapping
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a6-map && \
        echo ${1} begin `date` && \
        bash /root/.bashrc && \
         STAR --runThreadN 8 --alignEndsType EndToEnd \
         --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 \
         --genomeDir /home/obs/${obs_reference_genomeFile_star} \
         --readFilesIn /home/sfs/${JobName}/a5-rmrRNA/nonrRNA/nocontam_${1} \
         --outFileNamePrefix /home/sfs/${JobName}/a6-map/${1} \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts
         
      vars_iter:
        - '${fastq_files}'

# copy
  commands:
      - 'obsutil config -i=${AccessKey} -k=${SecretKey} -e=${endpoint}&& obsutil mkdir -p ${obs_location}/output/${JobName}/  && obsutil cp -r -f /home/sfs/${JobName}/ ${obs_location}/output/${JobName}/ && rm -rf /home/sfs/${JobName} && echo Check sfs && ls -al /home/sfs && ls -al /home/obs/output'
```

## RNA-seq small 
```sh
# cutadapter
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a2-cutadapter && \
        echo ${1} begin `date` /root/miniconda3/bin/cutadapt -m 18 \
        --match-read-wildcards -a ${adapter} \
        -o /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
        /home/obs/${obs_data_path}/${1} && \
        /root/miniconda3/bin/cutadapt -m 18 \
        --match-read-wildcards -a ${adapter} \
        -o /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
         /home/obs/${obs_data_path}/${1}
      vars_iter:
        - '${fastq_files}'


# fastx
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a3-filter && \
        echo ${1} begin `date` && \
        /home/test/bin/fastq_quality_filter \
        -Q33 -v -q 25 -p 75 \
        -i /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
        -o /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq 
      vars_iter:
        - '${fastq_files}'

# fastQC
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a4-qc && \
        echo ${1} begin `date` && \
        /home/test/FastQC/fastqc \
        /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq \
        -o /home/sfs/${JobName}/a4-qc
      vars_iter:
        - '${fastq_files}'

# remove rRNA
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a5-rmrRNA && \
        mkdir -p /home/sfs/${JobName}/a5-rmrRNA/nonrRNA && \
        echo ${1} begin `date` && \
        bash /root/.bashrc && \
        /home/test/bowtie-1.2.3-linux-x86_64/bowtie \
        -n 0 -norc --best -l 15 -p 8 \
         --un=/home/sfs/${JobName}/a5-rmrRNA/nonrRNA/nocontam_${1} /home/obs/${obs_reference_rRNA_bowtie} \
         -q /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq \
         /home/sfs/${JobName}/a5-rmrRNA/${1}.alin > \
         /home/sfs/${JobName}/a5-rmrRNA/${1}.err && \
         rm -rf /home/sfs/${JobName}/a5-rmrRNA/${1}.alin 
         
      vars_iter:
        - '${fastq_files}'

# bowtie mature
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a7-mature && \
        echo ${1} begin `date` && \
        bash /root/.bashrc && \
        /home/test/bowtie-1.2.3-linux-x86_64/bowtie \
        -n 0 -norc --best -l 15 -p 8 \
        /home/obs/${obs_reference_mature_bowtie} \
         -q /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq \
         /home/sfs/${JobName}/a7-mature/${1}.alin > \
         /home/sfs/${JobName}/a7-mature/${1}.err

         
      vars_iter:
        - '${fastq_files}'

# bowtie hairpin
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a8-hairpin && \
        echo ${1} begin `date` && \
        bash /root/.bashrc && \
        /home/test/bowtie-1.2.3-linux-x86_64/bowtie \
        -n 0 -norc --best -l 15 -p 8 \
         /home/obs/${obs_reference_hairpin_bowtie} \
         -q /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq \
         /home/sfs/${JobName}/a8-hairpin/${1}.alin > \
         /home/sfs/${JobName}/a8-hairpin/${1}.err
         
      vars_iter:
        - '${fastq_files}'

#  STAR
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a6-map && \
        echo ${1} begin `date` && \
        bash /root/.bashrc && \
         STAR --runThreadN 8 --alignEndsType EndToEnd \
         --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 \
         --genomeDir /home/obs/${obs_reference_genomeFile_star} \
         --readFilesIn /home/sfs/${JobName}/a5-rmrRNA/nonrRNA/nocontam_${1} \
         --outFileNamePrefix /home/sfs/${JobName}/a6-map/${1} \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts
         
      vars_iter:
        - '${fastq_files}'

# copy
    commands:
      - 'obsutil config -i=${AccessKey} -k=${SecretKey} -e=${endpoint}&& obsutil mkdir -p ${obs_location}/output/${JobName}/  && obsutil cp -r -f /home/sfs/${JobName}/ ${obs_location}/output/${JobName}/ && rm -rf /home/sfs/${JobName} && echo Check sfs && ls -al /home/sfs && ls -al /home/obs/output'

```

## RNA-seq ploy A
```sh
# cutadapter
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a2-cutadapter && \
        echo ${1} begin `date` /root/miniconda3/bin/cutadapt -m 18 \
        --match-read-wildcards -a ${adapter} \
        -o /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
        /home/obs/${obs_data_path}/${1} && \
        /root/miniconda3/bin/cutadapt -m 18 \
        --match-read-wildcards -a ${adapter} \
        -o /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
         /home/obs/${obs_data_path}/${1}
      vars_iter:
        - '${fastq_files}'

# fastx
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a3-filter && \
        echo ${1} begin `date` && \
        /home/test/bin/fastq_quality_filter \
        -Q33 -v -q 25 -p 75 \
        -i /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
        -o /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq 
      vars_iter:
        - '${fastq_files}'

# fastqc
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a4-qc && \
        echo ${1} begin `date` && \
        /home/test/FastQC/fastqc \
        /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq \
        -o /home/sfs/${JobName}/a4-qc
      vars_iter:
        - '${fastq_files}'

# remove rRNA
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a5-rmrRNA && \
        mkdir -p /home/sfs/${JobName}/a5-rmrRNA/nonrRNA && \
        echo ${1} begin `date` && \
        bash /root/.bashrc && \
        /home/test/bowtie-1.2.3-linux-x86_64/bowtie \
        -n 0 -norc --best -l 15 -p 8 \
         --un=/home/sfs/${JobName}/a5-rmrRNA/nonrRNA/nocontam_${1} /home/obs/${obs_reference_rRNA_bowtie} \
         -q /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq \
         /home/sfs/${JobName}/a5-rmrRNA/${1}.alin > \
         /home/sfs/${JobName}/a5-rmrRNA/${1}.err && \
         rm -rf /home/sfs/${JobName}/a5-rmrRNA/${1}.alin 
         
      vars_iter:
        - '${fastq_files}'

# STAR
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a6-map && \
        echo ${1} begin `date` && \
        bash /root/.bashrc && \
         STAR --runThreadN 8 --alignEndsType EndToEnd \
         --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 \
         --genomeDir /home/obs/${obs_reference_genomeFile_star} \
         --readFilesIn /home/sfs/${JobName}/a5-rmrRNA/nonrRNA/nocontam_${1} \
         --outFileNamePrefix /home/sfs/${JobName}/a6-map/${1} \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts
         
      vars_iter:
        - '${fastq_files}'

# copy
    commands:
      - 'obsutil config -i=${AccessKey} -k=${SecretKey} -e=${endpoint}&& obsutil mkdir -p ${obs_location}/output/${JobName}/  && obsutil cp -r -f /home/sfs/${JobName}/ ${obs_location}/output/${JobName}/ && rm -rf /home/sfs/${JobName} && echo Check sfs && ls -al /home/sfs && ls -al /home/obs/output'

```

## Ribo-seq Xtail

```sh
#cutadapter
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a2-cutadapter && \
        ls /home/obs/${obs_data_path}/ && \
        echo ${1} begin `date` /root/miniconda3/bin/cutadapt -m 18 \
        --match-read-wildcards -a ${adapter} \
        -o /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
        /home/sfs/${JobName}/a1-fastq/${1} && \
        /root/miniconda3/bin/cutadapt -m 18 \
        --match-read-wildcards -a ${adapter} \
        -o /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
         /home/obs/${obs_data_path}/${1}
    vars_iter:
        - '${fastq_files}'

# fastq_quality_filter
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a3-filter && \
        echo ${1} begin `date` && \
        /home/test/bin/fastq_quality_filter \
        -Q33 -v -q 25 -p 75 \
        -i /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
        -o /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq 
      vars_iter:
        - '${fastq_files}'

# fastQC
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a4-qc && \
        echo ${1} begin `date` && \
        /home/test/FastQC/fastqc \
        /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq \
        -o /home/sfs/${JobName}/a4-qc
      vars_iter:
        - '${fastq_files}'

#remove rRNA
            commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a5-rmrRNA && \
        mkdir -p /home/sfs/${JobName}/a5-rmrRNA/nonrRNA && \
        echo ${1} begin `date` && \
        bash /root/.bashrc && \
        /home/test/bowtie-1.2.3-linux-x86_64/bowtie \
        -n 0 -norc --best -l 15 -p 8 \
         --un=/home/sfs/${JobName}/a5-rmrRNA/nonrRNA/nocontam_${1} /home/obs/${obs_reference_rRNA_bowtie} \
         -q /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq \
         /home/sfs/${JobName}/a5-rmrRNA/${1}.alin > \
         /home/sfs/${JobName}/a5-rmrRNA/${1}.err && \
         rm -rf /home/sfs/${JobName}/a5-rmrRNA/${1}.alin 
         
      vars_iter:
        - '${fastq_files}'

#STAR
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a6-map && \
        echo ${1} begin `date` && \
        bash /root/.bashrc && \
         STAR --runThreadN 8 --alignEndsType EndToEnd \
         --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 \
         --genomeDir /home/obs/${obs_reference_genomeFile_star} \
         --readFilesIn /home/sfs/${JobName}/a5-rmrRNA/nonrRNA/nocontam_${1} \
         --outFileNamePrefix /home/sfs/${JobName}/a6-map/${1} \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts
         
      vars_iter:
        - '${fastq_files}'

# htseq
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/ && \
        mkdir -p /home/sfs/${JobName}/a7-htcount && \
        ls /root/miniconda2/bin/ && \
        ls /home/sfs/${JobName}/a6-map/ &&\
        echo ${1} begin `date` && \
        echo /root/miniconda2/bin/python \
        /root/miniconda2/bin/RPF_count_CDS.py \
        /home/sfs/${JobName}/a6-map/${1}Aligned.sortedByCoord.out.bam \
        /home/obs/${obs_reference_gtf} > /home/sfs/${JobName}/a7-htcount/${1}.count && \
        /root/miniconda2/bin/python \
        /root/miniconda2/bin/RPF_count_CDS.py \
        /home/sfs/${JobName}/a6-map/${1}Aligned.sortedByCoord.out.bam \
        /home/obs/${obs_reference_gtf} > /home/sfs/${JobName}/a7-htcount/${1}.count
      vars_iter:
        - '${fastq_files}'

# merge
commands:
      - ' ls -alh /home/sfs/${JobName}/a7-htcount && echo bash  /root/miniconda2/bin/merge.sh -a ${fastq_files_name} -b ${fastq_files_label} -c /home/sfs/${JobName}/a7-htcount  &&  bash  /root/miniconda2/bin/merge.sh -a ${fastq_files_name} -b ${fastq_files_label} -c /home/sfs/${JobName}/a7-htcount '

# xtail
        commands:
      - ' mkdir -p /home/sfs/${JobName}/a8-xtail && Rscript /home/test/xtail.r /home/sfs/${JobName}/a7-htcount/merge.counter ${xtail_ribo_vector} ${xtail_rna_vector} ${xtail_label} /home/sfs/${JobName}/a8-xtail '

# cp
commands:
      - 'obsutil config -i=${AccessKey} -k=${SecretKey} -e=${endpoint}&& obsutil mkdir -p ${obs_location}/output/${JobName}/  && obsutil cp -r -f /home/sfs/${JobName}/ ${obs_location}/output/${JobName}/ && rm -rf /home/sfs/${JobName} && echo Check sfs && ls -al /home/sfs && ls -al /home/obs/output'
   
```

## Ribo-seq Ribocode

```sh
#cutadapter
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a2-cutadapter && \
        ls /home/obs/${obs_data_path}/ && \
        echo ${1} begin `date` /root/miniconda3/bin/cutadapt -m 18 \
        --match-read-wildcards -a ${adapter} \
        -o /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
        /home/sfs/${JobName}/a1-fastq/${1} && \
        /root/miniconda3/bin/cutadapt -m 18 \
        --match-read-wildcards -a ${adapter} \
        -o /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
         /home/obs/${obs_data_path}/${1}
    vars_iter:
        - '${fastq_files}'

# fastq_quality_filter
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a3-filter && \
        echo ${1} begin `date` && \
        /home/test/bin/fastq_quality_filter \
        -Q33 -v -q 25 -p 75 \
        -i /home/sfs/${JobName}/a2-cutadapter/${1}_trimmed.fastq \
        -o /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq 
      vars_iter:
        - '${fastq_files}'

# fastQC
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a4-qc && \
        echo ${1} begin `date` && \
        /home/test/FastQC/fastqc \
        /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq \
        -o /home/sfs/${JobName}/a4-qc
      vars_iter:
        - '${fastq_files}'

#remove rRNA
            commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a5-rmrRNA && \
        mkdir -p /home/sfs/${JobName}/a5-rmrRNA/nonrRNA && \
        echo ${1} begin `date` && \
        bash /root/.bashrc && \
        /home/test/bowtie-1.2.3-linux-x86_64/bowtie \
        -n 0 -norc --best -l 15 -p 8 \
         --un=/home/sfs/${JobName}/a5-rmrRNA/nonrRNA/nocontam_${1} /home/obs/${obs_reference_rRNA_bowtie} \
         -q /home/sfs/${JobName}/a3-filter/${1}_trimmedQfilter.fastq \
         /home/sfs/${JobName}/a5-rmrRNA/${1}.alin > \
         /home/sfs/${JobName}/a5-rmrRNA/${1}.err && \
         rm -rf /home/sfs/${JobName}/a5-rmrRNA/${1}.alin 
         
      vars_iter:
        - '${fastq_files}'

#STAR
    commands_iter:
      command: |
        mkdir -p /home/sfs/${JobName}/a6-map && \
        echo ${1} begin `date` && \
        bash /root/.bashrc && \
         STAR --runThreadN 8 --alignEndsType EndToEnd \
         --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 \
         --genomeDir /home/obs/${obs_reference_genomeFile_star} \
         --readFilesIn /home/sfs/${JobName}/a5-rmrRNA/nonrRNA/nocontam_${1} \
         --outFileNamePrefix /home/sfs/${JobName}/a6-map/${1} \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts
         
      vars_iter:
        - '${fastq_files}'


# merge
commands:
      - ' ls -alh /home/sfs/${JobName}/a7-htcount && echo bash  /root/miniconda2/bin/merge.sh -a ${fastq_files_name} -b ${fastq_files_label} -c /home/sfs/${JobName}/a7-htcount  &&  bash  /root/miniconda2/bin/merge.sh -a ${fastq_files_name} -b ${fastq_files_label} -c /home/sfs/${JobName}/a7-htcount '

# xtail
        commands:
      - ' mkdir -p /home/sfs/${JobName}/a8-xtail && Rscript /home/test/xtail.r /home/sfs/${JobName}/a7-htcount/merge.counter ${xtail_ribo_vector} ${xtail_rna_vector} ${xtail_label} /home/sfs/${JobName}/a8-xtail '

# cp
commands:
      - 'obsutil config -i=${AccessKey} -k=${SecretKey} -e=${endpoint}&& obsutil mkdir -p ${obs_location}/output/${JobName}/  && obsutil cp -r -f /home/sfs/${JobName}/ ${obs_location}/output/${JobName}/ && rm -rf /home/sfs/${JobName} && echo Check sfs && ls -al /home/sfs && ls -al /home/obs/output'
   
```

### test on hub docker.
 
```sh

```

## 文件挂载

```sh
volumes:
  volumes-4ndk:
    mount_path: '/home/sfs'
    mount_from:
      pvc: '${GCS_SFS_PVC}'
  genobs:
    mount_path: '/home/obs'
    mount_from:
      pvc: '${GCS_DATA_PVC}'
```



## 

```sh
mkdir -p /home/sfs/a5-rmrRNA && \ mkdir -p /home/sfs/a5-rmrRNA/nonrRNA && \ echo SRR3498212.fq begin `date` && \ bash /root/.bashrc && \ /home/test/bowtie-1.2.3-linux-x86_64/bowtie \ -n 0 -norc --best -l 15 -p 8 \ --un=/home/sfs/a5-rmrRNA/nonrRNA/nocontam_SRR3498212.fq /home/obs/arabidopsis/huawei_file/refrence/tair_rRNA_bowtie_index/tair.rRNA.fa \ -q /home/sfs/a3-filter/SRR3498212.fq_trimmedQfilter.fastq \ /home/sfs/a5-rmrRNA/SRR3498212.fq.alin > \ /home/sfs/a5-rmrRNA/SRR3498212.fq.err && \ rm -rf /home/sfs/a5-rmrRNA/SRR3498212.fq.alin

obsutil config -i=5ULAGR0CWKBAEDV57Y6P -k=gvroYZE9uUmp3igpEPAEQRfuQzUjcVQn9kBoHz02 -e=https://obs.cn-north-4.myhuaweicloud.com&& obsutil mkdir -p obs://hw-gcs-logo-cn-north-4-06a54be3938010610f01c00da675d700/output/arabidopsis-smallrnaseq/ && obsutil cp -r -f /home/sfs/arabidopsis-smallrnaseq/ obs://hw-gcs-logo-cn-north-4-06a54be3938010610f01c00da675d700/output/arabidopsis-smallrnaseq/ && rm -rf /home/sfs/arabidopsis-smallrnaseq && echo Check sfs && ls -al /home/sfs && ls -al /home/obs/output
```


## bash with parameters

**merge.sh**

```sh
#!/bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -a FileName -b Lables -c Directory"
   echo -e "\t-a Description of what is Filename,like file1.counter;file2.counter;file3.counter;file4.counter"
   echo -e "\t-b Description of what is Label, like dark1;dark2;dark3;dark4"
   echo -e "\t-c Description of what is Diectory path:/Users/xugang/Downloads"
   echo -e '\t ./merge.sh -a "SRR3498212.fq;SRR3498213.fq;SRR966479.fq;SRR966480.fq" -b "control1;control2;case1;case2" -c "./"'
   echo "chmod 755 ./merge.sh";
   exit 1 # Exit script after printing help
 }

 while getopts "a:b:c:" opt
 do
    case "$opt" in
      a ) parameterA="$OPTARG" ;;
        b ) parameterB="$OPTARG" ;;
      c ) parameterC="$OPTARG" ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
   done

   # Print helpFunction in case parameters are empty
   if [ -z "$parameterA" ] || [ -z "$parameterB" ] || [ -z "$parameterC" ]
   then
      echo "Some or all of the parameters are empty";
     helpFunction
 fi

 # Begin script in case all parameters are correct
 echo "$parameterA"
 echo "$parameterB"
 echo "$parameterC"

 cd $parameterC

 IFS=';' read -ra name <<< "$parameterA"
 IFS=';' read -ra label <<< "$parameterB"

 for i in "${name[@]}";
 do
 echo $i;
 done

 merge_file(){
 head='gene'
 for i in ${label[@]};
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
  mv merge.counter2 heat.counter
  }

merge_file

```

**test script:**

```sh
sh merge.sh -a "SRR3498212.fq;SRR3498213.fq;SRR966479.fq;SRR966480.fq" -b "control1;control2;case1;case2" -c "./"

sh ./myscript -a "riboseq_heat;RNAseq_heat;riboseq_control_heat2;RNAseq_control_heat2;riboseq_control_heat1;RNAseq_control_heat1" -b "riboseq_heat;RNAseq_heat;riboseq_control_heat2;RNAseq_control_heat2;riboseq_control_heat1;RNAseq_control_heat1" -c "./"
String A
String B
String C

$ ./myscript -a "SRR966479.fq;" -c "String C" -b "String B"
String A
String B
String C

$ ./myscript -a "String A" -c "String C" -f "Non-existent parameter"
./myscript: illegal option -- f

Usage: ./myscript -a parameterA -b parameterB -c parameterC
    -a Description of what is parameterA
    -b Description of what is parameterB
    -c Description of what is parameterC

$ ./myscript -a "String A" -c "String C"
Some or all of the parameters are empty

Usage: ./myscript -a parameterA -b parameterB -c parameterC
    -a Description of what is parameterA
    -b Description of what is parameterB
    -c Description of what is parameterC


```
```r
args <- commandArgs(trailingOnly = TRUE)
# inputdata
filename=args[1]
# ribovector
ribovector=as.integer(unlist(strsplit(args[2],",")))
# rnavector
rnavector=as.integer(unlist(strsplit(args[3],",")))
# label mean
label=unlist(strsplit(args[4],","))
# output
output=args[5]
print(filename)
print(ribovector)
print(rnavector)
print(label)
print(paste(output,"/df",sep=""))


library(xtail)
lbxd=read.table(filename,header=T,row.name=1)
mrna=lbxd[,ribovector]
rpf=lbxd[,rnavector]
condition=label
test.results=xtail(mrna,rpf,condition,bins=1000,threads=2)
summary(test.results)

#
test.tab=resultsTable(test.results);
head(test.tab,5)

write.table(test.tab,paste(output,"/control_case_results.txt",sep=""),quote=F,sep="\t");

# Visualization
pdf(paste(output,"/control_caseFC.pdf",sep=""),width=6,height=4,paper='special')
lbxfc=plotFCs(test.results)
dev.off()
write.table(lbxfc$resultsTable,paste(output,"/control_casefc_results.txt",sep=""),quote=F,sep="\t");

pdf(paste(output,"control_caseRs.pdf",sep=""),width=6,height=4,paper='special')
lbxrs=plotRs(test.results)
dev.off()
write.table(lbxrs$resultsTable,paste(output,"/control_casers_results.txt",sep=""),quote=F,sep="\t");

pdf(paste(output,"/control_casevolcano.pdf",sep=""),width=6,height=4,paper='special')
volcanoPlot(test.results)
dev.off()
```

```sh
# run
Rscript xtail.r merge2.counter 1,2 3,4 control,case ./output/
```


## 对字符串进行分割成数组

```sh
string='["SRR3498212.fq","SRR966479.fq"]'
string=`sed "s/\"//g" <<<"$string"`
string=`sed "s/\[//g" <<<"$string"`
string=`sed "s/\]//g" <<<"$string"`
IFS=', ' read -r -a array <<< "$string"

for element in "${array[@]}"
do
    echo "$element"
done
```
