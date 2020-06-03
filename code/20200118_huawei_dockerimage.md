# 华为Docker image

## Content
* [上传镜像](#上传镜像)


## 上传镜像

```sh
docker save gangxu/base_ubuntu:4.0 | gzip > gangxu_base_ubuntu.tar.gz

docker login -u cn-north-4@Y3NHYJC8KGOGABXQMM9H -p 70216640613c345678fdb439ce901fe4cb83546ea2b429abc13ee23a60913fa0 swr.cn-north-4.myhuaweicloud.com
docker login -u cn-north-4@IES26GBWN2NH8G4ESXEZ -p a88225c20c59a275ee2083f03ea23be53185f93a4f06760978ff03814af3dbe3 swr.cn-north-4.myhuaweicloud.com

$ sudo docker tag [{镜像名称}:{版本名称}] swr.cn-north-4.myhuaweicloud.com/{组织名称}/{镜像名称}:{版本名称}
$ sudo docker tag gangxu/base_ubuntu:4.0 swr.cn-north-4.myhuaweicloud.com/gangxu/base:1.0
$ sudo docker push swr.cn-north-4.myhuaweicloud.com/{组织名称}/{镜像名称}:{版本名称}
$ sudo docker push swr.cn-north-4.myhuaweicloud.com/gangxu/base:1.0
docker push swr.cn-north-4.myhuaweicloud.com/gangxu/ribocode_ribominer:1.0
docker push swr.cn-north-4.myhuaweicloud.com/gangxu/htseq:1.0
docker tag swr.cn-north-4.myhuaweicloud.com/gangxu/htseq:1.0 gangxu/htseq:1.0

docker push gangxu/htseq:1.0

yanglab
62783319d226yang

```

```sh
docker run -dt --name base -v /home/xugang/app:/home/app gangxu/base_ubuntu:4.0
docker exec -it base bash
exit
docker stop base
docker rm base
```

```sh
docker ps|grep base
[xugang@hub app]$ docker ps|grep base
d881f9fdee83        gangxu/base_ubuntu:4.0     "/bin/bash"              About a minute ago   Up About a minute                                                                    base

```
fastqc
```sh
docker run -dt --name base -v /home/xugang/app:/home/app gangxu/base_ubuntu:4.0
docker exec -it base bash
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
apt-get update
apt-get install libcam-pdf-perl
apt-get install default-jdk
docker commit d881f9fdee83 swr.cn-north-4.myhuaweicloud.com/gangxu/fastqc:1.0
docker push swr.cn-north-4.myhuaweicloud.com/gangxu/fastqc:1.0
docker stop base
docker rm base

docker run -dt --name fastqc -v /home/xugang/singularity_image/huawei_file:/home/sfs swr.cn-north-4.myhuaweicloud.com/gangxu/fastqc:1.0
docker exec -it fastqc bash
/home/test/FastQC/fastqc
exit
docker stop fastqc
docker rm fastqc
```
xtail
```sh
docker tag xug15/xtail:latest swr.cn-north-4.myhuaweicloud.com/gangxu/xtail:1.0
docker push swr.cn-north-4.myhuaweicloud.com/gangxu/xtail:1.0
docker run -dt --name xtail -v /home/xugang/singularity_image/huawei_file:/home/sfs swr.cn-north-4.myhuaweicloud.com/gangxu/cutadapter:1.0
docker exec -it xtail bash
docker stop xtail
docker rm xtail
```
```sh
sudo singularity build --sandbox /home/xugang/singularity_image/xtail docker://gangxu/xtail:latest

```

bowtie
```sh
docker run -dt --name base -v /home/xugang/app:/home/app gangxu/base_ubuntu:4.0
docker exec -it base bash
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.3/bowtie-1.2.3-linux-x86_64.zip
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
apt-get update
apt-get install libcam-pdf-perl
apt-get update
  153  apt-get install python
  154  apt-get --reinstall install python-minimal
  155  apt-get --reinstall install python-minimal
exit
docker commit a9a6063a86fd swr.cn-north-4.myhuaweicloud.com/gangxu/bowtie12:1.1
docker push swr.cn-north-4.myhuaweicloud.com/gangxu/bowtie12:1.1
docker run -dt --name bowtie -v /home/xugang/singularity_image/huawei_file:/home/sfs swr.cn-north-4.myhuaweicloud.com/gangxu/bowtie12:1.1
/home/test/bowtie-1.2.3-linux-x86_64
/home/test/bowtie2-2.3.5.1-linux-x86_64
docker stop bowtie
docker rm bowtie

```
## HTSeq.
```sh
docker run -dt --name bowtie -v /home/xugang/singularity_image/huawei_file:/home/sfs swr.cn-north-4.myhuaweicloud.com/gangxu/bowtie12:1.1
docker exec -it bowtie bash
/root/miniconda3/bin/htseq-count

wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
cd /root/miniconda2/bin
./pip install numpy
./pip install matplotlib
./pip install pysam
./pip install HTSeq
./pip install 
rm -rf /home/test/*
exit
docker ps
docker commit 1f6dea975073  swr.cn-north-4.myhuaweicloud.com/gangxu/htseq:1.0
docker stop bowtie
docker rm bowtie

docker run -dt --name htseq -v /home/xugang/singularity_image/huawei_file:/home/sfs swr.cn-north-4.myhuaweicloud.com/gangxu/htseq:1.0
cp merge.sh /root/miniconda2/bin
chmod 755 /root/miniconda2/bin/merge.sh
# the file merge.sh is below.
docker commit 5beb5d631ab2 swr.cn-north-4.myhuaweicloud.com/gangxu/htseq:2.0
docker push swr.cn-north-4.myhuaweicloud.com/gangxu/htseq:2.0

docker stop htseq
docker rm htseq
docker run -dt --name htseq -v /home/xugang/singularity_image/huawei_file:/home/sfs gangxu/htseq:1.1

other verstion
cd /home/test/miniconda2/bin
./pip install numpy
./pip install matplotlib
./pip install pysam
./pip install HTSeq
exit
docker ps
docker commit ce6a6eeece26  gangxu/htseq:1.1
docker push gangxu/htseq:1.1
docker stop bowtie
docker rm bowtie
docker push  gangxu/htseq:1.1
sudo singularity build --sandbox /home/xugang/singularity_image/htseq docker://gangxu/htseq:1.1

sudo singularity shell -w  /home/xugang/singularity_image/htseq 
singularity exec -B /home/xugang/singularity_image/huawei_file:/home/sfs /home/xugang/singularity_image/htseq /home/test/miniconda2/bin/htseq-count
sudo singularity exec /home/xugang/singularity_image/htseq /home/test/miniconda2/bin/htseq-count


singularity exec /WORK/teaching/project/singularity_images/htseq /home/test/miniconda2/bin/htseq-count
```

## Ribocode
```sh
docker run -dt --name bowtie -v /home/xugang/singularity_image/huawei_file:/home/sfs swr.cn-north-4.myhuaweicloud.com/gangxu/bowtie12:1.1
docker exec -it bowtie bash

pip install ribocode
pip install RiboMiner
exit
docker commit ce6a6eeece26  swr.cn-north-4.myhuaweicloud.com/gangxu/ribocode_ribominer:1.0
```


cutadapter
```sh
# use bowtie images
conda install -c bioconda cutadapt
rm -rf /home/test/bowtie*
docker commit 4790821a127b swr.cn-north-4.myhuaweicloud.com/gangxu/cutadapter:1.0
docker push swr.cn-north-4.myhuaweicloud.com/gangxu/cutadapter:1.0
docker stop base
docker rm base

docker run -dt --name cutadapter -v /home/xugang/singularity_image/huawei_file:/home/sfs swr.cn-north-4.myhuaweicloud.com/gangxu/cutadapter:1.0

docker stop cutadapter
docker rm cutadapter
```

cutadapter test:

```sh
docker run -dt --name cutadapt -v /home/xugang/singularity_image/huawei_file:/home/sfs swr.cn-north-4.myhuaweicloud.com/gangxu/cutadapter:1.0

docker exec -it cutadapt bash
/root/miniconda3/bin/cutadapt -m 18 --match-read-wildcards -a CTGTAGGCACCATCAAT -o /home/sfs/a2-cutadapter/SRR3498212.fq_trimmed.fastq /home/sfs/a1-fastq/SRR3498212.fq 

docker stop cutadapt
docker rm cutadapt

```

fastx_toolkit
```sh
docker run -dt --name base -v /home/xugang/app:/home/app gangxu/base_ubuntu:4.0
docker exec -it base bash
wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
exit
docker commit 2c9d0b2e56e6 swr.cn-north-4.myhuaweicloud.com/gangxu/fastx_toolkit:1.0
docker push swr.cn-north-4.myhuaweicloud.com/gangxu/fastx_toolkit:1.0
docker stop base
docker rm base

docker run -dt --name fastx -v /home/xugang/singularity_image:/home/sfs swr.cn-north-4.myhuaweicloud.com/gangxu/fastx_toolkit:1.0
docker exec -it fastx bash
/home/test/bin/fastq_quality_filter
exit

```
bedtools
```sh
docker run -dt --name base  /home/xugang/app:/home/app gangxu/base_ubuntu:4.0
docker exec -it base bash
apt-get update
sudo apt-get install libboost-all-dev
sudo apt-get install libbz2-dev
apt-get install liblzma-dev
$ wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
$ tar -zxvf bedtools-2.29.1.tar.gz
$ cd bedtools2
$ make
exit
docker ps|grep base
docker commit ebe6a9bfece9 swr.cn-north-4.myhuaweicloud.com/gangxu/bedtools:1.0
docker push swr.cn-north-4.myhuaweicloud.com/gangxu/bedtools:1.0
docker stop base
docker rm base
docker run -dt --name bedtools -v /home/xugang/singularity_image:/home/sfs swr.cn-north-4.myhuaweicloud.com/gangxu/bedtools:1.0
docker exec -it bedtools bash



```
samtools
```sh
docker run -dt --name base -v /home/xugang/app:/home/app gangxu/base_ubuntu:4.0
docker exec -it base bash
wget https://sourceforge.net/projects/samtools/files/samtools/1.10/samtools-1.10.tar.bz2
wget https://sourceforge.net/projects/samtools/files/samtools/1.10/bcftools-1.10.tar.bz2
apt-get update
apt-get install libncurses5-dev
apt-get install libbz2-dev
apt-get install liblzma-dev
./configure
make
make install
exit
docker ps|grep base
docker commit fc5552888cdc swr.cn-north-4.myhuaweicloud.com/gangxu/samtools2bcftools:1.0
docker push swr.cn-north-4.myhuaweicloud.com/gangxu/samtools2bcftools:1.0
docker stop base
docker rm base
docker run -dt --name samtools -v /home/xugang/singularity_image:/home/sfs swr.cn-north-4.myhuaweicloud.com/gangxu/samtools2bcftools:1.0
samtools exec -it samtools bash
```




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

obsutil config -i=5ULAGR0CWKBAEDV57Y6P -k=gvroYZE9uUmp3igpEPAEQRfuQzUjcVQn9kBoHz02 -e=https://obs.cn-north-4.myhuaweicloud.com && obsutil cp -r -f -u obs://gene-container-xugang/gcs/huawei_file/a1-fastq/ /home/sfs && obsutil cp -r -f -u obs://gene-container-xugang/gcs/huawei_file/refrence/ /home/sfs && ls /home/sfs

https://support.huaweicloud.com/tr-gcs/gcs_tr_04_0002.html

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
    description: '将数据拷贝到制定目录下/home/sfs'
    description: 将原始的下机数据去除接头
    description: 去除低质量的数据    
    description: 对数据进行质量评估
    description: 去除核糖体的reads
    description: 比对到基因组上
    description: 将结果拷贝回obsvolumn中

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

```sh
mkdir -p /home/sfs/a5-rmrRNA && \ mkdir -p /home/sfs/a5-rmrRNA/nonrRNA && \ echo SRR3498212.fq begin `date` && \ bash /root/.bashrc && \ /home/test/bowtie-1.2.3-linux-x86_64/bowtie \ -n 0 -norc --best -l 15 -p 8 \ --un=/home/sfs/a5-rmrRNA/nonrRNA/nocontam_SRR3498212.fq /home/obs/arabidopsis/huawei_file/refrence/tair_rRNA_bowtie_index/tair.rRNA.fa \ -q /home/sfs/a3-filter/SRR3498212.fq_trimmedQfilter.fastq \ /home/sfs/a5-rmrRNA/SRR3498212.fq.alin > \ /home/sfs/a5-rmrRNA/SRR3498212.fq.err && \ rm -rf /home/sfs/a5-rmrRNA/SRR3498212.fq.alin

obsutil config -i=5ULAGR0CWKBAEDV57Y6P -k=gvroYZE9uUmp3igpEPAEQRfuQzUjcVQn9kBoHz02 -e=https://obs.cn-north-4.myhuaweicloud.com&& obsutil mkdir -p obs://hw-gcs-logo-cn-north-4-06a54be3938010610f01c00da675d700/output/arabidopsis-smallrnaseq/ && obsutil cp -r -f /home/sfs/arabidopsis-smallrnaseq/ obs://hw-gcs-logo-cn-north-4-06a54be3938010610f01c00da675d700/output/arabidopsis-smallrnaseq/ && rm -rf /home/sfs/arabidopsis-smallrnaseq && echo Check sfs && ls -al /home/sfs && ls -al /home/obs/output
```


```sh
docker save swr.cn-north-4.myhuaweicloud.com/gangxu/bowtie12:1.1 | gzip > /home/xugang/singularity_image/huawei_file/images/bowtie12.tar.gz
docker save swr.cn-north-4.myhuaweicloud.com/gangxu/samtools2bcftools:1.0 | gzip > /home/xugang/singularity_image/huawei_file/images/samtools2bcftools.tar.gz
docker save swr.cn-north-4.myhuaweicloud.com/gangxu/bedtools:1.0 | gzip > /home/xugang/singularity_image/huawei_file/images/bedtools.tar.gz
docker save swr.cn-north-4.myhuaweicloud.com/gangxu/fastx_toolkit:1.0 | gzip > /home/xugang/singularity_image/huawei_file/images/fastx_toolkit.tar.gz
docker save swr.cn-north-4.myhuaweicloud.com/gangxu/cutadapter:1.0 | gzip > /home/xugang/singularity_image/huawei_file/images/cutadapter.tar.gz
docker save swr.cn-north-4.myhuaweicloud.com/gangxu/fastqc:1.0 | gzip > /home/xugang/singularity_image/huawei_file/images/fastqc.tar.gz
docker save swr.cn-north-4.myhuaweicloud.com/gangxu/base:1.0 | gzip > /home/xugang/singularity_image/huawei_file/images/base.tar.gz
```
## bash with parameters.

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

