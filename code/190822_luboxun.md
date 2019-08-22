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
## 

##

