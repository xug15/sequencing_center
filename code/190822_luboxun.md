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

##

