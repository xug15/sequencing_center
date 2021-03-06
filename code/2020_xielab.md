# Data:

# Exprement design
|Sample ID | Descript|seeding|Data type|
|-|-|-|-|
|x-ten143-COL-O-DMSO-R_FKDL202558030-1a-33_H7MMGCCX2_L1|DMSO|COL|RPF|
|x-ten143-COL-O-DMSO-T_FKDL202558030-1a-25_H7MMGCCX2_L1|DMSO|COL|RNA|
|x-ten143-COL-O-SL-R_FKDL202558030-1a-34_H7MMGCCX2_L1|SL|COL|RPF|
|x-ten143-COL-O-SL-T_FKDL202558030-1a-26_H7MMGCCX2_L1|SL|COL|RNA|
|x-ten143-D14-DMSO-R_FKDL202558030-1a-35_H7MMGCCX2_L1|DMSO|D14|RPF|
|x-ten143-D14-DMSO-T_FKDL202558030-1a-27_H7MMGCCX2_L1|DMSO|D14|RNA|
|x-ten143-D14-SL-R_FKDL202558030-1a-36_H7MMGCCX2_L1|SL|D14|RPF|
|x-ten143-D14-SL-T_FKDL202558030-1a-28_H7MMGCCX2_L1|SL|D14|RNA|

# Fastqc
```sh
#!/bin/bash
#SBATCH -J o1.qc
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%j.out
#SBATH --error=%j.err
export SINGULARITY_BINDPATH='/WORK,/Share'
# Get the software
export PATH=/WORK/teaching/bin:$PATH

data_dir=/WORK/teaching/project/20200205_xielab/x-ten143-roboseq
name=(x-ten143-CHO-1-R_FKDL202558030-1a-29_H7MMGCCX2_L1 x-ten143-CHO-1-T_FKDL202558030-1a-21_H7MMGCCX2_L1 x-ten143-CHO-2-R_FKDL202558030-1a-30_H7MMGCCX2_L1 x-ten143-CHO-2-T_FKDL202558030-1a-22_H7MMGCCX2_L1 x-ten143-COL-O-DMSO-R_FKDL202558030-1a-33_H7MMGCCX2_L1 x-ten143-COL-O-DMSO-T_FKDL202558030-1a-25_H7MMGCCX2_L1 x-ten143-COL-O-SL-R_FKDL202558030-1a-34_H7MMGCCX2_L1 x-ten143-COL-O-SL-T_FKDL202558030-1a-26_H7MMGCCX2_L1 x-ten143-D14-DMSO-R_FKDL202558030-1a-35_H7MMGCCX2_L1 x-ten143-D14-DMSO-T_FKDL202558030-1a-27_H7MMGCCX2_L1 x-ten143-D14-SL-R_FKDL202558030-1a-36_H7MMGCCX2_L1 x-ten143-D14-SL-T_FKDL202558030-1a-28_H7MMGCCX2_L1 x-ten143-FLY-B-R_FKDL202558030-1a-31_H7MMGCCX2_L1 x-ten143-FLY-B-T_FKDL202558030-1a-23_H7MMGCCX2_L1 x-ten143-FLY-WT-R_FKDL202558030-1a-32_H7MMGCCX2_L1 x-ten143-FLY-WT-T_FKDL202558030-1a-24_H7MMGCCX2_L1)

out_dir=/WORK/teaching/project/20200205_xielab/x-ten143-roboseq/a1-fq

[ -d ${out_dir} ] || mkdir -p ${out_dir}

for i in ${name[@]};do 
echo ${i};
echo "fastqc -t 8 ${data_dir}/${i}/${i}_1.fq.gz ${data_dir}/${i}/${i}_2.fq.gz -o ${out_dir}";
fastqc -t 8 ${data_dir}/${i}/${i}_1.fq.gz ${data_dir}/${i}/${i}_2.fq.gz -o ${out_dir};
done

```

# cutadapter

## a2.cutadaptRPF.sh

```sh
#!/bin/bash
#SBATCH -J o1.qc
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%j.out
#SBATH --error=%j.err
export SINGULARITY_BINDPATH='/WORK,/Share'
# Get the software
export PATH=/WORK/teaching/bin:$PATH

data_dir=/WORK/teaching/project/20200205_xielab/x-ten143-roboseq
out_dir=/WORK/teaching/project/20200205_xielab/x-ten143-roboseq//a2-cutadapter

rpf=(x-ten143-CHO-1-R_FKDL202558030-1a-29_H7MMGCCX2_L1 x-ten143-CHO-2-R_FKDL202558030-1a-30_H7MMGCCX2_L1 x-ten143-COL-O-DMSO-R_FKDL202558030-1a-33_H7MMGCCX2_L1 x-ten143-COL-O-SL-R_FKDL202558030-1a-34_H7MMGCCX2_L1 x-ten143-D14-DMSO-R_FKDL202558030-1a-35_H7MMGCCX2_L1 x-ten143-D14-SL-R_FKDL202558030-1a-36_H7MMGCCX2_L1 x-ten143-FLY-B-R_FKDL202558030-1a-31_H7MMGCCX2_L1 x-ten143-FLY-WT-R_FKDL202558030-1a-32_H7MMGCCX2_L1)
adapt=CTGTAGGCACCATCAAT

[ -d ${out_dir} ] || mkdir -p ${out_dir}

for i in ${rpf[@]}
do
        echo " cutadapt -m 25 -M 35 --match-read-wildcards -a $adapt -o ${out_dir}/${i}_trimmed.fastq ${data_dir}/${i}/${i}_1.fq.gz > ${out_dir}/${i}_trimmed.log ";
        cutadapt -m 25 -M 35 --match-read-wildcards -a $adapt -o ${out_dir}/${i}_trimmed.fastq ${data_dir}/${i}/${i}_1.fq.gz > ${out_dir}/${i}_trimmed.log;
        done

```

## a3.cutadaptTotal.sh

```sh
#!/bin/bash
#SBATCH -J o1.qc
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%j.out
#SBATH --error=%j.err
export SINGULARITY_BINDPATH='/WORK,/Share'
# Get the software
export PATH=/WORK/teaching/bin:$PATH

data_dir=/WORK/teaching/project/20200205_xielab/x-ten143-roboseq
out_dir=/WORK/teaching/project/20200205_xielab/x-ten143-roboseq/a2-cutadapter

rpf=(x-ten143-CHO-1-T_FKDL202558030-1a-21_H7MMGCCX2_L1 x-ten143-CHO-2-T_FKDL202558030-1a-22_H7MMGCCX2_L1 x-ten143-COL-O-DMSO-T_FKDL202558030-1a-25_H7MMGCCX2_L1 x-ten143-COL-O-SL-T_FKDL202558030-1a-26_H7MMGCCX2_L1 x-ten143-D14-DMSO-T_FKDL202558030-1a-27_H7MMGCCX2_L1 x-ten143-D14-SL-T_FKDL202558030-1a-28_H7MMGCCX2_L1 x-ten143-FLY-B-T_FKDL202558030-1a-23_H7MMGCCX2_L1 x-ten143-FLY-WT-T_FKDL202558030-1a-24_H7MMGCCX2_L1)
adapt=CTGTAGGCACCATCAAT

[ -d ${out_dir} ] || mkdir -p ${out_dir}

for i in ${rpf[@]}
do
        echo " cutadapt -m 17 --match-read-wildcards -a $adapt -o ${out_dir}/${i}_trimmed.fastq ${data_dir}/${i}/${i}_1.fq.gz > ${out_dir}/${i}_trimmed.log ";
        cutadapt -m 17 --match-read-wildcards -a $adapt -o ${out_dir}/${i}_trimmed.fastq ${data_dir}/${i}/${i}_1.fq.gz > ${out_dir}/${i}_trimmed.log 
done
```

```sh

for i in `ls /WORK/teaching/project/20200205_xielab/x-ten143-roboseq/a2-cutadapter|grep log$|grep R`;
do echo $i;
grep -A 6 Summary $i;
done;

for i in `ls /WORK/teaching/project/20200205_xielab/x-ten143-roboseq/a2-cutadapter|grep log$|grep T`;
do echo $i;
grep -A 6 Summary $i;
done;

```

x-ten143-FLY-WT-R_FKDL202558030-1a-32_H7MMGCCX2_L1

# Fastqc after qc
a4.fc.after.adapter.sh
```sh
#!/bin/bash
#SBATCH -J o1.qc
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%j.out
#SBATH --error=%j.err
export SINGULARITY_BINDPATH='/WORK,/Share'
# Get the software
export PATH=/WORK/teaching/bin:$PATH

data_dir=/WORK/teaching/project/20200205_xielab/x-ten143-roboseq/a2-cutadapter
name=(x-ten143-CHO-1-R_FKDL202558030-1a-29_H7MMGCCX2_L1 x-ten143-CHO-1-T_FKDL202558030-1a-21_H7MMGCCX2_L1 x-ten143-CHO-2-R_FKDL202558030-1a-30_H7MMGCCX2_L1 x-ten143-CHO-2-T_FKDL202558030-1a-22_H7MMGCCX2_L1 x-ten143-COL-O-DMSO-R_FKDL202558030-1a-33_H7MMGCCX2_L1 x-ten143-COL-O-DMSO-T_FKDL202558030-1a-25_H7MMGCCX2_L1 x-ten143-COL-O-SL-R_FKDL202558030-1a-34_H7MMGCCX2_L1 x-ten143-COL-O-SL-T_FKDL202558030-1a-26_H7MMGCCX2_L1 x-ten143-D14-DMSO-R_FKDL202558030-1a-35_H7MMGCCX2_L1 x-ten143-D14-DMSO-T_FKDL202558030-1a-27_H7MMGCCX2_L1 x-ten143-D14-SL-R_FKDL202558030-1a-36_H7MMGCCX2_L1 x-ten143-D14-SL-T_FKDL202558030-1a-28_H7MMGCCX2_L1 x-ten143-FLY-B-R_FKDL202558030-1a-31_H7MMGCCX2_L1 x-ten143-FLY-B-T_FKDL202558030-1a-23_H7MMGCCX2_L1 x-ten143-FLY-WT-R_FKDL202558030-1a-32_H7MMGCCX2_L1 x-ten143-FLY-WT-T_FKDL202558030-1a-24_H7MMGCCX2_L1)

out_dir=/WORK/teaching/project/20200205_xielab/x-ten143-roboseq/a3-fq

[ -d ${out_dir} ] || mkdir -p ${out_dir}

for i in ${name[@]};do 
echo ${i};
echo "fastqc -t 8 ${data_dir}/${i}_trimmed.fastq -o ${out_dir}";
fastqc -t 8 ${data_dir}/${i}_trimmed.fastq -o ${out_dir};
done


```


