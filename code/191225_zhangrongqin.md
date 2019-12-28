# zhangqingrong
## 1.remove adapter

**a1.remove_adapter.sh**

```sh
#!/bin/bash
#SBATCH -J adapter
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --output=%j.out
#SBATCH --error=%j.err

source /WORK/Samples/singularity.sh

adapt=AGATCGGAAGAGCACACGTCTGAACTC
dir=/WORK/teaching/project/20191225zqr/miRNA/
i=clean.fq

echo "cutadapt -m 17 --match-read-wildcards -a $adapt -o ${dir}${i}_trimmed.fastq ${dir}${i}  > ${dir}${i}_trimmed.log"
cutadapt -m 17 --match-read-wildcards -a $adapt -o ${dir}${i}_trimmed.fastq ${dir}${i}  > ${dir}${i}_trimmed.log

```

## 2. Build the STAR index in genome 2017.

**a2.build_star_index.sh**
```sh

grep ^PIN Pinctada.gene.final.chr.gff| sort -k 1,1 -k 4,4n >Pinctada.sort.gff
sd
```

**a3.build_star_index.sh**

```sh
#!/bin/bash
#SBATCH -J adapter
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err
source /WORK/Samples/singularity.sh


fasta=/WORK/teaching/project/20191225zqr/genome2017/Pinctada.chromosome.fasta
gtf=/WORK/teaching/project/20191225zqr/genome2017/Pinctada.sort.gff
genomeout=/WORK/teaching/project/20191225zqr/genome2017/Pinctada_star

mkdir $genomeout
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir ${genomeout} --genomeFastaFiles ${fasta} --sjdbGTFfile ${gtf} >log.txt 2>&1
```

**submit**

```sh
sbatch a3.build_star_index.sh

```

## 2. Mapping miRNA seq into 2017 genome.

```sh

```

## 3.miRNA

Download from http://www.mirbase.org/ftp.shtml

**a3.u2t.pl**

```pl
open DATA, "<$ARGV[0]";
open OUT,  ">$ARGV[0].t.fa";
while(<DATA>)
{
if($_=~/^>/)
    {
        print OUT $_;
    }else{
        $_=~s/U/T/g;
        print OUT $_;

    }
}
close DATA;
close OUT;

```

```sh
mkdir hairpin_bowtie_index
mature_bowtie_index
```

**a4.bowtie.sh**

```sh
#!/bin/bash
#SBATCH -J bowtie
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err
source /WORK/Samples/singularity.sh


out=/WORK/teaching/project/20191225zqr/miRNA_mapping
fasta=/WORK/teaching/project/20191225zqr/miRNA/clean.fq_trimmed.fastq
bowtie_index=/WORK/teaching/project/20191225zqr/miRNA_reference/mature_bowtie_index/mature
name=miRNA

mkdir -p $out

echo -e "bowtie -n 0 -norc --best -l 15 -p 8 --un=$out/unmap_${name}.fastq $bowtie_index -q $fasta $out/${name}.sam > $out/miRNA_mapping_log.err 2>&1"
bowtie -n 0 -norc --best -l 15 -p 8 --un=$out/unmap_${name}.fastq $bowtie_index -q $fasta $out/${name}.sam > $out/miRNA_mapping_log.err 2>&1

```

**a6.read2fa.pl**

```pl
open DATA, "<$ARGV[0]";
open OUT,  ">$ARGV[0].t.fa";
while(<DATA>)
{
$seq=<DATA>;
<DATA>;
<DATA>;
print OUT $seq;
}
close DATA;
close OUT;

```

**a7.new_miRNA.sh**

```sh
#!/bin/bash
#SBATCH -J bowtie
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err
source /WORK/Samples/singularity.sh


sort /WORK/teaching/project/20191225zqr/miRNA_mapping/unmap_miRNA.fastq.t.fa | uniq -c > /WORK/teaching/project/20191225zqr/miRNA_mapping/new.txt
```

```sh
perl -pi -e 's/^\s+//g' /WORK/teaching/project/20191225zqr/miRNA_mapping/new.txt;

```

**a8.blank2tab.pl**

```pl

open DATA, "<new.txt";
open OUT,  ">new.miRNA.sum.txt";
while(<DATA>)
{
chomp;
@data=split / /,$_;
print OUT "$data[1]\t$data[0]\n";
}
close DATA;
close OUT;
```

```sh
sort -k 2,2nr new.miRNA.sum.txt > new.miRNA.sum.sort.txt

```


## Tophat map to genome

**a5.tophat2017.sh**

```sh
#!/bin/bash
#SBATCH -J tophat
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err
export PATH=/WORK/teaching/bin:$PATH

out=/WORK/teaching/project/20191225zqr/genome_2017_mapping
fasta=/WORK/teaching/project/20191225zqr/miRNA/clean.fq_trimmed.fastq
bowtie_index=/WORK/teaching/project/20191225zqr/genome2017/Pictada_bowtie/Pinctada.chromosome.fasta
name=miRNA
gtf=/WORK/teaching/project/20191225zqr/genome2017/Pinctada.sort.gff

mkdir $out
mkdir $out/tophat
 
echo -e "tophat -p 8 -o ${out}/tophat -g 10 --bowtie1 --read-realign-edit-dist 0  -G $gtf --no-novel-juncs --segment-length=15 $bowtie_index $fasta > $out/${name}.log"
tophat -p 8 -o ${out}/tophat -g 10 --bowtie1 --read-realign-edit-dist 0  -G $gtf --no-novel-juncs --segment-length=15 $bowtie_index $fasta > $out/${name}.log

```


## . Build the STAR index in genome 2016.

```sh

```

## . Mapping miRNA seq into 2016 genome.

```sh

```


