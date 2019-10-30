# zhangyonghui
/Share/home/tiangeng/project_result/sgRNA_count/project111_20190815_zhangyonghui_ylj/YLJ-6


## id, index, name

```sh
HGLibA_00001,GTCGCTGAGCTCCGATTCGA,A1BG
HGLibA_00002,ACCTGTAGTTGCCGGCGTGC,A1BG
HGLibA_00003,CGTCAGCGTCACATTGGCCA,A1BG
HGLibA_00004,CGCGCACTGGTCCAGCGCAC,A1CF
HGLibA_00005,CCAAGCTATATCCTGTGCGC,A1CF
HGLibA_00006,AAGTTGCTTGATTGCATTCT,A1CF
HGLibA_00007,CGCTTCTTAAATTCTTGGGT,A2M
HGLibA_00008,TCACAGCGAAGGCGACACAG,A2M
HGLibA_00009,CAAACTCCTTCATCCAAGTC,A2M
HGLibA_00010,AAATTTCCCCTCCGTTCAGA,A2ML1
```

## a1.csv2fa.pl

```pl
open D, "<geckov2.csv";
open O, ">geckov2.fa";

while(<D>)
{
    chomp;
    @data=split/,/,$_;
    print O ">$data[0]\n$data[1]\n";
}
```

## Build bowtie index

```sh
cd /Share/home/tiangeng/reference_genome/zhangyonghui

bowtie-build geckov2.fa geckov2
```

```py

import sys
import re
f=open(sys.argv[1],'r')
line=f.readline()
barcode=[]
while line:
    #remove the \n
    line=line.strip("\n")
    print(line)
    #split barcode information into array.
    data=re.split(r'\s+',line)
    if len(data)!=3:
        data=re.split(r'\t+',line)
    # add the data into the barcode array
    barcode.append(data)
    # continue the line and give it new value.
    line=f.readline()

f.close()

print (barcode)

f=open('YLJ-6_FKDL190751693-1a_2.clean.fq','r')
line=f.readline()
reverse_seq=[]
n=0
while line:
    name=line
    seq=f.readline()
    seq_index=seq[0:20]
    print(seq_index)
    reverse_seq.append(seq_index)
    name2=f.readline()
    quali=f.readline()
    line=f.readline()
    n+=1
    if n>100:
        break;

f.close()
print (reverse_seq)
f=open('YLJ-6_FKDL190751693-1a_1.clean.fq','r')
line=f.readline()
n=0
while line:

   name=line
   seq=f.readline()
   seq_index=seq[0:20]
   name2=f.readline()
   quali=f.readline()
   print(seq_index)
   for i in barcode:
       #print(i[0])
       test=re.search(i[0],seq)
       if test :
           print(i)
           print(seq_index)
   line=f.readline()
   n+=1
   if n > 100:
       break
f.close()
```

```sh
GAGCGCTA    GATAGACA    H0_4.7
TAGCGAGT    ATGCCTAA    H8_6.22
TAGCGAGT    TCCGTCTA    H8_7.1
TAGCGAGT    ATCCTGTA    H9_7.1
TAGCGAGT    GATAGACA    H9_7.15
TAGCGAGT    GCTAACGA    H9-HMBPP_7.15
GTAGCTCC    ATGCCTAA    H9_7.19
GTAGCTCC    TCCGTCTA    H10_7.5
GTAGCTCC    ATCCTGTA    H10_7.19
GTAGCTCC    GATAGACA    H11_7.19
GTAGCTCC    GCTAACGA    M0_4.5
AGGCTCCG    ATGCCTAA    M11_6.22
AGGCTCCG    TCCGTCTA    M11_7.1
AGGCTCCG    ATCCTGTA    M12_7.1
AGGCTCCG    GATAGACA    M12_7.15
AGGCTCCG    GCTAACGA    M12-HMBPP_7.15
GAGCGCTA    ATGCCTAA    M13_7.5
GAGCGCTA    TCCGTCTA    M13_7.15
GAGCGCTA    ATCCTGTA    M14_7.15
```

## cut adapter

**reverse and compliment**
```pl
my $origin_seq=<STDIN>;
my $revcomp = reverse $origin_seq;
$revcomp =~ tr/ATGCatgc/TACGtacg/;

print "$revcomp\n";
```
a1.cu.sh
```sh

cutadapt -a ADAPT1 -g ADAPT2 -o out1.fastq in1.fastq

                     
```
|parameter|mean|
|:-|:-|
|-a ADAPTER, --adapter=ADAPTER |equence of an adapter that was ligated to the 3\' end.The adapter itself and anything that follows is trimmed. If the adapter sequence ends with the '$' character, the adapter is anchored to the end of the read and only found if it is a suffix of the read. |
|-g ADAPTER, --front=ADAPTER |Sequence of an adapter that was ligated to the 5\' end. If the adapter sequence starts with the character \'^\', the adapter is \'anchored\'. An anchored adapter must appear in its entirety at the 5\' end of the read (it is a prefix of the read). A non-anchored adapter may appear partially at the 5\' end, or it may occur within the read. If it is found within a read, the sequence preceding the adapter is also trimmed. In all cases, the adapter itself is trimmed. |

## remove adapter

```sh
mkdir b2-removeadapter
```

```sh
f1=GCGAGTTCTTGTGGAAAGGACGAAACACCG
f2=CTCCTTCTTGTGGAAAGGACGAAACACCG
f4=AGGCTCCGTCTTGTGGAAAGGACGAAACACCG
f7=GAGCGCTATCTTGTGGAAAGGACGAAACACCG
f8=GCTCAGTCTCTTGTGGAAAGGACGAAACACCG

r1=ACAGTGCAGGGGAAAAGAATAGTAGATTA
r4=ACAGTGCAGGGGAAAAGAATAGTAGATAC
r7=ACAGTGCAGGGGAAAAGAATAGT
r8=ACAGTGCAGGGGAAAAGAATAGTAGAGTCGTTAGC
r10=ACAGTGCAGGGGAAAAGAATAGTAGATAGACGGA
#-a 3' -g 5'
echo "cutadapt -g $f7 -a $r7 -o H0_4.7_out.fastq ../split_file/H0_4.7.fq"
cutadapt -g $f7 -a $r7 -o H0_4.7_out.fastq ../split_file/H0_4.7.fq  >H0_4.7.log 2>&1
echo "cutadapt -g $f1 -a $r1 -o H8_6.22_out.fastq ../split_file/H8_6.22.fq"
cutadapt -g $f1 -a $r1 -o H8_6.22_out.fastq ../split_file/H8_6.22.fq> H8_6.22.log 2>&1
echo "cutadapt -g $f1 -a $r10 -o H8_7.1_out.fastq ../split_file/H8_7.1.fq"
cutadapt -g $f1 -a $r10 -o H8_7.1_out.fastq ../split_file/H8_7.1.fq >H8_7.1.log 2>&1
echo "cutadapt -g $f1 -a $r4 -o H9_7.1_out.fastq ../split_file/H9_7.1.fq"
cutadapt -g $f1 -a $r4 -o H9_7.1_out.fastq ../split_file/H9_7.1.fq>H9_7.1.log 2>&1
echo "cutadapt -g $f1 -a $r7 -o H9_7.15_out.fastq ../split_file/H9_7.15.fq"
cutadapt -g $f1 -a $r7 -o H9_7.15_out.fastq ../split_file/H9_7.15.fq >H9_7.15.log 2>&1
echo "cutadapt -g $f1 -a $r8 -o H9_HMBPP_7.15H9_out.fastq ../split_file/H9-HMBPP_7.15.fq"
cutadapt -g $f1 -a $r8 -o H9_HMBPP_7.15H9_out.fastq ../split_file/H9-HMBPP_7.15.fq>H9-HMBPP_7.15.log 2>&1
echo "cutadapt -g $f2 -a $r1 -o H9_7.19_out.fastq ../split_file/H9_7.19.fq"
cutadapt -g $f2 -a $r1 -o H9_7.19_out.fastq ../split_file/H9_7.19.fq>H9_7.19.log 2>&1
echo "cutadapt -g $f2 -a $r10 -o H10_7.5_out.fastq ../split_file/H10_7.5.fq"
cutadapt -g $f2 -a $r10 -o H10_7.5_out.fastq ../split_file/H10_7.5.fq>H10_7.5.log 2>&1
echo "cutadapt -g $f2 -a $r4 -o H10_7.19_out.fastq ../split_file/H10_7.19.fq"
cutadapt -g $f2 -a $r4 -o H10_7.19_out.fastq ../split_file/H10_7.19.fq >H10_7.19.log 2>&1
echo "cutadapt -g $f2 -a $r7 -o H11_7.19_out.fastq ../split_file/H11_7.19.fq"
cutadapt -g $f2 -a $r7 -o H11_7.19_out.fastq ../split_file/H11_7.19.fq>H11_7.19.log 2>&1
echo "cutadapt -g $f2 -a $r8 -o M0_4.5_out.fastq ../split_file/M0_4.5.fq"
cutadapt -g $f2 -a $r8 -o M0_4.5_out.fastq ../split_file/M0_4.5.fq>M0_4.5.log 2>&1
echo "cutadapt -g $f4 -a $r1 -o M11_6.22_out.fastq ../split_file/M11_6.22.fq"
cutadapt -g $f4 -a $r1 -o M11_6.22_out.fastq ../split_file/M11_6.22.fq>M11_6.22.log 2>&1
```

## cut the read to 20bp

a1.run.sh

```sh

mkdir b3-cutlength

```


```sh
for i in `ls ../b2-removeadapter|grep out.fastq`;do
echo $i;
echo "perl a2.20.pl ../b2-removeadapter/$i $i";
perl /Share/home/tiangeng/apps/a2.20.pl ../b2-removeadapter/$i $i ;
done
```

a2.20.pl
```pl
open D, "<$ARGV[0]";
open O, ">$ARGV[1].20.fa";

while(<D>)
{
    $name=$_;
    $seq=<D>;
    <D>;
    $qua=<D>;
    $seq_s = substr($seq, 0, 20); 
    $qua_s = substr($qua, 0, 20); 
    print O "$name";
    print O "$seq_s\n";
    print O "+\n";
    print O "$qua_s\n";
}
close D;
close O;
```


## bowtie map

```sh
mkdir b4-bowtie
```

```sh
bowtieindex=/Share/home/tiangeng/reference_genome/zhangyonghui/geckov2

for i in `ls ../b3-cutlength|grep .20.fa`;do
echo $i;

echo "bowtie -n 0 -norc --best -l 15 -p 8 --un=nocontam_${i}.fastq $bowtieindex -q ../b3-cutlength/${i} ${i}.alin > ${i}.err 2>&1";
bowtie -n 0 -norc --best -l 15 -p 8 --un=nocontam_${i}.fastq $bowtieindex -q ../b3-cutlength/${i} ${i}.alin > ${i}.err 2>&1

done
```

a2.merge.sh

```sh 

name=(H0_4.7_out.fastq.20.fa H10_7.19_out.fastq.20.fa H10_7.5_out.fastq.20.fa H11_7.19_out.fastq.20.fa H8_6.22_out.fastq.20.fa H8_7.1_out.fastq.20.fa H9_7.15_out.fastq.20.fa H9_7.19_out.fastq.20.fa H9_7.1_out.fastq.20.fa H9_HMBPP_7.15H9_out.fastq.20.fa M0_4.5_out.fastq.20.fa M11_6.22_out.fastq.20.fa M11_7.1_out.fastq.20.fa M12_7.15_out.fastq.20.fa M12_7.1_out.fastq.20.fa M12_HMBPP_7.15_out.fastq.20.fa M13_7.15_out.fastq.20.fa M13_7.5_out.fastq.20.fa M14_7.15_out.fastq.20.fa)
head='Iterm'
for i in ${name[@]};
do
 head+="\t${i}";
done
echo -e $head >merge.counter;

for i in ${name[@]};
do
echo ${i}.err;
grep -A 2 'reads processed' ${i}.err > ${i}.err.tmp;
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
echo -e "Iterm\nTotal\nalin\nnon">name.txt;
paste -d $'\t' name.txt summary.txt > summary2.txt
mv summary2.txt summary.txt
rm name.txt merge2.tmp
```


