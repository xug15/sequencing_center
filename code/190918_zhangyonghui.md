# zhangyonghui
/Share/home/tiangeng/project_result/sgRNA_count/project111_20190815_zhangyonghui_ylj/YLJ-6

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

cut adapter
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

id, index, name

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

Build bowtie index

```sh
bowtie-build geckov2.fa geckov2
```


