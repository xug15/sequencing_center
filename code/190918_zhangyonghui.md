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