usage = \
"""
statistic the distribution of reads' length from .fastq file
Usage: python fq_len_stat.py filename.fastq output.png
"""
import matplotlib
matplotlib.use('Agg')
from sys import argv, exit
import numpy as np
import pylab as plb

if len(argv) != 3:
        print usage
        exit()

length_stat = {}
length_list = []
noreads = 0
length_stat[0] = 0
with open(argv[1]) as fqFile:
        for no,line in enumerate(fqFile):
                line = line.strip()
                no += 1
                if no % 4 == 2:
                        length = len(line)
                        if length == 0:
                                noreads += 1
                                length_stat[0] += 1
                                continue
                        else:
                                length_list.append(length)
                                if length in length_stat:
                                        length_stat[length] += 1
                                else:
                                        length_stat[length] = 1

length_array = np.array(length_list)
#print the information
print "total reads number:\t%i" % len(length_list)
print "Zero length number:\t%i" % noreads
print "mean of the length:\t%f" % length_array.mean()
print "std of the length:\t%f" % length_array.std()

#plot
labels = sorted(length_stat.keys())
value_list = [length_stat[i] for i in labels]
plb.bar(labels,value_list,width=0.5)
plb.xlabel("length of the reads")
#plb.xlim([0, 50])
#plb.ylim([0, 100000])
plb.ylabel("NO. of the reads")
plb.title(argv[2])
plb.savefig(argv[2])
plb.close()