import HTSeq
import os
import sys
import argparse
import random
import string
import glob
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import pylab as plb

###used for statistic the unique mapped reads
##unique mapped reads number of gene's exon and not

#change: bug while read the gtf,


#usage : python mappedNumber.py bamFile gffFile outfile.txt
def invert_strand( iv ):
        iv2 = iv.copy()
        if iv2.strand == "+":
                iv2.strand = "-"
        elif iv2.strand == "-":
                iv2.strand = "+"
        else:
                raise ValueError, "Illegal strand"
        return iv2

#bamFile = "PC003NR_rescue.sorted.bam"
bamFile = sys.argv[1]
gffFile = sys.argv[2]
#gffFile= "/data1/DataBase/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/xiao_rebuild/genes_filter.gtf"
feature_type = "exon"
id_attribute = "gene_id"
features = HTSeq.GenomicArrayOfSets( "auto", stranded = "yes" )
geneFeatures = HTSeq.GenomicArrayOfSets( "auto", stranded = "yes" )
geneRange = {}
#read gtffile

gff = HTSeq.GFF_Reader(gffFile)
i = 0
for f in gff:
        #if f.attr["gene_biotype"] ==  'rRNA' or f.attr["gene_biotype"] == 'Mt_tRNA': continue
        if f.type == feature_type:
                feature_id = f.attr[id_attribute]
                features[f.iv] += feature_id
                if feature_id not in geneRange:
                        geneRange[feature_id] = [f.iv.chrom,0,0,f.iv.strand]
                if geneRange[feature_id][1] != 0:
                        geneRange[feature_id][1] = min(geneRange[feature_id][1],f.iv.start)
                else:
                        geneRange[feature_id][1] = f.iv.start
                geneRange[feature_id][2] = max(geneRange[feature_id][2],f.iv.end)

        i += 1
        if i % 100000 == 0:
                sys.stderr.write("%d GFF lines processed.\r" % i)

##geneFeatures
for g,v in geneRange.iteritems():
        chrom,s,e,strand = v
        tmp_iv = HTSeq.GenomicInterval(chrom,s,e,strand)
        geneFeatures[tmp_iv] += g

bam = HTSeq.BAM_Reader(bamFile)
#bamOut = HTSeq.BAM_Writer.from_BAM_Reader("empty.bam",bam)
uniqueGene = 0
uniqueNotGene = 0
uniqueIntron = 0
uniqueAmbiguous = 0

# read length statistic
RNA_len = Counter()
DNA_len = Counter()
Intron_len = Counter()

i=0
for r in bam:
        i += 1
        if i % 1000000==0: sys.stderr.write("%d BAM alignments records processed. \r" % i)


        if not r.aligned:
                continue

        if r.optional_field("NH") == 1:
                r.read.seq = r.read_as_aligned.seq
                r.read.qual = r.read_as_aligned.qual
                iv_seq = ( co.ref_iv for co in r.cigar if co.type == "M" and co.size > 0 )
                fs = set()
                intron_fs = set()
                for iv in iv_seq:
                        for iv2, fs2 in features[iv].steps():
                                fs = fs.union(fs2)
                        for iv3,fs3 in geneFeatures[iv].steps():
                                intron_fs = intron_fs.union(fs3)

                if fs is None or len(fs) == 0:
                        if intron_fs is None or len(intron_fs) == 0:
                                uniqueNotGene += 1
                                DNA_len[len(r.read.seq)] += 1
                        else:
                                uniqueIntron += 1
                                Intron_len[len(r.read.seq)] += 1

                elif len(fs) == 1:
                        uniqueGene += 1
                        RNA_len[len(r.read.seq)] += 1
                else:
                        uniqueAmbiguous += 1
                        RNA_len[len(r.read.seq)] += 1


with open(sys.argv[3],"w") as fout:
        fout.write("unique mapped reads of RNA: %i\n" % uniqueGene)
        fout.write("unique mapped reads of DNA: %i\n" % uniqueNotGene)
        fout.write("unique mapped reads of Intron: %i\n" % uniqueIntron)
        fout.write("unique mapped ambiguous reads of RNA: %i\n" % uniqueAmbiguous)

#plot
def myhistPlot(mydict,name):
        labels = sorted(mydict.keys())
        value_list = [mydict[i] for i in labels]
        plb.bar(labels,value_list,width=0.5)
        plb.xlabel("length of the reads")
        plb.ylabel("NO. of the reads")
        plb.title(name)
        plb.savefig(name + ".png")
        plb.close()
myhistPlot(RNA_len,sys.argv[4]+"_RNA")
myhistPlot(DNA_len,sys.argv[4]+"_DNA")
myhistPlot(Intron_len,sys.argv[4]+"_Intron")