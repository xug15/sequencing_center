# normalization

```R
# normolizaton

## read data.
counts=read.table('/Users/xugang/Documents/c-pycharm/git/sequencing_center/result/result_191014_guoruixin/merge.counter',header=T,row.names=1)
head(counts)
#normolize by cpm.
counts_cpm=cpm(counts)
head(counts_cpm)
#read table
write.csv(counts_cpm,'/Users/xugang/Documents/c-pycharm/git/sequencing_center/result/result_191014_guoruixin/cpm.csv')

```

