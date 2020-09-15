
# DESeq 的用法。

## 目录

• DESeqDataSet - build the dataset, see tximeta & tximport packages for preparing input 
• DESeq - perform differential analysis  
• results - build a results table    
• lfcShrink - estimate shrunken LFC (posterior estimates) using apeglm & ashr pakges • vst - apply variance stabilizing transformation, e.g. for PCA or sample clustering  
• Plots, e.g.: plotPCA, plotMA, plotCounts

```r
data_origin=read.table('/home/test/merge.tsv',header=T,sep="\t",row.names=1)
data_matrix=data.matrix(data_origin)
cond=factor(rep(c('A','B'),3))
dds <- DESeqDataSetFromMatrix(data_matrix, DataFrame(cond), ~ cond)
dds <- DESeq(dds)

res <- results(dds, contrast=c("cond","B","A"))

write.table(res, file='/home/test/diff.tsv', quote=FALSE, sep='\t')
```
  
## 将数据导入

```R
cnts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10) 
cond <- factor(rep(1:2, each=5))

```

## 差异表达分析
```R
# see vignette for suggestions on generating
# count tables from RNA-Seq data
cnts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10) 
cond <- factor(rep(1:2, each=5))
# object construction
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
# standard analysis dds <- DESeq(dds) res <- results(dds)

dds <- DESeq(dds)

res <- results(dds, contrast=c("cond","2","1"))

```

## 结果导出为表格
```R
res <- results(dds, contrast=c("cond","2","1"))
```


## plot PCA


```R

# using rlog transformed data:
dds <- makeExampleDESeqDataSet(betaSD=1) rld <- rlog(dds)
plotPCA(rld)

cnts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10) 
cond <- factor(rep(1:2, each=5))
# object construction
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)


# also possible to perform custom transformation:
dds <- estimateSizeFactors(dds)
# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                               colData=colData(dds))
# the call to DESeqTransform() is needed to
# trigger our plotPCA method. 
plotPCA( DESeqTransform( se ) )


```

## plot MA


```R
dds <- makeExampleDESeqDataSet() dds <- DESeq(dds)
plotMA(dds)
res <- results(dds)
plotMA(res)


```

## plot Counts


```R

dds <- makeExampleDESeqDataSet() 
plotCounts(dds, "gene1")

```
