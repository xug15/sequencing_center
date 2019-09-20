# Differential Expression with Limma-Voom

## 0. Intro

limma is an R package that was originally developed for differential expression (DE) analysis of microarray data.

voom is a function in the limma package that modifies RNA-Seq data for use with limma.

Together they allow fast, flexible, and powerful analyses of RNA-Seq data. Limma-voom is our tool of choice for DE analyses because it:

Allows for incredibly flexible model specification (you can include multiple categorical and continuous variables, allowing incorporation of almost any kind of metadata)

Based on simulation studies, maintains the false discovery rate at or below the nominal rate, unlike some other packages

Empirical Bayes smoothing of gene-wise standard deviations provides increased power.

## 1. Setup
Input data for this example is on the course github page.

First, install the edgeR package if not already installed (which installs limma as a dependency)
```R
# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")
```
Load the edgeR package (which loads limma as a dependency)
```R
library(edgeR)

## Loading required package: limma
```
Read in the counts table
```R
counts <- read.delim("all_counts.txt", row.names = 1)
head(counts)
```
```R
##           C61 C62 C63 C64 C91  C92 C93 C94 I561 I562 I563 I564 I591 I592
## AT1G01010 341 371 275 419 400  542 377 372  677  522  455  508  821  466
## AT1G01020 164  94 176 155 200  183 166 115  172  157  122  152  189  171
## AT1G03987   0   0   0   0   0    0   0   0    0    0    0    0    0    0
## AT1G01030  20  34  40  27  28   36  22  40   20    7   57   38   25   10
## AT1G01040 738 487 610 690 945 1033 836 908  857  821  770  751  848  607
## AT1G03993   1   0   0   0   0    0   0   1    3    0    1    1    1    0
##           I593 I594 I861 I862 I863 I864 I891 I892 I893 I894
## AT1G01010  691  500  157  473  459  228  590  491  565  496
## AT1G01020  163  185   46  162  119   53  172  212  169  157
## AT1G03987    0    0    0    0    0    0    0    0    0    0
## AT1G01030   17   26   49   17   24   48   27   28   47   32
## AT1G01040  871  756  361  618  641  439  783  692  768  625
## AT1G03993    0    0    0    1    0    1    3    2    1    1
```

OR read the file directly from the github page:
```R
counts <- read.delim("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2018-June-RNA-Seq-Workshop/master/thursday/all_counts.txt")
head(counts)
```
```R
##           C61 C62 C63 C64 C91 C92 C93 C94 I561 I562 I563 I564 I591 I592
## AT1G01010 289 317 225 343 325 449 310 299  563  438  380  407  678  386
## AT1G01020 127  78 142 130 156 146 144  95  138  129   99  118  154  140
## AT1G03987   0   0   0   0   0   0   0   0    0    0    0    0    0    0
## AT1G01030  17  25  32  24  22  25  21  35   18    6   46   33   19    8
## AT1G01040 605 415 506 565 762 854 658 753  704  692  641  601  704  508
## AT1G03993   1   1   0   0   0   0   1   1    3    0    1    1    1    0
##           I593 I594 I861 I862 I863 I864 I891 I892 I893 I894
## AT1G01010  567  421  130  411  382  190  501  390  480  407
## AT1G01020  127  154   35  132   97   46  139  175  137  123
## AT1G03987    0    0    0    0    0    0    0    0    0    0
## AT1G01030   14   21   37   16   19   38   24   18   37   23
## AT1G01040  733  614  297  521  542  381  651  573  650  550
## AT1G03993    0    0    0    1    0    1    3    1    1    1
```
Create DGEList object
```R
d0 <- DGEList(counts)
```

## 2. Preprocessing
Calculate normalization factors
```R
d0 <- calcNormFactors(d0)
d0
```
```R
## An object of class "DGEList"
## $counts
##           C61 C62 C63 C64 C91 C92 C93 C94 I561 I562 I563 I564 I591 I592
## AT1G01010 289 317 225 343 325 449 310 299  563  438  380  407  678  386
## AT1G01020 127  78 142 130 156 146 144  95  138  129   99  118  154  140
## AT1G03987   0   0   0   0   0   0   0   0    0    0    0    0    0    0
## AT1G01030  17  25  32  24  22  25  21  35   18    6   46   33   19    8
## AT1G01040 605 415 506 565 762 854 658 753  704  692  641  601  704  508
##           I593 I594 I861 I862 I863 I864 I891 I892 I893 I894
## AT1G01010  567  421  130  411  382  190  501  390  480  407
## AT1G01020  127  154   35  132   97   46  139  175  137  123
## AT1G03987    0    0    0    0    0    0    0    0    0    0
## AT1G01030   14   21   37   16   19   38   24   18   37   23
## AT1G01040  733  614  297  521  542  381  651  573  650  550
## 34257 more rows ...
## 
## $samples
##     group lib.size norm.factors
## C61     1 10502901    1.0404937
## C62     1  9423745    0.9719425
## C63     1  9437115    1.0401597
## C64     1 10415490    1.0360806
## C91     1 10169158    1.0468854
## 19 more rows ...
```
Note: calcNormFactors doesn’t normalize the data, it just calculates normalization factors for use downstream.

Filter low-expressed genes
```R
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left
## [1] 21080    24
```
“Low-expressed” is subjective and depends on the dataset.

Derive experiment information from the sample names

Our experiment has two factors, cultivar (“C”, “I5”, or “I8”) and time (6 or 9)

The sample names are the cultivar, followed by the time, followed by the replicate
```R
snames <- colnames(counts) # Sample names
snames
```
```R
##  [1] "C61"  "C62"  "C63"  "C64"  "C91"  "C92"  "C93"  "C94"  "I561" "I562"
## [11] "I563" "I564" "I591" "I592" "I593" "I594" "I861" "I862" "I863" "I864"
## [21] "I891" "I892" "I893" "I894"
```
```R
cultivar <- substr(snames, 1, nchar(snames) - 2) 
time <- substr(snames, nchar(snames) - 1, nchar(snames) - 1)
cultivar
```
```R
##  [1] "C"  "C"  "C"  "C"  "C"  "C"  "C"  "C"  "I5" "I5" "I5" "I5" "I5" "I5"
## [15] "I5" "I5" "I8" "I8" "I8" "I8" "I8" "I8" "I8" "I8"
```
```R
time
```
```R
##  [1] "6" "6" "6" "6" "9" "9" "9" "9" "6" "6" "6" "6" "9" "9" "9" "9" "6"
## [18] "6" "6" "6" "9" "9" "9" "9"
```
Create a new variable “group” that combines cultivar and time
```R
group <- interaction(cultivar, time)
group
```
```R
##  [1] C.6  C.6  C.6  C.6  C.9  C.9  C.9  C.9  I5.6 I5.6 I5.6 I5.6 I5.9 I5.9
## [15] I5.9 I5.9 I8.6 I8.6 I8.6 I8.6 I8.9 I8.9 I8.9 I8.9
## Levels: C.6 I5.6 I8.6 C.9 I5.9 I8.9
```
Note: you can also enter group information manually, or read it in from an external file. If you do this, it is very,very,very important that you make sure the metadata is in the same order as the column names of the counts table.

Multidimensional scaling (MDS) plot
```R
plotMDS(d, col = as.numeric(group))
```

## 3. Voom transformation and calculation of variance weights
Specify the model to be fitted. We do this before using voom since voom uses variances of the model residuals (observed - fitted)
```R
mm <- model.matrix(~0 + group)
```
The above specifies a model where each coefficient corresponds to a group mean

Voom
```R
y <- voom(d, mm, plot = T)
```

What is voom doing?

Counts are transformed to log2 counts per million reads (CPM), where “per million reads” is defined based on the normalization factors we calculated earlier
A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot above)
The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs.
More details at https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29

The above is a “good” voom plot. If your voom plot looks like the below, you might want to filter more:
```R
tmp <- voom(d0, mm, plot = T)
```

## 4. Fitting linear models in limma
lmFit fits a linear model using weighted least squares for each gene:
```R
fit <- lmFit(y, mm)
head(coef(fit))
```
```R
##           groupC.6 groupI5.6 groupI8.6 groupC.9 groupI5.9 groupI8.9
## AT1G01010 4.837410 5.3738532  5.065354 5.043214 5.5240004  5.363809
## AT1G01020 3.530869 3.4993460  3.212860 3.689622 3.7209961  3.736297
## AT1G01030 1.250817 0.9293832  1.559242 1.285596 0.4831707  1.215591
## AT1G01040 5.676015 5.9469878  5.778889 6.182374 5.8641107  5.815498
## AT1G01050 6.598712 6.5013631  6.463936 6.619239 6.7532886  6.711772
## AT1G01060 7.807988 7.4624783  7.390741 8.966047 8.2706387  8.376129
```
Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models:

Specify which groups to compare:

Comparison between times 6 and 9 for cultivar I5
```R
contr <- makeContrasts(groupI5.9 - groupI5.6, levels = colnames(coef(fit)))
contr
```
```R
##            Contrasts
## Levels      groupI5.9 - groupI5.6
##   groupC.6                      0
##   groupI5.6                    -1
##   groupI8.6                     0
##   groupC.9                      0
##   groupI5.9                     1
##   groupI8.9                     0
```
Estimate contrast for each gene
```R
tmp <- contrasts.fit(fit, contr)
```
Empirical Bayes smoothing of standard errors (shrinks standard errors that are much larger or smaller than those from other genes towards the average standard error) (see https://www.degruyter.com/doi/10.2202/1544-6115.1027)
```R
tmp <- eBayes(tmp)
```
What genes are most differentially expressed?
```R
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
```
```R
##                logFC  AveExpr          t      P.Value    adj.P.Val
## AT5G37260  3.0480867 6.964609  24.154998 1.193341e-16 2.515563e-12
## AT3G02990  1.6484885 3.304597  13.605334 8.573977e-12 9.036971e-08
## AT2G29500 -5.0342224 5.525802 -12.062120 7.937586e-11 5.577477e-07
## AT3G24520  1.8715741 5.882965  11.710003 1.360008e-10 7.167244e-07
## AT3G46230 -6.7068648 4.544494 -11.189910 3.081581e-10 1.299195e-06
## AT5G65630  1.0468903 7.553775  10.657591 7.329270e-10 2.575017e-06
## AT3G24100 -1.2759156 5.988295 -10.237187 1.485118e-09 4.021633e-06
## AT1G53540 -7.1847583 4.316965 -10.160183 1.693866e-09 4.021633e-06
## AT4G21320 -2.8641453 4.705136 -10.152258 1.717016e-09 4.021633e-06
## AT5G48570 -4.0078194 6.422341 -10.036624 2.094832e-09 4.415906e-06
## AT5G52882  0.8786555 7.649621   9.828722 3.007116e-09 5.398811e-06
## AT5G02810 -4.8872284 3.419336  -9.816290 3.073327e-09 5.398811e-06
## AT1G09140 -1.2469280 7.027052  -9.496095 5.419786e-09 8.788392e-06
## AT4G24780  1.9288845 6.979568   9.380927 6.666788e-09 1.003828e-05
## AT5G59570  1.5315276 3.005741   9.330689 7.300769e-09 1.026001e-05
## AT4G11880  2.7436628 2.919357   9.183289 9.547854e-09 1.245739e-05
## AT4G25500  1.0517676 6.143738   9.155503 1.004628e-08 1.245739e-05
## AT1G22770 -3.5882843 3.953036  -9.090404 1.132259e-08 1.326002e-05
## AT1G31230  1.1574470 6.114996   8.721535 2.252449e-08 2.499033e-05
## AT2G27820 -0.8709118 5.735808  -8.680724 2.433146e-08 2.564536e-05
##                   B
## AT5G37260 26.622628
## AT3G02990 16.006216
## AT2G29500 14.645646
## AT3G24520 14.385545
## AT3G46230 12.519533
## AT5G65630 12.802559
## AT3G24100 12.109394
## AT1G53540 10.864734
## AT4G21320 11.825499
## AT5G48570 11.752325
## AT5G52882 11.430101
## AT5G02810 10.430908
## AT1G09140 10.854555
## AT4G24780 10.651922
## AT5G59570 10.146392
## AT4G11880  9.781987
## AT4G25500 10.250838
## AT1G22770  9.958652
## AT1G31230  9.459672
## AT2G27820  9.385693
```
```R
logFC: log2 fold change of I5.9/I5.6
AveExpr: Average expression across all samples, in log2 CPM
t: logFC divided by its standard error
P.Value: Raw p-value (based on t) from test that logFC differs from 0
adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
B: log-odds that gene is DE (arguably less useful than the other columns)
AT5G37260 has higher expression at time 9 than at time 6 (logFC is positive). AT2G29500 has lower expression at time 9 than at time 6 (logFC is negative).
```
How many DE genes are there?
```R
length(which(top.table$adj.P.Val < 0.05))
## [1] 4680
```
Write top.table to a file
```R
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "time9_v_time6_I5.txt", row.names = F, sep = "\t", quote = F)
```
Let’s say we want to compare cultivars C and I5 at time 6. The only thing we have to change is the call to makeContrasts:
```R
contr <- makeContrasts(groupI5.6 - groupC.6, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
```
```R
##                logFC    AveExpr          t      P.Value    adj.P.Val
## AT4G12520 -9.3396792  0.5671202 -14.118616 4.276991e-12 9.015898e-08
## AT3G30720  5.6550000  3.5130064  11.628491 1.543229e-10 1.626564e-06
## AT5G26270  2.3722487  4.4417048   9.728356 3.587045e-09 2.520497e-05
## AT3G33528 -4.7012468 -1.6468149  -8.175366 6.441093e-08 3.394456e-04
## AT3G05955 -3.7995906 -1.6958534  -7.210032 4.545141e-07 1.916232e-03
## AT4G28100 -0.8099448  4.5746806  -7.041358 6.477323e-07 1.986777e-03
## AT1G64795 -4.7147903 -1.1051093  -7.032657 6.597456e-07 1.986777e-03
## AT3G25730  1.3872163  3.9538977   6.242201 3.650963e-06 9.620287e-03
## AT5G05480 -0.4869486  4.6067884  -5.647317 1.392471e-05 3.005506e-02
## AT2G18193  1.0100187  3.8571257   5.626785 1.459359e-05 3.005506e-02
## AT4G01870  1.6373400  5.6170176   5.595305 1.568338e-05 3.005506e-02
## AT2G14878 -0.5162040  6.4733722  -5.531829 1.814030e-05 3.186647e-02
## AT2G06995 -3.2115972 -2.3585729  -5.439207 2.244863e-05 3.640132e-02
## AT4G15248 -1.6999189  2.1452779  -5.392221 2.501939e-05 3.673680e-02
## AT1G62280 -1.7551000  2.7349000  -5.373243 2.614099e-05 3.673680e-02
## AT3G28740  2.2032568  5.5413479   5.284877 3.207698e-05 4.226142e-02
## AT5G14370  0.8031693  3.7422520   5.174352 4.147430e-05 4.748115e-02
## AT2G25737  0.8732454  3.6798171   5.163677 4.251876e-05 4.748115e-02
## AT2G37760  0.7703328  5.4139388   5.141986 4.472409e-05 4.748115e-02
## AT2G30400  1.0087276  2.0253964   5.138887 4.504853e-05 4.748115e-02
##                     B
## AT4G12520  6.42577720
## AT3G30720  9.54391447
## AT5G26270 11.00209007
## AT3G33528  1.56894746
## AT3G05955  2.19224634
## AT4G28100  6.18609294
## AT1G64795  2.33667151
## AT3G25730  4.53222502
## AT5G05480  3.26213083
## AT2G18193  3.23531092
## AT4G01870  3.11677860
## AT2G14878  2.94274140
## AT2G06995  0.02144583
## AT4G15248  2.49861111
## AT1G62280  2.67518742
## AT3G28740  2.42959865
## AT5G14370  2.25381821
## AT2G25737  2.23281223
## AT2G37760  2.08993315
## AT2G30400  2.09909628
```
```R
length(which(top.table$adj.P.Val < 0.05)) # number of DE genes
## [1] 20
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "I5_v_C_time6.txt", row.names = F, sep = "\t", quote = F)
```
What if we refit our model as a two-factor model (rather than using the group variable)?

Create new model matrix:
```R
mm <- model.matrix(~cultivar*time)
```
We are specifying that model includes effects for cultivar, time, and the cultivar-time interaction (which allows the differences between cultivars to differ across time)
```R
colnames(mm)
## [1] "(Intercept)"      "cultivarI5"       "cultivarI8"      
## [4] "time9"            "cultivarI5:time9" "cultivarI8:time9"
y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
head(coef(fit))
```
```R
##           (Intercept)  cultivarI5 cultivarI8      time9 cultivarI5:time9
## AT1G01010    4.837410  0.53644370  0.2279446 0.20580445      -0.05565729
## AT1G01020    3.530869 -0.03152318 -0.3180096 0.15875297       0.06289715
## AT1G01030    1.250817 -0.32143420  0.3084243 0.03477863      -0.48099113
## AT1G01040    5.676015  0.27097286  0.1028739 0.50635951      -0.58923660
## AT1G01050    6.598712 -0.09734846 -0.1347759 0.02052702       0.23139851
## AT1G01060    7.807988 -0.34550979 -0.4172467 1.15805850      -0.34989810
##           cultivarI8:time9
## AT1G01010       0.09265044
## AT1G01020       0.36468449
## AT1G01030      -0.37842909
## AT1G01040      -0.46975071
## AT1G01050       0.22730960
## AT1G01060      -0.17267051
```
The coefficient cultivarI5 represents the difference in mean expression between cultivar I5 and the reference cultivar (cultivar C), for time 6 (the reference level for time)
The coefficient time9 represents the difference in mean expression between time 9 and time 6, for cultivar C
The coefficient cultivarI5:time9 is the difference between times 9 and 6 of the differences between cultivars I5 and C (interaction effect)
Let’s estimate the difference between cultivars I5 and C at time 6
```R
tmp <- contrasts.fit(fit, coef = 2) # Directly test second coefficient
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
```
```R
##                logFC    AveExpr          t      P.Value    adj.P.Val
## AT4G12520 -9.3396792  0.5671202 -14.118616 4.276991e-12 9.015898e-08
## AT3G30720  5.6550000  3.5130064  11.628491 1.543229e-10 1.626564e-06
## AT5G26270  2.3722487  4.4417048   9.728356 3.587045e-09 2.520497e-05
## AT3G33528 -4.7012468 -1.6468149  -8.175366 6.441093e-08 3.394456e-04
## AT3G05955 -3.7995906 -1.6958534  -7.210032 4.545141e-07 1.916232e-03
## AT4G28100 -0.8099448  4.5746806  -7.041358 6.477323e-07 1.986777e-03
## AT1G64795 -4.7147903 -1.1051093  -7.032657 6.597456e-07 1.986777e-03
## AT3G25730  1.3872163  3.9538977   6.242201 3.650963e-06 9.620287e-03
## AT5G05480 -0.4869486  4.6067884  -5.647317 1.392471e-05 3.005506e-02
## AT2G18193  1.0100187  3.8571257   5.626785 1.459359e-05 3.005506e-02
## AT4G01870  1.6373400  5.6170176   5.595305 1.568338e-05 3.005506e-02
## AT2G14878 -0.5162040  6.4733722  -5.531829 1.814030e-05 3.186647e-02
## AT2G06995 -3.2115972 -2.3585729  -5.439207 2.244863e-05 3.640132e-02
## AT4G15248 -1.6999189  2.1452779  -5.392221 2.501939e-05 3.673680e-02
## AT1G62280 -1.7551000  2.7349000  -5.373243 2.614099e-05 3.673680e-02
## AT3G28740  2.2032568  5.5413479   5.284877 3.207698e-05 4.226142e-02
## AT5G14370  0.8031693  3.7422520   5.174352 4.147430e-05 4.748115e-02
## AT2G25737  0.8732454  3.6798171   5.163677 4.251876e-05 4.748115e-02
## AT2G37760  0.7703328  5.4139388   5.141986 4.472409e-05 4.748115e-02
## AT2G30400  1.0087276  2.0253964   5.138887 4.504853e-05 4.748115e-02
##                     B
## AT4G12520  6.42577720
## AT3G30720  9.54391447
## AT5G26270 11.00209007
## AT3G33528  1.56894746
## AT3G05955  2.19224634
## AT4G28100  6.18609294
## AT1G64795  2.33667151
## AT3G25730  4.53222502
## AT5G05480  3.26213083
## AT2G18193  3.23531092
## AT4G01870  3.11677860
## AT2G14878  2.94274140
## AT2G06995  0.02144583
## AT4G15248  2.49861111
## AT1G62280  2.67518742
## AT3G28740  2.42959865
## AT5G14370  2.25381821
## AT2G25737  2.23281223
## AT2G37760  2.08993315
## AT2G30400  2.09909628
```
```R
length(which(top.table$adj.P.Val < 0.05)) # number of DE genes
## [1] 20
```
We get the same results as with the model where each coefficient corresponded to a group mean. In essence, these are the same model, so use whichever is most convenient for what you are estimating.

The interaction effects cultivarI5:time9 and cultivarI8:time9 are easier to estimate and test in this setup
```R
head(coef(fit))
```
```R
##           (Intercept)  cultivarI5 cultivarI8      time9 cultivarI5:time9
## AT1G01010    4.837410  0.53644370  0.2279446 0.20580445      -0.05565729
## AT1G01020    3.530869 -0.03152318 -0.3180096 0.15875297       0.06289715
## AT1G01030    1.250817 -0.32143420  0.3084243 0.03477863      -0.48099113
## AT1G01040    5.676015  0.27097286  0.1028739 0.50635951      -0.58923660
## AT1G01050    6.598712 -0.09734846 -0.1347759 0.02052702       0.23139851
## AT1G01060    7.807988 -0.34550979 -0.4172467 1.15805850      -0.34989810
##           cultivarI8:time9
## AT1G01010       0.09265044
## AT1G01020       0.36468449
## AT1G01030      -0.37842909
## AT1G01040      -0.46975071
## AT1G01050       0.22730960
## AT1G01060      -0.17267051
tmp <- contrasts.fit(fit, coef = 5) # Test cultivarI5:time9
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
```
```R
##                logFC  AveExpr         t      P.Value    adj.P.Val
## AT5G38200  5.2681384 5.459342 12.169571 6.750884e-11 1.423086e-06
## AT5G06860  3.9331835 6.449356  7.667986 1.770961e-07 1.306381e-03
## AT1G08430  7.3873178 5.946190  7.640080 1.874141e-07 1.306381e-03
## AT1G76990 -0.9683284 8.651885 -7.503011 2.478901e-07 1.306381e-03
## AT5G18170  1.4739946 6.811722  7.145965 5.197419e-07 2.191232e-03
## AT5G24030  1.5597102 6.422377  6.737394 1.238209e-06 3.821427e-03
## AT1G62280  3.3975934 2.734900  6.652058 1.488488e-06 3.821427e-03
## AT2G41380  3.7331437 5.866261  6.615737 1.610271e-06 3.821427e-03
## AT2G45360 -4.8763311 1.093918 -6.592316 1.694190e-06 3.821427e-03
## AT2G46790  3.2686459 3.970742  6.530515 1.937883e-06 3.821427e-03
## AT5G15950  3.2415341 5.516374  6.517393 1.994103e-06 3.821427e-03
## AT5G67520 -1.8307730 3.290640 -6.362042 2.802299e-06 4.649936e-03
## AT1G01650 -1.1996420 5.524074 -6.351573 2.867608e-06 4.649936e-03
## AT5G15330 -1.3784824 4.115181 -6.268550 3.444132e-06 5.065418e-03
## AT1G48100  1.9296680 2.376553  6.247994 3.604425e-06 5.065418e-03
## AT4G19160 -2.3207711 7.235367 -6.179190 4.198820e-06 5.356807e-03
## AT1G05890 -0.6529356 6.455076 -6.160149 4.380439e-06 5.356807e-03
## AT5G56010  1.2941382 8.452257  6.130758 4.676715e-06 5.356807e-03
## AT4G37870  0.8880398 8.434194  6.116456 4.828242e-06 5.356807e-03
## AT1G58110 -0.7654977 5.850866 -6.004302 6.205156e-06 6.540234e-03
##                   B
## AT5G38200 13.058923
## AT5G06860  7.282910
## AT1G08430  6.287578
## AT1G76990  7.111815
## AT5G18170  6.408931
## AT5G24030  5.581114
## AT1G62280  4.943495
## AT2G41380  5.221262
## AT2G45360  3.503402
## AT2G46790  4.956204
## AT5G15950  5.114242
## AT5G67520  4.615689
## AT1G01650  4.778464
## AT5G15330  4.553861
## AT1G48100  3.973489
## AT4G19160  4.412505
## AT1G05890  4.372387
## AT5G56010  4.309468
## AT4G37870  4.278843
## AT1G58110  4.041773
```
```R
length(which(top.table$adj.P.Val < 0.05)) 
## [1] 882
```
The log fold change here is the difference between cultivarI5 and cultivar C in the log fold changes between times 9 and 6. It is ALSO the difference between times 9 and 6 in the log fold changes between cultivarI5 and cultivar C.

A gene for which this interaction effect is significant is one for which the effect of time differs between cultivars, and for which the effect of cultivar differs between times.

## 5. More complicated models
Specifying a different model is simply a matter of changing the calls to model.matrix (and possibly to contrasts.fit).

Let’s say we have information on the RNA extraction batch:
```R
batch <- factor(rep(rep(1:2, each = 2), 6))
batch
##  [1] 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2
## Levels: 1 2
```
To adjust for batch in the analysis, add batch to the end of the call to model matrix. Everything else about the code stays the same:
```R
mm <- model.matrix(~0 + group + batch)
y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
contr <- makeContrasts(groupI5.6 - groupC.6, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
```
```R
##                logFC    AveExpr          t      P.Value    adj.P.Val
## AT4G12520 -9.1849480  0.5671202 -13.894622 1.254284e-11 2.644031e-07
## AT3G30720  5.6553502  3.5130064  11.626230 2.902705e-10 3.059451e-06
## AT5G26270  2.3756632  4.4417048  10.532157 1.568494e-09 1.102129e-05
## AT3G33528 -4.6930275 -1.6468149  -8.114771 1.054174e-07 5.555497e-04
## AT1G64795 -4.6966105 -1.1051093  -7.585967 2.920305e-07 1.231200e-03
## AT4G28100 -0.8142969  4.5746806  -7.332425 4.823635e-07 1.694704e-03
## AT3G05955 -3.7806772 -1.6958534  -7.158687 6.837456e-07 2.059051e-03
## AT1G62280 -1.7706304  2.7349000  -6.183451 5.223128e-06 1.278365e-02
## AT2G18193  1.0171979  3.8571257   6.155724 5.544006e-06 1.278365e-02
## AT3G25730  1.3861828  3.9538977   6.114089 6.064349e-06 1.278365e-02
## AT4G34138  1.2583243  6.2839069   5.874446 1.020572e-05 1.955787e-02
## AT1G67890  0.4262907  5.1117436   5.759619 1.312886e-05 2.306302e-02
## AT3G63430 -1.2020149  3.7073650  -5.617984 1.794992e-05 2.838999e-02
## AT5G05480 -0.4876253  4.6067884  -5.595805 1.885483e-05 2.838999e-02
## AT4G15248 -1.7089812  2.1452779  -5.492499 2.372592e-05 3.079319e-02
## AT2G14878 -0.5161993  6.4733722  -5.474049 2.472296e-05 3.079319e-02
## AT2G06995 -3.1599964 -2.3585729  -5.445134 2.637233e-05 3.079319e-02
## AT4G01870  1.6376531  5.6170176   5.444359 2.641803e-05 3.079319e-02
## AT1G15380 -0.6846986  7.6756786  -5.421227 2.782070e-05 3.079319e-02
## AT3G45980 -0.3280034  7.2242831  -5.399375 2.921555e-05 3.079319e-02
##                    B
## AT4G12520  5.6334654
## AT3G30720  8.9070420
## AT5G26270 11.6801080
## AT3G33528  1.1444487
## AT1G64795  2.4478920
## AT4G28100  6.4756663
## AT3G05955  1.7078155
## AT1G62280  4.1194174
## AT2G18193  4.1644124
## AT3G25730  4.0704267
## AT4G34138  3.5467936
## AT1G67890  3.3320337
## AT3G63430  3.0358961
## AT5G05480  2.9985093
## AT4G15248  2.5246096
## AT2G14878  2.6674208
## AT2G06995 -0.2708462
## AT4G01870  2.6397139
## AT1G15380  2.5436455
## AT3G45980  2.4978614
```
```R
length(which(top.table$adj.P.Val < 0.05))
## [1] 27
```
What if we want to adjust for a continuous variable like RIN score:
```R
# Generate example RIN data
set.seed(99)
RIN <- rnorm(n = 24, mean = 7.5, sd = 1)
RIN
```
```R
##  [1] 7.713963 7.979658 7.587829 7.943859 7.137162 7.622674 6.636155
##  [8] 7.989624 7.135883 6.205758 6.754231 8.421550 8.250054 4.991446
## [15] 4.459066 7.500266 7.105981 5.754972 7.998631 7.770954 8.598922
## [22] 8.252513 7.440583 7.155431
```
Model adjusting for RIN score
```R
mm <- model.matrix(~0 + group + RIN)
y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
contr <- makeContrasts(groupI5.6 - groupC.6, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
```
```R
##                logFC    AveExpr          t      P.Value    adj.P.Val
## AT3G30720  5.6296669  3.5130064  18.264852 8.228016e-14 1.078373e-09
## AT4G12520 -9.2914672  0.5671202 -18.052498 1.023125e-13 1.078373e-09
## AT5G26270  2.4136835  4.4417048  10.642693 1.300777e-09 9.140129e-06
## AT3G33528 -4.9109883 -1.6468149  -9.609516 7.125167e-09 3.754963e-05
## AT1G64795 -4.7654967 -1.1051093  -7.254218 5.606683e-07 2.363777e-03
## AT4G28100 -0.8023348  4.5746806  -6.994622 9.486720e-07 3.333001e-03
## AT3G25730  1.4185849  3.9538977   6.546417 2.402778e-06 7.235795e-03
## AT4G15248 -1.7173926  2.1452779  -6.178114 5.259531e-06 1.385886e-02
## AT3G28740  2.3364878  5.5413479   5.900899 9.592908e-06 2.114492e-02
## AT3G05955 -3.6893043 -1.6958534  -5.880481 1.003080e-05 2.114492e-02
## AT4G01870  1.6906107  5.6170176   5.821660 1.141042e-05 2.126536e-02
## AT1G60750  1.8497644  3.8830018   5.794730 1.210552e-05 2.126536e-02
## AT1G17180  1.5714528  4.2264956   5.523743 2.205086e-05 3.575631e-02
## AT2G18193  1.0022545  3.8571257   5.472489 2.472195e-05 3.722420e-02
## AT5G07870  0.8857028  5.0986202   5.401616 2.896986e-05 3.817990e-02
## AT3G27940  1.4469855 -0.4298164   5.401474 2.897905e-05 3.817990e-02
## AT4G34135  1.0569416  6.3882644   5.238191 4.183945e-05 4.792741e-02
## AT5G48010  2.1658813  5.3481653   5.225912 4.301563e-05 4.792741e-02
## AT5G05480 -0.4793068  4.6067884  -5.191084 4.653849e-05 4.792741e-02
## AT4G34138  1.3466558  6.2839069   5.184933 4.719053e-05 4.792741e-02
##                   B
## AT3G30720 17.859502
## AT4G12520 11.272023
## AT5G26270 12.140983
## AT3G33528  3.415722
## AT1G64795  3.228510
## AT4G28100  5.862248
## AT3G25730  4.969968
## AT4G15248  4.106498
## AT3G28740  3.577535
## AT3G05955  0.636844
## AT4G01870  3.399840
## AT1G60750  3.421066
## AT1G17180  2.815864
## AT2G18193  2.734741
## AT5G07870  2.518925
## AT3G27940  1.914480
## AT4G34135  2.102290
## AT5G48010  2.121253
## AT5G05480  2.095733
## AT4G34138  1.986345
```
```R
length(which(top.table$adj.P.Val < 0.05))
## [1] 28
```
What if we want to look at the correlation of gene expression with a continuous variable like pH?
```R
# Generate example pH data
set.seed(99)
pH <- rnorm(n = 24, mean = 8, sd = 1.5)
pH
```
```R
##  [1] 8.320944 8.719487 8.131743 8.665788 7.455743 8.184011 6.704232
##  [8] 8.734436 7.453825 6.058637 6.881346 9.382326 9.125082 4.237169
## [15] 3.438599 8.000399 7.408972 5.382459 8.747947 8.406431 9.648382
## [22] 9.128770 7.910875 7.483147
```
Specify model matrix:
```R
mm <- model.matrix(~pH)
head(mm)
```
```R
##   (Intercept)       pH
## 1           1 8.320944
## 2           1 8.719487
## 3           1 8.131743
## 4           1 8.665788
## 5           1 7.455743
## 6           1 8.184011
y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
tmp <- contrasts.fit(fit, coef = 2) # test "pH" coefficient
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
```
```R
##                 logFC     AveExpr         t      P.Value adj.P.Val
## AT3G63220  0.10930007  3.72231850  3.984813 0.0005224914 0.9979999
## AT4G22790 -0.18156791  2.62688789 -3.924496 0.0006095427 0.9979999
## AT3G25855 -0.18661169  0.39590825 -3.888223 0.0006686278 0.9979999
## AT5G48060 -0.14380178  3.08645012 -3.729192 0.0010015260 0.9979999
## AT1G60030 -0.18144388  1.90130680 -3.688498 0.0011101065 0.9979999
## AT5G50290 -0.19527796  0.87276139 -3.683545 0.0011240846 0.9979999
## AT1G69588 -0.18176455  1.07793155 -3.430115 0.0021225677 0.9979999
## AT1G17250  0.18143048  0.91048376  3.423565 0.0021574364 0.9979999
## AT5G66350 -0.18602484  0.59789247 -3.376519 0.0024247787 0.9979999
## AT5G49320 -0.17105888  1.05400394 -3.292867 0.0029815543 0.9979999
## AT3G17730 -0.16801297  1.56449303 -3.282393 0.0030594390 0.9979999
## AT5G61420 -0.18139376  5.05638487 -3.280444 0.0030741424 0.9979999
## AT4G05430 -0.20982687 -0.33505919 -3.272344 0.0031360082 0.9979999
## AT1G51460 -0.26599783 -0.19515578 -3.249322 0.0033184675 0.9979999
## AT1G59730 -0.26312342  2.24935746 -3.235196 0.0034354577 0.9979999
## AT2G38920 -0.25174501  0.02998871 -3.228661 0.0034909211 0.9979999
## AT4G11655  0.30473405 -0.14566023  3.217560 0.0035871102 0.9979999
## AT1G63020 -0.11884059  3.38076346 -3.210234 0.0036519841 0.9979999
## AT2G16400 -0.09043316  4.49774646 -3.198184 0.0037611636 0.9979999
## AT1G17240  0.22499668  1.07480792  3.181368 0.0039188025 0.9979999
##                    B
## AT3G63220 -0.3281633
## AT4G22790 -0.5717697
## AT3G25855 -1.5833958
## AT5G48060 -0.8575577
## AT1G60030 -1.1988239
## AT5G50290 -1.5961103
## AT1G69588 -1.8900360
## AT1G17250 -2.2539850
## AT5G66350 -2.1461291
## AT5G49320 -2.1077044
## AT3G17730 -1.9798846
## AT5G61420 -1.6446842
## AT4G05430 -2.6782487
## AT1G51460 -2.5644227
## AT1G59730 -1.8823323
## AT2G38920 -2.5302070
## AT4G11655 -3.0281427
## AT1G63020 -1.8275637
## AT2G16400 -1.8149709
## AT1G17240 -2.5144550
```
```R
length(which(top.table$adj.P.Val < 0.05))
## [1] 0
```
In this case, limma is fitting a linear regression model, which here is a straight line fit, with the slope and intercept defined by the model coefficients:
```R
AT1G60030 <- y$E["AT1G60030",]
plot(AT1G60030 ~ pH, ylim = c(0, 3.5))
intercept <- coef(fit)["AT1G60030", "(Intercept)"]
slope <- coef(fit)["AT1G60030", "pH"]
abline(a = intercept, b = slope)
```

In this example, the log fold change logFC is the slope of the line, or the change in gene expression (on the log2 CPM scale) for each unit increase in pH.

Here, a logFC of -0.19 means a 0.19 log2 CPM decrease in gene expression for each unit increase in pH, or a 14% decrease on the CPM scale (2^0.19 = 1.14).

## 6. A bit more on linear models
Limma fits a linear model to each gene.

Linear models include analysis of variance (ANOVA) models, linear regression, and any model of the form

Y=β0+β1X1+β2X2+⋯+βpXp+ϵ
The covariates X can be:

a continuous variable (pH, RIN score, age, weight, temperature, etc.)
Dummy variables coding a categorical covariate (like cultivar, time, and group)
The β’s are unknown parameters to be estimated.

In limma, the β’s are the log fold changes.

The error (residual) term ϵ is assumed to be normally distributed with a variance that is constant across the range of the data.

Normally distributed means the residuals come from a distribution that looks like this: 

The log2 transformation that voom applies to the counts makes the data “normal enough”, but doesn’t completely stabilize the variance:
```R
tmp <- voom(d, mm, plot = T)
```

The log2 counts per million are more variable at lower expression levels. The variance weights calculated by voom address this situation.
```R
sessionInfo()
```
```R
## R version 3.5.0 (2018-04-23)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 16299)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] edgeR_3.22.2 limma_3.36.1
## 
## loaded via a namespace (and not attached):
##  [1] locfit_1.5-9.1  Rcpp_0.12.17    lattice_0.20-35 digest_0.6.15  
##  [5] rprojroot_1.3-2 grid_3.5.0      backports_1.1.2 magrittr_1.5   
##  [9] evaluate_0.10.1 stringi_1.1.7   rmarkdown_1.10  tools_3.5.0    
## [13] stringr_1.3.1   yaml_2.1.19     compiler_3.5.0  htmltools_0.3.6
## [17] knitr_1.20
```