library(dgof)
dir='output/'

## load read each position
NAI_N3_1_unmod_time =read.table(paste(dir,"/a7-mod-sum/Tetra_NAI_N3_1.unmod.tsv",sep=''),sep="\t",header=F)
denature_unmod_time =read.table(paste(dir,"/a7-mod-sum/Tetra_NAI_N3_denatured.unmod.tsv",sep=''),sep="\t",header=F)

# rename colomn 
colnames(denature_unmod_time)=c('read_name',seq(0,439,by=1))
colnames(NAI_N3_1_unmod_time)=c('read_name',seq(0,439,by=1))

### function for Z-score standardization
zscore_nor=function(x){
  (x-mean(x))/sd(x)
}

#zscore
shiftt2= data.frame(pos=c(),DF=c(),pv=c())
for (b in 2:418){

  #print(b)
  denature=zscore_nor(denature_unmod_time[,b][!is.na(denature_unmod_time[,b])])
  nai_n3=zscore_nor(NAI_N3_1_unmod_time[,b][!is.na(NAI_N3_1_unmod_time[,b])])
  #nai_n3
  #nai_n3 vs denature.
  s2=ks.test(denature,nai_n3,alternative = c( "greater"))
  # get info unmod vs denature
  #
  div2=s2['statistic'][[1]][[1]]
  dip2=s2['p.value'][[1]]
  shiftt2=rbind(shiftt2,c(b-2,div2,dip2))
}

colnames(shiftt2)=c("position","df","pvalue")
## set sig
shiftt2$sig=0
shiftt2$sig[shiftt2$df>0.15 & shiftt2$pvalue<0.01]=1
## set df and -log10 pvalue
shiftt2$dfpv=shiftt2$df*(-log10(shiftt2$pvalue+10^-270))

## smooth dfpv
shiftt2$dfpv_sm=filter(shiftt2$dfpv,filter=c(rep(1/4,4)))

#smooth
shiftt2$lin<-filter(shiftt2$df,filter=c(rep(1/4,4)))

#trans data to ksdata.
ksdata2=shiftt2

#save file
write.table(ksdata2,paste(dir,'/a10-ks/nai_n3_ks.tsv',sep=''),sep="\t")



