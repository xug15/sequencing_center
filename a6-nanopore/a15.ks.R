##
## load function
transfor=function(denature_time,denature_name){
denature_time_t=t(denature_time)
denature_name_t=t(denature_name)

denature_t=cbind(denature_name_t,denature_time_t)

coln=dim(denature_t)[2]-1
colnames(denature_t)=c('ref',1:coln)
denature_t=denature_t[-1,]
denature_t[denature_t[,1]=="TRUE",1]="T"
return(denature_t)
}
### function for Z-score standardization
zscore_nor=function(x){
  (x-mean(x))/sd(x)
}

#load file with Rdata.
library(dgof)
load("a16.denature_time.RData")
args <- commandArgs(trailingOnly = TRUE)
#tab
filename1=args[1]
#filename1='ENST00000361624.tab.tsv'
#ref
filename2=args[2]
#filename2='ENST00000361624.fa.tsv'
#
rna_time=read.table(filename1)
rna_name=read.table(filename2)
rna_t=transfor(rna_time,rna_name)


#
position=dim(rna_t)[1]
totalnum=dim(rna_t)[2]
#
shiftt2= data.frame(pos=c(),DF=c(),pv=c())
#loop 
for (b in 1:position){
#b=1
content=rna_t[b,1]
rna_pos_time=rna_t[b,2:totalnum]
rna_pos_time2=zscore_nor(as.numeric(rna_pos_time[!is.na(rna_pos_time)]))
#print(b)
#print(rna_pos_time2)
if(sum(is.na(rna_pos_time2)) < 1 && sum(is.nan(rna_pos_time2)) < 1 ){
if(content=='A'){
s2=ks.test(timeA,rna_pos_time2,alternative = c( "greater"))
}else if(content=='T'){
s2=ks.test(timeT,rna_pos_time2,alternative = c( "greater"))
}else if(content=='C'){
s2=ks.test(timeC,rna_pos_time2,alternative = c( "greater"))
}else if(content=='G'){
s2=ks.test(timeG,rna_pos_time2,alternative = c( "greater"))
 }
}
  # get info unmod vs denature
  #
if(exists('s2')){
   div2=s2['statistic'][[1]][[1]]
   dip2=s2['p.value'][[1]]
   shiftt2=rbind(shiftt2,c(b,div2,dip2))
 }else{
 shiftt2=rbind(shiftt2,c(b,0,1))

 }
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
write.table(ksdata2,paste(filename1,'ks.tsv',sep=''),sep="\t")
