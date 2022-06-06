# metagenome 




```r
#
install.packages('vegan')
library(ggplot2)
library(vegan)
library(reshape)
dataplot=data.frame(id=c('CK1','CK2','CK3','LAB1','LAB2','LAB3','M1','M2','M3','LAB+M1','LAB+M2','LAB+M3'),
  type=c('CK','CK','CK','LAB','LAB','LAB','M','M','M','LAB+M','LAB+M','LAB+M'),
                    Shannon=c(3.0886,2.0661,3.1925,0.6484,1.7249,0.8461,2.1726,2.1694,1.9401,0.438,0.4275,0.4977))

dataplot$type=factor(dataplot$type,levels=unique(dataplot$type))
head(dataplot)



ggplot(dataplot, aes(x=type, y=Shannon, fill=type)) +
  geom_boxplot(fill=c("#B8D9E9", "#5E91BF", "#98B682","#72B161"),color="black")+
  scale_fill_manual(values=c("black", "grey", "white"))+
  theme(axis.text=element_text(size=22,face="bold"),
                       axis.title=element_text(size=24,face="bold"),
                       legend.title=element_text(size=24,face="bold"), 
                       legend.text=element_text(size=22,face="bold"),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
        )

#1.检验数据是否符合正态分布
shapiro.test(dataplot$Shannon)
## w 接近1，p大于0.05,符合正态分布。

#2.单因素方差分析
fit=aov(dataplot$Shannon~dataplot$type)
summary(fit)
##pr 小于0.01，有差异。

#3.方差齐性分析
bartlett.test(Shannon~type,data=dataplot)
## p 大于0.01，不存在离群点。

#4.多重比较
TukeyHSD(fit)

# 
'docker pull r-base:3.4.0
docker run -dt -v /home/xugang:/home/xugang --name=r.3.4 r-base:3.4.0
docker exec -it r.3.4 bash
R
install.packages("vegan")
install.packages("ape")
install.packages("ggplot2")
library(vegan)
library(ape)
library(ggplot2)

df<-varespec
rownames(df)<-paste0("site",1:24)
#首先是将数据集赋值给新的变量，并以site1-24对新的数据集的行进行命名
df<-varespec
rownames(df)<-paste0("site",1:24)
计算距离
bray_dist<-vegdist(df,method = "bray")
使用ape这个包中的pcoa()函数做PCoA分析
library(ape)
df.pcoa<-pcoa(bray_dist,correction = "cailliez")
df.pcoa$vectors能够获得用于画图的数据
df.pcoa$values可以获得坐标轴上显示的百分比
最后用ggplot2来画这个图
df.plot<-data.frame(df.pcoa$vectors)
head(df.plot)
library(ggplot2)
x_label<-round(df.pcoa$values$Rel_corr_eig[1]*100,2)
y_label<-round(df.pcoa$values$Rel_corr_eig[2]*100,2)
x_label
y_label
ggplot(data=df.plot,aes(x=Axis.1,y=Axis.2))+
geom_point()+
theme_bw()+
theme(panel.grid = element_blank())+
geom_vline(xintercept = 0,lty="dashed")+
geom_hline(yintercept = 0,lty="dashed")+
labs(x=paste0("PCoA1 ",x_label,"%"),
y=paste0("PCoA2 ",y_label,"%"))
通过上图我们可以看到这些样地大体上可以分为两组，如果自己手头有样地的分组数据就可以看看这个结果是不是和自己的分组数据一致。

下面人为的给他分个组，然后添加一个表示分组的椭圆
df.plot$group<-ifelse(df.plot$Axis.1<0,"AAA","BBB")
ggplot(data=df.plot,aes(x=Axis.1,y=Axis.2,
color=group,shape=group))+
geom_point(size=5)+
theme_bw()+
theme(panel.grid = element_blank())+
geom_vline(xintercept = 0,lty="dashed")+
geom_hline(yintercept = 0,lty="dashed")+
labs(x=paste0("PCoA1 ",x_label,"%"),
y=paste0("PCoA2 ",y_label,"%"))+
stat_ellipse(data=df.plot,
geom = "polygon",
aes(fill=group),
alpha=0.3)+
scale_fill_manual(values = c("#e31a1c","#1f78b4"))
df=tread.table("/Users/xugang/Desktop/sequencing_center_desktop/yanglab/简历/文章用微生物数据/2.abundance.txt",header=T)
rownames(df)=df[,1]
df2=t(df[,-1])
bray_dist<-vegdist(df2,method = "bray")
library(ape)
df.pcoa<-pcoa(bray_dist,correction = "cailliez")
df.pcoa$vectors能够获得用于画图的数据
df.pcoa$values $Rel_corr_eig 可以获得坐标轴上显示的百分比
save(df.pcoa,file="/Users/xugang/Desktop/sequencing_center_desktop/yanglab/简历/文章用微生物数据/2.pcoa.RData")
load(df.pcoa,file="/Users/xugang/Desktop/sequencing_center_desktop/yanglab/简历/文章用微生物数据/2.pcoa.RData")
'
load(file="/Users/xugang/Desktop/sequencing_center_desktop/yanglab/简历/文章用微生物数据/2.pcoa.RData")
df.pcoa
df.plot<-data.frame(df.pcoa$vectors)
head(df.plot)
library(ggplot2)
#percentage 1,2 
x_label<-round(df.pcoa$values$Rel_corr_eig[1]*100,2)
y_label<-round(df.pcoa$values$Rel_corr_eig[2]*100,2)
x_label
y_label
ggplot(data=df.plot,aes(x=Axis.1,y=Axis.2))+
  geom_point(color=c("#B8D9E9","#B8D9E9","#B8D9E9", "#5E91BF","#5E91BF","#5E91BF", "#98B682","#98B682","#98B682","#72B161","#72B161","#72B161"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PCoA1 ",x_label,"%"),
       y=paste0("PCoA2 ",y_label,"%"))


df.plot$group=type=c('CK','CK','CK','LAB','LAB','LAB','M','M','M','LAB+M','LAB+M','LAB+M')
head(df.plot)
ggplot(data=df.plot,aes(x=Axis.1,y=Axis.2,
                        color=group))+
  geom_point(size=3,color=c("#B8D9E9","#B8D9E9","#B8D9E9", "#5E91BF","#5E91BF","#5E91BF", "#98B682","#98B682","#98B682","#72B161","#72B161","#72B161"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PCoA1 ",x_label,"%"),
       y=paste0("PCoA2 ",y_label,"%"))


pcoabox=data.frame(value=df.plot$Axis.1
                   ,type=c('CK','CK','CK','LAB','LAB','LAB','M','M','M','LAB+M','LAB+M','LAB+M'))
pcoabox$type=factor(pcoabox$type,level=unique(pcoabox$type))

ggplot(pcoabox, aes(x=type, y=value, fill=type)) +
  geom_boxplot(fill=c("#B8D9E9", "#5E91BF", "#98B682","#72B161"),color="black")+
  scale_fill_manual(values=c("black", "grey", "white"))+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=4,face="bold"),
        legend.title=element_text(size=4,face="bold"), 
        legend.text=element_text(size=2,face="bold"),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
  )

#
pcoabox=data.frame(value=df.plot$Axis.2
                   ,type=c('CK','CK','CK','LAB','LAB','LAB','M','M','M','LAB+M','LAB+M','LAB+M'))
pcoabox$type=factor(pcoabox$type,level=unique(pcoabox$type))

ggplot(pcoabox, aes(x=type, y=value, fill=type)) +
  geom_boxplot(fill=c("#B8D9E9", "#5E91BF", "#98B682","#72B161"),color="black")+
  scale_fill_manual(values=c("black", "grey", "white"))+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=4,face="bold"),
        legend.title=element_text(size=4,face="bold"), 
        legend.text=element_text(size=2,face="bold"),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
  )
#

stock=read.table("/Users/xugang/Desktop/sequencing_center_desktop/yanglab/简历/文章用微生物数据/stackplot.txt",header=T,row.names = 1)
#row.names(stock)=stock[,1]
stock2=stock
head(stock2)
stock2$CK=apply(stock2[,1:3],1,FUN=mean)
stock2$LAB=apply(stock2[,4:6],1,FUN=mean)
stock2$M=apply(stock2[,7:9],1,FUN=mean)
stock2$LAB_M=apply(stock2[,10:12],1,FUN=mean)
head(stock2)
#import data and modify it
phylum=stock2[,13:16]
phylum.ave <- apply(phylum, 1, FUN=mean)

phylum.2 <- cbind(phylum, phylum.ave)[order(-phylum.ave),]

#sum(phylum.2$phylum.ave[1:4])  #check the summary of top 4 taxonomies
phylum.2 <- subset(phylum.2, select=-phylum.ave)
phylum.3 <- rbind(phylum.2, others=apply(phylum.2, 2, function(x){1-sum(x)}))

#apply(phylum.2, 2, function(x){sum(x)})  #check the summary of each column
phylum.3 <- cbind(PhylumID=row.names(phylum.3), phylum.3)

#convert the data format
library(reshape2)
phylum.gg <- melt(phylum.3, id.vars="PhylumID", variable.name="SampleID", value.name="Abundance")

#sum(phypum.gg$Abundance) #check
#plot bar chart
library(ggplot2)
phylum.gg$color=1
unique(phylum.gg$PhylumID)
phylum.gg[phylum.gg$PhylumID=='Lactiplantibacillus_pentosus',]$color='#49AA79'
phylum.gg[phylum.gg$PhylumID=='Enterobacter_hormaechei',]$color='#CE8844'
phylum.gg[phylum.gg$PhylumID=='Enterococcus_mundtii',]$color='#8BA156'
phylum.gg[phylum.gg$PhylumID=='Lactiplantibacillus_coryniformis',]$color='#CB7D7D'
phylum.gg[phylum.gg$PhylumID=='Peptoniphilus_stercorisuis',]$color='#4697AB'
phylum.gg[phylum.gg$PhylumID=='Enterococcus_faecalis',]$color='#A8514A'
phylum.gg[phylum.gg$PhylumID=='Lactiplantibacillus_acidipiscis',]$color='#6A6599'
phylum.gg[phylum.gg$PhylumID=='Enterococcus_flavescens',]$color='#4B78A5'
phylum.gg[phylum.gg$PhylumID=='uncultured_bacterium_g_Garciella',]$color='#3CBBBB'
phylum.gg[phylum.gg$PhylumID=='uncultured_bacterium_g_Anaerosporobacter',]$color='#F2716D'
phylum.gg[phylum.gg$PhylumID=='Weissella_paramesenteroides',]$color='#DC7A47'
#Enterococcus_flavescens
phylum.gg[phylum.gg$PhylumID=='others',]$color='#5A8CBF'
phylum.gg
ggplot(phylum.gg, aes(SampleID, Abundance, fill=PhylumID)) +
  geom_bar(stat="identity",fill=phylum.gg$color) +
  guides(fill=guide_legend(reverse=F)) +
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=14,face="bold"),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text=element_text(size=12,face="bold"))

#
species2=stock
titlev='Lactiplantibacillus_pentosus'
titlev='Enterobacter_hormaechei'
titlev='Enterococcus_mundtii'
titlev='Lactiplantibacillus_coryniformis'

species3=species2[titlev,]

species4=melt(species3)
species4$type=c('CK','CK','CK','LAB','LAB','LAB','M','M','M','LAB+M','LAB+M','LAB+M')
species4$type=factor(species4$type,unique(species4$type))
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
species5=data_summary(species4,varname = 'value',groupnames = c("type"))
head(species5)
ggplot(species5, aes(x=type, y=value, fill=type)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9))+
  labs(title=titlev)+
  scale_fill_manual(values=c("#B8D9E9", "#5E91BF", "#98B682","#72B161"))+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        legend.title=element_text(size=24,face="bold"),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text=element_text(size=22,face="bold"))

# environment and bacteria abundence.
abundance=read.table("/Users/xugang/Desktop/sequencing_center_desktop/yanglab/简历/文章用微生物数据/2.abundance.txt",header=T,row.names = 1)
  
environment=read.table("/Users/xugang/Desktop/sequencing_center_desktop/yanglab/简历/文章用微生物数据/3.environment.txt",header=T,row.names = 1)
head(environment)
head(abundance)
cor.test(mtcars$mpg, mtcars$hp,
         method = "spearman")
#
head(environment)
head(abundance)
#
dim(environment)
dim(abundance)

cortable <- setNames(data.frame(matrix(ncol = dim(environment)[2], nrow = dim(abundance)[1])), colnames(environment))
rownames(cortable)=rownames(abundance)
pvtable <- setNames(data.frame(matrix(ncol = dim(environment)[2], nrow = dim(abundance)[1])), colnames(environment))
rownames(pvtable)=rownames(abundance)



for (i in row.names(abundance) ){
  print(i)
  print(abundance[i,])
  typeof(abundance[i,])
  for(j in colnames(environment)){
    print(j)
    #print(environment[,j])
    corres=cor.test(unlist(abundance[i,]), environment[,j],method = "spearman")
    cortable[i,j]=corres$estimate
    pvtable[i,j]=corres$p.value
    
  }
}

cortable
pvtable
cortable2 <- melt(as.matrix(cortable))
pvtable2 <- melt(as.matrix(pvtable))

colnames(cortable2) <- c("y", "x", "value")
colnames(pvtable2) <- c("y", "x", "value")
cortable2$y=factor(cortable2$y,levels = unique(cortable2$y))
cortable2$x=factor(cortable2$x,levels = unique(cortable2$x))
pvtable2$y=factor(pvtable2$y,levels = unique(pvtable2$y))
pvtable2$x=factor(pvtable2$x,levels = unique(pvtable2$x))
ggplot(cortable2, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#39B5B5",
                       mid = "#FFFFFF",
                       high = "#F16864") +
  coord_fixed()+theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1),
                      axis.text=element_text(size=10,face="bold"),
                      axis.title=element_text(size=14,face="bold"),
                      legend.title=element_text(size=14,face="bold"), 
                      legend.text=element_text(size=12,face="bold"))

pvtable3=pvtable2
pvtable3$value[pvtable3$value>0.05]=1
ggplot(pvtable3, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#F16864",
                       mid = "#FFFFFF",
                       high = "#39B5B5") +
  geom_text(aes(label = round(value, 4))) +
  coord_fixed()+theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=1),)


# install.packages("reshape")
library(reshape)

# Data 
set.seed(8)
m <- matrix(round(rnorm(200), 2), 10, 10)
colnames(m) <- paste("Col", 1:10)
rownames(m) <- paste("Row", 1:10)

# Transform the matrix in long format
df <- melt(m)
colnames(df) <- c("x", "y", "value")



```

```r

library('ggplot2')
dir1='/Users/xugang/Desktop/sequencing_center_desktop/yanglab/metagenome/bactaria/'
#dir2='/Users/xugang/Desktop/sequencing_center_desktop/2020/rna_structure/a2-human/structure_analysis/'

meta_data=read.table(paste(dir1,'b1.meta.tsv',sep=''),header = T)

head(meta_data)

#RS vs BS
dim(meta_data)[1]


resultb=data.frame(df1=c(0),pv1=c(0),df2=c(0),pv2=c(0))
colnames(resultb)=c("df1","pv1","df2","pv2")
for ( ii in 1:dim(meta_data)[1])
  #for ( ii in 1:2)
{
  group1=c(meta_data[ii,]$B_C_8,meta_data[ii,]$B_C_10,meta_data[ii,]$B_C_1,meta_data[ii,]$B_C_3,meta_data[ii,]$B_C_4,meta_data[ii,]$B_C_6,meta_data[ii,]$B_H_3,meta_data[ii,]$B_H_10,meta_data[ii,]$B_S_1_10,meta_data[ii,]$B_T_4,meta_data[ii,]$B_T_2,meta_data[ii,]$B_H_7,meta_data[ii,]$B_H_6,meta_data[ii,]$B_H_1,meta_data[ii,]$B_C_9,meta_data[ii,]$B_H_4,meta_data[ii,]$B_H_8,meta_data[ii,]$B_H_9,meta_data[ii,]$B_H_5,meta_data[ii,]$B_C_7,meta_data[ii,]$B_C_2,meta_data[ii,]$B_T_7,meta_data[ii,]$B_T_5,meta_data[ii,]$B_T_1,meta_data[ii,]$B_C_5,meta_data[ii,]$B_T_6,meta_data[ii,]$B_S_1_6,meta_data[ii,]$B_T_8,meta_data[ii,]$B_T_3,meta_data[ii,]$B_T_10,meta_data[ii,]$B_T_9,meta_data[ii,]$B_S_1_4,meta_data[ii,]$B_S_2_3,meta_data[ii,]$B_H_2,meta_data[ii,]$B_S_1_3,meta_data[ii,]$B_S_2_2,meta_data[ii,]$B_S_2_4,meta_data[ii,]$B_S_2_5,meta_data[ii,]$B_S_2_8,meta_data[ii,]$B_S_1_2,meta_data[ii,]$B_S_2_7,meta_data[ii,]$B_S_1_8,meta_data[ii,]$B_S_2_1,meta_data[ii,]$B_S_2_9,meta_data[ii,]$B_S_2_6,meta_data[ii,]$B_S_1_7,meta_data[ii,]$B_S_1_5,meta_data[ii,]$B_S_1_1,meta_data[ii,]$B_S_2_10,meta_data[ii,]$B_S_1_9)
  #group2=c(meta_data[ii,]$RH_C_4,meta_data[ii,]$RH_C_8,meta_data[ii,]$RH_C_5,meta_data[ii,]$RH_C_2,meta_data[ii,]$RH_C_1,meta_data[ii,]$RH_C_6,meta_data[ii,]$RH_C_7,meta_data[ii,]$RH_C_9,meta_data[ii,]$RH_S_1_6,meta_data[ii,]$RH_C_10,meta_data[ii,]$RH_S_1_1,meta_data[ii,]$RH_T_6,meta_data[ii,]$RH_H_4,meta_data[ii,]$RH_H_7,meta_data[ii,]$RH_H_2,meta_data[ii,]$RH_T_5,meta_data[ii,]$RH_H_5,meta_data[ii,]$RH_H_8,meta_data[ii,]$RH_H_9,meta_data[ii,]$RH_H_3,meta_data[ii,]$RH_S_2_8,meta_data[ii,]$RH_T_7,meta_data[ii,]$RH_H_10,meta_data[ii,]$RH_T_8,meta_data[ii,]$RH_T_1,meta_data[ii,]$RH_T_3,meta_data[ii,]$RH_C_3,meta_data[ii,]$RH_T_4,meta_data[ii,]$RH_T_10,meta_data[ii,]$RH_S_1_10,meta_data[ii,]$RH_S_1_3,meta_data[ii,]$RH_S_1_2,meta_data[ii,]$RH_T_2,meta_data[ii,]$RH_T_9,meta_data[ii,]$RH_S_2_2,meta_data[ii,]$RH_S_1_4,meta_data[ii,]$RH_H_1,meta_data[ii,]$RH_S_2_9,meta_data[ii,]$RH_S_2_5,meta_data[ii,]$RH_S_2_6,meta_data[ii,]$RH_S_2_4,meta_data[ii,]$RH_S_1_7,meta_data[ii,]$RH_S_2_10,meta_data[ii,]$RH_S_2_3,meta_data[ii,]$RH_H_6,meta_data[ii,]$RH_S_2_7,meta_data[ii,]$RH_S_1_9,meta_data[ii,]$RH_S_2_1,meta_data[ii,]$RH_S_1_8,meta_data[ii,]$RH_S_1_5)
  #group2=c(meta_data[ii,]$R_C_5,meta_data[ii,]$R_S_2_8,meta_data[ii,]$R_H_5,meta_data[ii,]$R_C_10,meta_data[ii,]$R_C_1,meta_data[ii,]$R_C_8,meta_data[ii,]$R_H_10,meta_data[ii,]$R_S_2_1,meta_data[ii,]$R_H_3,meta_data[ii,]$R_S_1_8,meta_data[ii,]$R_T_4,meta_data[ii,]$R_H_9,meta_data[ii,]$R_C_9,meta_data[ii,]$R_S_1_10,meta_data[ii,]$R_T_3,meta_data[ii,]$R_C_6,meta_data[ii,]$R_C_7,meta_data[ii,]$R_S_1_2,meta_data[ii,]$R_S_2_10,meta_data[ii,]$R_C_3,meta_data[ii,]$R_T_8,meta_data[ii,]$R_T_5,meta_data[ii,]$R_H_2,meta_data[ii,]$R_T_2,meta_data[ii,]$R_S_1_1,meta_data[ii,]$R_H_8,meta_data[ii,]$R_S_2_6,meta_data[ii,]$R_S_1_6,meta_data[ii,]$R_S_2_3,meta_data[ii,]$R_H_1,meta_data[ii,]$R_T_1,meta_data[ii,]$R_C_4,meta_data[ii,]$R_S_2_9,meta_data[ii,]$R_C_2,meta_data[ii,]$R_S_2_7,meta_data[ii,]$R_T_7,meta_data[ii,]$R_T_10,meta_data[ii,]$R_H_6,meta_data[ii,]$R_T_6,meta_data[ii,]$R_S_2_2,meta_data[ii,]$R_S_1_4,meta_data[ii,]$R_S_2_5,meta_data[ii,]$R_S_1_3,meta_data[ii,]$R_H_4,meta_data[ii,]$R_S_1_5,meta_data[ii,]$R_T_9,meta_data[ii,]$R_H_7,meta_data[ii,]$R_S_1_9,meta_data[ii,]$R_S_2_4,meta_data[ii,]$R_S_1_7)
  #group2=c(meta_data[ii,]$N_D_8,meta_data[ii,]$N_H_3,meta_data[ii,]$N_T_7,meta_data[ii,]$N_H_9,meta_data[ii,]$N_D_9,meta_data[ii,]$N_H_10,meta_data[ii,]$N_H_4,meta_data[ii,]$N_C_2,meta_data[ii,]$N_C_7,meta_data[ii,]$N_D_1,meta_data[ii,]$N_T_6,meta_data[ii,]$N_D_5,meta_data[ii,]$N_C_4,meta_data[ii,]$N_H_6,meta_data[ii,]$N_C_8,meta_data[ii,]$N_D_2,meta_data[ii,]$N_C_6,meta_data[ii,]$N_H_7,meta_data[ii,]$N_C_3,meta_data[ii,]$N_C_9,meta_data[ii,]$N_T_5,meta_data[ii,]$N_T_10,meta_data[ii,]$N_H_8,meta_data[ii,]$N_H_2,meta_data[ii,]$N_T_3,meta_data[ii,]$N_H_5,meta_data[ii,]$N_D_7,meta_data[ii,]$N_C_1,meta_data[ii,]$N_C_10,meta_data[ii,]$N_T_9,meta_data[ii,]$N_D_10,meta_data[ii,]$N_D_4,meta_data[ii,]$N_Y_1,meta_data[ii,]$N_D_6,meta_data[ii,]$N_D_3,meta_data[ii,]$N_H_1,meta_data[ii,]$N_T_1,meta_data[ii,]$N_T_8,meta_data[ii,]$N_T_4,meta_data[ii,]$N_Y_6,meta_data[ii,]$N_Y_9,meta_data[ii,]$N_T_2,meta_data[ii,]$N_Y_4,meta_data[ii,]$N_Y_8,meta_data[ii,]$N_Y_2,meta_data[ii,]$N_Y_5,meta_data[ii,]$N_Y_3,meta_data[ii,]$N_Y_7,meta_data[ii,]$N_Y_10,meta_data[ii,]$N_C_5)
  group2=c(meta_data[ii,]$F_C_7,meta_data[ii,]$F_H_10,meta_data[ii,]$F_T_4,meta_data[ii,]$F_H_5,meta_data[ii,]$F_T_9,meta_data[ii,]$F_T_1,meta_data[ii,]$F_H_4,meta_data[ii,]$F_H_8,meta_data[ii,]$F_C_4,meta_data[ii,]$F_C_10,meta_data[ii,]$F_T_2,meta_data[ii,]$F_H_3,meta_data[ii,]$F_C_5,meta_data[ii,]$F_H_2,meta_data[ii,]$F_H_7,meta_data[ii,]$F_T_5,meta_data[ii,]$F_D_2,meta_data[ii,]$F_T_7,meta_data[ii,]$F_C_2,meta_data[ii,]$F_H_1,meta_data[ii,]$F_T_8,meta_data[ii,]$F_T_10,meta_data[ii,]$F_C_9,meta_data[ii,]$F_C_6,meta_data[ii,]$F_T_6,meta_data[ii,]$F_C_8,meta_data[ii,]$F_C_1,meta_data[ii,]$F_C_3,meta_data[ii,]$F_D_4,meta_data[ii,]$F_D_1,meta_data[ii,]$F_T_3,meta_data[ii,]$F_D_9,meta_data[ii,]$F_D_8,meta_data[ii,]$F_D_5,meta_data[ii,]$F_D_10,meta_data[ii,]$F_Y_3,meta_data[ii,]$F_D_3,meta_data[ii,]$F_Y_7,meta_data[ii,]$F_H_9,meta_data[ii,]$F_Y_1,meta_data[ii,]$F_Y_8,meta_data[ii,]$F_D_7,meta_data[ii,]$F_Y_2,meta_data[ii,]$F_Y_5,meta_data[ii,]$F_D_6,meta_data[ii,]$F_Y_10,meta_data[ii,]$F_H_6,meta_data[ii,]$F_Y_6,meta_data[ii,]$F_Y_4,meta_data[ii,]$F_Y_9)
  
  wilt=wilcox.test(group1,group2) 
  wilt2=wilcox.test(group2,group1) 
  result=c(wilt$statistic,wilt$p.value,wilt2$statistic,wilt2$p.value)
  #print(result)
  resultb=rbind(resultb,result)
}

resultb2=resultb[-1,]

write.csv(resultb2,paste(dir1,"bs_vs_rhi.csv",sep=''))
write.csv(resultb2,paste(dir1,"bs_vs_rootendosphere.csv",sep=''))
write.csv(resultb2,paste(dir1,"bs_vs_leefendosphere.csv",sep=''))
write.csv(resultb2,paste(dir1,"bs_vs_phylloplane.csv",sep=''))

Rhizo_vs_soil_deplet=dim(na.omit(resultb2[resultb2$pv1<0.05 & resultb2$df1>1 & resultb2$df1 > resultb2$df2,]))[1]
Rhizo_vs_soil_enrich=dim(na.omit(resultb2[resultb2$pv2<0.05 & resultb2$df2>1 & resultb2$df2 > resultb2$df1,]))[1]
Rootendo_vs_soil_deplet=dim(na.omit(resultb2[resultb2$pv1<0.05 & resultb2$df1>1 & resultb2$df1 > resultb2$df2,]))[1]
Rootendo_vs_soil_enrich=dim(na.omit(resultb2[resultb2$pv2<0.05 & resultb2$df2>1 & resultb2$df2 > resultb2$df1,]))[1]
leefendo_vs_soil_deplet=dim(na.omit(resultb2[resultb2$pv1<0.05 & resultb2$df1>1 & resultb2$df1 > resultb2$df2,]))[1]
leefendo_vs_soil_enrich=dim(na.omit(resultb2[resultb2$pv2<0.05 & resultb2$df2>1 & resultb2$df2 > resultb2$df1,]))[1]
phylloplane_vs_soil_deplet=dim(na.omit(resultb2[resultb2$pv1<0.05 & resultb2$df1>1 & resultb2$df1 > resultb2$df2,]))[1]
phylloplane_vs_soil_enrich=dim(na.omit(resultb2[resultb2$pv2<0.05 & resultb2$df2>1 & resultb2$df2 > resultb2$df1,]))[1]


num=c(Rhizo_vs_soil_enrich,Rhizo_vs_soil_deplet,Rootendo_vs_soil_enrich,Rootendo_vs_soil_deplet,leefendo_vs_soil_enrich,leefendo_vs_soil_deplet,phylloplane_vs_soil_enrich,phylloplane_vs_soil_deplet)
location=c('Rhizopl.vs Soil','Rhizopl.vs Soil','Root Endosp. vs Soil','Root Endosp. vs Soil',"Leef Endosp. vs Soil","Leef Endosp. vs Soil","Phyllop. vs Soil","Phyllop. vs Soil")
group=c("Enriched","Depleted","Enriched","Depleted","Enriched","Depleted","Enriched","Depleted")

bsdata=data.frame(num,location,group)



ggplot(bsdata, aes(x=location, y=num, fill=group)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=num), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(plot.title = element_text( size = 24),
        axis.text.x = element_text( size = 20),
        axis.text.y = element_text( size = 20),  
        axis.title.x = element_text( size = 20),
        axis.title.y = element_text( size = 20))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                       panel.background = element_blank(), axis.line = element_line(colour = "black"))
#################################################################################################
#
#
#################################################################################################
dir1='/Users/xugang/Desktop/sequencing_center_desktop/yanglab/metagenome/'
#dir1='/Users/xugang/Desktop/sequencing_center_desktop/yanglab/metagenome/bactaria/'
#dir2='/Users/xugang/Desktop/sequencing_center_desktop/2020/rna_structure/a2-human/structure_analysis/'
#set average degree table.
avgtable=data.frame(name=character(),node=integer(),positive_edge=integer(),negative_edge=integer(),avg_degree=double(),edge=integer())

meta_data=read.table(paste(dir1,'b1.meta.tsv',sep=''),header = T)
head(meta_data)
#B_
#Bulk Soil
#nametitle='Bulk Soil'

#RH_
#Rhizophsere Soil
#nametitle='Rhizophsere Soil'

#R_
#Root Endosphere
#nametitle='Root Endosphere'

#F_
#Phylloplane
#nametitle='Phylloplane'

#N_
#Leaf Endosphere
#nametitle='Leaf Endosphere'



nametitle='Leaf Endosphere'

## get bulk soil
select_col_name=c('p')
for  (i in colnames(meta_data) ){
  #print(i)
  if( grepl('^N_',i) ){
    print(i)
    select_col_name=append(select_col_name,i)
  }
}
#head(meta_data)
#select_col_name
meta_select=meta_data[,select_col_name]
rownames(meta_select)=meta_data$name
#head(meta_select)
#rowSums(meta_select[,-1])

meta_select=meta_select[rowSums(meta_select[,-1])>0,]
meta_select=meta_select[meta_select[,2]>0,]
dim(meta_select)
groupv=c()
sizev=c()
idv=c()
id=0

#meta_select[1,1][1]

for  (i in rownames(meta_select) ){
  #
  id=id+1
  idv=append(idv,id)
  sizev=append(sizev,1)
  #print(meta_select[i,c('p')])
  groupv=append(groupv,as.character(meta_select[i,c('p')]))
}
#groupv
#sizev
#idv
nodes_nkk=data.frame(id=idv,label=idv,group=groupv,value=sizev)
dim(nodes_nkk)
cort2=round(cor(t(meta_select[,-1]),method="spearman"), 2)
#dim(cort2)
#####################################################
calu_cor_test=function(){
  cor_matrix=t(meta_select[,-1])
  spcnum=dim(cor_matrix)[2]
  cor_matrix_calcute_co=data.frame(matrix(0,nrow=spcnum,ncol=spcnum))
  #cor_matrix_calcute_pv=data.frame(matrix(0,nrow=spcnum,ncol=spcnum))
  row.names(cor_matrix_calcute_co)=colnames(cor_matrix)
  colnames(cor_matrix_calcute_co)=colnames(cor_matrix)
  #row.names(cor_matrix_calcute_pv)=colnames(cor_matrix)
  
  for (i in 1:spcnum){
    for(j in 1:spcnum){
      ab=cor.test(cor_matrix[,i],cor_matrix[,j],method="spearm")
      #print(ab$p.value)
      #print(ab$estimate)
      if(ab$p.value < 0.05){
        cor_matrix_calcute_co[i,j]=ab$estimate
      }else{
        cor_matrix_calcute_co[i,j]=0
      }
      #
      #cor_matrix_calcute_pv[i,j]=ab$p.value
    }
  }
  
  cort2=cor_matrix_calcute_co
  return(cor_matrix_calcute_co)
}

#cor_matrix_calcute_pv
##########################################


#head(cort2[1:35,1:15])

sourcev=c()
targetv=c()
colorv=c()
sizev=c()
postivev=0
negativev=0
################################################

################################################
for (r in 1:nrow(cort2) ){
  #get number
  nlength=0
  for(c in 1:ncol(cort2)){
    if(r != c & cort2[r,c] >0.6 ){
      nlength=nlength+1

    }
  }
  for(c in 1:ncol(cort2)){
    if(r != c & cort2[r,c] < -0.6 ){
      nlength=nlength+1

    }
  }
  
  sizev=append(sizev,nlength)
  for(c in 1:ncol(cort2)){
    if(r != c & cort2[r,c] >0.6 & r >= c){
      
      sourcev=append(sourcev,r)
      targetv=append(targetv,c)
      colorv=append(colorv,'red')
      postivev=postivev+1
    }
    if(r != c & cort2[r,c] < -0.6 & r >= c){
      
      sourcev=append(sourcev,r)
      targetv=append(targetv,c)
      colorv=append(colorv,'green')
      negativev=negativev+1
    }
  }
}

#sizev
valuev=rowSums(meta_select[,-1])
#nodes_nkk=data.frame(id=idv,label=idv,group=groupv,value=valuev*100)
#edges_nkk=data.frame(from=sourcev, to=targetv,color=colorv)
#edges_nkk[edges_nkk$color=='green',]
nodes_out=data.frame(id=idv,label=idv,group=groupv,value=valuev*100)
nodes_out$color='#C0C0C0'
#levels(factor(nodes_out$group))
#######################################################
#fungi
######################################################
fungi_color=function(nodes_out){
  nodes_out$color[nodes_out$group=='Ascomycota']='#FF0000'
  nodes_out$color[nodes_out$group=='Basidiomycota']='#FFC000'
  nodes_out$color[nodes_out$group=='Glomeromycota']='#92D050'
  nodes_out$color[nodes_out$group=='Mortierellomycota']='#00B050'
  nodes_out$color[nodes_out$group=='Chytridiomycota']='#00B0F0'
  nodes_out$color[nodes_out$group=='Kickxellomycota']='#0070C0'
  nodes_out$color[nodes_out$group=='Zoopagomycota']='#DE89FF'
  nodes_out$color[nodes_out$group=='Aphelidiomycota']='#FF5583'
  nodes_out$color[nodes_out$group=='Rozellomycota']='#D4B2B1'
  nodes_out$color[nodes_out$group=='unclassified']='#C0C0C0'
  return(nodes_out)
}
nodes_out=fungi_color(nodes_out)
bacteria_color=function(nodes_out){
  nodes_out$color[nodes_out$group=='Proteobacteria']='#FF0000'
  nodes_out$color[nodes_out$group=='Firmicutes']='#FFC000'
  nodes_out$color[nodes_out$group=='Actinobacteriota']='#92D050'
  nodes_out$color[nodes_out$group=='Chloroflexi']='#00B050'
  nodes_out$color[nodes_out$group=='Bacteroidota']='#00B0F0'
  nodes_out$color[nodes_out$group=='Myxococcota']='#0070C0'
  nodes_out$color[nodes_out$group=='Bdellovibrionota']='#DE89FF'
  nodes_out$color[nodes_out$group=='Acidobacteriota']='#FF5583'
  nodes_out$color[nodes_out$group=='Verrucomicrobiota']='#D4B2B1'
  nodes_out$color[nodes_out$group=='unclassified']='#C0C0C0'
  nodes_out$color[nodes_out$group=='Gemmatimonadota']='#DE70A2'
  nodes_out$color[nodes_out$group=='Planctomycetota']='#36976A'
  nodes_out$color[nodes_out$group=='Patescibacteria']='#0687C3'
  nodes_out$color[nodes_out$group=='SAR324_cladeMarine_group_B']='#B13080'
  nodes_out$color[nodes_out$group=='Desulfobacterota']='#6F488F'
  nodes_out$color[nodes_out$group=='Elusimicrobiota']='#73B35D'
  nodes_out$color[nodes_out$group=='Abditibacteriota']='#918D20'
  nodes_out$color[nodes_out$group=='Entotheonellaeota']='#36976A'
  #nodes_out$color[nodes_out$group=='Armatimonadota']=''
  #nodes_out$color[nodes_out$group=='Deinococcota']=''
  return(nodes_out)
}
nodes_out=bacteria_color(nodes_out)

edges_out=data.frame(Source=sourcev,
                     Target=targetv,color=colorv)

nodes_out2=data.frame(name=character(),label=character(),group=character(),value=character(),color=character())

for (i in unique(c(edges_out$Source,edges_out$Target))){
  #print(i)
  nodes_out2=rbind(nodes_out2,nodes_out[nodes_out$id==i,])
} 
nodes_out=nodes_out2
write.csv(nodes_out,paste(dir1,'out.nodes_',nametitle,'_nkk.csv',sep=''),row.names = F)
write.csv(edges_out,paste(dir1,'out.edges_',nametitle,'_nkk.csv',sep=''),row.names = F)


#c(edges_out$Source,edges_out$Target)

#
avg=dim(edges_out)[1]/dim(nodes_out)[1]*2
head(avgtable)
rv=data.frame(name=nametitle,node=dim(nodes_out)[1],postive_edge=postivev,negative_edge=negativev,avg_degree=avg,edge=dim(edges_out)[1])
avgtable=rbind(avgtable,rv)
avgtable


#write.csv(avgtable,paste(dir1,'out.degree.csv',sep=''),row.names = F)
##########################################################
library(igraph)
net <- graph_from_data_frame(d=edges_out, vertices=nodes_out, directed=F) 
ave_degree=mean(degree(net, mode="in"))
ave_degree_vector=degree(net, mode="in")

#Inverse of the node’s average geodesic distance to others in the network.
closeness_centr=closeness(net,mode='all', weights=NA,normalized=T) 
#closeness_centr=closeness(net,mode='in', weights=NA,normalized=T) 

#closeness_centr=closeness(net,mode='out', weights=NA,normalized=T) 
#closeness_centr=centr_eigen(net, directed=F, normalized=T)$vector
#max(closeness_centr)


#length(ave_degree_vector)
#length(closeness_centr)
#nodes_out
closeness_degree=data.frame(closeness=closeness_centr,degree=ave_degree_vector,color=nodes_out$color)
head(closeness_degree)
#closeness_degree2=closeness_degree[closeness_degree$degree>60,]
#plot(closeness_degree)
hub_nodes_num=dim(closeness_degree[closeness_degree$closeness>0.04 & closeness_degree$degree>60,])[1]
print(hub_nodes_num)

#
#closeness_degree=closeness_degree5

#nametitle='Bulk Soil'
#nametitle='Rhizophsere Soil'
#nametitle='Root Endosphere'
#nametitle='Phylloplane'
#nametitle='Leaf Endosphere'

ggplot(closeness_degree,aes(x=closeness,y=degree) )+
geom_point(color=closeness_degree$color)+scale_x_continuous(limits = c(0, 0.25))+
  scale_y_continuous(limits = c(0, 150))+
  ggtitle(paste(  "degree and closeness", nametitle, sep=' ' )) + xlab("Closeness centrality")+
  theme(plot.title = element_text( size = 24),
        axis.text.x = element_text( size = 20),
        axis.text.y = element_text( size = 20),  
        axis.title.x = element_text( size = 20),
        axis.title.y = element_text( size = 20))+
        geom_hline(yintercept=60, linetype="dashed")+
        geom_vline(xintercept = 0.04,linetype="dashed")
pdf(paste(dir1,'out.closeness_',nametitle,'.pdf',sep=''))
ggplot(closeness_degree,aes(x=closeness,y=degree) )+
  geom_point(color=closeness_degree$color)+
  scale_x_continuous(limits = c(0, 0.25))+
  scale_y_continuous(limits = c(0, 150))+
  ggtitle(paste(  "degree and closeness", nametitle, sep=' ' )) + xlab("Closeness centrality")+
  theme(plot.title = element_text( size = 24),
        axis.text.x = element_text( size = 20),
        axis.text.y = element_text( size = 20),  
        axis.title.x = element_text( size = 20),
        axis.title.y = element_text( size = 20))+
  geom_hline(yintercept=60, linetype="dashed")+
  geom_vline(xintercept = 0.04,linetype="dashed")
dev.off()
write.csv(closeness_degree,paste(dir1,'out.closeness_degree_',nametitle,'_nkk.csv',sep=''),row.names = F)

#####################################################
#calculate average distance_table
ave_path_distance=mean_distance(net, directed=F)
#
print(ave_path_distance)
avg_cluster_coefficient=transitivity(net)
print(avg_cluster_coefficient)
ceb <- cluster_edge_betweenness(net,modularity = TRUE) 
modularity(ceb) # how modular the graph partitioning is
####################################################



pdf(paste(dir1,'out.edges_',nametitle,'.pdf',sep=''))
ggplot(closeness_degree,aes(x=closeness,y=degree) )+
  geom_point(color=closeness_degree$color)+
  scale_x_continuous(limits = c(0, 0.25))+
  scale_y_continuous(limits = c(0, 150))+
  ggtitle(paste(  "degree and closeness", nametitle, sep=' ' )) + xlab("Closeness centrality")+
  theme(plot.title = element_text( size = 24),
        axis.text.x = element_text( size = 20),
        axis.text.y = element_text( size = 20),  
        axis.title.x = element_text( size = 20),
        axis.title.y = element_text( size = 20))+
  geom_hline(yintercept=60, linetype="dashed")+
  geom_vline(xintercept = 0.04,linetype="dashed")
dev.off()


####

closeness_degree1=closeness_degree
closeness_degree2=closeness_degree
closeness_degree3=closeness_degree
closeness_degree4=closeness_degree
closeness_degree5=closeness_degree

head(closeness_degree1)
head(closeness_degree2)
head(closeness_degree3)
head(closeness_degree4)
head(closeness_degree5)
write.csv(closeness_degree1,paste(dir1,'out.closeness_degree1.csv',sep=''),row.names = F)
write.csv(closeness_degree2,paste(dir1,'out.closeness_degree2.csv',sep=''),row.names = F)
write.csv(closeness_degree3,paste(dir1,'out.closeness_degree3.csv',sep=''),row.names = F)
write.csv(closeness_degree4,paste(dir1,'out.closeness_degree4.csv',sep=''),row.names = F)
write.csv(closeness_degree5,paste(dir1,'out.closeness_degree5.csv',sep=''),row.names = F)
#dir1='/Users/xugang/Desktop/sequencing_center_desktop/yanglab/metagenome/bactaria/'
#dir1='/Users/xugang/Desktop/sequencing_center_desktop/yanglab/metagenome/'
closeness_degree1=read.csv(paste(dir1,'out.closeness_degree1.csv',sep=''))
closeness_degree2=read.csv(paste(dir1,'out.closeness_degree2.csv',sep=''))
closeness_degree3=read.csv(paste(dir1,'out.closeness_degree3.csv',sep=''))
closeness_degree4=read.csv(paste(dir1,'out.closeness_degree4.csv',sep=''))
closeness_degree5=read.csv(paste(dir1,'out.closeness_degree5.csv',sep=''))

###############################################################


###############################################################
max(closeness_degree1$closeness)
max(closeness_degree2$closeness)
max(closeness_degree3$closeness)
max(closeness_degree4$closeness)
max(closeness_degree5$closeness)
#####################################################
#calculate average distance_table
ave_path_distance=mean_distance(net, directed=F)
#
transitivity(net)
ceb <- cluster_edge_betweenness(net,modularity = TRUE) 
modularity(ceb) # how modular the graph partitioning is
####################################################

dendPlot(ceb, mode="hclust")
#plot(ceb, net) 

#Let’s examine the community detection igraph object:
class(ceb)
length(ceb)     # number of communities
membership(ceb) # community membership for each node
modularity(ceb) # how modular the graph partitioning is
crossing(ceb, net)   # boolean vector: TRUE for edges across communities
(dir1)
#####################################################

plot(closeness_degree1[,1:2])
plotdata=closeness_degree5
nametitle2='Bulk Soil'
nametitle2='Rhizophsere Soil'
nametitle2='Root Endosphere'
nametitle2='Phylloplane'
nametitle2='Leaf Endosphere'
pdf(paste(dir1,'out.edges_',nametitle2,'.pdf',sep=''))
ggplot(plotdata,aes(x=closeness,y=degree) )+
  geom_point(color=plotdata$color)+scale_x_continuous(limits = c(0, 0.25))+
  scale_y_continuous(limits = c(0, 200))+
  ggtitle(paste(  "degree and closeness", nametitle2, sep=' ' )) + xlab("Closeness centrality")+
  theme(plot.title = element_text( size = 24),
        axis.text.x = element_text( size = 20),
        axis.text.y = element_text( size = 20),  
        axis.title.x = element_text( size = 20),
        axis.title.y = element_text( size = 20))
dev.off()
#####################################################
cor_matrix=t(meta_select[,-1])
spcnum=dim(cor_matrix)[2]
cor_matrix_calcute_co=data.frame(matrix(0,nrow=spcnum,ncol=spcnum))
cor_matrix_calcute_pv=data.frame(matrix(0,nrow=spcnum,ncol=spcnum))
row.names(cor_matrix_calcute_co)=colnames(cor_matrix)
row.names(cor_matrix_calcute_pv)=colnames(cor_matrix)
head(cor_matrix_calcute)
for (i in 1:spcnum){
  for(j in 1:spcnum){
    ab=cor.test(cor_matrix[,i],cor_matrix[,j],method="spearm")
    #print(ab$p.value)
    #print(ab$estimate)
    if(ab$p.value < 0.05){
      cor_matrix_calcute_co[i,j]=ab$estimate
    }else{
      cor_matrix_calcute_co[i,j]=0
    }
    #
    #cor_matrix_calcute_pv[i,j]=ab$p.value
  }
}
cor_matrix_calcute_co
#cor_matrix_calcute_pv
##########################################
### Example 1: Violent crime rates by US state
hc <- hclust(dist(USArrests), "ave")
plot(hc)
plot(hc, hang = -1)
###########################################
library(igraph)
## An in-star
nodes <- read.csv("/Users/xugang/Downloads/netscix2016/Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("/Users/xugang/Downloads/netscix2016/Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)
head(nodes)
head(links)

#Examine the data:
head(nodes)
head(links)
nrow(nodes); length(unique(nodes$id))
nrow(links); nrow(unique(links[,c("from", "to")]))

# collapse same links
links <- aggregate(links[,3], links[,-3], sum)
links <- links[order(links$from, links$to),]
colnames(links)[4] <- "weight"
rownames(links) <- NULL

#3.2 DATASET 2: matrix
#Two-mode or bipartite graphs have two different types of actors and links that go across, but not within each type. Our second media example is a network of that kind, examining links between news sources and their consumers.
nodes2 <- read.csv("/Users/xugang/Downloads/netscix2016/Dataset2-Media-User-Example-NODES.csv", header=T, as.is=T)
links2 <- read.csv("/Users/xugang/Downloads/netscix2016/Dataset2-Media-User-Example-EDGES.csv", header=T, row.names=1)

#Examine the data:
head(nodes2)
head(links2)

#We can see that links2 is an adjacency matrix for a two-mode network:
links2 <- as.matrix(links2)
dim(links2)
dim(nodes2)

#4. Turning networks into igraph objects
#We start by converting the raw data to an igraph network object. Here we use igraph’s graph.data.frame function, which takes two data frames: d and vertices.

#d describes the edges of the network. Its first two columns are the IDs of the source and the target node for each edge. The following columns are edge attributes (weight, type, label, or anything else). 
#vertices starts with a column of node IDs. Any following columns are interpreted as node attributes.
#4.1 Dataset 1
library(igraph)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net <- graph_from_data_frame(d=edges_out, vertices=nodes_out, directed=F) 
#nodes_out
#edges_out

class(net)
net
###We also have easy access to nodes, edges, and their attributes with:
E(net)       # The edges of the "net" object
V(net)       # The vertices of the "net" object
E(net)$type  # Edge attribute "type"
V(net)$media # Vertex attribute "media"

#Now that we have our igraph network object, let’s make a first attempt to plot it.
plot(net, edge.arrow.size=.4,vertex.label=NA)

#That doesn’t look very good. Let’s start fixing things by removing the loops in the graph.
net <- simplify(net, remove.multiple = F, remove.loops = T) 

#Or data frames describing nodes and edges:
as_data_frame(net, what="edges")
as_data_frame(net, what="vertices")

#4.2 Dataset 2
#As we have seen above, this time the edges of the network are in a matrix format. 
#We can read those into a graph object using graph_from_incidence_matrix(). 
#In igraph, bipartite networks have a node attribute called type that is FALSE (or 0) for vertices in one mode and TRUE (or 1) for those in the other mode.
head(nodes2)
head(links2)

net2 <- graph_from_incidence_matrix(links2)
table(V(net2)$type)

# 5. Plotting networks with igraph
plot(net, edge.arrow.size=.4, edge.curved=.1)


plot(net, edge.arrow.size=.2, edge.curved=0,
     vertex.color="orange", vertex.frame.color="#555555",
     vertex.label=V(net)$media, vertex.label.color="black",
     vertex.label.cex=.7) 
#The second way to set attributes is to add them to the igraph object.
#Let’s say we want to color our network nodes based on type of media, and size them based on audience size (larger audience -> larger node). 
#We will also change the width of the edges based on their weight.
# Generate colors based on media type:
colrs <- c("gray50", "tomato", "gold")
V(net)$color <- colrs[V(net)$media.type]

# Set node size based on audience size:
V(net)$size <- V(net)$audience.size*0.7

# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(net)$label.color <- "black"
V(net)$label <- NA

# Set edge width based on weight:
E(net)$width <- E(net)$weight/6

#change arrow size and edge color:

E(net)$arrow.size <- .2
E(net)$edge.color <- "gray80"
E(net)$width <- 1+E(net)$weight/12
plot(net, edge.color="orange", vertex.color="gray50") 

#It helps to add a legend explaining the meaning of the colors we used:
plot(net) 
legend(x=-1.5, y=-1.1, c("Newspaper","Television", "Online News"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)


# 6.5 Node degrees
deg <- degree(net, mode="all")
plot(net, vertex.size=deg*3)
hist(deg, breaks=1:vcount(net)-1, main="Histogram of node degree")

#6.7 Centrality & centralization
#Degree (number of ties)
mean(degree(net, mode="in"))
degree(net, mode="in")
centr_degree(net, mode="in", normalized=T)

#Closeness (centrality based on distance to others in the graph)
#Inverse of the node’s average geodesic distance to others in the network.
closeness(net, mode="all", weights=NA,normalized=T)
max(closeness(net, mode="all", weights=NA,normalized=T))
centr_clo(net, mode="all", normalized=T) 
max(centr_clo(net, mode="all", normalized=T)$res )


#Eigenvector (centrality proportional to the sum of connection centralities)
#Values of the first eigenvector of the graph matrix.
eigen_centrality(net, directed=F, weights=NA)
max(eigen_centrality(net, directed=F, weights=NA)$res)
centr_eigen(net, directed=F, normalized=T)
max(centr_eigen(net, directed=F, normalized=T)$vector)
#Betweenness (centrality based on a broker position connecting others)
#Number of geodesics that pass through the node or the edge.
betweenness(net, directed=F, weights=NA)
max(betweenness(net, directed=F, weights=NA))
edge_betweenness(net, directed=F, weights=NA)
max(edge_betweenness(net, directed=F, weights=NA))
centr_betw(net, directed=F, normalized=T)
max(centr_betw(net, directed=F, normalized=T))
# 6.8 Hubs and authorities
#The hubs and authorities algorithm developed by Jon Kleinberg was initially used to examine web pages. Hubs were expected to contain catalogs with a large number of outgoing links; while authorities would get many incoming links from hubs, presumably because of their high-quality relevant information.
hs <- hub_score(net, weights=NA)$vector
as <- authority_score(net, weights=NA)$vector

par(mfrow=c(1,2))
plot(net, vertex.size=hs*50, main="Hubs")
plot(net, vertex.size=as*30, main="Authorities")
dev.off()

#7. Distances and paths
#Average path length: the mean of the shortest distance between each pair of nodes in the network (in both directions for directed graphs).
mean_distance(net, directed=F)
## [1] 2.058824
mean_distance(net, directed=T)
## [1] 2.742188

#We can also find the length of all shortest paths in the graph:
distances(net) # with edge weights
distances(net, weights=NA) # ignore weights

#We can extract the distances to a node or set of nodes we are interested in. Here we will get the distance of every media from the New York Times.
dist.from.NYT <- distances(net, v=V(net)[media=="NY Times"], to=V(net), weights=NA)
# Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.NYT)+1)
col <- col[dist.from.NYT+1]
plot(net, vertex.color=col, vertex.label=dist.from.NYT, edge.arrow.size=.6, 
     vertex.label.color="white")

#We can also find the shortest path between specific nodes. Say here between MSNBC and the New York Post:
news.path <- shortest_paths(net, 
                              from = V(net)[media=="MSNBC"], 
                              to  = V(net)[media=="New York Post"],
                              output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(net))
ecol[unlist(news.path$epath)] <- "orange"

# Generate edge width variable to plot the path:
ew <- rep(2, ecount(net))
ew[unlist(news.path$epath)] <- 4

# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(net))
vcol[unlist(news.path$vpath)] <- "gold"

plot(net, vertex.color=vcol, edge.color=ecol, 
     edge.width=ew, edge.arrow.mode=0)

#8. Subgroups and communities
net.sym <- as.undirected(net, mode= "collapse",
                         edge.attr.comb=list(weight="sum", "ignore"))
## 8.1 Cliques
## Find cliques (complete subgraphs of an undirected graph)
cliques(net.sym) # list of cliques       
sapply(cliques(net.sym), length) # clique sizes
largest_cliques(net.sym) # cliques with max number of nodes
vcol <- rep("grey80", vcount(net.sym))
vcol[unlist(largest_cliques(net.sym))] <- "gold"
plot(as.undirected(net.sym), vertex.label=V(net.sym)$name, vertex.color=vcol)

#8.2 Community detection
##A number of algorithms aim to detect groups that consist of densely connected nodes with fewer connections across groups.

#Community detection based on edge betweenness (Newman-Girvan)
#High-betweenness edges are removed sequentially (recalculating at each step) and the best partitioning of the network is selected.
ceb <- cluster_edge_betweenness(net,modularity = TRUE) 
dendPlot(ceb, mode="hclust")
#plot(ceb, net) 

#Let’s examine the community detection igraph object:
class(ceb)
length(ceb)     # number of communities
membership(ceb) # community membership for each node
modularity(ceb) # how modular the graph partitioning is
crossing(ceb, net)   # boolean vector: TRUE for edges across communities
#High modularity for a partitioning reflects dense connections within communities and sparse connections across communities.
#Community detection based on based on propagating labels
#Assigns node labels, randomizes, than replaces each vertex’s label with the label that appears most frequently among neighbors. 
#Those steps are repeated until each vertex has the most common label of its neighbors.

clp <- cluster_label_prop(net)
plot(clp, net)

#Community detection based on greedy optimization of modularity
cfg <- cluster_fast_greedy(as.undirected(net))
plot(cfg, as.undirected(net))

#We can also plot the communities without relying on their built-in plot:
V(net)$community <- cfg$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
plot(net, vertex.color=colrs[V(net)$community])

#8.3 K-core decomposition
#The k-core is the maximal subgraph in which every node has degree of at least k. The result here gives the coreness of each vertex in the network. A node has coreness D if it belongs to a D-core but not to (D+1)-core.
kc <- coreness(net, mode="all")
plot(net, vertex.size=kc*6, vertex.label=kc, vertex.color=colrs[kc])

#9. Assortativity and Homophily
#Homophily: the tendency of nodes to connect to others who are similar on some variable.
#assortativity_nominal() is for categorical variables (labels)
#assortativity() is for ordinal and above variables
#assortativity_degree() checks assortativity in node degrees
assortativity_nominal(net, V(net)$media.type, directed=F)
## [1] 0.1715568
assortativity(net, V(net)$audience.size, directed=F)
## [1] -0.1102857
assortativity_degree(net, directed=F)
## [1] -0.009551146





```





