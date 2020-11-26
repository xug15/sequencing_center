args <- commandArgs(trailingOnly = TRUE)
filename=args[1]
orfstart=as.integer(args[2])
orfend=as.integer(args[3])
name=args[4]
data=read.table(filename,sep="\t",header=T)
pdf(paste(filename,'.pdf',sep=''),width=6,height=4,paper='special')
mp <- barplot(data$number+1,col= rainbow(3),border = NA,main=name,log="y")
axis(1,at=mp,labels=data$position)
dev.off()

pdf(paste(filename,'orf.pdf',sep=''),width=6,height=4,paper='special')
mp <- barplot(data$number[orfstart:orfend]+1,col= rainbow(3),border = NA,main=name,log="y")
axis(1,at=mp,labels=data$position[orfstart:orfend])
dev.off()


