exp=read.csv("exp.csv",head=TRUE,row.name=1)
library(reshape2)
melt(meta)
library(limma)
library(tidyr)
library(bioMart)
ensembl<-  useMart("ensembl", host="grch37.ensembl.org",dataset="hsapiens_gene_ensembl")
ltime=getBM(attributes=c("hgnc_symbol"),
filters = "entrezgene", values = names(trend_time), mart=ensembl)

MEAN=apply(exp[,1:12],1,mean)
data=foreach(i=l,.combine=rbind) %dopar% {
x=which(exp$Ensembl==i)
y=MEAN[x]
y=names(y)[which(y==max(y))]
exp[y,]
}

rownames(data)=data$Ensembl
symbol=data$GeneSymbol
data=data[,1:12]
rm(exp)
dim(data)
design=model.matrix(~0+meta$point+meta$sample)
library(limma)
fit=lmFit(data[,as.character(meta[,1])],design)
fit <- eBayes(fit)

X20=topTable(fit,coef="pointX2",n=22086)
X20_sig=X20[which(X20[,5]<0.05&abs(X20[,1])>1),]
X20$sig=ifelse(X20[,5]<0.05&abs(X20[,1])>1,1,0)
X20$dir=ifelse(X20[,1]>0,1,-1)
X20$sig=apply(X20[,7:8],1,prod)

draw=X20[,c(1,4,7)]

draw[,3]=as.character(draw[,3])
library(ggplot2)
colnames(draw)=c("logFC","P.Value","trend")
p=ggplot(draw,aes(logFC,-log10(P.Value),color=trend))
p=p+geom_point()
p=p+scale_color_manual(breaks=c("-1","0","1"),values=c("#FB61D7","grey","#00C094"),labels=c("down","non","up"))
p=p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.y = element_text(size=14),plot.title = element_text(size=10),axis.line = element_line(colour = "black"))
p=p+geom_hline(yintercept=3.5)
p=p+geom_vline(xintercept=1)
p=p+geom_vline(xintercept=-1)
p

tiff("volcano_X21.tiff",units="in",height=6,width=4,res = 300)
p
dev.off()
gene_hm=unique(c(gene_hm,rownames(X20_sig)))

colnames(draw)=c("A1","A2","A0","B1","B2","B0","C0","D0","C1","C2","D1","D2")
X20DifGene_allSample=pheatmap(draw,scale="row")

library("org.Hs.eg.db")
library(clusterProfiler)
names(trend_2VS0)[which(trend_2VS0==1)]


egotime_down <- enrichGO(gene         =names(trend_time)[which(trend_time==1)],
                OrgDb         = org.Hs.eg.db,
                
                ont           = "BP",
                pAdjustMethod = "BH",
              
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
egotime_up=data.frame(egotime_down)[1:30,c(2,6,9)]

draw=egotime_down[c(1,2,15,21,30),]
library(stringr)
draw$Description = str_wrap(draw$Description, width = 25)
colnames(draw)=c("term","p.adj")
draw$term=factor(draw$term,levels=draw$term[6:1])
p=ggplot(draw,aes(term,-log10(p.adj)))
p=p+geom_bar(stat = "identity",fill="#00C094")
p=p+coord_flip()
p=p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.title= element_blank(),
panel.background = element_blank(),axis.text.y = element_text(size=12),plot.title = element_text(size=10),axis.line = element_line(colour = "black"))
p=p+geom_hline(yintercept=1.3)
p
tiff("egoX20_up.tif", units="in",height=2.6,width=4.2,res = 300)
p
dev.off()



