#load data and save, get mean expression for each gene level and protein level for chimp and human for both data sets
r.143<-read.table(file="TableS3.txt",header=TRUE,sep="\t")
r.159<-read.table(file="TableS4.txt",header=TRUE,sep="\t")
names.143<-rownames(r.143)
names.159<-rownames(r.159)

par(mfrow=c(2,2))
chimp.143.rna<-apply(r.143[,3:8],1,mean,na.rm=T);chimp.143.rna
human.143.rna<-apply(r.143[,9:14],1,mean,na.rm=T);human.143.rna
plot(chimp.143.rna,human.143.rna,xlab="human mRNA expression",ylab=" human mRNA expression",main="Mean mRNA expression for chimp and human in 143 genes",cex.lab=.75,cex.main=.75,pch=20,col="blueviolet")

chimp.143.pr<-apply(r.143[,15:17],1,mean,na.rm=T);chimp.143.pr
human.143.pr<-apply(r.143[,18:20],1,mean,na.rm=T);human.143.pr
plot(chimp.143.pr,human.143.pr,xlab="chimp protein expression",ylab="human protein expression",main="Mean protein expression for chimp and human in 143 genes",cex.lab=.75,cex.main=.75,pch=20,col="blueviolet")

chimp.159.rna<-apply(r.159[,3:8],1,mean,na.rm=T);chimp.159.rna
human.159.rna<-apply(r.159[,9:14],1,mean,na.rm=T);human.159.rna
plot(chimp.159.rna,human.159.rna,xlab="chimp mRNA expression",ylab="human mRNA expression",main="Mean mRNA expression for chimp and human in 159 genes",cex.lab=.75,cex.main=.75,pch=20,col="darkgoldenrod2")

chimp.159.pr<-apply(r.159[,15:17],1,mean,na.rm=T);chimp.159.pr
human.159.pr<-apply(r.159[,18:20],1,mean,na.rm=T);human.159.pr
plot(chimp.159.pr,human.159.pr,xlab="chimp protein expression",ylab="human protein expression",main="Mean protein expression for chimp and human in 159 genes",cex.lab=.75,cex.main=.75,pch=20,col="darkgoldenrod2")

#Question2
#Use cor to get pearson correlation coefficient between mean mrna/prot expression per gene, per organism,per sample
cor(chimp.143.rna,human.143.rna)
cor(chimp.143.pr,human.143.pr)
cor(chimp.159.rna,human.159.rna)
cor(chimp.159.pr,human.159.pr)
#both pearson correlation values for the relationship between chimp mRNA expression and human mRNA expression show a strong positive relationship exists between the two with values that were both positive and nearly 1 for both 143 and 159 samples. On the other hand, the pearson correlation demonstrates almost no relationship between chimp and human protein expression
#with both tests resulting in values that were close to zero.  If any relationship exists, it is both a weak and negative relationship.

#Question3
#fit 4 linear models to the previous means. make dfs first for analysis. plot graphs with regression line
df.143<-data.frame(chimp.rna=chimp.143.rna,human.rna=human.143.rna,chimp.pr=chimp.143.pr,human.pr=human.143.pr)
df.159<-data.frame(chimp.rna=chimp.159.rna,human.rna=human.159.rna,chimp.pr=chimp.159.pr,human.pr=human.159.pr)

rna.model.143<-lm(df.143$human.rna~df.143$chimp.rna,df.143);rna.model.143
hist(resid(rna.model.143))
shapiro.test(resid(rna.model.143))$p.value
plot(rna.model.143)
anova(rna.model.143)
plot(df.143$chimp.rna,df.143$human.rna,xlab="chimp mRNA expression",ylab="human mRNA expression",main="Mean mRNA expression for chimp and human in 143 genes",cex.lab=.75,cex.main=.75,col="blueviolet")
points(df.143$chimp.rna[!is.na(df.143$human.rna)],predict(rna.model.143),col="darkturquoise",pch=20)

pr.model.143<-lm(df.143$human.pr~df.143$chimp.pr,df.143);pr.model.143
hist(resid(pr.model.143))
shapiro.test(resid(pr.model.143))$p.value
plot(pr.model.143)
anova(pr.model.143)
plot(df.143$chimp.pr,df.143$human.pr,xlab="chimp protein expression",ylab="human protein expression",main="Mean protein expression for chimp and human in 143 genes",cex.lab=.75,cex.main=.75,col="blueviolet")
points(df.143$chimp.pr[!is.na(df.143$human.pr)],predict(pr.model.143),col="darkturquoise",pch=20)

rna.model.159<-lm(df.159$human.rna~df.159$chimp.rna,df.159);rna.model.159
hist(resid(rna.model.159))
shapiro.test(resid(rna.model.159))$p.value
plot(rna.model.159)
anova(rna.model.159)
plot(df.159$chimp.rna,df.159$human.rna,xlab="chimp mRNA expression",ylab="human mRNA expression",main="Mean mRNA expression for chimp and human in 159 genes",cex.lab=.75,cex.main=.75,col="darkgoldenrod2",pch=20)
points(df.159$chimp.rna[!is.na(df.159$human.rna)],predict(rna.model.159),col="purple",pch=20)

pr.model.159<-lm(df.159$human.pr~df.159$chimp.pr,df.159);pr.model.159
hist(resid(pr.model.159))
shapiro.test(resid(pr.model.159))$p.value
plot(pr.model.159)
anova(pr.model.159)
plot(df.159$chimp.pr,df.159$human.pr,xlab="chimp protein expression",ylab="human protein expression",main="Mean protein expression for chimp and human in 159 genes",cex.lab=.75,cex.main=.75,col="darkgoldenrod2",pch=20)
points(df.159$chimp.pr[!is.na(df.159$human.pr)],predict(pr.model.159),col="purple",pch=20)

#Question 4 Plot mean protein vs mean mRNA  for chimp and for human
par(mfrow=c(2,2))
plot(df.143$chimp.rna,df.143$chimp.pr,xlab="chimp mRNA expression",ylab="chimp protein expression",main="Relationship between mean protein and mRNA expression for chimp across 143 genes",cex.lab=.75,cex.main=.75,col="dodgerblue3",pch=20)
plot(df.143$human.rna,df.143$human.pr,xlab="human mRNA expression",ylab="human protein expression",main="Relationship between mean protein and mRNA expression for human across 143 genes",cex.lab=.75,cex.main=.75,col="dodgerblue3",pch=20)
plot(df.159$chimp.rna,df.159$chimp.pr,xlab="chimp mRNA expression",ylab="chimp protein expression",main="Relationship between mean protein and mRNA expression for chimp across 159 genes",cex.lab=.75,cex.main=.75,col="darkorange",pch=20)
plot(df.159$human.rna,df.159$human.pr,xlab="human mRNA expression",ylab="human protein expression",main="Relationship between mean protein and mRNA expression for human across 159 genes",cex.lab=.75,cex.main=.75,col="darkorange",pch=20)


#Question 5 
#find differences in mean rna values between species and mean protein values, plot 
par(mfrow=c(2,2))
#My result 143
rna.diff.143<-human.143.rna-chimp.143.rna
pr.diff.143<-human.143.pr-chimp.143.pr
plot(rna.diff.143,pr.diff.143,xlab="mRNA expression",ylab="protein expression",main="Relationship of 143 mRNA and protein expression differences",col="turquoise4",cex.lab=.75,cex.main=.55,pch=1)

rna.diff.159<-human.159.rna-chimp.159.rna
pr.diff.159<-human.159.pr-chimp.159.pr
plot(rna.diff.159,pr.diff.159,xlab="mRNA expression",ylab="protein expression",main="Relationship of 159 mRNA and protein expression differences",col="turquoise4",cex.lab=.75,cex.main=.55,pch=1)

#Their result 143
Fu.diff.143<-human.143.rna-chimp.143.rna
Fu.pr.diff.143<-chimp.143.pr-human.143.pr
plot(Fu.diff.143,Fu.pr.diff.143,xlab="mRNA expression",ylab="protein expression",main="Relationship for 143 found in Fu et al.,2007",col="turquoise4",cex.lab=.75,cex.main=.55,pch=1)

#Question 6
#find correlation coeff for relationships between differences in rna/protein, for both data sets.
cor(rna.diff.143,pr.diff.143)
cor(rna.diff.159,pr.diff.159)

#Question 7 - Linear model for Question 5
par(mfrow=c(2,2))
lm.143<-lm(pr.diff.143~rna.diff.143)
hist(resid(lm.143))
shapiro.test(resid(lm.143))$p.value
plot(lm.143)
anova(lm.143)
plot(rna.diff.143,pr.diff.143,xlab="mRNA expression",ylab="protein expression",main="Relationship of 143 mRNA and protein expression differences",col="turquoise4",cex.lab=.75,cex.main=.55,pch=1)
points(rna.diff.143[!is.na(pr.diff.143)],predict(lm.143),col="hotpink",pch=20)

lm.159<-lm(pr.diff.159~rna.diff.159)
hist(resid(lm.159))
shapiro.test(resid(lm.159))$p.value
plot(lm.159)
anova(lm.159)
plot(rna.diff.159,pr.diff.159,xlab="mRNA expression",ylab="protein expression",main="Relationship of 159 mRNA and protein expression differences",col="turquoise4",cex.lab=.75,cex.main=.55,pch=1)
points(rna.diff.159[!is.na(pr.diff.159)],predict(lm.159),col="hotpink",pch=20)

#Question 8 find p value via ttest for each gene prot express and mRNA express for equal variances
m.143.eq<-matrix(0,nrow=nrow(r.143),ncol=4)
colnames(m.143.eq)<-c("rna.pval","pr.pval","r.anova","p.anova")
rownames(m.143.eq)<-rownames(r.143)
m.159.eq<-matrix(0,nrow=nrow(r.159),ncol=4)
colnames(m.159.eq)<-c("rna.pval","pr.pval","r.anova","p.anova")
rownames(m.159.eq)<-rownames(r.159)
count.143.rna<-0
count.143.pr<-0
count.159.rna<-0
count.159.pr<-0
for (i in 1:nrow(r.143)){
  p.rna.143<-t.test(r.143[i,3:8],r.143[i,9:14],var.equal=T,conf.level=.99)
  m.143.eq[i,1]<-p.rna.143$p.value
  if (p.rna.143$p.value<0.01){count.143.rna<-count.143.rna+1}
  p.pr.143<-t.test(r.143[i,15:17],r.143[i,18:20],var.equal=T,conf.level=.99)
  m.143.eq[i,2]<-p.pr.143$p.value
  if (p.pr.143$p.value<.01){count.143.pr<-count.143.pr+1}
}

for (i in 1:nrow(r.159)){
  p.rna.159<-t.test(r.159[i,3:8],r.159[i,9:14],var.equal=T,conf.level=.99)
  m.159.eq[i,1]<-p.rna.159$p.value
  if (p.rna.159$p.value<0.01){count.159.rna<-count.159.rna+1}
  p.pr.159<-t.test(r.159[i,15:17],r.159[i,18:20],var.equal=T,conf.level=.99)
  m.159.eq[i,2]<-p.pr.159$p.value
  if (p.pr.159$p.value<.01){count.159.pr<-count.159.pr+1}
  }
count.143.rna
count.143.pr
count.159.rna
count.143.pr


par(mfrow=c(2,2))
hist(m.143.eq[,"rna.pval"],col="orange",xlab="mRNA p-values",main="t-test p-value distribution for mRNA")
hist(m.143.eq[,"pr.pval"],col="orange",xlab="protein p-values",main="t-test p-value distribution for protein")
hist(m.159.eq[,"rna.pval"],col="orange",xlab="mRNA p-values",main="t-test p-value distribution for mRNA")
hist(m.159.eq[,"pr.pval"],col="orange",xlab="protein p-values",main="t-test p-value distribution for protein")
plot(m.143.eq[,"rna.pval"],m.143.eq[,"pr.pval"],xlab="mRNA p-values",ylab="protein p-values",main="relationship between mRNA and protein p-values for sample 143",cex.main=.75,cex.lab=.75)
plot(m.159.eq[,"rna.pval"],m.159.eq[,"pr.pval"],xlab="mRNA p-values",ylab="protein p-values",main="relationship between mRNA and protein p-values for sample 159",cex.main=.75,cex.lab=.75)

#Question 9 unequal variances(default t.test)
m.143<-matrix(0,nrow=nrow(r.143),ncol=4)
colnames(m.143)<-c("rna.pval","pr.pval","ANOVA rna","ANOVA pr")
rownames(m.143)<-rownames(r.143)
m.159<-matrix(0,nrow=nrow(r.159),ncol=4)
colnames(m.159)<-c("rna.pval","pr.pval","ANOVA rna","ANOVA pr")
rownames(m.159)<-rownames(r.159)
for (i in 1:nrow(r.143)){
  p.rna.143<-t.test(r.143[i,3:8],r.143[i,9:14],conf.level=.95)
  m.143[i,1]<-p.rna.143$p.value
  p.pr.143<-t.test(r.143[i,15:17],r.143[i,18:20],conf.level=.95)
  if (p.rna.143$p.value<0.01){count.143.rna<-count.143.rna+1}
  m.143[i,2]<-p.pr.143$p.value
  if (p.pr.143$p.value<.01){count.143.pr<-count.143.pr+1}
}
for (i in 1:nrow(r.159)){
  p.rna.159<-t.test(r.159[i,3:8],r.159[i,9:14],conf.level=.99)
  m.159[i,1]<-p.rna.159$p.value
  if (p.rna.159$p.value<0.01){count.159.rna<-count.159.rna+1}
  p.pr.159<-t.test(r.159[i,15:17],r.159[i,18:20],conf.level=.99)
  m.159[i,2]<-p.pr.159$p.value
  if (p.pr.159$p.value<.01){count.159.pr<-count.159.pr+1}
}
par(mfrow=c(2,2))
hist(m.143[,"rna.pval"],col="blue",xlab="mRNA p-values",main="t-test p-value distribution for mRNA")
hist(m.143[,"pr.pval"],col="blue",xlab="protein p-values",main="t-test p-value distribution for protein")
hist(m.159[,"rna.pval"],col="blue",xlab="mRNA p-values",main="t-test p-value distribution for mRNA")
hist(m.159[,"pr.pval"],col="blue",xlab="protein p-values",main="t-test p-value distribution for protein")
plot(m.143[,"rna.pval"],m.143[,"pr.pval"],xlab="mRNA p-values",ylab="protein p-values",main="relationship between mRNA and protein p-values for sample 143")
plot(m.159[,"rna.pval"],m.159[,"pr.pval"],xlab="mRNA p-values",ylab="protein p-values",main="relationship between mRNA and protein p-values for sample 159")

count.143.rna
count.143.pr
count.159.rna
count.159.pr
#9-Linear models 
chimp.mRNAs <- paste("C",1:6,sep="")
human.mRNAs <- paste("H",1:6,sep="")
chimp.prots <- paste("c",1:3,sep="")
human.prots <- paste("h",1:3,sep="") 

anova.hc.mRNA <- function(x){
  x.tmp <- x[c(human.mRNAs,chimp.mRNAs)]
  s.tmp <- factor(c(rep("H",length(human.mRNAs)),rep("C",length(chimp.mRNAs))))
  anova(lm(x.tmp~s.tmp))[1,"Pr(>F)"] # get the ANOVA p-value
}

anova.hc.prot <- function(x){
  x.tmp <- x[c(human.prots,chimp.prots)]
  s.tmp <- factor(c(rep("H",length(human.prots)),rep("C",length(chimp.prots))))
  anova(lm(x.tmp~s.tmp))[1,"Pr(>F)"]
}
mRNA.143<-apply(r.143,1,anova.hc.mRNA)
m.143.eq[,3]<-mRNA.143
mRNA.159<-apply(r.159,1,anova.hc.mRNA)
m.159.eq[,3]<-mRNA.159
prot.143<-apply(r.143,1,anova.hc.prot)
m.143.eq[,4]<-prot.143
prot.159<-apply(r.159,1,anova.hc.prot)
m.159.eq[,4]<-prot.159
m.143.eq
m.159.eq
plot(m.143.eq[,"rna.pval"],m.143.eq[,"r.anova"])
plot(m.143.eq[,"pr.pval"],m.143.eq[,"p.anova"])
plot(m.159.eq[,"rna.pval"],m.159.eq[,"r.anova"])
plot(m.159.eq[,"pr.pval"],m.159.eq[,"p.anova"])
#10 & #11 - Find common genes and scatterplot mean rna and protein values for human across sets and chimp across sets 
seq_id<-c(rownames(r.143)[c(rownames(r.143)) %in% c(rownames(r.159))])
m.mu<-matrix(0,nrow=length(seq_id),ncol=12)
rownames(m.mu)<-seq_id
colnames(m.mu)<-c("h.r.143","h.r.159","h.r.pval","h.p.143","h.p.159","h.p.pval","c.r.143","c.r.159","c.r.pval","c.p.143","c.p.159","c.p.pval")

for (id in 1:length(seq_id)){
  h.mu.143<-apply(r.143[seq_id[id],9:14],1,mean)
  m.mu[id,1]<-h.mu.143
  h.mu.159<-apply(r.159[seq_id[id],9:14],1,mean)
  m.mu[id,2]<-h.mu.159
  m.mu[id,3]<-t.test(r.143[seq_id[id],9:14],r.159[seq_id[id],9:14])$p.value
  h.mu.143<-apply(r.143[seq_id[id],18:20],1,mean)
  m.mu[id,4]<-h.mu.143
  h.mu.159<-apply(r.159[seq_id[id],18:20],1,mean)
  m.mu[id,5]<-h.mu.159
  m.mu[id,6]<-t.test(r.143[seq_id[id],18:20],r.159[seq_id[id],18:20])$p.value
  c.mu.143<-apply(r.143[seq_id[id],3:8],1,mean)
  m.mu[id,7]<-c.mu.143
  c.mu.159<-apply(r.159[seq_id[id],3:8],1,mean)
  m.mu[id,8]<-c.mu.159
  m.mu[id,9]<-t.test(r.143[seq_id[id],3:8],r.159[seq_id[id],3:8])$p.value
  c.mu.143<-apply(r.143[seq_id[id],15:17],1,mean)
  m.mu[id,10]<-c.mu.143
  c.mu.159<-apply(r.159[seq_id[id],15:17],1,mean)
  m.mu[id,11]<-c.mu.159
  m.mu[id,12]<-t.test(r.143[seq_id[id],15:14],r.159[seq_id[id],15:17])$p.value
}
plot(m.mu[,1],m.mu[,2],xlab="human 159 mRNA",ylab="human 143 mRNA",main="relationship of human mRNA levels in two samples")
plot(m.mu[,4],m.mu[,5],xlab="human 159 protein",ylab="human 143 protein",main="relationship of human protein levels in two samples")
plot(m.mu[,7],m.mu[,8],xlab="chimp 159 mRNA",ylab="chimp 143 mRNA",main="relationship of chimp mRNA  levels in two samples")
plot(m.mu[,10],m.mu[,11],xlab="chimp 159 protein",ylab="chimp 143 protein",main="relationship of chimp protein levels in two samples")
hist(m.mu[,3],xlab="human mRNA p-values",main="human mRNA p-value distribution")
hist(m.mu[,6],xlab="human protein p-values",main="human protein p-value distribution")
hist(m.mu[,9],xlab="chimp mRNA p-values",main="chimp mRNA p-value distribution")
hist(m.mu[,12],xlab="chimp protein p-values",main="chimp protein p-value distribution")
#this data illustrates a near perfect, if not perfect, linear relationship between sample 143 human mRNA expression levels and sample 159 human mRNA expression levels for the 96 identified shared genes.  The same pattern applies to chimpanzee as well.For each of the 96 genes, the scatter plots imply that if gene x is expressed highly in sample 159, it is also expressed as highly in sample 143.
#for protein expression, the relationship is less clear.  There appears to be a potential negative linear relationship between sample 143 human protein expression levels and sample 159 human protein expression levels for the 96 identified shared genes.Again,the same pattern applies to chimpanzee protein expression between samples.This indicates that if gene x is expressed highly in sample 159, that gene is actually expressed at a lower level in sample 143.
#p-value distributions for human mRNA expression between samples 143 and 159 confirm the visual inspection of the scatter plot in that the p-values are all equal to 1.The same distribution of p-values is seen in chimpanzee as well. This confirms that the two samples are not significantly different from each other and would both be drawn from the same population.  It is curious that the means for both samples for both chimpanzee and human were identical with no variation in expression levels.It makes one wonder if the samples were actually the same sample repeated twice. 
#p-value distributions for human protein expression between samples 143 and 159 show varying levels of p-values though mostly non-significant ones. The same pattern holds true in chimpanzee.This means that the two samples are not significantly different from each other and would be both drawn from the same population.

