#BiocManager::install("TCGAbiolinks")
#BiocManager::install("EDASeq")
browseVignettes('TCGAbiolinks')
library(TCGAbiolinks); library(EDASeq)
projects<-getGDCprojects()
projects[,c("id", "project_id", "released", "tumor")]

#check the database for sample sizes for each transcriptome and expression study. Need most data possible for case + control for modeling.
id <- projects$id
smpls <- list()
for(i in 1:length(id)){
  temp <- NULL
  query_Target <- NULL
  tryCatch({
    query_Target <- GDCquery(project = id[i], 
                             data.category = "Transcriptome Profiling", 
                             data.type = "Gene Expression Quantification",
                             workflow.type = "STAR - Counts")
  }, error= function(e){cat("ERROR", id[i], ":" ,conditionMessage(e), "\n")})
  if(!is.null(query_Target)){
    samplesDown_Target <- getResults(query_Target)#, cols = c("cases"))
    temp[[1]] <- table(samplesDown_Target$sample_type)
    names(temp) <- id[i]
  } else {
    temp[[1]] <- NA
    names(temp) <- id[i]
  }
  smpls <- c(smpls, temp)
}

smpls

#now can narrow down to a few based on number of cases/controls
posIDs <- c("CPTAC-3", "REBC-THYR", "TCGA-COAD", "TCGA-BRCA")
projects[projects$id %in% posIDs, ]
smpls[names(smpls) %in% posIDs[3:4]]

#LOAD DATA OF INTEREST
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")#tells type of genomic data available

#### Downloading and prepare TARGET CASE ####
TargetSamples <- GDCquery(project = "TCGA-BRCA", 
                          data.category = "Transcriptome Profiling", 
                          data.type = "Gene Expression Quantification",
                          workflow.type = "STAR - Counts")
#### obtain case information ####
CaseInfo <- getResults(TargetSamples)#, cols = c("cases"))
head(CaseInfo)

#### subset samples so that there is an equal number of cancer and control samples ####
dataPrimary_Target <- TCGAquery_SampleTypes(barcode = CaseInfo$cases, typesample = "TP") # primary tumor
dataNormal_Target <- TCGAquery_SampleTypes(barcode = CaseInfo$cases, typesample = "NT") # normal tissue
dataPrimary_Target <- dataPrimary_Target[1:113]
dataNormal_Target <- dataNormal_Target[1:113]
#### downloaded samples of interest ####
TargetSamples <- GDCquery(project = "TCGA-BRCA",
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "STAR - Counts",
                          barcode = c(dataPrimary_Target, dataNormal_Target))
#### Download the data (Note: Depending on your computer, you may not have enough RAM to process this amount of data. If this happens, please subset the data
#### to include 50 cancer and 50 normal tissue samples)
GDCdownload(TargetSamples) # will download 226 files and about 960 MB of data

###PREPARE DATA

data<-GDCprepare(TargetSamples)
brca_m<-assay(data,'unstranded')
colData(data)
rowData(data)
table(rowData(data)$gene_type)
SECoding<-data[rowData(data)$gene_type=="protein_coding",]
#below:returns data from specified rows in Summarized Experiment Object
dataPrep_Coding <- TCGAanalyze_Preprocessing(object = SECoding, cor.cut = 0.6,  datatype = "fpkm_unstrand")                      
#look for batch effects. Are quantiles similar?
boxplot(dataPrep_Coding, outline = FALSE)

#next two functions for quantile normalization
dataNorm_Coding <- TCGAanalyze_Normalization(tabDF = dataPrep_Coding, geneInfo = geneInfoHT, method = "geneLength")
dataFilt_Coding <- TCGAanalyze_Filtering(tabDF = dataPrep_Coding, method = "quantile", qnt.cut =  0.25)   
#prev. result = matrix of quantile normalized expr. data for remaining prot. coding genes
boxplot(dataNorm_Coding, outline = FALSE)

#Look at differential expression between normal and tumor
DEGsCoding <- TCGAanalyze_DEA(mat1 = dataFilt_Coding[,dataNormal_Target],
                              mat2 = dataFilt_Coding[,dataPrimary_Target],
                              pipeline="limma",
                              Cond1type = "Normal" ,
                              Cond2type = "Tumor",
                              fdr.cut = 0.01 ,
                              logFC.cut = 1,
                              method = "glmLRT", ClinicalDF = data.frame())

#DEGsCoding holds the diff.exprs. genes only (FDR pval<=.01,logFC>=1)
head(DEGsCoding)

#set-up dataframes; samples and genes with Diff Expr.only, bind clincal contin.vars for correlation
top.de.genes<-dataFilt_Coding[rownames(DEGsCoding),]
t.de.genes<-t(top.de.genes)
de.genes<-log2(top.de.genes +1)
t.log.genes<-t(de.genes)
clin<-as.data.frame(colData(data))
genomic_idx <- match(rownames(clin), rownames(t.de.genes))
clin.order<-clin[order(genomic_idx),]
age.diag<-clin.order$age_at_diagnosis
top.de.age<-data.frame(cbind(t.log.genes),age.diag)
dim(top.de.age)

tumor.stat<-clin.order$shortLetterCode
tumor.stat.1<- factor(tumor.stat, levels = c("NT", "TP"), labels = c(0, 1))
top.de.tumor<-cbind.data.frame(t.log.genes,tumor.stat.1, stringsAsFactors=TRUE)
dim(top.de.tumor)


##GENE-CLINICAL CORRELATIONS
#correlation with age_at_diagnosis
cor.age.p<-apply(top.de.age[,1:7534],2,function(x){cor.test(top.de.age$age.diag,x)$p.value})
cor.age.cor<-apply(top.de.age[,1:7534],2,function(x){cor.test(top.de.age$age.diag,x)$estimate})
cor.age.df<-data.frame('p.val'=cor.age.p,'cor'=cor.age.cor)
rownames(cor.age.df)<-colnames(top.de.age[1:7534])
cor.df.2<-cor.age.df[order(cor.age.df$cor),]
hist(cor.df.2$cor,main="Distribution of correlations between age at diagnosis and gene expression",cex.main=.65,col='seagreen',xlab='Correlation coefficients')
#although there were significant correlations between age and gene expression, the correlations
#are weak (.25 or -.25) and are likely not confounding

#are there any correlations between gene expression and the presence/absence of a tumor
#is there any relationship between most significant gene and race
#boxplot of gene expression 1 (lowest p-value and race
x<-subset(clin.order,race!="not reported")
race.rows<-rownames(x)
x$race <- ifelse(x$race %in% c("asian", "black or african american"), "other", x$race)
y<-subset(de.genes[,race.rows])
boxplot(y[1,]~x$race,xlab="Race",ylab="Expression of ENSG00000158966.16",col="gold2",cex.lab=.5)
t.test(y[1,]~x$race)



##CLIN-CLIN CORRELATIONS----
#expected correlation
age.cor<-cor(clin$days_to_birth,clin$year_of_birth,'pairwise.complete.obs')
#is there a correlation between year of birth and days to death
age.d.cor<-cor(clin$year_of_birth,clin$days_to_death,'pairwise.complete.obs')
#is there a correlation between age at diagnosis and days to death?
ind.d.cor<-cor(clin$age_at_diagnosis,clin$days_to_death,'pairwise.complete.obs')
#is there a correlation between days to collection and days to death?
col.d.cor<-cor(clin$days_to_collection,clin$days_to_death,'pairwise.complete.obs')
#is there a correlation between age at diagnosis and vital status?
vs<-as.numeric(as.factor(clin$vital_status))
age.r.cor<-cor.test(clin$age_at_diagnosis,vs);age.r.cor

#gene-gene correlation
ge.cor<-cor(t.de.genes)
ge.cor[1:3,1:3]
hist(ge.cor, main="Distribution of correlation coefficients",col="lightblue",xlab="correlation coefficient")
# Define threshold
threshold <- 0.8
# Create an empty list to store correlated gene pairs
correlated_gene_pairs <- list()
# Loop through the correlation matrix and identify correlated gene pairs
for (i in 1:(ncol(ge.cor) - 1)) {
  for (j in (i + 1):ncol(ge.cor)) {
    if (ge.cor[i, j] > threshold || ge.cor[i, j] < -threshold) {
      correlated_gene_pairs <- append(correlated_gene_pairs, list(c(rownames(ge.cor)[i], rownames(ge.cor)[j],ge.cor[i,j])))
    }
  }
}
pr.ge.cor<-matrix(nrow=length(correlated_gene_pairs),ncol=3)
for(i in 1:length(correlated_gene_pairs)){
  pr.ge.cor[i,1]<-correlated_gene_pairs[[i]][[1]]
  pr.ge.cor[i,2]<-correlated_gene_pairs[[i]][[2]]
  pr.ge.cor[i,3]<-correlated_gene_pairs[[i]][[3]]
}
pr.cor.uniq<-as.data.frame(unique(pr.ge.cor));dim(pr.cor.uniq)
neg.cor.genes <- pr.cor.uniq[order(pr.cor.uniq$V3),];neg.cor.genes[1:10,1:3]
pos.cor.genes <- pr.cor.uniq[order(pr.cor.uniq$V3,decreasing=TRUE),];pos.cor.genes[1:10,1:3]
#there are 39,980 unique gene pair correlations >.80 and 1 < -.80

###-----VARIANCE
#get variance of every DE gene
ge.var<-apply(top.de.genes,1,var)
ge.var.sort <- order(ge.var, decreasing = TRUE);rownames(top.de.genes[ge.var.sort[1:10],]);rownames(top.de.genes[ge.var.sort[7424:7434],])
hist(ge.var,main="Distribution of Gene Expression Variances",col="lavender",xlab="Log(Variance)")
#The genes with the most variance are likely the ones with the most influence/correlation to an outcome
#I might use variance to decide which genes to include in the model by looking at the distribution of the top 100,1000 and 5000 genes
#However in deciding on a model, I wouldn't want to include too many of the highly correlated genes I discovered earlier.  
#Forward/backward selection could help with this. 
de.ge.var<-apply(de.genes,1,var)
de.ge.var.sort<-order(de.ge.var,decreasing=TRUE);rownames(de.genes[ge.var.sort[1:10],]);rownames(de.genes[ge.var.sort[7424:7434],])
hist(log(de.ge.var),main="Distribution of Gene Expression Variances",col="lavender",xlab="Log(Variance)")
#####PCA----------
#Try PCA on gene expression, see what we get:
pca.genes<-prcomp(t.de.genes,scale=TRUE)
plot(pca.genes$x[,1],pca.genes$x[,2],xlab="PC1",ylab="PC2",pch=20,main="PCA of Gene Expression of BRCA Tissue Samples")
plot(pca.genes$x[,1],pca.genes$x[,3],xlab="PC1",ylab="PC3",pch=20,main="PCA of Gene Expression of BRCA Tissue Samples")
plot(pca.genes$x[,2],pca.genes$x[,3],xlab="PC2",ylab="PC3",pch=20,main="PCA of Gene Expression of BRCA Tissue Samples")

#try coloring points by tumor stats (tumor/normal)
cols<-c("dodgerblue3","orange2")#assigns a color to every tissue type
point.cols<-cols[as.factor(tumor.stat.1)]#assigns a color to every sample based on the color assigned to that tissue type
point.cols

plot(pca.genes$x[,1],pca.genes$x[,2],xlab="PC1",ylab="PC2",pch=20,main="PCA of Gene Expression of BRCA Tissue Samples",col=point.cols)
legend("topright",
       legend = levels(as.factor(clin$shortLetterCode)),
       pch = 20,
       col = cols,
       cex=.75)
plot(pca.genes$x[,1],pca.genes$x[,3],xlab="PC1",ylab="PC3",pch=20,main="PCA of Gene Expression of BRCA Tissue Samples",col=point.cols)
legend("topright",
       legend = levels(as.factor(clin$shortLetterCode)),
       pch = 20,
       col = cols,
       cex=.75)
plot(pca.genes$x[,2],pca.genes$x[,3],xlab="PC2",ylab="PC3",pch=20,main="PCA of Gene Expression of BRCA Tissue Samples",col=point.cols)
legend("topright",
       legend = levels(as.factor(clin$shortLetterCode)),
       pch = 20,
       col = cols,
       cex=.75)
#look in 3D to see if there is anything else to see here for 3PCs together
library(scatterplot3d)
pca.genes.3d=scatterplot3d(pca.genes$x[,1:3],pch=20,color=point.cols,cex.symbol=1)
genes.coords = pca.genes.3d$xyz.convert(pca.genes$x[,1:3])
legend("topright", legend = levels(as.factor(clin$shortLetterCode)),
       col = cols , pch = 16,xpd=TRUE,cex=.75)

#try clusters again but with top 200 variance genes this time
ge.var.sort.1 <- order(ge.var, decreasing = TRUE)[1:200]
top.de.var<-t.de.genes[,ge.var.sort.1]
pca.genes.var<-prcomp(top.de.var,scale=TRUE)

plot(pca.genes.var$x[,1],pca.genes.var$x[,2],xlab="PC1",ylab="PC2",pch=20,main="PCA of Gene Expression of BRCA Tissue Samples",col=point.cols)
legend("topright",
       legend = levels(as.factor(clin$shortLetterCode)),
       pch = 20,
       col = cols,
       cex=.75)
plot(pca.genes.var$x[,1],pca.genes.var$x[,3],xlab="PC1",ylab="PC3",pch=20,main="PCA of Gene Expression of BRCA Tissue Samples",col=point.cols)
legend("topright",
       legend = levels(as.factor(clin$shortLetterCode)),
       pch = 20,
       col = cols,
       cex=.75)
plot(pca.genes.var$x[,2],pca.genes.var$x[,3],xlab="PC2",ylab="PC3",pch=20,main="PCA of Gene Expression of BRCA Tissue Samples",col=point.cols)
legend("topright",
       legend = levels(as.factor(clin$shortLetterCode)),
       pch = 20,
       col = cols,
       cex=.75)
#look in 3D to see if there is anything else to see here for 3PCs together
library(scatterplot3d)
pca.genesvar.3d=scatterplot3d(pca.genes.var$x[,1:3],pch=20,color=point.cols,cex.symbol=1)
genes.coords = pca.genesvar.3d$xyz.convert(pca.genes.var$x[,1:3])
legend("topright", legend = levels(as.factor(clin$shortLetterCode)),
       col = cols , pch = 16,xpd=TRUE,cex=.75)

#Clustering is a little more apparent even, with tighter clustering, when 
#only top 200 variance genes are used.

#Variable selection-----------
#Using forward selection to avoid getting only the genes most correlated with the outcome and with each other
#generate p-values for gene expression and tumor status to look at relationship between each gene and tumor or not
hist(DEGsCoding$adj.P.Val,col='lightblue',xlab='p-values',main='Distribution of p-values')
top.p.genes<-as.data.frame(t(de.genes[order(DEGsCoding$adj.P.Val),]))
top.p.genes[1:3,1:3]

#get pseudo correlations from p-values(transform p-values to correlation coefficient scale)
C=function(p){
  1/(1+exp(log10(p)+1))}
corrs.to.outcome<-sapply(DEGsCoding$adj.P.Val,C)

#get correlations between genes
ge.cor<-cor(top.p.genes,use="pairwise")
ge.cor[1:3,1:3]

#now select the genes
N=c(2,3,5,10)
a = 0.5 
term1 = corrs.to.outcome
for (i in 1:length(N)){
  selected = which.max(corrs.to.outcome)
  for ( j in 1:(N[i]-1) ) {
    term2 = apply(abs(ge.cor[selected,,drop=F]),2,sum)
    score = term1 - a*term2
    for ( s in order(score,decreasing=T) ) {
      if ( ! s %in% selected ) {
        selected = c(selected,s)
        break
      }
    }
    if (i==1){
      selected2<-as.numeric(selected)}
    else if (i==2){
      selected3<-as.numeric(selected)}
    else if (i==3){
      selected5<-as.numeric(selected)}
    else{
      selected10<-as.numeric(selected)} 
  }    
}
two.sel<-colnames(top.p.genes[selected2])
three.sel<-colnames(top.p.genes[selected3])
five.sel<-colnames(top.p.genes[selected5])
ten.sel<-colnames(top.p.genes[selected10])
ge.cor[selected,selected]
###MODELS------------------------
#cross-validate top 10 and top 20 forward selected genes in LR and SVM
#LR----
do.check=function(data,model) {
  if ( is.null(data) ) { stop("New data must be specified") }
  if ( is.null(model) ) {
    stop("The model has not been trained yet!")}}

predictor.LR = list(
  model = NULL,
  train = function(f,data,...) {
    
    predictor.LR$model <<- glm(f,data, 
                               family="binomial",na.action="na.exclude",...)
  },
  predict=function(newdata=NULL) {
    do.check(newdata,predictor.LR$model)
    as.numeric(
      predict(predictor.LR$model,newdata,type="response") > 0.5
    )
  }
)

#bring in prediction metric assessment for LR
assess.prediction=function(truth,predicted) {
  predicted = predicted[ ! is.na(truth) ]
  truth = truth[ ! is.na(truth) ]
  truth = truth[ ! is.na(predicted) ]
  predicted = predicted[ ! is.na(predicted) ]
  cat("Total cases that are not NA: ",length(truth),"\n",sep="")
  cat("Correct predictions (accuracy): ",sum(truth==predicted), "(",signif(sum(truth==predicted)*100/length(truth),3),"%)\n",sep="")
  TP = sum(truth==1 & predicted==1)
  TN = sum(truth==0 & predicted==0)
  FP = sum(truth==0 & predicted==1)
  FN = sum(truth==1 & predicted==0)
  P = TP+FN 
  N = FP+TN 
  cat("TPR (sensitivity)=TP/P: ", signif(100*TP/P,3),"%\n",sep="")
  cat("TNR (specificity)=TN/N: ", signif(100*TN/N,3),"%\n",sep="")
  cat("PPV (precision)=TP/(TP+FP): ", signif(100*TP/(TP+FP),3),"%\n",sep="")
  cat("FDR (false discovery)=1-PPV: ", signif(100*FP/(TP+FP),3),"%\n",sep="")
  cat("FPR =FP/N=1-TNR: ", signif(100*FP/N,3),"%\n",sep="")
}

#cross-validation for LR
cross.validate=function(predictor,formula,data=NULL,
                        method="random", N=100, n.out=5,...) {
  if ( is.null(data) ) { stop("data must be specified") }
  f.str = deparse(formula)
  dependent.name = sub("\\s*~.*","",f.str)
  if ( ! dependent.name %in% names(data) ) {
    dependent.data = get(dependent.name,envir=environment(formula))
    data=cbind(dependent.data,data)
    names(data)[1] = dependent.name
  } else {
    ind = match(dependent.name,names(data))
    data = cbind( data[,ind,drop=F],data[,-ind,drop=F] )
  }
  
  truth = data[,dependent.name] 
  truth = truth[0] 
  prediction=numeric()
  for ( i in 1:N ) {
    leave.out = sample(nrow(data),size=n.out)
    training.data = data[-leave.out,,drop=F]
    test.data = data[leave.out,,drop=F]
    predictor$train(formula , data=training.data,...)
    pred=predictor$predict(test.data[,-1,drop=F])
    truth[ (length(truth)+1):(length(truth)+n.out) ] =
      test.data[,dependent.name]
    prediction = c(prediction, pred)
  }
  list(truth=truth,prediction=prediction)
}



cv.LR.2=cross.validate(predictor.LR, tumor.stat.1 ~ ., top.p.genes[two.sel])
assess.prediction(cv.LR.2$truth,cv.LR.2$prediction)

cv.LR.3=cross.validate(predictor.LR, tumor.stat.1 ~ . , top.p.genes[three.sel])
assess.prediction(cv.LR.3$truth,cv.LR.3$prediction)

cv.LR.5=cross.validate(predictor.LR, tumor.stat.1 ~ . , top.p.genes[five.sel])
assess.prediction(cv.LR.5$truth,cv.LR.5$prediction)



#NEURAL NETWORK----------------------------
library(neuralnet)
expand.formula = function(f,data=NULL) {
  f.str = deparse(f)
  if ( grepl("\\.$",f.str)[1] ) {
    if ( is.null(data) ) { stop("Shortcut formula ~. requires a dataframe") }
  } else {
    return(f)
  }
  dependent.name = sub("\\s*~.*","",f.str)
  n = names(data) 
  n = n[ n != dependent.name ]
  rhs = paste(n,collapse=" + ")
  f.str = sub("\\.$",rhs,f.str)
  f = as.formula(f.str,env=environment(f))
  return(f)
}

do.check=function(data,model) {
  if ( is.null(data) ) { stop("New data must be specified") }
  if ( is.null(model) ) {
    stop("The model has not been trained yet!")}}

predictor.NN = list(
  model = NULL,
  train = function(f,data,...) {
    
    predictor.NN$model <<- {
      f = expand.formula(f,data)
      f.str = deparse(f)
      f.str = paste(f.str,collapse="")
      dependent.name = sub("\\s*~.*","",f.str)
      if ( ! dependent.name %in% names(data) ) {
        dependent.data = get(dependent.name,envir=environment(f))
        data=cbind(dependent.data,data)
        names(data)[1] = dependent.name
      }
      has.na = apply(data,1,function(x) { any(is.na(x)) } )
      data = data[! has.na,,drop=F]
      if ( is.factor( data[,dependent.name] ) ) {
        data[,dependent.name] = as.numeric(as.vector(data[,dependent.name]))
      }
      neuralnet(f,data,...)
    }
  },
  predict=function(newdata=NULL) {
    do.check(newdata,predictor.NN$model)
    as.numeric(
      predict(predictor.NN$model,newdata,type="response") > 0.5
    )
  }
)

#bring in prediction metric assessment for NN
assess.prediction=function(truth,predicted) {
  predicted = predicted[ ! is.na(truth) ]
  truth = truth[ ! is.na(truth) ]
  truth = truth[ ! is.na(predicted) ]
  predicted = predicted[ ! is.na(predicted) ]
  cat("Total cases that are not NA: ",length(truth),"\n",sep="")
  cat("Correct predictions (accuracy): ",sum(truth==predicted), "(",signif(sum(truth==predicted)*100/length(truth),3),"%)\n",sep="")
  TP = sum(truth==1 & predicted==1)
  TN = sum(truth==0 & predicted==0)
  FP = sum(truth==0 & predicted==1)
  FN = sum(truth==1 & predicted==0)
  P = TP+FN 
  N = FP+TN 
  cat("TPR (sensitivity)=TP/P: ", signif(100*TP/P,3),"%\n",sep="")
  cat("TNR (specificity)=TN/N: ", signif(100*TN/N,3),"%\n",sep="")
  cat("PPV (precision)=TP/(TP+FP): ", signif(100*TP/(TP+FP),3),"%\n",sep="")
  cat("FDR (false discovery)=1-PPV: ", signif(100*FP/(TP+FP),3),"%\n",sep="")
  cat("FPR =FP/N=1-TNR: ", signif(100*FP/N,3),"%\n",sep="")
}

#cross-validation for NN
cross.validate=function(predictor,formula,data=NULL,
                        method="random", N=100, n.out=5,...) {
  if ( is.null(data) ) { stop("data must be specified") }
  f.str = deparse(formula)
  dependent.name = sub("\\s*~.*","",f.str)
  if ( ! dependent.name %in% names(data) ) {
    dependent.data = get(dependent.name,envir=environment(formula))
    data=cbind(dependent.data,data)
    names(data)[1] = dependent.name
  } else {
    ind = match(dependent.name,names(data))
    data = cbind( data[,ind,drop=F],data[,-ind,drop=F] )
  }
  
  truth = data[,dependent.name] 
  truth = truth[0] 
  prediction=numeric()
  for ( i in 1:N ) {
    leave.out = sample(nrow(data),size=n.out)
    training.data = data[-leave.out,,drop=F]
    test.data = data[leave.out,,drop=F]
    predictor$train(formula , data=training.data,...)
    pred=predictor$predict(test.data[,-1,drop=F])
    truth[ (length(truth)+1):(length(truth)+n.out) ] =
      test.data[,dependent.name]
    prediction = c(prediction, pred)
  }
  list(truth=truth,prediction=prediction)
}

df=as.data.frame(t(de.genes)[,two.sel])
cv.NN.2=cross.validate(predictor.NN, tumor.stat.1 ~ . , data=df[,1:2,drop=F],hidden=0,linear.output=FALSE)
assess.prediction(cv.NN.2$truth,cv.NN.2$prediction)

df=as.data.frame(t(de.genes)[,three.sel])
cv.NN.3=cross.validate(predictor.NN, tumor.stat.1 ~ . , data=df[,1:3,drop=F],hidden=4,linear.output=FALSE)
assess.prediction(cv.NN.3$truth,cv.NN.3$prediction)

df=as.data.frame(t(de.genes)[,five.sel])
cv.NN.5=cross.validate(predictor.NN, tumor.stat.1 ~ . , data=df[,1:5,drop=F],hidden=5,linear.output=FALSE)
assess.prediction(cv.NN.5$truth,cv.NN.5$prediction)



 



#Random Forest--------------------------------------
library(tree)
library(randomForest)


predictor.RF = list(
  model = NULL,
  train = function(f,data,...) {
    predictor.RF$model <<- 
      randomForest(f,data,na.action="na.exclude",importance=TRUE,
                   mtry=1)
     
  },
  predict=function(newdata=NULL) {
    do.check(newdata,predictor.RF$model)
    as.numeric(
      predict(predictor.RF$model,newdata) 
    )
  }
)

#bring in prediction metric assessment for RF
assess.prediction=function(truth,predicted) {
  x<-sum(predictor.RF$model$confusion)
  y<-sum(predictor.RF$model$confusion[1,1],predictor.RF$model$confusion[2,2])
  cat("Total cases that are not NA: ",x,"\n",sep="")
  cat("Correct predictions (accuracy): ",y/x*100,"%\n",sep="")
  TP = predictor.RF$model$confusion[2,2]
  TN = predictor.RF$model$confusion[1,1]
  FP = predictor.RF$model$confusion[1,2]
  FN = predictor.RF$model$confusion[2,1]
  P = TP+FN # total number of positives in the truth data
  N = FP+TN # total number of negatives
  cat("TPR (sensitivity)=TP/P: ", signif(100*TP/P,3),"%\n",sep="")
  cat("TNR (specificity)=TN/N: ", signif(100*TN/N,3),"%\n",sep="")
  cat("PPV (precision)=TP/(TP+FP): ", signif(100*TP/(TP+FP),3),"%\n",sep="")
  cat("FDR (false discovery)=1-PPV: ", signif(100*FP/(TP+FP),3),"%\n",sep="")
}

#cross-validation for RF
cross.validate=function(predictor,formula,data=NULL,
                        method="random", N=100, n.out=5,...) {
  if ( is.null(data) ) { stop("data must be specified") }
  f.str = deparse(formula)
  dependent.name = sub("\\s*~.*","",f.str)
  if ( ! dependent.name %in% names(data) ) {
    dependent.data = get(dependent.name,envir=environment(formula))
    data=cbind(dependent.data,data)
    names(data)[1] = dependent.name
  } else {
    ind = match(dependent.name,names(data))
    data = cbind( data[,ind,drop=F],data[,-ind,drop=F] )
  }
  
  truth = data[,dependent.name] 
  truth = truth[0] 
  prediction=numeric()
  for ( i in 1:N ) {
    leave.out = sample(nrow(data),size=n.out)
    training.data = data[-leave.out,,drop=F]
    test.data = data[leave.out,,drop=F]
    predictor$train(formula , data=training.data,...)
    pred=predictor$predict(test.data[,-1,drop=F])
    truth[ (length(truth)+1):(length(truth)+n.out) ] =
      test.data[,dependent.name]
    prediction = c(prediction, pred)
  }
  list(truth=truth,prediction=prediction)
}
df=as.data.frame(t(de.genes)[,two.sel])
cv.RF.2=cross.validate(predictor.RF, tumor.stat.1 ~ . , data=df[,1:2])
assess.prediction(cv.RF.2$truth,cv.RF.2$prediction)

df=as.data.frame(t(de.genes)[,three.sel])
cv.RF.3=cross.validate(predictor.RF, tumor.stat.1 ~ . , data=df[,1:3])
assess.prediction(cv.RF.3$truth,cv.RF.3$prediction)

df=as.data.frame(t(de.genes)[,five.sel])
cv.RF.5=cross.validate(predictor.RF, tumor.stat.1 ~ . , data=df[,1:5])
assess.prediction(cv.RF.5$truth,cv.RF.5$prediction)



#SVM----------------------
library(e1071)
#tune kernel;gamma
predictor.svm = list(
  model = NULL,
  train = function(f,data,g,...) {
    predictor.svm$model <<- svm(f,data,kernel="linear",cost=1)
  },
  predict=function(newdata=NULL) {
    # check that we got data and that the model was already trained:
    do.check(newdata,predictor.svm$model)
    # here we wrap the specifics of predicting with LR: the original
    # predict would return the probabilities, so we have to apply the
    # decision cutoff and convert to numeric in order to get predicted
    # classes. All these gory details are hidden in this wrapper:
    as.numeric(as.character(predict(predictor.svm$model,newdata))
    )
  }
)
assess.prediction=function(truth,predicted) {
  # check for missing values (we are going to compute metrics
  # on non-missing values of course)
  predicted = predicted[ ! is.na(truth) ]
  truth = truth[ ! is.na(truth) ]
  truth = truth[ ! is.na(predicted) ]
  predicted = predicted[ ! is.na(predicted) ]
  cat("Total cases that are not NA: ",length(truth),"\n",sep="")
  # overall accuracy of the test: how many cases (both positive and
  # negative) we got right:
  cat("Correct predictions (accuracy): ",sum(truth==predicted), "(",signif(sum(truth==predicted)*100/length(truth),3),"%)\n",sep="")
  # how predictions align against known training/testing outcomes:
  # TP/FP= true/false positives, TN/FN=true/false negatives:
  TP = sum(truth==1 & predicted==1)
  TN = sum(truth==0 & predicted==0)
  FP = sum(truth==0 & predicted==1)
  FN = sum(truth==1 & predicted==0)
  P = TP+FN # total number of positives in the truth data
  N = FP+TN # total number of negatives
  cat("TPR (sensitivity)=TP/P: ", signif(100*TP/P,3),"%\n",sep="")
  cat("TNR (specificity)=TN/N: ", signif(100*TN/N,3),"%\n",sep="")
  cat("PPV (precision)=TP/(TP+FP): ", signif(100*TP/(TP+FP),3),"%\n",sep="")
  cat("FDR (false discovery)=1-PPV: ", signif(100*FP/(TP+FP),3),"%\n",sep="")
  cat("FPR =FP/N=1-TNR: ", signif(100*FP/N,3),"%\n",sep="")
}
cross.validate=function(predictor,formula,data=NULL,
                        method="random", N=1000, n.out=5,...) {
  if ( is.null(data) ) { stop("data must be specified") }
  f.str = deparse(formula)
  dependent.name = sub("\\s*~.*","",f.str)
  if ( ! dependent.name %in% names(data) ) {
    dependent.data = get(dependent.name,envir=environment(formula))
    data=cbind(dependent.data,data)
    names(data)[1] = dependent.name
  } else {
    ind = match(dependent.name,names(data))
    data = cbind( data[,ind,drop=F],data[,-ind,drop=F] )
  }
  truth = data[,dependent.name] 
  truth = truth[0] 
  prediction=numeric()
  for ( i in 1:N ) {
    leave.out = sample(nrow(data),size=n.out)
    training.data = data[-leave.out,,drop=F]
    test.data = data[leave.out,,drop=F]
    predictor$train(formula , data=training.data,...)
    pred=predictor$predict(test.data[,-1,drop=F])
    truth[ (length(truth)+1):(length(truth)+n.out) ] =
      test.data[,dependent.name]
    prediction = c(prediction, pred)
  }
  list(truth=truth,prediction=prediction)
}

cv.SVM.2=cross.validate(predictor.svm, tumor.stat.1 ~ . , top.p.genes[two.sel])
assess.prediction(cv.SVM.2$truth,cv.SVM.2$prediction)

cv.SVM.3=cross.validate(predictor.svm, tumor.stat.1 ~ . , top.p.genes[three.sel])
assess.prediction(cv.SVM.3$truth,cv.SVM.3$prediction)

cv.svm.5=cross.validate(predictor.svm, tumor.stat.1 ~ .,top.p.genes[five.sel] )
assess.prediction(cv.svm.5$truth,cv.svm.5$prediction)


#Final Model Choices
cv.LR.2=cross.validate(predictor.LR, tumor.stat.1 ~ . , top.p.genes[two.sel])
assess.prediction(cv.LR.2$truth,cv.LR.2$prediction)

df=as.data.frame(top.p.genes[,three.sel])
cv.NN.3=cross.validate(predictor.NN, tumor.stat.1 ~ . , data=df[,1:3,drop=F],hidden=4,linear.output=FALSE)
assess.prediction(cv.NN.3$truth,cv.NN.3$prediction)

cv.RF.3=cross.validate(predictor.RF, tumor.stat.1 ~ . , data=df[,1:3,drop=F])
assess.prediction(cv.RF.3$truth,cv.RF.3$prediction)

cv.svm.5=cross.validate(predictor.svm, tumor.stat.1 ~ .,top.p.genes[five.sel] )
assess.prediction(cv.svm.5$truth,cv.svm.5$prediction)




#---------OTHER INVESTIGATIONS
#boxplot of gene expression 1 (lowest p-value and outcome)
boxplot(top.de.genes[1,]~tumor.stat.1)
#boxplot of gene expression 1 (lowest p-value and race, white is significantly higher expression)
x<-subset(clin.order,race!="not reported")
y<-subset(top.de.genes[,colnames(top.de.genes)%in%rownames(x)])
boxplot(y[1,]~x$race)
anova.race <- aov(y[1,] ~ x$race)

#is there a difference in outcome based on race?  
sig.gene1<-t.log.genes[race.rows,1]
sig.gene2<-t.log.genes[race.rows,2]
sig.gene3<-t.log.genes[race.rows,3]
t.s.1<-top.de.tumor[race.rows,]$tumor.stat.1
df1=cbind.data.frame(sig.gene1,sig.gene2,sig.gene3,race=x$race,t.s.1,stringsAsFactors=TRUE)
t.race.table<-table(df1$t.s.1,df1$race)
fisher.test(t.race.table)

#model of just most significant gene and outcome
ge.model <- glm(t.s.1 ~ sig.gene1, data = df1, family = binomial)
summary(ge.model)

#is there an additive effect of gene and race?
ge.race.model <- glm(t.s.1 ~ race + sig.gene, data = df1, family = binomial)
summary(ge.race.model)

#is there an interaction effect of gene and race?
ge.race.itx <- glm(t.s.1 ~ race * sig.gene, data = df1, family = binomial)
summary(ge.race.itx)

#x-val of most sig. gene
cv.LR.1=cross.validate(predictor.LR, t.s.1 ~ sig.gene2, data=df1)
assess.prediction(cv.LR.1$truth,cv.LR.1$prediction)


class_0_data <- top.de.tumor[top.de.tumor$tumor.stat.1 == 0, ]
class_1_data <- top.de.tumor[top.de.tumor$tumor.stat.1 == 1, ]

# Create a histogram for class 0
hist(class_0_data$ENSG00000136158.12, col = "blue", main = "Histogram of Input Variable by Class",
     xlab = "Input Variable", ylab = "Frequency", xlim = c(min(top.de.tumor$ENSG00000136158.12), max(top.de.tumor$ENSG00000136158.12)))

# Add a histogram for class 1 with a different color
hist(class_1_data$ENSG00000136158.12, col = "red", add = TRUE)

# Add a legend to differentiate the classes
legend("topright", legend = c("Class 0", "Class 1"), fill = c("blue", "red"))


race1<-as.factor(x$race)
levels(race1)<-c(0,1)
ge.race.itx <- glm(t.s.1 ~ race1 + sig.gene, data = df1, family = binomial)
summary(ge.race.itx)
ge.race.itx1 <- glm(t.s.1 ~ race1 * sig.gene, data = df1, family = binomial)
summary(ge.race.itx1)
plot(ge.race.itx)


race.tumor<-subset(top.de.tumor[race.rows,])
top.p.g.r<-cbind(top.p.genes[race.rows,two.sel],race1,tissue=race.tumor$tumor.stat.1)

cv.LR.2=cross.validate(predictor.LR, tissue ~.,data=top.p.g.r)
assess.prediction(cv.LR.2$truth,cv.LR.2$prediction)

# Gene annotation of the 5 genes in the models
library(Homo.sapiens)
five.sel.1 <- sub("\\..*", "", five.sel)
select(org.Hs.eg.db,keys=five.sel.1,keytype="ENSEMBL",
       columns=c("ENSEMBL","SYMBOL","GENENAME"))

#Gene annotation of the 3 differentially expressed genes most closely related to the outcome
top.3.genes<-colnames(t.log.genes[,1:3])
top.3.genes.1<-sub("\\..*", "", top.3.genes)
select(org.Hs.eg.db,keys=top.3.genes.1,keytype="ENSEMBL",
       columns=c("ENSEMBL","SYMBOL","GENENAME"))

summary(glm(tumor.stat.1~t.log.genes[,2], family=binomial))
summary(glm(tumor.stat.1~t.log.genes[,1]+t.log.genes[,2], family=binomial))

