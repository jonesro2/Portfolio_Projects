library(GEOquery)
library(preprocessCore)
library(samr)
library(impute)
library(shiny)
library(GSA)
library(dplyr)
library(stringr)
library(factoextra)
library(cluster)
library(ComplexHeatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(hta20transcriptcluster.db)
library("fgsea")
library(WGCNA)
library("tidyverse")
library("gridExtra")


#download data of interest
my_id <- "GSE143754"
gse <- getGEO(my_id)
str(gse)
#look at the data features
#View(exprs(gse[[1]]))
#View(pData(gse[[1]]))
#View(fData(gse[[1]]))
#check for normalization; values fall between 0 and 16, so likely to be log2 transformed already
summary(exprs(gse[[1]]))
#boxplot to check normalization further; looks good
boxplot(exprs(gse[[1]]),outline=FALSE)

#save sample and gene data to variables
samps<-pData(gse[[1]])
exp_mat<-exprs(gse[[1]])

#subset samples dataframe to exclude the PDAC samples
smpls <- samps %>%
  dplyr::select("title","disease state:ch1")%>%
  dplyr::rename(patient=title,group="disease state:ch1")%>%
  dplyr::filter(group!="Tumor")

#subset gene df to only samples remaining in sample df
cp_dat<-exp_mat[,colnames(exp_mat)%in%rownames(smpls)]

#get gene symbols from features data into gene expr df
#annotation for limma -FDR corrected genes
features<-fData(gse[[1]])
mapping <- AnnotationDbi::mapIds(
  hta20transcriptcluster.db,
  keys = rownames(features),
  column = 'SYMBOL',
  keytype = 'PROBEID',multiVals='first')

mapping2 <- AnnotationDbi::mapIds(
  hta20transcriptcluster.db,
  keys = rownames(features),
  column = 'ENSEMBL',
  keytype = 'PROBEID',multiVals='first')


cp_dat_sym<-as.data.frame(cbind(cp_dat,sym=mapping,ID=mapping2))

#Prepare data with named genes from data set only
sym.na<-is.na(cp_dat_sym$sym)
cp_named<-cp_dat_sym[!sym.na,]
cp_named$ID<-NULL
unique(cp_named$sym)

#use data.table package for so many genes...need to get only unique genes for GSEA
library(data.table)

# Convert data.frame to data.table
setDT(cp_named)

# Specify columns to be averaged
columns <- names(cp_named)[1:15]  
cp_named[, (columns) := lapply(.SD, as.numeric), .SDcols = columns]

# Group by 'sym' column and calculate mean for other columns
cp_dat_avg <- cp_named[, lapply(.SD, mean, na.rm = TRUE), by = sym, .SDcols = columns]
cp_dat_avg<-as.data.frame(cp_dat_avg)
rownames(cp_dat_avg)<-cp_dat_avg$sym
cp_dat_nsym <- cp_dat_avg[, -which(names(cp_dat_avg) %in% "sym")]

library(limma)
design<-model.matrix(~0+factor(smpls$group))
colnames(design) <- c("Normal","CP")
fit<-(lmFit(cp_dat_nsym,design))
head(fit$coefficients)
contrasts <- makeContrasts(CP - Normal, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
#get pvalues <.05 and see how many genes there are
sig.p.vals<-fit2$p.value<.05

#get adjusted pvalues <.05 and see how many genes there are
pval.adj<-p.adjust(fit2$p.value,method="BH")
sig.p.adj<-pval.adj<.05
dim(cp_dat_nsym[sig.p.adj,])

#save gene rankings and subset to just ones with adjust p<.05
all.genes<-topTable(fit2,number=23922)
adj_limma_genes<-all.genes[all.genes$adj.P.Val<.05,]
adj_limma_genes2<-adj_limma_genes[abs(adj_limma_genes$logFC)>1,]


#just limma -FDR corrected genes
lfc<-adj_limma_genes2$logFC
adj.p<-adj_limma_genes2$adj.P.Val
DEGs<-cp_dat_nsym[rownames(adj_limma_genes2),]
DEGS_sym<-cp_dat_avg[rownames(adj_limma_genes2),]
DEGS_sym<-data.frame(DEGS_sym, LFC = lfc, adjp=adj.p)
top_5_genes <- DEGS_sym %>%
  arrange(desc(abs(lfc))) %>%
  head(5)

#hierarchical clustering genes
d.gene<-dist(DEGs)
avg.gene.clust<-hclust(d.gene,method="average")
ward.gene.clust<-hclust(d.gene,method="ward.D2")
plot(ward.gene.clust,main="Ward",labels=FALSE)
abline(h = 40, col = "green", lwd = 2)
h.clusters<-cutree(ward.gene.clust,h=40)
gne.clusters<-as.data.frame(h.clusters)
table(gne.clusters)
colnames(DEGs)<-smpls$group
Heatmap(DEGs, show_row_names = FALSE,clustering_method_rows="ward.D2",clustering_method_columns = "ward.D2",row_split = h.clusters,heatmap_legend_param = list(at = sort(unique(h.clusters)), labels = sort(unique(h.clusters)),column_labels=smpls$group))
#after multiple clusterings, a height of 40 results in 3 clusters that all separate based on tissue type
#get genes in specific clusters for further analysis

gc_1<-rownames(gne.clusters)[gne.clusters$h.clusters==1]
gc_1.df<-DEGs[rownames(DEGs)%in%gc_1,]
gc_2<-rownames(gne.clusters)[gne.clusters$h.clusters==2]
gc_2.df<-DEGs[rownames(DEGs)%in%gc_2,]
gc_3<-rownames(gne.clusters)[gne.clusters$h.clusters==3]
gc_3.df<-DEGs[rownames(DEGs)%in%gc_3,]


Heatmap(gc_1.df,show_row_names=TRUE,clustering_method_rows="ward.D2",clustering_method_columns = "ward.D2",row_names_gp = gpar(fontsize =4))#clusters by samples
Heatmap(gc_2.df,show_row_names=TRUE,clustering_method_rows="ward.D2",clustering_method_columns = "ward.D2",row_names_gp = gpar(fontsize =4))
Heatmap(gc_3.df,show_row_names=TRUE,clustering_method_rows="ward.D2",clustering_method_columns = "ward.D2",row_names_gp = gpar(fontsize =4))#clusters on samples

#Look at ontologies for Hierarchical 
gc_1sym.df<-DEGS_sym[rownames(DEGS_sym)%in%gc_1,]
gc1_up<-gc_1sym.df[gc_1sym.df$LFC>0,] #27 are upreg.
gc1_dn<-gc_1sym.df[gc_1sym.df$LFC<0,] #100 are downreg.

gc_2sym.df<-DEGS_sym[rownames(DEGS_sym)%in%gc_2,]
gc2_up<-gc_2sym.df[gc_2sym.df$LFC>0,] #11 are upreg.
gc2_dn<-gc_2sym.df[gc_2sym.df$LFC<0,] #20 are downreg.

gc_3sym.df<-DEGS_sym[rownames(DEGS_sym)%in%gc_3,]
gc3_up<-gc_3sym.df[gc_3sym.df$LFC>0,] #56 are upreg.
gc3_dn<-gc_3sym.df[gc_3sym.df$LFC<0,] #8 are downreg.

library(topGO)

#define the gene universe and set levels for upregulated genes of interest in cluster 1
all.gene.ids<-cp_dat_avg$sym
names.dn1<-gc1_dn$sym
all.genes.dn1<-factor(as.integer( all.gene.ids %in% names.dn1 ) )
names(all.genes.dn1) = all.gene.ids

names.up1<-gc1_up$sym
all.genes.up1<-factor(as.integer( all.gene.ids %in% names.up1 ) )
names(all.genes.up1) = all.gene.ids
godata.up1=new("topGOdata",ontology="MF",allGenes=all.genes.up1,annot=annFUN.org,
               mapping="org.Hs.eg.db",ID="SYMBOL")
resultFisher.up1 <- runTest(godata.up1, algorithm = "classic",
                            statistic = "fisher")
GenTable(godata.up1,classicFisher=resultFisher.up1,
         ranksOf="classicFisher",topNodes=20)

godata.up1=new("topGOdata",ontology="BP",allGenes=all.genes.up1,annot=annFUN.org,
               mapping="org.Hs.eg.db",ID="SYMBOL")
resultFisher.up1 <- runTest(godata.up1, algorithm = "classic",
                            statistic = "fisher")
GenTable(godata.up1,classicFisher=resultFisher.up1,
         ranksOf="classicFisher",topNodes=20)

#make GO object, fisher test, get table 
godata.dn=new("topGOdata",ontology="MF",allGenes=all.genes.dn1,annot=annFUN.org,
              mapping="org.Hs.eg.db",ID="SYMBOL")
resultFisher.dn <- runTest(godata.dn, algorithm = "classic",
                           statistic = "fisher")
GenTable(godata.dn,classicFisher=resultFisher.dn,
         ranksOf="classicFisher",topNodes=20)

godata.dn=new("topGOdata",ontology="BP",allGenes=all.genes.dn1,annot=annFUN.org,
              mapping="org.Hs.eg.db",ID="SYMBOL")
resultFisher.dn <- runTest(godata.dn, algorithm = "classic",
                           statistic = "fisher")
GenTable(godata.dn,classicFisher=resultFisher.dn,
         ranksOf="classicFisher",topNodes=20)
#define the gene universe and set levels for upregulated genes of interest in cluster 4
all.gene.ids<-cp_dat_avg$sym
names.dn2<-gc2_dn$sym
all.genes.dn2<-factor(as.integer( all.gene.ids %in% names.dn2 ) )
names(all.genes.dn2) = all.gene.ids

names.up2<-gc2_up$sym
all.genes.up2<-factor(as.integer( all.gene.ids %in% names.up2 ) )
names(all.genes.up2) = all.gene.ids

godata.dn2=new("topGOdata",ontology="MF",allGenes=all.genes.dn2,annot=annFUN.org,
              mapping="org.Hs.eg.db",ID="SYMBOL")
resultFisher.dn2 <- runTest(godata.dn2, algorithm = "classic",
                           statistic = "fisher")
GenTable(godata.dn2,classicFisher=resultFisher.dn2,
         ranksOf="classicFisher",topNodes=20)

godata.up2=new("topGOdata",ontology="MF",allGenes=all.genes.up2,annot=annFUN.org,
               mapping="org.Hs.eg.db",ID="SYMBOL")
resultFisher.up2 <- runTest(godata.up2, algorithm = "classic",
                            statistic = "fisher")
GenTable(godata.up2,classicFisher=resultFisher.up2,
         ranksOf="classicFisher",topNodes=20)

#define the gene universe and set levels for upregulated genes of interest in cluster 3
all.gene.ids<-cp_dat_avg$sym
names.dn3<-gc3_dn$sym
all.genes.dn3<-factor(as.integer( all.gene.ids %in% names.dn3 ) )
names(all.genes.dn3) = all.gene.ids

godata.dn3=new("topGOdata",ontology="MF",allGenes=all.genes.dn3,annot=annFUN.org,
               mapping="org.Hs.eg.db",ID="SYMBOL")
resultFisher.dn3 <- runTest(godata.dn3, algorithm = "classic",
                            statistic = "fisher")
GenTable(godata.dn3,classicFisher=resultFisher.dn3,
         ranksOf="classicFisher",topNodes=20)

names.up3<-gc3_up$sym
all.genes.up3<-factor(as.integer( all.gene.ids %in% names.up3 ) )
names(all.genes.up3) = all.gene.ids

godata.up3=new("topGOdata",ontology="MF",allGenes=all.genes.up3,annot=annFUN.org,
               mapping="org.Hs.eg.db",ID="SYMBOL")
resultFisher.up3 <- runTest(godata.up3, algorithm = "classic",
                            statistic = "fisher")
GenTable(godata.up3,classicFisher=resultFisher.up3,
         ranksOf="classicFisher",topNodes=20)


#WGCNA
dataExpr<-t(DEGs)
softPower <- pickSoftThreshold(dataExpr, powerVector = seq(from = 1, to = 20, by = 1), networkType = "signed")
net = blockwiseModules(dataExpr, power = 7,
                       TOMType = "signed", minModuleSize = 10,
                       reassignThreshold =10, mergeCutHeight = 0.1,
                       numericLabels = FALSE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase="TOM", verbose=3, ds=6)

mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#save modules
genes=colnames(dataExpr)
moduleColors=labels2colors(net$colors)
mymodules=cbind(genes,moduleColors)
moduleColors=labels2colors(net$colors)
mymodules=cbind(genes,moduleColors)

#correlate with phenotype
phenotype=list(0,0,0,0,0,1,1,1,1,0,0,0,0,1,1) #if you have 12 arrays at timepoint 0, 1 and 7 your list would be: (0,0,0,0,1,1,1,7,7,7,7)
#get Eigengenes
MEs0=moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
module_order = names(MEs) %>% gsub("ME","", .)
moduleTraitCor = cor(MEs, phenotype, use = "p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor, 15)
MEs$type<-as.factor(unlist(phenotype))
mME = MEs %>%
  pivot_longer(-type) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )
mME %>% ggplot(., aes(x=type, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")


tom<-load("TOM-block.1.RData")
tomdata<-get(tom)
str(tomdata)

exportNetworkToCytoscape(tomdata, edgeFile="edge_info_simple1.tsv", nodeFile="node_info_simple1.tsv", weighted=TRUE, threshold=0.3
                         , nodeNames=m.mymodules$genes, altNodeNames=NULL, nodeAttr=NULL, includeColNames=TRUE)
m.mymodules<-data.frame(mymodules)
write.table(mymodules, file = "mymodules.txt", sep = "\t", quote = FALSE, row.names = FALSE)
