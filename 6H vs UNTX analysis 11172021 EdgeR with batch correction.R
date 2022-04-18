#https://rnabio.org/module-03-expression/0003/05/01/Batch-Correction/

library("sva") #Note this analysis requires sva (>= v3.36.0) which is only available for R (>= v4.x)
library("ggplot2")
library("gridExtra")
library("edgeR")
library("UpSetR")

setwd("~/Desktop/Berman Lab/RNA sequencing analysis HIV+METH/EdgeR analysis VC 11172021/Raw counts")
uncorrected_data <- read.delim("rawcountsun6h.txt", check.names=FALSE, stringsAsFactors=FALSE) ###loaded already mapped counts
head(uncorrected_data)

###sample_names
names(uncorrected_data) = c("ID", "VC87.6H.MA","VC87.UNTX","VC88.6H.MA","VC88.UNTX","VC89.6H.MA","VC89.UNTX")
sample_names = names(uncorrected_data)[2:length(names(uncorrected_data))]
sample_names

#review data structure
head(uncorrected_data)
dim(uncorrected_data)

#define conditions, and batch i.e each donor is considered as a batch
treatment = c("6H", "UN","6H", "UN","6H", "UN")
batch = c("1", "1", "2", "2", "3", "3")

#calculate principal components for the uncorrected data
pca_uncorrected_obj = prcomp(uncorrected_data[,sample_names])

#pull PCA values out of the PCA object
pca_uncorrected = as.data.frame(pca_uncorrected_obj[2]$rotation)

#assign labels to the data frame
pca_uncorrected[,"treatment"] = treatment
pca_uncorrected[,"batch"] = batch

#plot the PCA
#create a classic 2-dimension PCA plot (first two principal components) with conditions and batches indicated
cols <- c("6H" = "#481567FF", "UN" = "#1F968BFF")
p1 = ggplot(data=pca_uncorrected, aes(x=PC1, y=PC2, color=treatment, shape=batch))
p1 = p1 + geom_point(size=10)
p1 = p1 + stat_ellipse(type="norm", linetype=2)
p1 = p1 + labs(title="PCA, RNA-seq counts for 6H vs UN (uncorrected data)", color="Treatment", shape="Batch")
p1 = p1 + scale_colour_manual(values = cols)
p1 = p1 + theme(text=element_text(size=20), #change font size of all text
          axis.text=element_text(size=20), #change font size of axis text
          axis.title=element_text(size=20), #change font size of axis titles
          plot.title=element_text(size=20), #change font size of plot title
          legend.text=element_text(size=20), #change font size of legend text
          legend.title=element_text(size=20)) 

p1 ##can see that each donor is clustering together

#perform the batch correction
groups = sapply(as.character(treatment), switch, "6H" = 1, "UN" = 2, USE.NAMES = F)
batches = sapply(as.character(batch), switch, "1" = 1, "2" = 2, "3" = 3,USE.NAMES = F)
corrected_data = ComBat_seq(counts = as.matrix(uncorrected_data[,sample_names]), batch = batches, group = groups)
str(corrected_data)

corrected_data = cbind(uncorrected_data[,("ID")], corrected_data)
colnames(corrected_data) <- c("ID", "VC87.6H.MA","VC87.UNTX","VC88.6H.MA","VC88.UNTX","VC89.6H.MA","VC89.UNTX")

#compare dimensions of corrected and uncorrected data sets
dim(uncorrected_data)
str(uncorrected_data)
dim(corrected_data)
str(corrected_data)

#visually compare values of corrected and uncorrected data sets
head(uncorrected_data)
head(corrected_data)

#convert the corrected data to dataframe
corrected_data <- as.data.frame((corrected_data))
str(corrected_data)

###convert values from characters to integers
corrected_data$VC87.UNTX <- as.integer((corrected_data$VC87.UNTX))
corrected_data$VC87.6H.MA <- as.integer((corrected_data$VC87.6H.MA))
corrected_data$VC88.UNTX <- as.integer((corrected_data$VC88.UNTX))
corrected_data$VC88.6H.MA <- as.integer((corrected_data$VC88.6H.MA))
corrected_data$VC89.UNTX <- as.integer((corrected_data$VC89.UNTX))
corrected_data$VC89.6H.MA <- as.integer((corrected_data$VC89.6H.MA))

str(corrected_data)

#calculate principal components for the corrected data
pca_corrected_obj = prcomp(corrected_data[,sample_names])

#pull PCA values out of the PCA object
pca_corrected = as.data.frame(pca_corrected_obj[2]$rotation)

#assign labels to the data frame
pca_corrected[,"treatment"] = treatment
pca_corrected[,"batch"] = batch

#as above, create a PCA plot for comparison to the uncorrected data
cols <- c("6H" = "#481567FF", "UN" = "#1F968BFF")
p2 = ggplot(data=pca_corrected, aes(x=PC1, y=PC2, color=treatment, shape=batch))
p2 = p2 + geom_point(size=10)
p2 = p2 + stat_ellipse(type="norm", linetype=2)
p2 = p2 + labs(title="PCA, RNA-seq counts for 6H METH vs UNTREATED (batch corrected data)", color="Treatment", shape="Batch")
p2 = p2 + scale_colour_manual(values = cols)
p2 = p2 + theme(text=element_text(size=20), #change font size of all text
                axis.text=element_text(size=20), #change font size of axis text
                axis.title=element_text(size=20), #change font size of axis titles
                plot.title=element_text(size=20), #change font size of plot title
                legend.text=element_text(size=20), #change font size of legend text
                legend.title=element_text(size=20)) 
p2 

#pdf(file="6H vs UNTREATED Uncorrected-vs-BatchCorrected-PCA.pdf")
#grid.arrange(p1, p2, nrow = 2)
dev.off()

rm(pca_corrected,pca_corrected_obj,pca_uncorrected,pca_uncorrected_obj,cols,p1,p2,uncorrected_data)
rm(batch,batches,groups,sample_names,treatment)

###now lets proceed with DEG analysis using edgeR, 
corrected_data 
head(corrected_data)

data <- corrected_data
head(data)


library(edgeR)
y <- DGEList(counts=data[,2:7], genes=data[,1])
dim(y)

###annotation of genes
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

idfound <- y$genes$genes %in% mappedRkeys(org.Hs.egENSEMBL)
y <- y[idfound,]
dim(y)

##add Entrez GeneIDSto annotation
egENSEMBL <- toTable(org.Hs.egENSEMBL)
head(egENSEMBL)

m <- match(y$genes$genes, egENSEMBL$ensembl_id)
y$genes$EntrezGene <- egENSEMBL$gene_id[m]

##use EntrezGeneIDs to update gene symbols
egSYMBOL <- toTable(org.Hs.egSYMBOL)
head(egSYMBOL)

m <- match(y$genes$EntrezGene, egSYMBOL$gene_id)
y$genes$Symbol <- egSYMBOL$symbol[m]
head(y$genes)

###filtering and normalizing
##choose transcripts with highest count, keep ohe transcript per gene symbol
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$Symbol)
y <- y[!d,]
nrow(y)

###filter transcripts with low count size
keep <- rowSums(cpm(y)>1) > 3
y <- y[keep,]
dim(y)

###recompute library size
y$samples$lib.size <- colSums(y$counts)
y$samples$lib.size

##use Entrez Gene IDs as row names
rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene
y$genes$EntrezGene <- NULL

###TMM normalization 
##accounts for differences in library sizes
y <- calcNormFactors(y)
y$samples

#View a summary of the normalized counts
summary(y)

##data exploration
plotMDS(y) ##dim one separates based on treatment, dim 2 separates based on samples

#design matrix
Donor <- factor(c("VC87","VC87","VC88","VC88","VC89","VC89"))
Treatment <- factor(c("6H","UN","6H", "UN","6H","UN"))
data.frame(Sample=colnames(y),Donor,Treatment)
Treatment <- relevel(Treatment, "UN")
data.frame(Sample=colnames(y),Donor,Treatment)
design <- model.matrix(~Donor+Treatment) ##appropriate for paired samples or those with batch effects
rownames(design) <- colnames(y)
design

# estimate the NB dispersion
#install.packages("statmod")
library('statmod')
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion ##sqrt of this gives coefficient of biological variation
sqrt(y$common.dispersion)
plotBCV(y)

###differential expression
fit <- glmFit(y, design)
fit

#conduct likelihood ratio test
lrt <- glmLRT(fit)
topTags(lrt)
colnames(lrt)

#top DE tags have tiny p-values and FDR values, as well as large fold changes.
## closer look at the counts-per-million in individual samples for the top genes
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]

#total number of differentially expressed genes at 5% FDR
summary(decideTests(lrt))


####save the lrt results and counts
names(lrt)
names(lrt$table)
topTags(lrt, n=Inf, adjust.method = "BH")
df <- topTags(lrt, n=Inf, adjust.method = "BH")
DEG <- as.data.frame(df$table)
detagsnames<- rownames(df$table)
countlist <- y$counts[detagsnames, ]
results <- cbind(countlist, 
                 # Get the actual values from the 'listData' slot
                 as.data.frame(df$table))
results[1,1:6]
####save the files
#write.table(topTags(lrt, n=Inf, adjust.method = "BH"), file="6HvsUNTX_DEGstatistics11172021.csv", sep=",", row.names=TRUE)
#write.table(countlist, file="countlistfor6HvsUNTX_DEG11172021.csv", sep=",", row.names=TRUE)
write.table(results, file="Combined counts and statistics  6HvsUNTX 11172021.csv", sep=",", row.names=TRUE)

#Plot log-fold change against log-counts per million, with DE genes highlighted
plotMD(lrt)
abline(h=c(-1, 1), col="blue")

plotMD(lrt)
abline(h=c(-0.25, 0.25), col="blue")

rm(d, detagsnames,Donor,idfound,keep,m,o,Treatment,design,egENSEMBL,egSYMBOL,fit)
