


library(dplyr)
library(limma)
library(edgeR)
readr::local_edition()
library(GEOquery)
library(DESeq2)

my_id <- "GSE35863"


## extract geo expression, fData, eData
Sys.setenv(VROOM_CONNECTION_SIZE = 256000000)
###extracting the file from GEO
expr <- getGEO(my_id)[[1]]

sampleInfo <- pData(expr)
edata <- exprs(expr) #to compare the edata matrix to the raw data that I'll obtain later
annot <- fData(expr) #annotation data

boxplot(edata)


#################################################################
####read in the data and convert to an ElistRaw object
baseDir <- './'

#import data
x <- read.table(paste0(baseDir, 'GSE32863_non-normalized.txt.gz'),
                header = TRUE, sep = '\t', stringsAsFactors = FALSE, skip = 0)
x <- x[order(x$ID_REF),] 


#extract pvalues
detectionpvalues <- x[,grep('Detection.Pval', colnames(x))]
rownames(detectionpvalues) <- x$ID_REF

#remove pvalues from dataframe
x <- x[,-grep('Detection.Pval', colnames(x))]

#extract probes, 1st column
probes <- x$ID_REF


#convert to matrix
x <- data.matrix(x[,2:ncol(x)])

#columns and rows names
rownames(x) <- probes
colnames(x) <- colnames(edata)
colnames(detectionpvalues) <- colnames(x)


#################################
################# data from soft file has been log2 normalised


####################################
#http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#starting-from-count-matrices

sampleInfo$source_name_ch1 <- gsub(" ", "_", sampleInfo$source_name_ch1)


dds <- DESeqDataSetFromMatrix(
  countData = round(x),
  colData = sampleInfo,
  design= ~source_name_ch1)



nrow(dds)


keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)


library(vsn)
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)

####Sample distances###

sampleDists <- dist(t(assay(vsd)))
sampleDists

library("pheatmap")
library("RColorBrewer")


sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(vsd$geo_accession, vsd$source_name_ch1, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

plotPCA(vsd, intgroup = "source_name_ch1")


#########differential expression
dds <- DESeq(dds)

#building results table
res <- results(dds)
res

mcols(res, use.names = TRUE)

summary(res)


res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)


topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("source_name_ch1"))


library("apeglm")
resultsNames(dds)
res <- lfcShrink(dds, coef="source_name_ch1_Lung_adenocarcinoma_vs_Adjacent_non.tumor_lung", type="apeglm")
plotMA(res, ylim = c(-5, 5))


hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")




anno <- select(annot,Symbol,Entrez_Gene_ID,Chromosome,Cytoband)
res$genes <- anno
head(res)


filter(res)

############################################log2 normalisation
x_log <- log2(x)

Control <- annot$Source=="ILMN_Controls"
NoSymbol <- annot$Symbol == ""
isexpr <- rowSums(detectionpvalues <= 0.05) >= 3
x_log <- x_log[!Control & !NoSymbol & isexpr, ]

sampleInfo <- select(sampleInfo, source_name_ch1,characteristics_ch1.1)

## Optionally, rename to more convenient column names


library(pheatmap)
## argument use="c" stops an error if there are any missing data points

corMatrix <- cor(x_log,use="c")
pheatmap(corMatrix)     
rownames(sampleInfo)
colnames(corMatrix)

rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo)    


library(ggrepel)
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(x_log))

## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=source_name_ch1,label=paste("batch", characteristics_ch1.1))) + geom_point() + geom_text_repel()


library(readr)
features <- select(annot,Symbol,Entrez_Gene_ID,Chromosome,Cytoband)
full_output <- cbind(features, x_log)
write_csv(full_output, path="gse_full_output.csv")

##differential expression########################################
library(limma)
design <- model.matrix(~0+sampleInfo$source_name_ch1)
design

## the column names are a bit ugly, so we will rename
colnames(design) <- c("Normal","Tumour")

summary(x_log)

## calculate median expression level
cutoff <- median(x_log)

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- x_log > cutoff

## Identify genes expressed in more than 2 samples
keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
gse <- x_log[keep,]

fit <- lmFit(gse, design)
head(fit$coefficients)


contrasts <- makeContrasts(Tumour - Normal, levels=design)

## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

fit2 <- contrasts.fit(fit, contrasts)


fit2 <- eBayes(fit2)

topTable(fit2)
decideTests(fit2)
table(decideTests(fit2))


## calculate relative array weights
aw <- arrayWeights(x_log,design)
aw

fit <- lmFit(x_log, design,
             weights = aw)
contrasts <- makeContrasts(Tumour - Normal, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend = TRUE)


anno <- select(anno,Symbol,Entrez_Gene_ID,Chromosome,Cytoband)
fit2$genes <- anno
topTable(fit2)

full_results <- topTable(fit2, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

## Make sure you have ggplot2 loaded
library(ggplot2)
ggplot(full_results,aes(x = logFC, y=B)) + geom_point()



## change according to your needs
p_cutoff <- 0.05
fc_cutoff <- 1

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()


filter(full_results,  grepl("SLC22", Symbol))
plotSA(fit2)


plotMD(fit2, column=1, status=dt[,1], main=colnames(fit2)[1], xlim=c(6,16))







