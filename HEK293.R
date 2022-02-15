# Introduction
#Install the folowing packages:

install.packages("BiocManager")
install.packages("forcats")
install.packages("stringr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readr")
install.packages("tidyr")
install.packages("survminer")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("org.Hs.eg.db")


# Importing the data
#The data from this experiment comprises the DNA methylation profiles of 28 adenocarcinomas of the lungs of never-smokers with paired adjacent nonmalignant lung tissue. We correlated differential methylation changes with gene expression changes from the same 28 samples.
#The function to download a GEO dataset is getGEO from the GEOquery package. You have to specify the ID of the dataset that you want. 

library(rmarkdown) # paged table
library(GEOquery) #geo query access
library(knitr) #tables
library(dplyr)
library(devtools)
library(ggplot2)
library(limma)
library(Glimma)
library(edgeR)
library(ggrepel)
Sys.setenv(VROOM_CONNECTION_SIZE = 25600000)
my_id <- "GSE75037"
gse <- getGEO(my_id)


#Some datasets on GEO may be derived from different microarray platforms. Therefore the object gse is a list of different datasets. You can find out how many were used by checking the length of the gse object. Usually there will only be one platform and the dataset we want to analyse will be the first object in the list (gse[[1]]).
length(gse)

# Extract the data
gse <- gse[[1]]
gse

# Exploratory analysis
#The `exprs` function can retrieve the expression values as a data frame; with one column per-sample and one row per-gene.

pdata= pData(gse) #sample information
edata= exprs(gse) #expression data
fdata = fData(gse) #gene annotation

## Inspect the clinical variables
#Data submitted to GEO contain sample labels assigned by the experimenters, and some information about the processing protocol. All these data can be extracted by the pData function.

paged_table(pdata)
## Tables for factor/character variables

#Tables are good for looking at factor or character variables, especially in phenotype data

kable(table(pdata[43]))

kable(table(pdata$characteristics_ch1.2, pdata$source_name_ch1))


## Look for missing values

# Use option useNA to include NA's in table
kable(table(pdata$characteristics_ch1,useNA="ifany"))



# is.na checks for NA values
kable(table(is.na(pdata$characteristics_ch1)))


# Check for other common missing names
sum(pdata$characteristics_ch1==" ")

# Check genomic data for NAs
sum(is.na(edata))


# Make the distribution of NA's by genes
gene_na = rowSums(is.na(edata))
kable(head(gene_na))

## Make sure dimensions match up

dim(fdata)
dim(pdata)
dim(edata)

## Look at overall distributions

#For visualisation and statistical analysis, we will inspect the data to discover what scale the data are presented in. The methods we will use assume the data are on a log2 scale; typically in the range of 0 to 16.

#The `summary` function can then be used to print the distributions.
## exprs get the expression levels as a data frame and get the distribution
kable(summary(exprs(gse)))

#A boxplot can also be generated to see if the data have been normalised. If so, the distributions of each sample should be highly similar.


boxplot(edata,outline=FALSE)


# Sample clustering and Principal Components Analysis 

#The function `cor` can calculate the correlation (on scale 0 - 1) in a pairwise fashion between all samples. This can be then visualised on a heatmap. Among the many options for creating heatmaps in R, the `pheatmap` library is one of the more popular ones. The only argument it requires is a matrix of numerical values (such as the correlation matrix).

##Let's pick just those columns that we might need for the analysis
sampleInfo <- select(pdata, 44,"characteristics_ch1.4")
## Optionally, rename to more convenient column names
sampleInfo <- dplyr::rename(sampleInfo, "group" = "histology:ch1", "patient"="characteristics_ch1.4")
paged_table(sampleInfo)


library(pheatmap)
## argument use="c" stops an error if there are any missing data points
corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix) 


#We can incorporate sample information onto the plot to try and understand the clustering. We have already created such a data frame previously (`pdata`). However, we need to take care that the rownames of these data match the columns of the correlation matrix.

rownames(sampleInfo)
colnames(corMatrix)
rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo)


#As PCA is an unsupervised method, the known sample groups are not taken into account. However, we can add labels when we plot the results. The `ggplot2` package is particularly convenient for this. The `ggrepel` package can be used to postion the text labels more cleverly so they can be read.

#It is important to *transpose* the expression matrix, otherwise R will try and compute PCA on the genes (instead of samples) and quickly run out of memory.

#### Principal Components Analysis

plotMDS(gse, labels=sampleInfo[,"group"],
        gene.selection="common")
#Plot principal components labeled by group
pca <- prcomp(t(edata))
## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste("Patient", patient))) + geom_point() + geom_text_repel()


### What happens if we spot a batch effect?

#Nothing at this stage. Provided the experimental design is sensible (i.e. representatives from all samples groups are present in each batch) we can correct for batch when we run the differential expression analysis.

### What happens if we detect outliers?

#If we suspect some samples are outliers we can remove them for further analysis.

# Differential Expression

#By far the most-popular package for performing differential expression is `limma`. The user-guide is extensive and covers the theory behind the analysis and many use-cases (Chapters 9 and 17 for single-channel data such as Illumina and Affymetrix)

#https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

#Crucially, we have to allocate the samples in our dataset to the sample groups of interest. A useful function is  `model.matrix`, which will create a *design matrix* from one of the columns in your `sampleInfo`. Here I choose `sampleInfo$group`.

#The design matrix is a matrix of `0` and `1`s; one row for each sample and one column for each sample group. A `1` in a particular row and column indicates that a given sample (the row) belongs to a given group (column).

library(limma)
design <- model.matrix(~0+sampleInfo$group)
kable(head(design))

#Column names are a bit ugly, so we will rename

#rename
colnames(design) <- c("adenocarcinoma","non_malignant")

#Count the number of samples modeled by each coefficient

colSums(design)

#In order to perform the differential analysis, we have to define the
#contrast that we are interested in. In our case we only have two groups
#and one contrast of interest.
contrasts <- makeContrasts(adenocarcinoma - non_malignant, levels=design)
contrasts


## Fit the model
#Using limma package
fit <- lmFit(exprs(gse), design)
head(fit$coefficients)

## Fit the contrasts
fit2 <- contrasts.fit(fit,contrasts = contrasts)

#calculate the t-statistics for the contrasts apply the empirical Bayes’ step to get our differential expression statistics and p-values
efit <- eBayes(fit2)

### Add genes IDs
anno <- fdata
paged_table(anno)
anno <- select(anno,Symbol, Entrez_Gene_ID,Chromosome,Cytoband)
efit$genes <- anno
topTable(efit)

#visualise
plotSA(efit, main="Final model: Mean-variance trend")

### Summarize results

#Quick look at differential expression levels
topTable(efit)

#If we want to know how many genes are differentially-expressed overall

results <- (decideTests(efit))
kable(summary(results))

#Significance is defined using an adjusted p-value cutoff that is set at 5% by default comparison between expression levels in adenocarcinoma and non_mallignant

#8622 genes are found to be down-regulated in adenocarcinoma relative to non_malignant
#10214 are up-regulated un adenocarcinoma relative to non_malignant

#Some studies require more than an adjusted p-value cut-off for a stricter definition on significance, one may require log-fold-changes (log-FCs) to be above a minimum value.
tfit <- treat(efit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
adeno_vs_non <- topTreat(tfit, coef=1, n=Inf)
head(adeno_vs_non)
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(2,16))


group=as.factor(sampleInfo$group)
##### GLIMMA ########3
#BiocManager::install("Glimma")
library(Glimma)
library(edgeR)
cpm <- cpm(edata)
lcpm <- cpm(edata, log=TRUE)
d <-glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
             side.main="Entrez_Gene_ID", counts=edata, groups=group, path = "..", launch=TRUE)
d


# save plot
library(htmlwidgets)
#htmlwidgets::saveWidget(d, "glimma-plot.html", selfcontained = T)


# Session Info
#devtools::session_info()