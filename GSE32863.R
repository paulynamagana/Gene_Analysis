#https://www.biostars.org/p/473877/#474950
#https://www.biostars.org/p/424944/
#https://www.biostars.org/p/403439/#403446
  


####download raw file
#url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE32863&format=file"
#dir.create("GSE32863")
#utils::download.file(url, destfile="GSE32863/GSE32863_RAW.tar", mode="wb") 

#BiocManager::install("illuminaio")
#library(illuminaio)
#utils::untar("./GSE32863/GSE32863_RAW.tar", exdir = "./GSE32863")  #untar
#bgx <- readBGX(file.path("./GSE32863/GPL6884_HumanWG-6_V3_0_R0_11282955_A.bgx.gz")) #read file


########
readr::local_edition(1) #some problem with geoquery, run this before getGEO
library(GEOquery)
library(edgeR)
my_id <- "GSE32863"

##expression supp File
gse <- getGEOSuppFiles(my_id)

## extract geo expression, fData, eData
Sys.setenv(VROOM_CONNECTION_SIZE = 256000000)
expr <- getGEO(my_id)[[1]]

sampleInfo <- pData(expr)
edata <- exprs(expr) #to compare the edata matrix to the raw data that I'll obtain later
fdata <- fData(expr)



#################

library(readr)
library(dplyr)

####read in the data and convert to an ElistRaw object
x <- read.table('./GSE32863/GSE32863_non-normalized.txt.gz',
  header = TRUE, sep = '\t', stringsAsFactors = FALSE, skip = 0)

#extract detection p-value columns
detectionpvalues <- x[,grep('Detection.Pval', colnames(x))]
x <- x[,-grep('Detection.Pval', colnames(x))]

# set rownames and tidy up final expression matrix
probes <- x$ID_REF
x <- data.matrix(x[,2:ncol(x)])
rownames(x) <- probes
colnames(x) <- sampleInfo$geo_accession


##check normalisation and scales
summary(edata)
summary(x)
summary(project.bgcorrect.norm.filt.mean$E)


# read in annotation and align it with the expression data
anno <- fdata

# create a custom EListRaw object
project <- new('EListRaw')
project@.Data[[1]] <- 'illumina'
project@.Data[[2]] <- sampleInfo
project@.Data[[3]] <- NULL
project@.Data[[4]] <- x
project@.Data[[5]] <- NULL
project$E <- x
project$targets <- sampleInfo
project$genes <- NULL
project$other$Detection <- detectionpvalues


# for BeadArrays, background correction and normalisation are handled by a single function: neqc()
# this is the same as per Agilent single colour arrays
#
# perform background correction on the fluorescent intensities
#   'normexp' is beneficial in that it doesn't result in negative values, meaning no data is lost
#   for ideal offset, see Gordon Smyth's answer, here: https://stat.ethz.ch/pipermail/bioconductor/2006-April/012554.html
# normalize the data with the 'quantile' method, to be consistent with RMA for Affymetrix arrays
project.bgcorrect.norm <- neqc(project, offset = 16)

# filter out control probes, those with no symbol, and those that failed
annot <- anno[which(anno$ID %in% rownames(project.bgcorrect.norm)),]
project.bgcorrect.norm <- project.bgcorrect.norm[which(rownames(project.bgcorrect.norm) %in% fdata$ID),]
annot <- annot[match(rownames(project.bgcorrect.norm), annot$ID),]
project.bgcorrect.norm@.Data[[3]] <- annot
project.bgcorrect.norm$genes <- annot
Control <- project.bgcorrect.norm$genes$Source=="ILMN_Controls"
NoSymbol <- project.bgcorrect.norm$genes$Symbol == ""
isexpr <- rowSums(project.bgcorrect.norm$other$Detection <= 0.05) >= 3
project.bgcorrect.norm.filt <- project.bgcorrect.norm[!Control & !NoSymbol & isexpr, ]
dim(project.bgcorrect.norm)
dim(project.bgcorrect.norm.filt)



# summarise across genes by mean
# ID is used to identify the replicates
project.bgcorrect.norm.filt.mean <- avereps(project.bgcorrect.norm.filt,
                                            ID = project.bgcorrect.norm.filt$genes$Symbol)
dim(project.bgcorrect.norm.filt.mean)


############
##check normalisation and scales
summary(edata)
summary(x)
summary(project.bgcorrect.norm.filt.mean$E)


boxplot(exprs(expr),outline=FALSE)
boxplot(x,outline=FALSE)
boxplot(project.bgcorrect.norm.filt.mean$E,outline=FALSE)


corMatrix <- cor(project.bgcorrect.norm.filt.mean$E,use="c")
pheatmap(corMatrix)     

## Print the rownames of the sample information and check it matches the correlation matrix
rownames(sampleInfo)
colnames(corMatrix)

pheatmap(corMatrix,
         annotation_col = sampleInfo$characteristics_ch1)



library(ggplot2)
library(ggrepel)
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(project.bgcorrect.norm.filt.mean$E))

## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=characteristics_ch1,label=paste("Patient", characteristics_ch1.1))) + geom_point() + geom_text_repel()








################



###########################
library(limma)

design<- model.matrix(~0 +sampleInfo$source_name_ch1)  
design

## the column names are a bit ugly, so we will rename
colnames(design) <- c("Non_tumor","Adenocarcinoma")

#fit the model to the data
#The result of which is to estimate the expression level in each of the groups that we specified.
fit <- lmFit(project.bgcorrect.norm.filt.mean$E, design)
fit
head(fit$coefficients)


#In order to perform the differential analysis, we have to define the contrast that we are interested in
#In our case we only have two groups and one contrast of interest.
#Multiple contrasts can be defined in the makeContrasts function.
contrasts<-makeContrasts(Adenocarcinoma-Non_tumor,levels = design)
fit2 <- contrasts.fit(fit, contrasts)

#apply the empirical Bayes step to get our differential expression statistics and p-values.
fit2 <- eBayes(fit2)
topTable(fit2)


#If we want to know how many genes are differentially-expressed overall
#we can use the decideTests function.
summary(decideTests(fit2))


#########################Coping with outliers
###A compromise, which has been shown to work well is to calculate weights to define the reliability of each sample.
## calculate relative array weights
aw <- arrayWeights(x,design)
aw

#The lmFit function can accept weights, and the rest of the code proceeds as above.
fit <- lmFit(project.bgcorrect.norm.filt.mean$E, design,
             weights = aw)
contrasts <- makeContrasts(Non_tumor-Adenocarcinoma,levels = design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)



##################################

##Further processing and visualisation of DE results
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



library(ggrepel)
p_cutoff <- 0.05
fc_cutoff <- 1
topN <- 20

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, ID,"")) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")


## Get the results for particular gene of interest
filter(full_results, ID == "SLC22A1")
filter(full_results, ID == "SLC22A4")
filter(full_results, ID == "SLC22A5")


## Get results for genes with SLC22 in the name
filter(full_results, grepl("SLC22", ID))


##save csv
library(readr)
write.csv(full_results, "./GSE32863/GSE32863_normalised.csv")







































########################
samples <- read.delim("./GSE32863/GSE32863_non-normalized.txt", stringsAsFactors = TRUE, row.names = 1)

#Sys.setenv(VROOM_CONNECTION_SIZE = 25600000)
expr <- getGEO(my_id)[[1]]

sampleInfo <- pData(expr)
edata <- exprs(expr)
fdata <- fData(expr)


#delete columns that contain pvalues
gse <- select(samples, -contains("Detection"))
summary(gse)

colnames(edata)
colnames(gse)

colnames(gse) <- colnames(edata)
table(colnames(gse)==sampleInfo$geo_accession)

y <- DGEList(gse)
names(y)

#crete groups for samples
group <- sampleInfo$source_name_ch1
group
#convert to factor
group <- factor(group)
#add the group information into the DGElist
y$samples$group <- group
y$samples

gender <- sampleInfo$characteristics_ch1.10
ethnicity <- sampleInfo$characteristics_ch1.11
smoker <- sampleInfo$characteristics_ch1.13

y$samples$smoker <- smoker
y$samples$ethnicity <- ethnicity
y$samples$gender <- gender
y

## source_name_ch1 and characteristics_ch1.1 seem to contain factors we might need for the analysis. Let's pick just those columns
sampleInfo <- select(sampleInfo, source_name_ch1,characteristics_ch1.1, characteristics_ch1.10)

## Optionally, rename to more convenient column names
sampleInfo <- rename(sampleInfo,group = source_name_ch1, batch=characteristics_ch1.1, gender=characteristics_ch1.10)



############################################
#Filtering lowly expressed genes
myCPM <- cpm(gse)
head(myCPM)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
# There are 48803 genes that have TRUEs in all 12 samples.
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
summary(keep)
# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(myCPM[,1],gse[,1])


y <- y[keep, keep.lib.sizes=FALSE]

#read per sample
y$samples$lib.size
#plot library sizes
barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
title("Barplot of library sizes")

# we can also adjust the labelling if we want
barplot(y$samples$lib.size/1e06, names=colnames(y), las=2, ann=FALSE, cex.names=0.75)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")


# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

#plotMDS
plotMDS(y)

col.cell <- c("purple","orange")[y$samples$group]
data.frame(sampleInfo$group,col.cell)

# Redo the MDS with cell type colouring
plotMDS(y,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("topleft",fill=c("purple","orange"),legend=levels(sampleInfo$group))
# Add a title
title("Cell type")





# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)

# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

library(RColorBrewer)
library(ggplot2)
## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[y$samples$group]

# Plot the heatmap
heatmap(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
# Save the heatmap
png(file="./GSE32863/High_var_genes.heatmap.png", width = 1200, height = 1100, res = 200)
heatmap(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell, scale="row")
dev.off()

# Apply normalisation to DGEList object
y <- calcNormFactors(y)
y$samples


#plot
par(mfrow=c(1,2))
plotMD(logcounts,column = 7)
abline(h=0,col="grey")
plotMD(logcounts,column = 11)
abline(h=0,col="grey")


##design
# Specify a design matrix without an intercept term
design <- model.matrix(~ 0 + group)
design
## Make the column names of the design matrix a bit nicer
colnames(design) <- c("Adjacent_non_tumor", "Adenocarcinoma")
design

par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)
v


par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")



# Fit the linear model
fit <- lmFit(v)
names(fit)

cont.matrix <- makeContrasts(NonvsAdeno=Adjacent_non_tumor - Adenocarcinoma,levels=design)
cont.matrix



fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)

summa.fit <- decideTests(fit.cont)
summary(summa.fit)


# We want to highlight the significant genes. We can get this from decideTests.

plotMD(fit.cont,coef=1,status=summa.fit[,"NonvsAdeno"], values = c(-1, 1), hl.col=c("blue","red"))

# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit.cont,coef=1,highlight=100,names=fit.cont$genes$SYMBOL, main="NonvsAdeno")
