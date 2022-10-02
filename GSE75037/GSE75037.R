library(readr)
library(dplyr)
library(limma)
library(edgeR)
library(GEOquery)
library(edgeR)
library(illuminaio)

my_id <- "GSE75037"

####download raw file
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75037/suppl/GSE75037_non_normalized.txt.gz"
dir.create(my_id)
utils::download.file(url, destfile="GSE75037/GSE75037_non_normalized.txt.gz", mode="wb") 

## extract geo expression, fData, eData
Sys.setenv(VROOM_CONNECTION_SIZE = 256000000)
expr <- getGEO(my_id)[[1]]

sampleInfo <- pData(expr)
edata <- exprs(expr) #to compare the edata matrix to the raw data that I'll obtain later
annot <- fData(expr) #annotation data

####read in the data and convert to an ElistRaw object
baseDir <- 'GSE75037/'

#import data
x <- read.table(paste0(baseDir, '/GSE75037_non_normalized.txt.gz'),
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
#verify rownames in x match those in annot
#dataframes should have same number of observations and rows
annot <- annot[which(annot$ID %in% rownames(x)),]

# create a custom EListRaw object
project <- new('EListRaw')
project@.Data[[1]] <- 'illumina'
project@.Data[[2]] <- sampleInfo
#project@.Data[[3]] <- annot
project@.Data[[3]] <- NULL
project@.Data[[4]] <- x
project@.Data[[5]] <- NULL
project$E <- x
project$targets <- sampleInfo
#project$genes <- annot
project$genes <- NULL
project$E <- x
project$other$Detection <- detectionpvalues

#boxplot
#The intensities vary from about 5 to 16 on the log2 scale:
boxplot(log2(project$E),range=0,ylab="log2 intensity")

#How many probes are truly expressed?
#The detection values contain p-values for testing whether each probe is more intense than the negative control probes. Small values are evidence that the probe corresponds to a truly expressed gene:
head(project$other$Detection)

project.bgcorrect.norm <- neqc(project, offset = 16)


# filter out control probes, those with no symbol, and those that failed
annot <- annot[which(annot$ID %in% rownames(project.bgcorrect.norm)),]
project.bgcorrect.norm <- project.bgcorrect.norm[which(rownames(project.bgcorrect.norm) %in% annot$ID),]
annot <- annot[match(rownames(project.bgcorrect.norm), annot$ID),]
project.bgcorrect.norm@.Data[[3]] <- annot
project.bgcorrect.norm$genes <- annot

#Filter out probes that are not expressed. We keep probes that are expressed in at least three arrays according to a detection p-values of 5%
expressed <- rowSums(project$other$Detection < 0.05) >= 3
dim(project.bgcorrect.norm)

project.bgcorrect.norm.filt <- project.bgcorrect.norm[expressed,]
dim(project.bgcorrect.norm.filt)

# summarise across genes by mean - use ID
project.bgcorrect.norm.filt.mean <- avereps(project.bgcorrect.norm.filt,
                                            ID = project.bgcorrect.norm.filt$genes$Symbol)

#MDS plot
plotMDS(project.bgcorrect.norm.filt.mean, labels=project.bgcorrect.norm.mean$targets$source_name_ch1)

ct <- factor(sampleInfo$source_name_ch1)
design<- model.matrix(~0 + ct)
colnames(design) <- c("adjacent_non_tumor", "adenocarcinoma") #rename columns

fit <- lmFit(project.bgcorrect.norm.filt.mean$E, design)

contrasts <- makeContrasts(adenocarcinoma - adjacent_non_tumor, levels=design)
contrasts

fit2 <- contrasts.fit(fit, contrasts)
efit <- eBayes(fit2)
summary(decideTests(fit2, method="global"))
topTable(efit, coef=1)

#volcano plot
library(ggplot2)
full_results <- topTable(efit, number=Inf)
full_results <- tibble::rownames_to_column(full_results, "ID")
ggplot(full_results, aes(x=logFC, y=B)) +
  geom_point()

####SA
plotSA(efit, main="Final model: Mean-variance trend")


tfit <- treat(efit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
adeno_vs_non <- topTreat(tfit, coef=1, n=Inf)
head(adeno_vs_non)


plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(2,16))


############################################




