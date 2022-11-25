


library(dplyr)
library(limma)
library(edgeR)
library(GEOquery)
library(DESeq2)
library(readr)

my_id <- "GSE32863"
bgxfile <- "GPL6884_HumanWG-6_V3_0_R0_11282955_A.bgx.gz"
targetsfile <- 'targets.txt'



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

# read in annotation
annot <- illuminaio::readBGX(bgxfile)$probes
annot <- annot[,which(colnames(annot) %in% c('Source','Symbol','Transcript','ILMN_Gene','RefSeq_ID',
                                             'Entrez_Gene_ID','Symbol','Protein_Product','Probe_Id','Probe_Type',
                                             'Probe_Start','Chromosome','Probe_Chr_Orientation','Probe_Coordinates',
                                             'Cytoband', 'Definition', 'Ontology_Component', 'Ontology_Process',
                                             'Ontology_Function', 'Synonyms'))]
annot <- annot[which(annot$Probe_Id %in% rownames(x)),]
annot <- annot[match(rownames(x), annot$Probe_Id),]

# update the target info
targetinfo <- readTargets(targetsfile, sep = '\t')
rownames(targetinfo) <- targetinfo$IDATfile
x <- x[,match(rownames(targetinfo), colnames(x))]
if (!all(colnames(x) == rownames(targetinfo)))
  stop('Target info is not aligned to expression data - they must be in the same order')




# create a custom EListRaw object
project <- new('EListRaw')
project@.Data[[1]] <- 'illumina'
project@.Data[[2]] <- targetinfo
#project@.Data[[3]] <- annot
project@.Data[[3]] <- NULL
project@.Data[[4]] <- x
project@.Data[[5]] <- NULL
project$E <- x
project$targets <- targetinfo
project$genes <- annot
project$other$Detection <- detectionpvalues



# normalize the data with the 'quantile' method, to be consistent with RMA for Affymetrix arrays
project.bgcorrect.norm <- neqc(project, offset = 16)




# filter out control probes, those with no symbol, and those that failed
annot <- annot[which(annot$Probe_Id %in% rownames(project.bgcorrect.norm)),]
project.bgcorrect.norm <- project.bgcorrect.norm[which(rownames(project.bgcorrect.norm) %in% annot$Probe_Id),]
annot <- annot[match(rownames(project.bgcorrect.norm), annot$Probe_Id),]
project.bgcorrect.norm@.Data[[3]] <- annot
project.bgcorrect.norm$genes <- annot



















Control <- project.bgcorrect.norm$genes$Source=="ILMN_Controls"
NoSymbol <- project.bgcorrect.norm$genes$Symbol == ""
isexpr <- rowSums(project.bgcorrect.norm$other$Detection <= 0.05) >= 2
project.bgcorrect.norm.filt <- project.bgcorrect.norm[!Control & !NoSymbol & isexpr, ]
dim(project.bgcorrect.norm)
dim(project.bgcorrect.norm.filt)















# remove annotation columns we no longer need
project.bgcorrect.norm.filt$genes <- project.bgcorrect.norm.filt$genes[,c(
  'Probe_Id',
  'Definition','Ontology_Component','Ontology_Process','Ontology_Function',
  'Chromosome','Probe_Coordinates','Cytoband','Probe_Chr_Orientation',
  'RefSeq_ID','Entrez_Gene_ID','Symbol')]
head(project.bgcorrect.norm.filt$genes)

# summarise across genes by mean
# ID is used to identify the replicates
project.bgcorrect.norm.filt.mean <- avereps(project.bgcorrect.norm.filt,
                                            ID = project.bgcorrect.norm.filt$genes$Symbol)
dim(project.bgcorrect.norm.filt.mean)



### diferential expression ####


design<- model.matrix(~0 + targetinfo$Group)
colnames(design)

## the column names are a bit ugly, so we will rename
colnames(design) <- c("Adjacent_non_tumor","Lung_Adenocarcinoma")

aw <- arrayWeights(project.bgcorrect.norm.filt.mean, design)

fit <- lmFit(project.bgcorrect.norm.filt.mean, design, weights= aw)


contrasts2<-makeContrasts(Adjacent_non_tumor-Lung_Adenocarcinoma,levels=design)

contr.fit2<-eBayes(contrasts.fit(fit,contrasts2))






volcanoplot(contr.fit2,main="non-tumor-Adenocarcinoma")

results <- decideTests(contr.fit2, method= "global")

vennDiagram(results)

topTable(contr.fit2,coef=1)


plotMD(contr.fit2, status= results)



















