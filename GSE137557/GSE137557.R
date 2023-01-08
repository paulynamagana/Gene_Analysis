#use to untar file from GEO
utils::untar("./GSE137557_RAW.tar", exdir = "./GSE137557")  #untar

dir.create("./plots")
save_dir <- "./plots/" #folder to save plots

############################ LOAD DATA ################################
library("affy")
library("GEOquery")
library("limma")
library("ggplot2")
library("pheatmap")
library("dplyr")
library("oligo") #load withj oligo
library(xps)

# specify the path on computer where the folder that contains the CEL-files is located
celpath = "./GSE137557/"

# import CEL files containing raw probe-level data into an R AffyBatch object
celFiles <- list.files("./GSE137557/", full.names=TRUE)
raw_data <- read.celfiles(celFiles)


## quality control on raw_data
exprs(raw_data)[1:5, 1:5]



##########------------------- DATA FROM GEO ##############################################
## extract geo expression, fData, eData
Sys.setenv(VROOM_CONNECTION_SIZE = 256000000)

###extracting the file from GEO
my_id <- "GSE137557" #id from GEO
#extract data
expr <- getGEO(my_id)[[1]]


# print data info from GEO
print("Abstract from project")
print(expr@experimentData@abstract) # print the abstract
exp_data <- expr@experimentData@abstract # get the abstract
fileConn<-file(paste0(my_id,"_EXPDATA.txt"))# save#
writeLines(exp_data, fileConn) #save file
close(fileConn)

 ###############PCA##############
gc() #clears garbage

# prepare data for PCA plot #
exp_raw <- log2(exprs(raw_data))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

data_raw <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                       Disease = pData(expr)$characteristics_ch1,
                       Phenotype = pData(expr)[36])



############# build GeneFeatureSet with raw data ####################
geo_sampleInfo <- pData(expr) #expression 
geo_annot <- fData(expr) #annotation data
geo_edata <- exprs(expr)
edata <- exprs(raw_data) #rawdata

rownames(edata) <- rownames(geo_edata)




