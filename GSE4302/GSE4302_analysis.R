#use to untar file from GEO
#utils::untar("./GSE4302_RAW.tar", exdir = "./GSE4302")  #untar

dir.create("./plots")
save_dir <- "./plots/" #folder to save plots

############################ LOAD DATA ################################
library("affy")
library("GEOquery")
library("limma")
library("ggplot2")
library(pheatmap)
library(dplyr)
# specify the path on your computer where the folder that contains the CEL-files is located
celpath = "./GSE4302/"

# import CEL files containing raw probe-level data into an R AffyBatch object
raw_data = ReadAffy(celfile.path=celpath)

## quality control on raw_data
exprs(raw_data)[1:5, 1:5]



########## DATA FROM GEO ##############################################
## extract geo expression, fData, eData
Sys.setenv(VROOM_CONNECTION_SIZE = 256000000)

###extracting the file from GEO
my_id <- "GSE4302" #id from GEO
#extract data
expr <- getGEO(my_id)[[1]]

# print data info from GEO
print("Abstract from project")
print(expr@experimentData@abstract)
# get the abstract
exp_data <- expr@experimentData@abstract
# save#
fileConn<-file("GSE32863_EXPDATA.txt")
writeLines(exp_data, fileConn)
close(fileConn)


### in case data from GEO is needed
#geo_sampleInfo <- pData(expr) #expression 
#geo_edata <- exprs(expr) #to compare the edata matrix to the raw data that I'll obtain later
#geo_annot <- fData(expr) #annotation data

#####-------------------------------------------------------------------------------
################################QUALITY CONTROL ON THE RAW DATA ###################
###### PCA plot ####
gc() #clears garbage

# prepare data for PCA plot #
exp_raw <- log2(exprs(raw_data))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

data_raw <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                       Disease = pData(expr)$characteristics_ch1,
                       Phenotype = pData(expr)[36])


# ggplot PCA and save plot#
png(paste0(save_dir, "pca_rawdata.png"), width=1200, height=850)
ggplot(data_raw, aes(PC1, PC2)) +
  geom_point(aes(colour = Disease, shape = sample.type.ch1), size=4) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


rm(data_raw, exp_raw, PCA_raw, percentVar, sd_ratio) # get rid of unused files from now on

### boxplot raw intensities
png(paste0(save_dir, "boxplot_rawdata.png"), width=1200, height=750)
oligo::boxplot(raw_data, target = "core", 
               main = "Boxplot of log2-intensitites for the raw data")
dev.off()


# save density raw_data
png(paste0(save_dir, "hist_rawdata.png"), width=800, height=750)
hist(raw_data, main = "histogram of RAW data")
dev.off()






#####################################Relative Log Expression (RLE) #######################
#perform RMA without normalisation
#change names
#non-normalize data
rma_non_normalised <- affy::rma(raw_data, normalize = FALSE)

###add phenodata and feature data to raw file
rma_non_normalised@phenoData <- expr@phenoData
rma_non_normalised@featureData <- expr@featureData

#extract edata and change colnames to match pdata
edata <- exprs(rma_non_normalised)
colnames(edata) <- colnames(exprs(expr))

#change colname for sample type
#sampleInfo <- pData(rma_non_normalised)
#sampleInfo["Sample_type"] <-sampleInfo[36] # change name




###Plotting the RLE

row_medians_assayData <- Biobase::rowMedians(as.matrix(edata))

RLE_data <- sweep(edata, 1, row_medians_assayData)

RLE_data <- as.data.frame(RLE_data)

RLE_data_gathered <- 
  tidyr::gather(RLE_data, patient_array, log2_expression_deviation)

##save boxplot RLE values
png(paste0(save_dir, "Boxplot for the RLE values.png"), width=800, height=750)
ggplot(RLE_data_gathered, aes(patient_array, log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2, 2)) + 
  theme(axis.text.x = element_text(colour = "aquamarine4", 
                                   angle = 60, size = 6.5, hjust = 1 ,
                                   face = "bold"))
dev.off()

write.exprs(rma_non_normalised, file="rma_non_normalised.txt")

#get rid of data
rm(edata, RLE_data, RLE_data_gathered, rma_non_normalised, row_medians_assayData)




#Note that the y-axis now displays for each microarray the deviation of expression
#intensity from the median expression of the respective single transcripts across arrays.

#Boxes with a larger extension therefore indicate an unusually high deviation from the
#median in a lot of transcripts, suggesting that these arrays are different from most of the others in some way.

#Boxes that are shifted in y-direction indicate a systematically higher or lower
#expression of the majority of transcripts in comparison to most of the other arrays.
#This could be caused by quality issues or batch effects.

#Therefore, if shape and median of a given box varies too much from the bulk, they
#should be inspected and potentially removed.

#By inspecting the boxplot in Figure, five arrays could be considered as outliers:  are negatively y-shifted.

#We will keep these samples in mind for heatmap cluster analysis later on in the
#workflow. Arrays that are confirmed to be outliers by heatmap analysis could be removed for subsequent analysis.




#############------------------------------------------------------------################
#################### RMA CALIBRATION OF THE DATA ################
#Now, we can apply the full RMA algorithm to our data in order to background-correct, normalize and summarize:
eset_norm <- rma(raw_data, target = "core")

## ADD PHENO AND FEATURE DATA
eset_norm@phenoData <- expr@phenoData
eset_norm@featureData <- expr@featureData

#extract edata and change colnames to match pdata
edata <- exprs(eset_norm)
colnames(edata) <- colnames(exprs(expr))

#extracta all data
pdata <- pData(eset_norm)
annot <- fData(eset_norm)


###delete not needed data
rm(raw_data)
################################QUALITY CONTROL ON NORM RAW DATA ###################
###### PCA plot ####
gc() #clears garbage

# prepare data for PCA plot #
PCA_raw <- prcomp(t(edata), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

data.frame <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                       Disease = pData(eset_norm)$characteristics_ch1,
                       Phenotype = pData(eset_norm)[36])



# ggplot PCA and save plot#
png(paste0(save_dir, "pca_normdata.png"), width=1200, height=850)
ggplot(data_raw, aes(PC1, PC2)) +
  geom_point(aes(colour = Disease, shape = sample.type.ch1), size=4) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


######################## Heatmap clustering analysis
library(pheatmap)

corMatrix <- cor(edata,use="c")

png(paste0(save_dir, "heatmap.png"), width=1000, height=750)
pheatmap(corMatrix)     
dev.off()


#############histogram

medians <- rowMedians(edata)

#plot
png(paste0(save_dir, "Histogram of the median intensities per gene.png"), width=1200, height=850)
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
dev.off()

#plot
png(paste0(save_dir, " Histogram of the median intensities per gene with manual intensity filtering threshold.png"), width=1200, height=850)
man_threshold <- 4
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
abline(v = man_threshold, col = "coral4", lwd = 2)
dev.off()



###############################FLTERING INTENSITY ###########################
########cutoff
no_of_samples <- table(pdata["sample type:ch1"])
no_of_samples 
samples_cutoff <- min(no_of_samples)

idx_man_threshold <- apply(edata, 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)
#FALSE  TRUE 
# 9428 45247 
edata_manfiltered <- subset(edata, idx_man_threshold)

#filter out those ids that are in annot
annot <- annot[which(rownames(annot) %in% rownames(edata_manfiltered)),]


##### REMOVE MULTIPLE MAPPINGS

annot <- annot %>%
  rename("SYMBOL" = "Gene Symbol",
         "Gene_Title" = "Gene Title",
         "PROBEID" = "ID")

anno_grouped <- group_by(annot, PROBEID)
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))

head(anno_summarized)

anno_filtered <- filter(anno_summarized, no_of_matches > 1)
head(anno_filtered)
probe_stats <- anno_filtered 
nrow(probe_stats) # NOTHING TO FILTER OUT


########### no need to filter them #

annot_select <- select(annot, columns = "PROBEID", "SYMBOL", "Gene_Title")
head(annot_select)

###feature data, filter out no symbols ################
#change empty cells for NA
annot_select <- annot_select %>% mutate_all(na_if,"")

table(is.na(annot_select$SYMBOL))
#FALSE  TRUE 
#39023  6224 

#get rid of NA in gene Symbol
annot_select <- subset(annot_select, !is.na(SYMBOL))

####subset rows that are found in annot_select after getting rid of all the empty SYMBOLs
edata_final <- edata_manfiltered[rownames(edata_manfiltered) %in% annot_select$columns,]
edata_final


#####  FIX PDATA COLUMN ####################
library(tidyverse)

pdata <- pdata %>% rename("Sample_type" = "sample type:ch1")
pdata$Sample_type <- gsub(" ", "_", pdata$Sample_type)

sample <- factor(pdata$Sample_type)
sample

levels(sample) <- c("Healthy_control", "Smoker", "Asthmatic_at_baseline", "Asthmatic_after_Flovent", "Asthmatic_after_Placebo")
levels(sample)

design <- model.matrix(~0 +sample)
colnames(design)



contrasts <- makeContrasts(
  AsthPlacebovsHealthy = sampleAsthmatic_after_Placebo-sampleHealthy_control , 
  AsthBaselinevsHealthy = sampleAsthmatic_at_baseline-sampleHealthy_control ,
  smokervsHealthy = sampleSmoker-sampleHealthy_control,
  AsthFloventvsHealthy = sampleAsthmatic_after_Flovent-sampleHealthy_control,
  AsthFloventvsAsthBaseline = sampleAsthmatic_after_Flovent- sampleAsthmatic_at_baseline,
  AsthPlacebovsAsthBaseline = sampleAsthmatic_after_Placebo-sampleAsthmatic_at_baseline,
  AsthFloventvsAsthPlacebo = sampleAsthmatic_after_Flovent-sampleAsthmatic_after_Placebo,
  levels=design)

contrasts


# create a custom ExpressionSet object
edata_final <- edata_final[match(rownames(pdata), colnames(edata_final)),]
edata_final <- edata_final[match(rownames(ann), rownames(edata_final)),]

pd <- new("AnnotatedDataFrame", data = pdata)
ann <- new("AnnotatedDataFrame", data=annot_select)
eset <- ExpressionSet(edata_final, phenoData = pd, featureData = ann)




#calculate weigths
aw <- arrayWeights(eset, design)
fit <- lmFit(eset, design, weights= aw)
fit

## plot weights ##
png(paste0(save_dir, "weights.png"), width=1000, height=750)
barplot(aw, xlab="Array", ylab="Weight", col="white", las=2)
abline(h=1, lwd=1, lty=2)
dev.off()

## histogram fit Amean
png(paste0(save_dir, "histogram_fit_mean.png"), width=1000, height=750)
hist(fit$Amean)
dev.off()

#### plot sa ####
png(paste0(save_dir, "plot_SA.png"), width=1000, height=750)
plotSA(fit, main="Final model: Mean-variance trend")
dev.off()


#Finally, apply the empirical Bayes’ step to get our differential expression statistics and p-values.
contr.fit <- eBayes(contrasts.fit(fit, contrasts), trend = TRUE)
topTable(contr.fit)


#decidetests and save venn diagram
results <- decideTests(contr.fit, method= "global", lfc=1)

result1 = topTable(contr.fit, n=Inf, coef = "AsthPlacebovsHealthy", p.value = 0.05)
result2 = topTable(contr.fit, n=Inf, coef = "AsthBaselinevsHealthy", p.value = 0.05)
result3 = topTable(contr.fit, n=Inf, coef = "smokervsHealthy", p.value = 0.05)
result4 = topTable(contr.fit, n=Inf, coef = "AsthFloventvsHealthy", p.value = 0.05)
result5 = topTable(contr.fit, n=Inf, coef = "AsthFloventvsAsthBaseline", p.value = 0.05)
result6 = topTable(contr.fit, n=Inf, coef = "AsthPlacebovsAsthBaseline", p.value = 0.05)
result7 = topTable(contr.fit, n=Inf, coef = "AsthFloventvsAsthPlacebo", p.value = 0.05)


###### plot hists ##
library(RColorBrewer)
png(paste0(save_dir, "hist_AsthPlacebovsHealthy.png"), width=1000, height=750)
hist(result1$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "hist_AsthPlacebovsHealthy", xlab = "p-values")
dev.off()


png(paste0(save_dir, "hist_AsthBaselinevsHealthy.png"), width=1000, height=750)
hist(result2$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "AsthBaselinevsHealthy", xlab = "p-values")
dev.off()

png(paste0(save_dir, "hist_smokervsHealthy.png"), width=1000, height=750)
hist(result3$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "smokervsHealthy", xlab = "p-values")
dev.off()

png(paste0(save_dir, "hist_AsthFloventvsHealthy.png"), width=1000, height=750)
hist(result4$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "AsthFloventvsHealthy", xlab = "p-values")
dev.off()

png(paste0(save_dir, "hist_AsthFloventvsAsthBaseline.png"), width=1000, height=750)
hist(result5$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "AsthFloventvsAsthBaseline", xlab = "p-values")
dev.off()

png(paste0(save_dir, "hist_AsthPlacebovsAsthBaseline.png"), width=1000, height=750)
hist(result6$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "AsthPlacebovsAsthBaseline", xlab = "p-values")
dev.off()

png(paste0(save_dir, "hist_AsthFloventvsAsthPlacebo.png"), width=1000, height=750)
hist(result7$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "AsthFloventvsAsthPlacebo", xlab = "p-values")
dev.off()


############
## plot mds ## 
png(paste0(save_dir, "MDS_fit.png"), width=1000, height=750)
plotMD(fit, coef=1,main="Mean-Difference Plot of fit, coef=1")
dev.off()


png(paste0(save_dir, "SA_fit.png"), width=1000, height=750)
plotSA(fit,main="Residual standard deviation versus average log expression for fit")
dev.off()


## plot MDs after ebayes
png(paste0(save_dir, "MDS_fit_after_ebayes.png"), width=1000, height=750)
plotMD(contr.fit, coef=1,main="Mean-Difference Plot of fit2 (after ebayes), coef=1")
abline(0,0,col="blue")
dev.off()

png(paste0(save_dir, "MDS_fit_after_ebayes_coeff2.png"), width=1000, height=750)
plotMD(contr.fit, coef=2,main="Mean-Difference Plot of fit2 (after ebayes), coef=2")
abline(0,0,col="blue")
dev.off()


png(paste0(save_dir, "MDS_fit_after_ebayes_coeff3.png"), width=1000, height=750)
plotMD(contr.fit, coef=3,main="Mean-Difference Plot of fit2 (after ebayes), coef=3")
abline(0,0,col="blue")
dev.off()


png(paste0(save_dir, "MDS_fit_after_ebayes_coeff4.png"), width=1000, height=750)
plotMD(contr.fit, coef=4,main="Mean-Difference Plot of fit2 (after ebayes), coef=4")
abline(0,0,col="blue")
dev.off()

png(paste0(save_dir, "MDS_fit_after_ebayes_coeff5.png"), width=1000, height=750)
plotMD(contr.fit, coef=5,main="Mean-Difference Plot of fit2 (after ebayes), coef=5")
abline(0,0,col="blue")
dev.off()

png(paste0(save_dir, "MDS_fit_after_ebayes_coeff6.png"), width=1000, height=750)
plotMD(contr.fit, coef=6,main="Mean-Difference Plot of fit2 (after ebayes), coef=6")
abline(0,0,col="blue")
dev.off()

png(paste0(save_dir, "MDS_fit_after_ebayes_coeff7.png"), width=1000, height=750)
plotMD(contr.fit, coef=7,main="Mean-Difference Plot of fit2 (after ebayes), coef=7")
abline(0,0,col="blue")
dev.off()




###### venn diagram #####





################summary #
print("Sumary expression")
summary(results)



######## save results table #############
full_results <- topTable(contr.fit, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

#save data table
library(readr)
filter(full_results, adj.P.Val < 0.05, abs(F) > 1) %>%
  write_csv(paste0(save_dir,my_id,path="_filtered_DE_results.csv"))


genes_interest <- c("SLC22A1", "SLC22A4", "SLC22A5")
filter(full_results, SYMBOL %in% genes_interest) %>%
  write_csv(paste0(save_dir,my_id,path="_filtered_genes_results.csv"))



####plot genes expression#######################

























