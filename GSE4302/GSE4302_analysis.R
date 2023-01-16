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
library(annotate)
library(hgu133plus2.db)



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


#rm(data_raw, exp_raw, PCA_raw, percentVar, sd_ratio) # get rid of unused files from now on

### boxplot raw intensities
png(paste0(save_dir, "boxplot_rawdata.png"), width=1200, height=750)
oligo::boxplot(raw_data, target = "core", 
               main = "Boxplot of log2-intensitites for the raw data")
dev.off()


# save density raw_data
png(paste0(save_dir, "hist_rawdata.png"), width=800, height=750)
hist(raw_data, main = "histogram of RAW data")
dev.off()



geo_pdata <-pData(expr)


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
#rm(edata, RLE_data, RLE_data_gathered, rma_non_normalised, row_medians_assayData)




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


## retrieve symbols from probeids
eset_norm_data<- exprs(eset_norm)



#extract the rownames for hg133plus2
ID <- rownames(eset_norm_data)
## extract SYMBOl names from the IDs from probes
gene.symbols <- getSYMBOL(ID, "hgu133plus2")
gene.symbols


#add gene.symbols as fdata for the eset_norm
fData(eset_norm) <- data.frame(Symbol=gene.symbols)

#extract the fdata for corroboration
pdata<- fData(eset_norm)
pdata$probe <- rownames(pdata)

gname <- as.data.frame(gene.symbols)
gname$probe <- rownames(gname)

gname <- gname[match(gname$probe, rownames(eset_norm_data)),]


###### later
eset_norm_data <- as.data.frame(eset_norm_data) #convert into a dataframe to add column for probes and symbols
##add column probes that is the same as rownames
eset_norm_data$probe <- row.names(eset_norm_data)

mydata <- merge(eset_norm_data, gname, by = "probe")

###filtered those with na in symbols
mydatafilt <- subset(mydata, !is.na(gene.symbols)) ###MIGHT HAVE TO ADJUST I CHANGED IT TO SYMBOL BEFORE


#### filterig those with no symbols
NoSymbol <- mydatafilt$gene.symbols == ""


##averge expression by symbol column
A <- limma::avereps(mydatafilt, mydatafilt$gene.symbols)


#filter so only those probes in a are in annot data
annot <- fData(expr) #extract annot 
names(annot)[names(annot) == "Gene Symbol"] <- "SYMBOL" #change column name to Symbol
annot <- annot[which(annot$ID %in% A),] #subset to get only those in A

rows_A <- annot$ID

#match order of annot to A 
#A <- as.data.frame(A)
#drop columns
A<- A[,-1]
A<- A[,-119]

pdata <- pData(expr)
colnames(A)<- pdata$geo_accession
rownames(A)<- rows_A


library(dplyr)
#A <- as.numeric(A)
A<- as.data.frame(A)

df2 <- mutate_all(A, function(x) as.numeric(as.character(x)))


pd <- new("AnnotatedDataFrame", data = pdata) #conver into datafrme
ann <- new("AnnotatedDataFrame", data = annot) #convert annot that has been filtered into dataframee


## check namnes
featureNames(ann)
rownames(df2)

identical(featureNames(ann), rownames(df2))

colnames(df2)
identical(featureNames(pd), colnames(df2)) #TRUE


norm_data <- new("ExpressionSet", exprs=as.matrix(df2), phenoData = pd, featureData = ann)

#### complete the expressionset
norm_data@experimentData <- expr@experimentData #copy experimentdata from the geo dataset




#extract edata and change colnames to match pdata
edata <- exprs(norm_data)

#extracta all data
pdata <- pData(norm_data)
annot <- fData(norm_data)


###delete not needed data
#rm(raw_data)
############################################QUALITY CONTROL ON NORM RAW DATA ###################
###### PCA plot ####
gc() #clears garbage

# prepare data for PCA plot #
PCA_raw <- prcomp(t(edata), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

data.frame <- data.frame(PC1 = PCA_raw$x[,1],
                         PC2 = PCA_raw$x[,2],
                       Disease = pData(norm_data)$characteristics_ch1,
                       Phenotype = pdata$`sample type:ch1`)



# ggplot PCA and save plot#
png(paste0(save_dir, "pca_normdata.png"), width=1200, height=850)
ggplot(data.frame, aes(PC1, PC2)) +
  geom_point(aes(color = Phenotype), size=4) +
  ggtitle(paste0(my_id, " PCA plot of the log-transformed raw expression data")) +
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
                 main = paste0(my_id," Histogram of the median intensities"), 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
dev.off()

#plot
png(paste0(save_dir, "Histogram of the median intensities per gene with manual intensity filtering threshold.png"), width=1200, height=850)
man_threshold <- 4
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = paste0(my_id," Histogram of the median intensities"), 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
abline(v = man_threshold, col = "coral4", lwd = 2)
dev.off()



###############################FILTERING INTENSITY ###########################
########cutoff
no_of_samples <- table(pdata["sample type:ch1"])
no_of_samples 

samples_cutoff <- min(no_of_samples)
samples_cutoff

idx_man_threshold <- apply(norm_data, 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})


table(idx_man_threshold)
#TRUE 
#20825 

#########filterin for intensity
manfiltered <- subset(norm_data, idx_man_threshold)



#extract all data
pdata <- pData(expr)


library(tidyverse)

names(pdata)[names(pdata) == "sample type:ch1"] <- "Sample_type"
pdata$Sample_type <- gsub(" ", "_", pdata$Sample_type)

sample <- factor(pdata$Sample_type)
sample

levels(sample) <- c("Healthy_control", "Smoker", "Asthmatic_at_baseline", "Asthmatic_after_Flovent", "Asthmatic_after_Placebo")
levels(sample)

design <- model.matrix(~0 +sample)
colnames(design)

colnames(design) <- c("Healthy_control", "Smoker", "Asthmatic_at_baseline", "Asthmatic_after_Flovent", "Asthmatic_after_Placebo")
design

contrasts <- makeContrasts(Asthmatic_at_baselinevshealthy = Asthmatic_at_baseline - Healthy_control,
  levels=design)

contrasts


#calculate weigths
aw <- arrayWeights(manfiltered, design)
fit <- lmFit(manfiltered, design, weights= aw)
fit





## plot weights ##
png(paste0(save_dir, "weights.png"), width=1000, height=750)
barplot(aw, xlab="Array", ylab="Weight", col="white", las=2, main = paste0(my_id, " weights"))
abline(h=1, lwd=1, lty=2)
dev.off()

## histogram fit Amean
png(paste0(save_dir, "histogram_fit_mean.png"), width=1000, height=750)
hist(fit$Amean, main = paste0(my_id, " Amean"))
dev.off()

#### plot sa ####
png(paste0(save_dir, "plot_SA.png"), width=1000, height=750)
plotSA(fit, main= paste0(my_id," Final model: Mean-variance trend"))
dev.off()


#Finally, apply the empirical Bayes’ step to get our differential expression statistics and p-values.
#### make contrasts and ebayes ####
contrasts <- makeContrasts(Asthmatic_at_baseline - Healthy_control, levels = design)
#Finally, apply the empirical Bayes’ step to get our differential expression statistics and p-values.
contr.fit <- eBayes(contrasts.fit(fit, contrasts), trend = TRUE)
contr.fit

topTable(contr.fit)
#decidetests and save venn diagram
results <- decideTests(contr.fit, method="global", lfc=1)
summary(results)

vennDiagram(results)

####save results
write.table(summary(results), "regulation.csv")

result1 = topTable(contr.fit, n=Inf, p.value = 0.05)


#result2 = topTable(contr.fit, n=Inf, coef = "AsthBaselinevsHealthy", p.value = 0.05)
#result3 = topTable(contr.fit, n=Inf, coef = "smokervsHealthy", p.value = 0.05)
#result4 = topTable(contr.fit, n=Inf, coef = "AsthFloventvsHealthy", p.value = 0.05)
#result5 = topTable(contr.fit, n=Inf, coef = "AsthFloventvsAsthBaseline", p.value = 0.05)
#result6 = topTable(contr.fit, n=Inf, coef = "AsthPlacebovsAsthBaseline", p.value = 0.05)
#result7 = topTable(contr.fit, n=Inf, coef = "AsthFloventvsAsthPlacebo", p.value = 0.05)


###### plot hists ##
library(RColorBrewer)
png(paste0(save_dir, "AsthBaselinevsHealthy.png"), width=1000, height=750)
hist(result1$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "hist_AsthPlacebovsHealthy", xlab = "p-values")
dev.off()


#png(paste0(save_dir, "hist_AsthBaselinevsHealthy.png"), width=1000, height=750)
#hist(result2$P.Value, col = brewer.pal(3, name = "Set2")[1],
  #   main = "AsthBaselinevsHealthy", xlab = "p-values")
#dev.off()

#png(paste0(save_dir, "hist_smokervsHealthy.png"), width=1000, height=750)
#hist(result3$P.Value, col = brewer.pal(3, name = "Set2")[1],
#     main = "smokervsHealthy", xlab = "p-values")
#dev.off()

#png(paste0(save_dir, "hist_AsthFloventvsHealthy.png"), width=1000, height=750)
#hist(result4$P.Value, col = brewer.pal(3, name = "Set2")[1],
  #   main = "AsthFloventvsHealthy", xlab = "p-values")
#dev.off()

#png(paste0(save_dir, "hist_AsthFloventvsAsthBaseline.png"), width=1000, height=750)
#hist(result5$P.Value, col = brewer.pal(3, name = "Set2")[1],
#     main = "AsthFloventvsAsthBaseline", xlab = "p-values")
#dev.off()

#png(paste0(save_dir, "hist_AsthPlacebovsAsthBaseline.png"), width=1000, height=750)
#hist(result6$P.Value, col = brewer.pal(3, name = "Set2")[1],
#     main = "AsthPlacebovsAsthBaseline", xlab = "p-values")
#dev.off()

#png(paste0(save_dir, "hist_AsthFloventvsAsthPlacebo.png"), width=1000, height=750)
#hist(result7$P.Value, col = brewer.pal(3, name = "Set2")[1],
#     main = "AsthFloventvsAsthPlacebo", xlab = "p-values")
#dev.off()


############
## plot mds ## 
png(paste0(save_dir, "MDS_fit.png"), width=1000, height=750)
plotMD(fit, coef=1, main= paste0(my_id, " Mean-Difference Plot of fit, coef=1"))
dev.off()


png(paste0(save_dir, "SA_fit.png"), width=1000, height=750)
plotSA(fit,main="Residual standard deviation versus average log expression for fit")
dev.off()


## plot MDs after ebayes
png(paste0(save_dir, "MDS_fit_after_ebayes.png"), width=1000, height=750)
plotMD(contr.fit, coef=1,main="Mean-Difference Plot of fit2 (after ebayes), coef=1")
abline(0,0,col="blue")
dev.off()

#png(paste0(save_dir, "MDS_fit_after_ebayes_coeff2.png"), width=1000, height=750)
#plotMD(contr.fit, coef=2,main="Mean-Difference Plot of fit2 (after ebayes), coef=2")
#abline(0,0,col="blue")
#dev.off()


#png(paste0(save_dir, "MDS_fit_after_ebayes_coeff3.png"), width=1000, height=750)
#plotMD(contr.fit, coef=3,main="Mean-Difference Plot of fit2 (after ebayes), coef=3")
#abline(0,0,col="blue")
#dev.off()


#png(paste0(save_dir, "MDS_fit_after_ebayes_coeff4.png"), width=1000, height=750)
#plotMD(contr.fit, coef=4,main="Mean-Difference Plot of fit2 (after ebayes), coef=4")
#abline(0,0,col="blue")
#dev.off()

#png(paste0(save_dir, "MDS_fit_after_ebayes_coeff5.png"), width=1000, height=750)
#plotMD(contr.fit, coef=5,main="Mean-Difference Plot of fit2 (after ebayes), coef=5")
#abline(0,0,col="blue")
#dev.off()

#png(paste0(save_dir, "MDS_fit_after_ebayes_coeff6.png"), width=1000, height=750)
#plotMD(contr.fit, coef=6,main="Mean-Difference Plot of fit2 (after ebayes), coef=6")
#abline(0,0,col="blue")
#dev.off()

#png(paste0(save_dir, "MDS_fit_after_ebayes_coeff7.png"), width=1000, height=750)
#plotMD(contr.fit, coef=7,main="Mean-Difference Plot of fit2 (after ebayes), coef=7")
#abline(0,0,col="blue")
#dev.off()




###### venn diagram #####





################summary #
print("Sumary expression")
summary(results)



######## save results table #############
full_results <- topTable(contr.fit, number=Inf)
full_results<- full_results[,-1]
full_results <- tibble::rownames_to_column(full_results,"ID")

#save data table
library(readr)
filter(full_results, adj.P.Val < 0.05, abs(F) > 1) %>%
  write_csv(paste0(save_dir,my_id,path="_filtered_DE_results.csv"))


genes_interest <- c("SLC22A1", "SLC22A4", "SLC22A5")
filter(full_results, SYMBOL %in% genes_interest) %>%
  write_csv(paste0(save_dir,my_id,path="_filtered_genes_results.csv"))



####plot genes expression#######################
#### plot genes expression ############
library(dplyr); library(tidyr); library(ggplot2); library(stringr)

edata_final <- exprs(manfiltered)
as_data <- as.data.frame(edata_final, SKIP=0)

annot_select <- fData(manfiltered)
as_data$probe <- rownames(as_data)
as_data$genes <- annot_select$SYMBOL

data_long <- gather(as_data, IDATfile, log_fold, GSM98141:GSM98258, factor_key=FALSE)
data_long


names(pdata)[names(pdata) == "geo_accession"] <- "IDATfile"
targetinfo <- pdata %>%
  dplyr::select(c("IDATfile", "sample type:ch1"))
data_long <- merge(data_long, targetinfo, by  = "IDATfile")
names(data_long)[names(data_long) == "sample type:ch1"] <- "Group"


data_long <- data_long %>%
  filter(stringr::str_detect(Group, "Asthmatic at baseline|Healthy"))

table(data_long$Group)



#save data
write.csv(data_long, "data_long_asthmatic_healthy.csv")


#function for plotting
plot_gene <- function(data, title){
  ggplot(data, aes(x= Group, log_fold)) +
    geom_boxplot(outlier.shape = NA, color= "black",fill= c("gray60", "gray33")) +
    geom_jitter(width=0.08, height = 0.5, size= 2.0, color= "black") +
    ggtitle(data$genes, subtitle = my_id) +# We'll make this a jitter plot
    ylab("Expression") + 
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size=24), 
          axis.title = element_text(size = 20), # font size of axis
          axis.text.x = element_text(size=16), #font size of x ticks
          axis.text.y = element_text(size=14),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+ # font size of y ticks
    scale_y_continuous(breaks = round(seq(min(data$log_fold), max(data$log_fold), by = 0.5),1))
  ggsave(title, width = 15,height = 25, units="cm")
}


data_long <- data_long %>%
  mutate(across("Group", str_replace, "Asthmatic_at_baseline", "Asthmatic")) %>%
  mutate(across("Group", str_replace, "Healthy_control", "Control"))

###filter for genes of interest
SLC22A1_table <- data_long %>%
  filter(genes == "SLC22A1")

SLC22A4_table <- data_long %>%
  filter(genes == "SLC22A4")

SLC22A5_table <- data_long %>%
  filter(genes == "SLC22A5")

## save plots ##
plot_gene(SLC22A1_table, paste0(save_dir,"SLC22A1_expression.png"))
plot_gene(SLC22A4_table, paste0(save_dir,"SLC22A4_expression.png"))
plot_gene(SLC22A5_table, paste0(save_dir,"SLC22A5_expression.png"))






#############################plot spider####

data_longer <- data_long %>%
  group_by(Group, genes) %>%
  summarise(mean_log_fold = mean(log_fold))


library(reshape2)

SLC22A <- data_longer %>%
  filter(stringr::str_detect(genes, "SLC22A"))

#plot spider#
max_min <- dcast(SLC22A, Group ~ genes)
max_min[max_min < 9] <- 9
max_min[nrow(max_min) + 1,] <- 0
max_min <- max_min[-1,-1]
#change rownames

SLC22A <- dcast(SLC22A, Group ~ genes)
SLC22A <- SLC22A[-1]
df <- rbind(max_min, SLC22A)

rownames(df)<- c("Max", "Min", "Adjacent_non_tumor", "Adenocarcinoma")


library(fmsb)

png(paste0(save_dir, "radar_chart.png"), width=1000, height=750)
radarchart(df)
dev.off()



create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}


png(paste0(save_dir, "radarchart_SLC22A.png"), width=1000, height=750)
op <- par(mar = c(1,2,2,2))  
create_beautiful_radarchart(df, caxislabels = c(0,6,7,8,9), color = c("#00AFBB", "#E7B800"))
legend(x="bottom", legend = rownames(df[-c(1,2),]), horiz = TRUE, bty = "n", pch=20, col = c("#00AFBB", "#E7B800"), text.col = "black", cex=1, pt.cex= 1.5)
par(op)
dev.off()




#############################
######################plot all SLC22 genes in boxplot
#table
SLC22 <- data_long %>%
  filter(stringr::str_detect(genes, "SLC22A"))



plot_all_genes <- function(data, title){
  ggplot(data, aes(x= genes ,log_fold, fill=Group))  +
    geom_boxplot(outlier.shape = NA, color= "black", position="dodge") +
    ggtitle(my_id, subtitle="cohort: Asthmatic samples(n=42) and Healthy control (n=28)") +
    ylab("Log2(counts+1)") + 
    xlab("") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size=20), 
          axis.title = element_text(size = 18), # font size of axis
          axis.text.x = element_text(size=16, angle = 90), #font size of x ticks
          axis.text.y = element_text(size=12))+ # font size of y ticks
    scale_y_continuous(breaks = round(seq(min(data$log_fold), max(data$log_fold), by = 0.5),1))+
    scale_fill_manual(values=c("gray60", "gray33"))
  ggsave(title, width = 35,height = 20, units="cm")
}


#plot and save
plot_all_genes(SLC22, paste0(save_dir,"SLC22_expression.png"))



stats = full_results[,c("SYMBOL","logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
write.csv(stats, "all_stats_proteins.csv")








#################new volcano
library(EnhancedVolcano)
#The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.


keyvals.colour <- ifelse(
  full_results$logFC < -1, 'royalblue',
  ifelse(full_results$logFC > 1, 'red',"grey"))
names(keyvals.colour)[keyvals.colour == 'red'] <- 'Up-regulated'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'Not-Significant'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'Down-regulated'



png(paste0(save_dir, "enhanced_volcano.png"), width=700, height=900)
EnhancedVolcano(full_results,
                lab = "",
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = my_id,
                subtitle = "Differential expression: Asthmatic samples(n=42) and Healthy control (n=28)",
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                colCustom = keyvals.colour,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'),
                col = c('grey', 'grey', 'grey', 'red3'),
                legendPosition = 'bottom',
                legendLabSize = 16,
                legendIconSize = 5.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'full',
                borderWidth = 1.0,
                borderColour = 'black')
dev.off()



