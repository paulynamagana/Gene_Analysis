#use to untar file from GEO
#utils::untar("./GSE67472_RAW.tar", exdir = "./GSE67472")  #untar

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

# specify the path on computer where the folder that contains the CEL-files is located
celpath = "./GSE67472/"


###############################create cdf file
#HGU133Plus2cdf<- make.cdf.package("GPL16311_HGU133Plus2_Hs_ENTREZG.cdf.gz", species = "Homo sapiens",
 #                                 compress = TRUE)

#library(HGU133Plus2cdf)

##open cdf package from the folder data

#load cel files
raw_affy<- ReadAffy(celfile.path = celpath)


library(hgu133plus2.db)
library(annotate)




### boxplot raw intensities
png(paste0(save_dir, "boxplot_rawdata.png"), width=1200, height=750)
boxplot(log2(raw_affy@assayData$exprs),range=0, ylab="log2 intensity",
        main= paste0(my_id, " Boxplot of log2-intensiyties for RAW data"))
dev.off()




#normalise usinf rma from affy
project.norm <- affy::rma(raw_affy)


#extract IDs to match 
eset_norm <- exprs(project.norm)

ID <- rownames(eset_norm) #extract rownames
head(ID) #check probes
gene.symbols <- getSYMBOL(ID, "hgu133plus2") #match the ID to the library
gene.symbols #print and check

#retrieve annotation from library (hgu133plus2.db)
columns(hgu133plus2.db)
annot <- AnnotationDbi::select(hgu133plus2.db, keys = ID, columns = c("SYMBOL", "GENENAME" ))

#copy gene symbols into the data
fData(project.norm) <- data.frame(Symbol=gene.symbols)



##plot boxplot again in log2 scale
png(paste0(save_dir, "boxplots_probes_after_norm.png"), width=1200, height=750)
boxplot(log2(eset_norm),range=0,ylab="log2 intensity", main = paste0(my_id, " Boxplot probes log2 after normlisation"))
dev.off()












########## DATA FROM GEO ##############################################
## extract geo expression, fData, eData
Sys.setenv(VROOM_CONNECTION_SIZE = 256000000)

###extracting the file from GEO
my_id <- "GSE67472" #id from GEO
#extract data
expr <- getGEO(my_id)[[1]]

# print data info from GEO
print("Abstract from project")
print(expr@experimentData@abstract)
# get the abstract
exp_data <- expr@experimentData@abstract
# save#
fileConn<-file(paste0(my_id, "_EXPDATA.txt"))
writeLines(exp_data, fileConn)
close(fileConn)


#extract the data for samples from GEO
geo_pdata <- pData(expr)
colnames(eset_norm) <- geo_pdata$geo_accession


fdata<- fData(project.norm)
fdata$probe <- rownames(fdata)

###match names
fdata <- fdata[match(fdata$probe, rownames(eset_norm)),]


###### later
eset_norm <- as.data.frame(eset_norm) #convert into a dataframe to add column for probes and symbols
##add column probes that is the same as rownames



pd <- new("AnnotatedDataFrame", data = geo_pdata) #conver into datafrme
ann <- new("AnnotatedDataFrame", data = fdata) #convert annot that has been filtered into dataframee


## check namnes for rows and columns
featureNames(ann)
rownames(eset_norm)

identical(featureNames(ann), rownames(eset_norm))

colnames(eset_norm)
identical(featureNames(pd), colnames(eset_norm)) #TRUE


norm_data <- new("ExpressionSet", exprs=as.matrix(eset_norm), phenoData = pd, featureData = ann)




####filtering based on intensity
medians <- rowMedians(Biobase::exprs(norm_data))

png(paste0(save_dir, "hist_medians_normdata.png"), width=1200, height=850)
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = paste0(my_id, " Histogram of the median intensities"), 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
dev.off()

######## define threshold
man_threshold <- 2

png(paste0(save_dir, "hist_medians_cutoff_normdata.png"), width=1200, height=850)
hist_res <- hist(medians, 100, col = "cornsilk", freq = FALSE, 
                 main = paste0(my_id, " Histogram of the median intensities"), 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
abline(v = man_threshold, col = "coral4", lwd = 2)
dev.off()


##########filter intensity
no_of_samples <- table(geo_pdata["disease state:ch1"])
no_of_samples 

samples_cutoff <- 20

idx_man_threshold <- apply(Biobase::exprs(norm_data), 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)

##removing multiple mapping with probes
norm_data_filtered <- subset(norm_data, idx_man_threshold)

#remove multiple mapping
annot_filtered <- subset(annot, !is.na(SYMBOL)) #get rid of empty
anno_grouped <- group_by(annot_filtered, PROBEID)
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized) 

#look at the probes with more than 1
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
head(anno_filtered)

probe_stats <- anno_filtered 
nrow(probe_stats)


ids_to_exlude <- (featureNames(norm_data_filtered) %in% probe_stats$PROBEID)
table(ids_to_exlude)


norm_data_filtered_final <- subset(norm_data_filtered, !ids_to_exlude)
validObject(norm_data_filtered_final)


# summarise across genes by mean - use ID
data <- exprs(norm_data_filtered_final)
symbols <- fData(norm_data_filtered_final)


library(dplyr)

data<- as.data.frame(data)
dim(data) #53140   105

data$SYMBOL <- symbols$Symbol

#average expression for symbols
A <- limma::avereps(data, data$SYMBOL)
dim(A) #20825   106
A<- A[,-106] #get rid of te symbol column
A<- as.data.frame(A)

A <- mutate_all(A, function(x) as.numeric(as.character(x)))


######################create filtered expressionset again
pd <- new("AnnotatedDataFrame", data = geo_pdata) #conver into datafrme

ann <- fData(norm_data_filtered_final)
ann <- ann[which(rownames(ann) %in% rownames(A)),] #subset to get only those in A

ann <- new("AnnotatedDataFrame", data = ann)


## check namnes for rows and columns
featureNames(ann)
rownames(A)

identical(featureNames(ann), rownames(A))

colnames(A)
identical(featureNames(pd), colnames(A)) #TRUE

norm_data_filtered_final_mean <- new("ExpressionSet", exprs=as.matrix(A), phenoData = pd, featureData = ann)





##############################################################
#####################linear models ##################################
pdata <- pData(norm_data_filtered_final_mean) #extact pdata from final eset


library(tidyverse)

#extract variables
names(pdata)[names(pdata) == "disease state:ch1"] <- "disease"
sample <- factor(pdata$disease)
sample
#create levels
levels(sample) <- c("healthy", "asthma")
levels(sample)

#create designs
design <- model.matrix(~0 +sample)
colnames(design)

#change colnames
colnames(design) <- c("healthy", "asthma")
design

#create contrasts
contrasts <- makeContrasts(asthma - healthy,
                           levels=design)
contrasts


#calculate weigths
aw <- arrayWeights(norm_data_filtered_final_mean, design)
fit <- lmFit(norm_data_filtered_final_mean, design, weights= aw)
fit


##########################visualize###


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








###################
#Finally, apply the empirical Bayes’ step to get our differential expression statistics and p-values.
#### make contrasts and ebayes ####
contrasts <- makeContrasts(asthma - healthy, levels = design)
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

results = topTable(contr.fit, n=Inf, p.value = 0.05)


###### plot hists ##
library(RColorBrewer)
png(paste0(save_dir, "AstvsHealthy.png"), width=1000, height=750)
hist(results$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = paste0(my_id, " hist_vsHealthy"), xlab = "p-values")
dev.off()






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






######## save results table #############
full_results <- topTable(contr.fit, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

#save data table
library(readr)
write_csv(full_results, paste0(save_dir,my_id,path="_filtered_DE_results.csv"))


genes_interest <- c("SLC22A1", "SLC22A4", "SLC22A5")
filter(full_results, Symbol %in% genes_interest) %>%
  write_csv(paste0(save_dir,my_id,path="_filtered_genes_results.csv"))



####plot genes expression#######################
#### plot genes expression ############
library(dplyr); library(tidyr); library(ggplot2); library(stringr)

edata_final <- exprs(norm_data_filtered_final_mean)
as_data <- as.data.frame(edata_final, SKIP=0)

annot_select <- fData(norm_data_filtered_final_mean)
as_data$probe <- rownames(as_data)
as_data$genes <- annot_select$Symbol

data_long <- gather(as_data, IDATfile, log_fold, GSM1647628:GSM1647732, factor_key=FALSE)
data_long


names(pdata)[names(pdata) == "geo_accession"] <- "IDATfile"
targetinfo <- pdata %>% dplyr::select(c("IDATfile", "disease"))
data_long <- merge(data_long, targetinfo, by  = "IDATfile")

#samples
table(data_long$disease)

#save data
write.csv(data_long, "data_long_asthmatic_healthy.csv")






#function for plotting
plot_gene <- function(data, title){
  ggplot(data, aes(x= disease, log_fold)) +
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





#############################
######################plot all SLC22 genes in boxplot
#table
SLC22 <- data_long %>%
  filter(stringr::str_detect(genes, "SLC22A"))



plot_all_genes <- function(data, title){
  ggplot(data, aes(x= genes ,log_fold, fill=disease))  +
    geom_boxplot(outlier.shape = NA, color= "black", position="dodge") +
    ggtitle(my_id, subtitle="cohort: mild to moderate asthma (n=62) and control subjects without asthma (n=43)") +
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



stats = full_results[,c("Symbol","logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
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
                subtitle = "Differential expression: mild to moderate asthma (n=62) and healthy controls (n=43)",
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







