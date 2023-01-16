
#-----------------------------building targetsfile from pdata
#targets <- sampleInfo[,which(colnames(sampleInfo) %in% c("title", "geo_accession", "source_name_ch1"))]
#rownames(targets) <- NULL
#targets$title <-sub("^","X",targets$title)
#write.csv(targets, file ="targets.csv")


########1 SET FOR SAVING PLOTS #################

save_dir <- "./data_plots/"
dir.create(save_dir)

###################1.1 LIBRARIES ##################

library(dplyr)
library(limma)
library(tinytex)
library(edgeR)
library(Biobase)
library(GEOquery)
require(RColorBrewer)
require(PCAtools)

my_id <- "GSE75037"
bgxfile <- "GPL6884_HumanWG-6_V3_0_R0_11282955_A.bgx.gz"
targetsfile <- 'targets.csv'


############ 1.2 LOAD DATA #################

# DATA FROM GEO #
## extract geo expression, fData, eData
Sys.setenv(VROOM_CONNECTION_SIZE = 256000000)
###extracting the file from GEO
expr <- getGEO(my_id)[[1]]

sampleInfo <- pData(expr)
edata <- exprs(expr) #to compare the edata matrix to the raw data that I'll obtain later
#annot <- fData(expr) #annotation data

#plot series from GEO
png(paste0(save_dir, "densities_from_GEOseries.png"), width=800, height=600)
boxplot(edata, main = paste0(my_id, " boxplot from expr data from GEO"))
dev.off()


# print data info from GEO
print("Abstract from project")
print(expr@experimentData@abstract)


# get the abstract
exp_data <- expr@experimentData@abstract
# save#
fileConn<-file(paste0(my_id,"_data_abstract.txt"))
writeLines(exp_data, fileConn)
close(fileConn)





################################## read in the raw data 
baseDir <- './'
bgxfile <- "GPL6884_HumanWG-6_V3_0_R0_11282955_A.bgx.gz"
targetsfile <- 'targets.csv'
file <- "GSE75037_non_normalized.txt.gz"

#read in data
x <- read.table(paste0(baseDir, file),
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





######################################## load the target info
targetinfo <- read.csv(targetsfile)



rownames(targetinfo) <- targetinfo$geo_accession
#columns and rows names
rownames(x) <- probes
colnames(x) <- colnames(edata)
colnames(detectionpvalues) <- colnames(x)


x <- x[,match(rownames(targetinfo), colnames(x))]

if (!all(colnames(x) == rownames(targetinfo))) stop('Target info is not aligned to expression data - they must be in the same order')



############################################### read in annotation

annot <- illuminaio::readBGX(bgxfile)$probes

annot <- annot[,which(colnames(annot) %in% c('Source','Transcript','ILMN_Gene','RefSeq_ID',
                                             'Entrez_Gene_ID','Symbol','Protein_Product','Probe_Id','Probe_Type',
                                             'Probe_Start','Chromosome','Probe_Chr_Orientation','Probe_Coordinates',
                                             'Cytoband', 'Definition', 'Ontology_Component', 'Ontology_Process',
                                             'Ontology_Function', 'Synonyms'))]

annot <- annot[which(annot$Probe_Id %in% rownames(x)),]

annot <- annot[match(rownames(x), annot$Probe_Id),]







###############1.3 create ElistRaw for neqc #####################

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
project$other$detection <- detectionpvalues







###########################2 QC PLOTS FOR ALL RAW DATA #############

### Dimensions of the project
print("dimensions of project")
dim(project) # 48803   166

write.table(dim(project), file = paste0(save_dir, "dimensions_project.csv"))

# plot densities raw data and save
png(paste0(save_dir, "densities_rawdata.png"), width=1200, height=850)
plotDensities(project, legend=FALSE, main = paste0(my_id ," Densities raw data"))
dev.off() # a function call to save the file


### boxplot raw intensities
png(paste0(save_dir, "boxplot_rawdata.png"), width=1200, height=750)
boxplot(log2(project$E),range=0, ylab="log2 intensity",
        main= paste0(my_id, " Boxplot of log2-intensiyties for RAW data"))
dev.off()





# prepare data for PCA plot #
exp_raw <- log2(project$E)
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

data_raw <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Sample = sampleInfo$`histology:ch1`,
                     Batch = sampleInfo$characteristics_ch1.2)


# ggplot PCA and save plot#
png(paste0(save_dir, "pca_rawdata.png"), width=1200, height=850)
ggplot(data_raw, aes(PC1, PC2)) +
  geom_point(aes(colour = Batch, shape = Sample), size=4) +
  ggtitle(my_id, subtitle = "PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5, size=25, face="bold"), plot.subtitle = element_text(size=20))
dev.off()








#################3 Background correction / normalisation ##################
project.bgcorrect.norm <- neqc(project, offset = 16)



# save hisgoram from data after neqc#
png(paste0(save_dir, "hist_data_after_neqc.png"), width=800, height=750)
hist(project.bgcorrect.norm$E, main = paste0(my_id, " histogram of NEQC data"))
dev.off()







### filter out control probes, those with no symbol, and those that failed #
annot <- annot[which(annot$Probe_Id %in% rownames(project.bgcorrect.norm)),]

project.bgcorrect.norm <- project.bgcorrect.norm[which(rownames(project.bgcorrect.norm) %in% annot$Probe_Id),]
annot <- annot[match(rownames(project.bgcorrect.norm), annot$Probe_Id),]

project.bgcorrect.norm@.Data[[3]] <- annot
project.bgcorrect.norm$genes <- annot

#print dimensions
print("dimensions of project after neqc")
dim(project.bgcorrect.norm) #48803   166

write.table(dim(project.bgcorrect.norm), file = paste0(save_dir, "dimensions_project_bgcorrect_norm.csv"))

# save densities from data after neqc#
png(paste0(save_dir, "densities_data_after_neqc.png"), width=800, height=750)
plotDensities(project.bgcorrect.norm, legend=FALSE, main = paste0(my_id, " Densities data after neqc"))
dev.off()



####
###
## plots for raw data
png(paste0(save_dir, "boxplots_probes_rawdata.png"), width=1200, height=750)
par(mfrow = c(2,1))

boxplot(log2(project$E[project$genes$Source == "ILMN_Controls", ]),
        range = 0, las = 2, xlab = "", ylab = expression(log[2](intensity)),
        main = paste0(my_id, " Control probes, RAW data"))

boxplot(log2(project$E),
        range = 0, las = 2, xlab = "", ylab = expression(log[2](intensity)),
        main = paste0(my_id, " Regular probes, RAW data"))
dev.off()
###




### plots for neqc data##
png(paste0(save_dir, "boxplots_probes_after_neqc.png"), width=1200, height=750)
par(mfrow = c(2,1))

boxplot(project.bgcorrect.norm$E[project$genes$Source == "ILMN_Controls", ],
        range = 0, las = 2, xlab = "", ylab = expression(log[2](intensity)),
        main = paste0(my_id, " Control probes, NEQC data"))

boxplot(project.bgcorrect.norm$E, range = 0, ylab = expression(log[2](intensity)),
        las = 2, xlab = "", main = paste0(my_id, " Regular probes, NEQC normalized"))
dev.off()
####
###


# print table with the different ILMN probes in the data
table(project.bgcorrect.norm$genes$Source)
##


##plot boxplot again in log2 scale
png(paste0(save_dir, "boxplots_probes_after_neqc_log2.png"), width=1200, height=750)
boxplot(log2(project.bgcorrect.norm$E),range=0,ylab="log2 intensity", main = paste0(my_id, " Boxplot probes log2 after neqc"))
dev.off()



###
## pheatmaps ##
library(pheatmap)

corMatrix <- cor(project.bgcorrect.norm$E,use="c")

png(paste0(save_dir, "heatmap_neqc.png"), width=1000, height=750)
pheatmap(corMatrix)     
dev.off()



####
###
groups <- targetinfo[4:4]

png(paste0(save_dir, "heatmap_neqc_groups.png"), width=1000, height=750)
pheatmap(corMatrix,
         annotation_col=groups)
dev.off()
###




##### PCA neqc data #################
exp_raw <- log2(project.bgcorrect.norm$E)
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Sample = sampleInfo$`histology:ch1`,
                     Batch = sampleInfo$characteristics_ch1.2)

png(paste0(save_dir, "PCA_neqc_data.png"), width=1000, height=750)
ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Batch, shape = Sample)) +
  ggtitle(paste0(my_id, " PCA plot of the log-transformed neqc expression data")) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()





####3.1 dealing with batch effects ##################
# PLOT MDS ####
png(paste0(save_dir, "MDS_neqc_data.png"), width=1000, height=750)
plotMDS(project.bgcorrect.norm$E, labels=targetinfo$source_name_ch1,
        main= paste0(my_id," MDS data after neqc"))
dev.off()
###





#### avg signal from neqc data ####
aveSignal <- rowMeans(project.bgcorrect.norm$E)
png(paste0(save_dir, "avg_signal_neqc_data.png"), width=400, height=600)
boxplot(aveSignal, main = paste0(my_id, " Average signal data after neqc"))
dev.off()





############ filtering based on probe annotation ##############

Control <- project.bgcorrect.norm$genes$Source=="ILMN_Controls"
NoSymbol <- project.bgcorrect.norm$genes$Symbol == ""
isexpr <- rowSums(project.bgcorrect.norm$other$detection <= 0.05) >= 3

table(Control) # FALSE  TRUE 
#                 48799     4 

table(NoSymbol) # FALSE  TRUE 
#                 35966 12837

table(isexpr)#    FALSE  TRUE 
#                 16432  32371 


## filter out probes
project.bgcorrect.norm.filt <- project.bgcorrect.norm[!Control & !NoSymbol & isexpr, ]
print("Dimensions of data after neqc and filtering")
dim(project.bgcorrect.norm.filt)     #27097      166

write.table(dim(project.bgcorrect.norm.filt), file = paste0(save_dir, "dimensions_project_bgcorrect_norm_filt.csv"))




## densities from neqc and filtered data ##
png(paste0(save_dir, "densities_data_after_neqc_filtering.png"), width=1000, height=750)
plotDensities(project.bgcorrect.norm.filt, legend=FALSE, main = paste0(my_id," Densities data after bg correction"))
dev.off()



### IQR plot ###
IQR <- apply(project.bgcorrect.norm.filt$E, 1, IQR, na.rm = TRUE)
topVar <- order(IQR, decreasing = TRUE)[1:500]
d <- dist(t(project.bgcorrect.norm.filt$E[topVar, ]))
png(paste0(save_dir, "cluster_dendogram_after_neqc_and_filtering.png"), width=1000, height=750)
plot(hclust(d), labels = project.bgcorrect.norm.filt$targets$source_name_ch1)
dev.off()

##
## heatmap
png(paste0(save_dir, "heatmap_after_neqc_and_filtering.png"), width=1000, height=750)
heatmap(project.bgcorrect.norm.filt$E[topVar, ])
dev.off()


## avg singal after neqc and filt
aveSignal <- rowMeans(project.bgcorrect.norm.filt$E)
png(paste0(save_dir, "avg_signal_neqc_filt.png"), width=400, height=600)
boxplot(aveSignal, main= paste0(my_id, " Avg Signal after normalization and filt"))
dev.off()




### prepare for differential expression
####
####
# remove annotation columns we no longer need
project.bgcorrect.norm.filt$genes <- project.bgcorrect.norm.filt$genes[,c(
  'Probe_Id',
  'Definition','Ontology_Component','Ontology_Process','Ontology_Function',
  'Chromosome','Probe_Coordinates','Cytoband','Probe_Chr_Orientation',
  'RefSeq_ID','Entrez_Gene_ID','Symbol')]


head(project.bgcorrect.norm.filt$genes)




#### summarise across genes by mean ####
# ID is used to identify the replicates
project.bgcorrect.norm.filt.mean <- avereps(project.bgcorrect.norm.filt,
                                            ID = project.bgcorrect.norm.filt$genes$Symbol)
print("Dimensions after neqc, filt and mean")
dim(project.bgcorrect.norm.filt.mean) # 20008    166

write.table(dim(project.bgcorrect.norm.filt.mean), file = paste0(save_dir, "dimensions_project_bgcorrect_norm_filt_mean.csv"))




#plot histograms samples
for (i in 1:6)
{
  name = paste(save_dir, "QC_neqc_filt_mean_histogram",i,".jpg",sep="")
  jpeg(name)
  hist(project.bgcorrect.norm.filt.mean$E[,i],lwd=2, ylab='Density',xlab='Log2 intensities', main= paste0(my_id, " ", project.bgcorrect.norm.filt.mean$targets$IDATfile[i]))
  dev.off()
}








# 4. DIFFERENTIAL EXPRESSION ######################
design<- model.matrix(~0 + targetinfo$source_name_ch1)
colnames(design)

## the column names are a bit ugly, so we will rename
colnames(design) <- c("Lung_cancer","non_malignant")

#calculate arrays
aw <- arrayWeights(project.bgcorrect.norm.filt.mean, design)
fit <- lmFit(project.bgcorrect.norm.filt.mean, design, weights= aw)

## plot weights ##
png(paste0(save_dir, "weights.png"), width=1000, height=750)
barplot(aw, xlab="Array", ylab="Weight", col="white", las=2)
abline(h=1, lwd=1, lty=2)
dev.off()



## histogram fit Amean
png(paste0(save_dir, "histogram_fit_mean.png"), width=800, height=600)
hist(fit$Amean, main = paste0(my_id, " histogram fit Amean"))
dev.off()

#### plot sa ####
png(paste0(save_dir, "plot_SA.png"), width=1000, height=750)
plotSA(fit, main= paste0(my_id, " Final model: Mean-variance trend"))
dev.off()





#### make contrasts and ebayes ####
contrasts <- makeContrasts(Lung_cancer - non_malignant, levels = design)
#Finally, apply the empirical Bayes’ step to get our differential expression statistics and p-values.
contr.fit <- eBayes(contrasts.fit(fit, contrasts), trend = TRUE)
topTable(contr.fit)
#decidetests and save venn diagram
results <- decideTests(contr.fit, method= "global", lfc=1)




## plot mds ## after ebayes
png(paste0(save_dir, "MDS_fit.png"), width=1000, height=750)
plotMD(fit, coef=1,main= paste0(my_id, " Mean-Difference Plot of fit, coef=1"))
dev.off()



png(paste0(save_dir, "SA_fit.png"), width=1000, height=750)
plotSA(fit, main= paste0(my_id, " Residual standard deviation versus average log expression for fit"))
dev.off()



## plot MDs after ebayes
png(paste0(save_dir, "MDS_fit_after_ebayes.png"), width=1000, height=750)
plotMD(contr.fit, coef=1, main= paste0(my_id, " Mean-Difference Plot of fit2 (after ebayes), coef=1"))
abline(0,0,col="blue")
dev.off()



png(paste0(save_dir, "SA_fit_after_ebayes.png"), width=1000, height=750)
plotSA(contr.fit,main= paste0(my_id, " Residual standard deviation versus average log expression for fit2 (after ebayes)"))
abline(0,0,col="blue")
dev.off()








## venn diagram #####
png(paste0(save_dir, "venn_diagram.png"), width=1000, height=750)
vennDiagram(results)
dev.off()

#summary #
print("Summary expression")
summary(results)

write.table(summary(results), "regulation.csv")




######## save results table #############
full_results <- topTable(contr.fit, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

#write csv full results
write.csv(full_results, paste0(save_dir ,"FULL_RESULTS.csv"))


#save data table
library(readr)
filter(full_results, adj.P.Val < 0.05, abs(logFC) > 1) %>%
  write_csv(paste0(save_dir,my_id,path="_filtered_DE_results.csv"))

genes_interest <- c("SLC22A1", "SLC22A4", "SLC22A5")
filter(full_results, Symbol %in% genes_interest) %>%
  write_csv(paste0(save_dir,my_id,path="_filtered_SLC22sgenes_results.csv"))



### volcano plot ####
png(paste0(save_dir, "volcanoplot.png"), width=1000, height=750)
ggplot(full_results, aes(x = logFC, y=B)) + geom_point() +
  ggtitle(paste0(my_id," Volcano Plot")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))# font size of y ticks)
dev.off()


###
###
###
## change according to your needs
p_cutoff <- 0.001
fc_cutoff <- 2

#volcanoplot with cutoffs
png(paste0(save_dir, "volcanoplot_cutoff.png"), width=600, height=800)
full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point() +
  ggtitle(paste0(my_id, " Volcano Plot")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))# font size of y ticks)
dev.off()





## plot SA ######################
png(paste0(save_dir, "plot_SA_contrfit.png"), width=1000, height=750)
plotSA(contr.fit, main= paste0(my_id," Final model: Mean-variance trend"))
dev.off()



## plot MD #######################
png(paste0(save_dir, "plot_MD.png"), width=1000, height=750)
plotMD(contr.fit, status=results)
dev.off()


print("AdjPval > 0.05")
length(which(full_results$adj.P.Val < 0.05)) #13752


#### plot genes expression ############
library(dplyr); library(tidyr); library(ggplot2); library(stringr)

#manipulate dataset
as_data <- as.data.frame(project.bgcorrect.norm.filt.mean$E, SKIP=0 )
as_data$genes <- rownames(as_data)
data_long <- gather(as_data, source_name_ch1, log_fold, GSM1941120:GSM1941285, factor_key=FALSE)
targetinfo <- targetinfo[3:4]
data_long <- merge(data_long, targetinfo, by.x = "source_name_ch1", by.y="geo_accession") 

colnames(data_long) <- c("Sample", "genes", "log_fold", "Group")

#save data
write.csv(data_long, paste0(save_dir, "data_long.csv"))


###################plot per gene

plot_gene <- function(data, title){
  ggplot(data, aes(x= Group, log_fold)) +
    geom_boxplot(outlier.shape = NA, color= "black", fill= c("gray60", "gray33")) +
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
  filter(genes == "SLC22A1") %>%
  write.csv("SLC22A1.csv")

SLC22A4_table <- data_long %>%
  filter(genes == "SLC22A4")  %>%
  write.csv("SLC22A4.csv")

SLC22A5_table <- data_long %>%
  filter(genes == "SLC22A5")  %>%
  write.csv("SLC22A5.csv")



## save plots ##
plot_gene(SLC22A1_table, paste0(save_dir,"SLC22A1_expression.png"))
plot_gene(SLC22A4_table, paste0(save_dir,"SLC22A4_expression.png"))
plot_gene(SLC22A5_table, paste0(save_dir,"SLC22A5_expression.png"))




#############################
###################### plot all SLC22 genes in boxplot
#table
SLC22 <- data_long %>%
  filter(stringr::str_detect(genes, "SLC22A")) %>%
  write.csv("SLC22A.csv")

plot_all_genes <- function(data, title){
  ggplot(data, aes(x= genes ,log_fold, fill=Group))  +
    geom_boxplot(outlier.shape = NA, color= "black", position="dodge") +
    ggtitle(my_id, subtitle="cohort: Adenocarcinomas(n=83) and matched adjacent non-malignant lung (n=83)") +
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
write.csv(stats, paste0(save_dir, "all_stats_proteins.csv"))




#################new volcano
library(EnhancedVolcano)
#The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.

keyvals.colour <- ifelse(
  full_results$logFC < -2, 'royalblue',
  ifelse(full_results$logFC > 2, 'red',"grey"))
names(keyvals.colour)[keyvals.colour == 'red'] <- 'Up-regulated'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'Not-Significant'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'Down-regulated'



png(paste0(save_dir, "enhanced_volcano.png"), width=700, height=900)
EnhancedVolcano(full_results,
                lab = "",
                x = 'logFC',
                y = 'adj.P.Val',
                FCcutoff = 2.0,
                title = my_id,
                subtitle = "Differential expression Adenocarcinomas(n=83) and matched adjacent non-malignant lung (n=83)",
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


