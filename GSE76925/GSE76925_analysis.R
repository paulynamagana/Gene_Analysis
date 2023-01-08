#utils::untar("GSE76925_RAW.tar", exdir = "./GSE76925")  #untar


#save plots
dir.create("./plots")
save_dir <- "./plots/"


## libraries
library(dplyr)
library(limma)
library(tinytex)
library(edgeR)
library(Biobase)
library(GEOquery)
require(RColorBrewer)
require(PCAtools)


baseDir <- './'
my_id <- "GSE76925"
bgxfile <- "./GSE76925/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz"
file <- "GSE76925_non-normalized.txt.gz"





########### 1.2 LOAD DATA #################

# DATA FROM GEO #
## extract geo expression, fData, eData
Sys.setenv(VROOM_CONNECTION_SIZE = 256000000)
###extracting the file from GEO
expr <- getGEO(my_id)[[1]]

#sampleInfo <- pData(expr)
GEO_edata <- exprs(expr) #to compare the edata matrix to the raw data that I'll obtain later
#annot <- fData(expr) #annotation data

#boxplot(GEO_edata)

# print data info from GEO
print("Abstract from project")
print(expr@experimentData@abstract)

# get the abstract
exp_data <- expr@experimentData@abstract
# save#
fileConn<-file("GSE76925_EXPerimentDATA.txt")
writeLines(exp_data, fileConn)
close(fileConn)






###################
#read in data
x <- read.table(paste0(baseDir, file),
                header = TRUE, sep = '\t', stringsAsFactors = FALSE, skip = 0)
#extract pvalues
detectionpvalues <- x[,grep('Detection.Pval', colnames(x))]
rownames(detectionpvalues) <- rownames(x)
#remove pvalues from dataframe
x <- x[,-grep('Detection.Pval', colnames(x))]

#convert to matrix
x <- data.matrix(x[,1:ncol(x)])







######################################## prepare targetinfo
sampleInfo <- pData(expr)

targetinfo <- data.frame(matrix(NA, nrow = 151, ncol=0))
rownames(targetinfo) <- sampleInfo$geo_accession
targetinfo = sampleInfo[, c("description.1", "characteristics_ch1", "age:ch1", "Sex:ch1", "tissue:ch1", "packyears:ch1", "ID:ch1", "copd:ch1", "bmi:ch1")]

#columns and rows names
targetinfo$description.1 #make sure this column matches the colnames in x
colnames(x)
#change columns "sample" for the GSM sample number
colnames(x) <- rownames(targetinfo)





################### read annotation
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
dim(project) # 47249   151


write.table(dim(project), file = paste0(save_dir, "dimensions_project.csv"))

# plot densities raw data and save
png(paste0(save_dir, "densities_rawdata.png"), width=1200, height=850)
plotDensities(project, legend=FALSE, main= "Densities raw data")
dev.off() # a function call to save the file


### boxplot raw intensities
png(paste0(save_dir, "boxplot_rawdata.png"), width=1200, height=750)
boxplot(log2(project$E),range=0, ylab="log2 intensity",
        main= "Boxplot of log2-intensiyties for RAW data")
dev.off()


# prepare data for PCA plot #
exp_raw <- log2(project$E)
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

data_raw <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                       Sample = targetinfo$`copd:ch1`,
                       Batch = targetinfo$`Sex:ch1`)

# ggplot PCA and save plot#
png(paste0(save_dir, "pca_rawdata.png"), width=1200, height=850)
ggplot(data_raw, aes(PC1, PC2)) +
  geom_point(aes(colour = Batch, shape = Sample), size=4) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()







################################
#################3 Background correction / normalisation ##################
project.bgcorrect.norm <- neqc(project, offset = 16)


# save hisgoram from data after neqc#
png(paste0(save_dir, "hist_data_after_neqc.png"), width=800, height=750)
hist(project.bgcorrect.norm$E, main = "histogram of NEQC data")
dev.off()


# plot densities after neqd data and save
png(paste0(save_dir, "densities_after_neqc.png"), width=1200, height=850)
plotDensities(project.bgcorrect.norm, legend=FALSE, main= "Densities after neqc data")
dev.off() # a function call to save the file






####################
########## filter out control probes, those with no symbol, and those that failed #
annot <- annot[which(annot$Probe_Id %in% rownames(project.bgcorrect.norm)),]
project.bgcorrect.norm <- project.bgcorrect.norm[which(rownames(project.bgcorrect.norm) %in% annot$Probe_Id),]

annot <- annot[match(rownames(project.bgcorrect.norm), annot$Probe_Id),]

project.bgcorrect.norm@.Data[[3]] <- annot
project.bgcorrect.norm$genes <- annot

#print dimensions
print("dimensions of prohect after neqc")
dim(project.bgcorrect.norm) #47249   151

write.table(dim(project.bgcorrect.norm), file = paste0(save_dir, "dimensions_project_bgcorrect_norm.csv"))

# save densities from data after neqc#
png(paste0(save_dir, "densities_data_after_neqc.png"), width=800, height=750)
plotDensities(project.bgcorrect.norm, legend=FALSE, main = "Densities data after neqc")
dev.off()





##############################
###
## plots for raw data
png(paste0(save_dir, "boxplots_probes_rawdata.png"), width=1200, height=750)
par(mfrow = c(2,1))

boxplot(log2(project$E[project$genes$Source == "ILMN_Controls", ]),
        range = 0, las = 2, xlab = "", ylab = expression(log[2](intensity)),
        main = "Control probes, RAW data")

boxplot(log2(project$E),
        range = 0, las = 2, xlab = "", ylab = expression(log[2](intensity)),
        main = "Regular probes, RAW data")
dev.off()
###


### plots for neqc data##
png(paste0(save_dir, "boxplots_probes_after_neqc.png"), width=1200, height=750)
par(mfrow = c(2,1))

boxplot(project.bgcorrect.norm$E[project$genes$Source == "ILMN_Controls", ],
        range = 0, las = 2, xlab = "", ylab = expression(log[2](intensity)),
        main = "Control probes, NEQC data")

boxplot(project.bgcorrect.norm$E, range = 0, ylab = expression(log[2](intensity)),
        las = 2, xlab = "", main = "Regular probes, NEQC normalized")
dev.off()
####

# print table with the different ILMN probes in the data
table(project.bgcorrect.norm$genes$Source)
##



## pheatmaps ##
library(pheatmap)

corMatrix <- cor(project.bgcorrect.norm$E,use="c")

png(paste0(save_dir, "heatmap_neqc.png"), width=1000, height=750)
pheatmap(corMatrix)     
dev.off()
####
###
groups <- targetinfo[, c("Sex:ch1", "copd:ch1")]


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
                     Sample = targetinfo$`copd:ch1`,
                     Batch = targetinfo$`Sex:ch1`)

png(paste0(save_dir, "PCA_neqc_data.png"), width=1000, height=750)
ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Batch, shape = Sample)) +
  ggtitle("PCA plot of the log-transformed neqc expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))
dev.off()



####3.1 dealing with batch effects ##################
# PLOT MDS ####
png(paste0(save_dir, "MDS_neqc_data.png"), width=1000, height=750)
plotMDS(project.bgcorrect.norm$E, labels=targetinfo$`copd:ch1`, main= "MDS data after neqc")
dev.off()
###



#### avg signal from neqc data ####
aveSignal <- rowMeans(project.bgcorrect.norm$E)
png(paste0(save_dir, "avg_signal_neqc_data.png"), width=1000, height=750)
boxplot(aveSignal, main = "Average signal data after neqc")
dev.off()






#####################################
######### filtering based on probe annotation ##############

Control <- project.bgcorrect.norm$genes$Source=="ILMN_Controls"
NoSymbol <- project.bgcorrect.norm$genes$Symbol == ""
isexpr <- rowSums(project.bgcorrect.norm$other$detection <= 0.05) >= 3

table(Control) # FALSE  TRUE 
#                 47156    93 

table(NoSymbol) # FALSE  TRUE 
#                 43984  3265 

table(isexpr)#    FALSE  TRUE 
#                 12835 34414


## filter out probes
project.bgcorrect.norm.filt <- project.bgcorrect.norm[!Control & !NoSymbol & isexpr, ]
print("Dimensions of data after neqc and filtering")
dim(project.bgcorrect.norm.filt)     #31665   151

write.table(dim(project.bgcorrect.norm.filt), file = paste0(save_dir, "dimensions_project_bgcorrect_norm_filt.csv"))


## densities from neqc and filtered data ##
png(paste0(save_dir, "densities_data_after_neqc_filtering.png"), width=1000, height=750)
plotDensities(project.bgcorrect.norm.filt, legend=FALSE, main = "Densities data after bg correction")
dev.off()


### IQR plot ###
IQR <- apply(project.bgcorrect.norm.filt$E, 1, IQR, na.rm = TRUE)
topVar <- order(IQR, decreasing = TRUE)[1:500]
d <- dist(t(project.bgcorrect.norm.filt$E[topVar, ]))
png(paste0(save_dir, "cluster_dendogram_after_neqc_and_filtering.png"), width=1000, height=750)
plot(hclust(d))
dev.off()

##
## heatmap
png(paste0(save_dir, "heatmap_after_neqc_and_filtering.png"), width=1000, height=750)
heatmap(project.bgcorrect.norm.filt$E[topVar, ])
dev.off()


## avg singal after neqc and filt
aveSignal <- rowMeans(project.bgcorrect.norm.filt$E)
png(paste0(save_dir, "avg_signal_neqc_filt.png"), width=1000, height=750)
boxplot(aveSignal)
dev.off()


########################################## ----------------
### prepare for differential expression
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
dim(project.bgcorrect.norm.filt.mean) # 23644   151

write.table(dim(project.bgcorrect.norm.filt.mean), file = paste0(save_dir, "dimensions_project_bgcorrect_norm_filt_mean.csv"))

#plot histograms samples
for (i in 1:6)
{
  name = paste(save_dir, "QC_neqc_filt_mean_histogram",i,".jpg",sep="")
  jpeg(name)
  hist(project.bgcorrect.norm.filt.mean$E[,i],lwd=2, ylab='Density',xlab='Log2 intensities',main=project.bgcorrect.norm.filt.mean$targets$IDATfile[i])
  dev.off()
}







################################------------------------------
############### 4. DIFFERENTIAL EXPRESSION ######################
targets <- factor(targetinfo$`copd:ch1`, levels= c("cont", "case")) 
design<- model.matrix(~0 +targets)
colnames(design) <- c("Control", "COPD")
design

#calculate arrays
aw <- arrayWeights(project.bgcorrect.norm.filt.mean, design)
fit <- lmFit(project.bgcorrect.norm.filt.mean, design, weights = aw)
fit

#calculate array weight
aw <- arrayWeights(project.bgcorrect.norm.filt.mean, design)
fit <- lmFit(project.bgcorrect.norm.filt.mean, design, weights= aw)


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




#### make contrasts and ebayes ####
contrasts <- makeContrasts(COPD- Control, levels = design)
#Finally, apply the empirical Bayes’ step to get our differential expression statistics and p-values.
contr.fit <- eBayes(contrasts.fit(fit, contrasts), trend = TRUE)
topTable(contr.fit)
#decidetests and save venn diagram
results <- decideTests(contr.fit, method= "global", lfc=1)



##########################################################
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


png(paste0(save_dir, "SA_fit_after_ebayes.png"), width=1000, height=750)
plotSA(contr.fit,main="Residual standard deviation versus average log expression for fit2 (after ebayes)")
abline(0,0,col="blue")
dev.off()








## venn diagram #####
png(paste0(save_dir, "venn_diagram.png"), width=1000, height=750)
vennDiagram(results)
dev.off()

#summary #
print("Sumary expression")
summary(results)







######## save results table #############
full_results <- topTable(contr.fit, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

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
  ggtitle("Volcano Plot") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()



## change according to your needs
p_cutoff <- 0.001
fc_cutoff <- 2



#volcanoplot with cutoffs
png(paste0(save_dir, "volcanoplot_cutoff.png"), width=1000, height=750)
full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()
ggtitle("Volcano Plot") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#
## plot SA ######################
png(paste0(save_dir, "plot_SA_contrfit.png"), width=1000, height=750)
plotSA(contr.fit, main="Final model: Mean-variance trend")
dev.off()

## plot MD #######################
png(paste0(save_dir, "plot_MD.png"), width=1000, height=750)
plotMD(contr.fit, status=results)
dev.off()


print("AdjPval > 0.05")
length(which(full_results$adj.P.Val < 0.05))





##########-------------
#### plot genes expression ############
library(dplyr); library(tidyr); library(ggplot2)

as_data <- as.data.frame(project.bgcorrect.norm.filt.mean$E, SKIP=0 )
as_data$genes <- rownames(as_data)


data_long <- gather(as_data, IDATfile, log_fold, GSM2040792:GSM2040942, factor_key=FALSE)
data_long

targetinfo$IDATfile <- rownames(targetinfo)
targetinfo <- targetinfo[,c("IDATfile", "copd:ch1", "description.1")]
data_long <- merge(data_long, targetinfo, by  = "IDATfile") 
data_long <- data_long %>% rename("Group" = "copd:ch1")









############################# new

plot_gene <- function(data, title){
  ggplot(data, aes(x= Group, log_fold)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width=0.1, height = 0.5, size= 2.0, color= "darkblue") +
    ggtitle(data$genes) +# We'll make this a jitter plot
    ylab("Expression") + 
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size=24), 
          axis.title = element_text(size = 20), # font size of axis
          axis.text.x = element_text(size=16), #font size of x ticks
          axis.text.y = element_text(size=14))+ # font size of y ticks
    scale_y_continuous(breaks = round(seq(min(data$log_fold), max(data$log_fold), by = 0.5),1))
  ggsave(title)
}


write.csv(data_long, file = "data_long.csv")



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


png(paste0(save_dir, "radarchart_SLC22A.png"), width=600, height=600)
op <- par(mar = c(1,2,2,2))  
create_beautiful_radarchart(df, caxislabels = c(0,6,7,8,9), color = c("#00AFBB", "#E7B800"))
legend(x="bottom", legend = rownames(df[-c(1,2),]), horiz = TRUE, bty = "n", pch=20, col = c("#00AFBB", "#E7B800"), text.col = "black", cex=1, pt.cex= 1.5)
par(op)
dev.off()













