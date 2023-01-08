dir.create("./plots")
########1 SET FOR SAVING PLOTS #################

save_dir <- "./plots/"


#gc() # free unused memory

#untar
utils::untar("./GSE63383_RAW.tar", exdir = "./GSE63383")  #untar

###################1.1 LIBRARIES ##################

library(dplyr)
library(limma)
library(tinytex)
library(edgeR)
library(Biobase)
library(GEOquery)
require(RColorBrewer)
require(PCAtools)
library(DESeq2)