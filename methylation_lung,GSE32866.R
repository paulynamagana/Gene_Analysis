library(rmarkdown) # paged table
library(GEOquery) #geo query access
library(knitr) #tables
library(dplyr)
library(devtools)
library(ggplot2)
library(limma)
library(Glimma)
library(edgeR)
library(ggrepel)
Sys.setenv(VROOM_CONNECTION_SIZE = 25600000)

#Genome-scale DNA methylation profiling of lung adenocarcinoma

my_id <- "GSE32866"
gse <- getGEO(my_id)