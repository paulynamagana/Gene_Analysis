#use to untar file from GEO
utils::untar("./GSE67472_RAW.tar", exdir = "./GSE67472")  #untar

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
library(BiocManager)
install("makecdfenv")
library(makecdfenv)

# specify the path on computer where the folder that contains the CEL-files is located
celpath = "./GSE67472/"


######################################
library(HGU133Plus2cdf)
library(affy)



HGU133Plus2cdf<- make.cdf.package("GPL16311_HGU133Plus2_Hs_ENTREZG.cdf.gz", species = "Homo sapiens",
                                  compress = TRUE)



raw <- read.celfiles(
  filenames = list.files('GSE67472/', pattern = '*CEL', full.names = TRUE),
  sampleNames = gsub('HGU133Plus2cdf', '', list.files('GSE67472/', pattern = '*CEL')))
dim(exprs(raw))

#open the file
names <- names(HGU133Plus2cdf)
dim(names)

library(hgu133plus2.db)
library(annotate)
library(limma)


gns2 <- select(hgu133plus2.db, keys= names,
              c("ENTREZID","SYMBOL","GENENAME"))











