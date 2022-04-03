
## Dependencies

library(httr)
library(jsonlite)
library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(rmarkdown)


## Import data
SLC22A4url <- "https://www.proteinatlas.org/ENSG00000197208.json"
SLC22A4_json <- jsonlite::fromJSON(SLC22A4url, "text")

## Arrange data
SLC22A4RAW <- SLC22A4_json[82:136]
SLC22A4RAW <- as.data.frame(SLC22A4RAW)
rownames(SLC22A4RAW) <- "nTPM"

## Separate the cell lines from string
colnames(SLC22A4RAW) <- gsub("Tissue.RNA...", "", names(SLC22A4RAW), fixed = TRUE)
colnames(SLC22A4RAW) <- gsub("..nTPM.", "", names(SLC22A4RAW), fixed = TRUE)
colnames(SLC22A4RAW) <- gsub(".", " ", names(SLC22A4RAW), fixed = TRUE)

#AS DATAFRAME
SLC22A4 <- as.data.frame(t(SLC22A4RAW))
SLC22A4$tissue <- rownames(SLC22A4)

SLC22A4 <-transform(SLC22A4, nTPM=as.numeric(nTPM))
SLC22A4$transporter <- "SLC22A4"
rownames(SLC22A4) <- NULL

#plot
ggplot(SLC22A4, aes(x=tissue, y=nTPM))+
  geom_bar(stat='identity')





###############SLC22A1
## Import data
SLC22A1url <- "https://www.proteinatlas.org/ENSG00000175003.json"
SLC22A1_json <- jsonlite::fromJSON(SLC22A1url, "text")

## Arrange data
SLC22A1RAW <- SLC22A1_json[82:136]
SLC22A1RAW <- as.data.frame(SLC22A1RAW)
rownames(SLC22A1RAW) <- "nTPM"

## Separate the cell lines from string
colnames(SLC22A1RAW) <- gsub("Tissue.RNA...", "", names(SLC22A1RAW), fixed = TRUE)
colnames(SLC22A1RAW) <- gsub("..nTPM.", "", names(SLC22A1RAW), fixed = TRUE)
colnames(SLC22A1RAW) <- gsub(".", " ", names(SLC22A1RAW), fixed = TRUE)

#AS DATAFRAME
SLC22A1 <- as.data.frame(t(SLC22A1RAW))
SLC22A1$tissue <- rownames(SLC22A1)
SLC22A1 <-transform(SLC22A1, nTPM=as.numeric(nTPM))
SLC22A1$transporter <- "SLC22A1"
rownames(SLC22A1) <- NULL



###############SLC22A5
## Import data
SLC22A5url <- "https://www.proteinatlas.org/ENSG00000197375.json"
SLC22A5_json <- jsonlite::fromJSON(SLC22A5url, "text")

## Arrange data
SLC22A5RAW <- SLC22A5_json[82:136]
SLC22A5RAW <- as.data.frame(SLC22A5RAW)
rownames(SLC22A5RAW) <- "nTPM"

## Separate the cell lines from string
colnames(SLC22A5RAW) <- gsub("Tissue.RNA...", "", names(SLC22A5RAW), fixed = TRUE)
colnames(SLC22A5RAW) <- gsub("..nTPM.", "", names(SLC22A5RAW), fixed = TRUE)
colnames(SLC22A5RAW) <- gsub(".", " ", names(SLC22A5RAW), fixed = TRUE)

#AS DATAFRAME
SLC22A5 <- as.data.frame(t(SLC22A5RAW))
SLC22A5$tissue <- rownames(SLC22A5)
SLC22A5 <-transform(SLC22A5, nTPM=as.numeric(nTPM))
SLC22A5$transporter <- "SLC22A5"
rownames(SLC22A5) <- NULL


###############SLC22A2
## Import data
SLC22A2url <- "https://www.proteinatlas.org/ENSG00000112499.json"
SLC22A2_json <- jsonlite::fromJSON(SLC22A2url, "text")

## Arrange data
SLC22A2RAW <- SLC22A2_json[82:136]
SLC22A2RAW <- as.data.frame(SLC22A2RAW)
rownames(SLC22A2RAW) <- "nTPM"

## Separate the cell lines from string
colnames(SLC22A2RAW) <- gsub("Tissue.RNA...", "", names(SLC22A2RAW), fixed = TRUE)
colnames(SLC22A2RAW) <- gsub("..nTPM.", "", names(SLC22A2RAW), fixed = TRUE)
colnames(SLC22A2RAW) <- gsub(".", " ", names(SLC22A2RAW), fixed = TRUE)

#AS DATAFRAME
SLC22A2 <- as.data.frame(t(SLC22A2RAW))
SLC22A2$tissue <- rownames(SLC22A2)
SLC22A2 <-transform(SLC22A2, nTPM=as.numeric(nTPM))
SLC22A2$transporter <- "SLC22A2"
rownames(SLC22A2) <- NULL


###############SLC22A3
## Import data
SLC22A3url <- "https://www.proteinatlas.org/ENSG00000146477.json"
SLC22A3_json <- jsonlite::fromJSON(SLC22A3url, "text")

## Arrange data
SLC22A3RAW <- SLC22A3_json[82:136]
SLC22A3RAW <- as.data.frame(SLC22A3RAW)
rownames(SLC22A3RAW) <- "nTPM"

## Separate the cell lines from string
colnames(SLC22A3RAW) <- gsub("Tissue.RNA...", "", names(SLC22A3RAW), fixed = TRUE)
colnames(SLC22A3RAW) <- gsub("..nTPM.", "", names(SLC22A3RAW), fixed = TRUE)
colnames(SLC22A3RAW) <- gsub(".", " ", names(SLC22A3RAW), fixed = TRUE)

#AS DATAFRAME
SLC22A3 <- as.data.frame(t(SLC22A3RAW))
SLC22A3$tissue <- rownames(SLC22A3)
SLC22A3 <-transform(SLC22A3, nTPM=as.numeric(nTPM))
SLC22A3$transporter <- "SLC22A3"
rownames(SLC22A3) <- NULL




###############SLC22A6
## Import data
SLC22A6url <- "https://www.proteinatlas.org/ENSG00000197901.json"
SLC22A6_json <- jsonlite::fromJSON(SLC22A6url, "text")

## Arrange data
SLC22A6RAW <- SLC22A6_json[82:136]
SLC22A6RAW <- as.data.frame(SLC22A6RAW)
rownames(SLC22A6RAW) <- "nTPM"

## Separate the cell lines from string
colnames(SLC22A6RAW) <- gsub("Tissue.RNA...", "", names(SLC22A6RAW), fixed = TRUE)
colnames(SLC22A6RAW) <- gsub("..nTPM.", "", names(SLC22A6RAW), fixed = TRUE)
colnames(SLC22A6RAW) <- gsub(".", " ", names(SLC22A6RAW), fixed = TRUE)

#AS DATAFRAME
SLC22A6 <- as.data.frame(t(SLC22A6RAW))
SLC22A6$tissue <- rownames(SLC22A6)
SLC22A6 <-transform(SLC22A6, nTPM=as.numeric(nTPM))
SLC22A6$transporter <- "SLC22A6"
rownames(SLC22A6) <- NULL











###############paste all
#paste
transporters <- rbind(SLC22A1, SLC22A4, SLC22A5, SLC22A2, SLC22A3)

transporters$log10nTPM <- log10(transporters$nTPM+1)

write.csv(transporters, file = "./data/transportershumantlas.csv")




#plot
ggplot(transporters, aes(x=tissue, y=log10nTPM, fill=transporter))+
  geom_bar(stat="identity", position="dodge") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle=90)) 


pdf(file = "./PDFs_outcome/expression_proteinatlas.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 12) 

p <- ggplot(transporters, aes(y=tissue, x=transporter, fill=log10nTPM)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="black") +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "SLC22 Abundance expression in human") +
  scale_y_discrete(limits = rev(levels(as.factor(transporters$tissue)))) +
  labs(fill='log TPM')  +
  theme_classic()
p

dev.off()














