library(illuminaio)
library(limma)

## read files
idatfiles <- dir(pattern="idat")
bgxfile <- dir(pattern="bgx")
x <- read.idat(idatfiles, bgxfile)

## extract geo expression, fData, eData
my_id <- "GSE75037"
Sys.setenv(VROOM_CONNECTION_SIZE = 256000000)
expr <- getGEO(my_id)[[1]]

#extract data from expr
sampleInfo <- pData(expr)
edata <- exprs(expr) #to compare the edata matrix to the raw data that I'll obtain later
annot <- fData(expr) #annotation data

#add sample info to the dataframe
x$targets <- sampleInfo

data2 <- neqc(x)
##The intensities vary from about 5 to 16 on the log2 scale:
boxplot(log2(data2$E),range=0, ylab="log2 intensity")

design <- model.matrix(~0 + x$targets$source_name_ch1)
colnames(design) <- c("cancer", "non_malignant")

# Apply the intensity values to lmFit. 
fit <- lmFit(data2, design)

# Create a contrast matrix. In this example, all combinations of contrasts can be set up as below. 
contrast.matrix <- makeContrasts("cancer-non_malignant",  levels=design)

# Apply this contrast matrix to the modeled data and compute statistics for the data.
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


summary(decideTests(fit2, method="global"))
topTable(fit2, coef=1)


####SA
plotSA(fit2, main="Final model: Mean-variance trend")
