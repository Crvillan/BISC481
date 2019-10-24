## Cecilia Villanueva Code

##Question 4- Using MLR to predict binding sites

## Install packages
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
# DNAshapeR
BiocManager::install("DNAshapeR")
# Caret
install.packages("caret")
install.packages ('data.table')

## Initialization
library(DNAshapeR)
library(caret)
workingPath <- "/Users/ceciliavillanueva/Desktop/BISC481/gcPBM/"

## Predict DNA shapes
fn_fasta <- paste0(workingPath, "Mad.txt.fa")
pred <- getShape(fn_fasta)

## Encode feature vectors
featureType <- c("1-mer", "1-shape")
featureVector <- encodeSeqShape(fn_fasta, pred, featureType)
head(featureVector)

## Build MLR model by using Caret
# Data preparation
fn_exp <- paste0(workingPath, "Mad.txt")
exp_data <- read.table(fn_exp)
df <- data.frame(affinity=exp_data$V2, featureVector)

# Arguments setting for Caret
trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)

# Prediction without L2-regularized
model <- train (affinity~ ., data = df, trControl=trainControl, 
                method = "lm", preProcess=NULL)
summary(model)

# Prediction with L2-regularized
model2 <- train(affinity~., data = df, trControl=trainControl, 
                method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model2
result <- model2$results$Rsquared[1]
result

###To generate the plot comparing the two methods we use the following code which takes into account that x-axis/values = "1-mer" and y-axis/values = "1-mer" + "1-shape"
## Install and initialize packages
install.packages("ggplot2")
install.packages("grid")
library(ggplot2)
library(grid)

## Theme
my.theme <- theme(
  plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm"),
  axis.text.x = element_text(colour="black", size=12),
  axis.text.y = element_text(colour="black", size=12),
  axis.title.x = element_text(colour="black", size=12),
  axis.title.y = element_text(colour="black", size=12),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(colour ="black"),
  axis.ticks = element_line(colour = "black")
)

#DATA TO INPUT
#Mad: 1-mer only R=0.854836, 1-mer+1-shape  R= 0.8639159
#Max: 1-mer only R=0.854836, 1-mer+1-shape R=0.8644361
#Myc: 1-mer only R=0.854836, 1-mer+1-shape R=0.854836
## Data preparation
one-mer_Rsquared <- c(0.854836,0.854836,0.854836)
one-mer+shape_Rsquared <- c(0.8639159,0.8644361, 0.854836)

## Ploting
ggplot() +
  geom_point(aes(x = one-mer_Rsquared, y = one-mer+shape_Rsquared), color = "red", size=1) +
  geom_abline(slope=1) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
  coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  my.theme

############################## New Problem
##Question 5- High-throughput in vivo data analysis
## Install packages
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
# DNAshapeR
BiocManager::install("DNAshapeR")
# Caret
install.packages("caret")

## Install packages
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("AnnotationHub")
BiocManager::install("GenomicRanges")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

## Initialization
library(DNAshapeR)
library(AnnotationHub)
library(BSgenome.Mmusculus.UCSC.mm10)

## Data retreival
seqLength <- 50 #500
sampleSize <- 2000 #42045
workingPath <- "/Users/ceciliavillanueva/Desktop/BISC481-master/CTCF/"

# Bound (ChIP-seq)
ah <- AnnotationHub()
unique(ah$dataprovider)
query(ah, "H3K4me3")
ctcfPeaks <- ah[["AH46100"]]
seqlevelsStyle(ctcfPeaks) <- "UCSC"
getFasta( GR = sample(ctcfPeaks, sampleSize), BSgenome = Mmusculus, width = seqLength, 
          filename = paste0(workingPath, "bound_500.fa"))

# Unbound (random regions w/o overlapping)
chrName <- names(Mmusculus)[1:22]
chrLength <- seqlengths(Mmusculus)[1:22]
randomGr <- GRanges()

while ( length(randomGr) < sampleSize )  {
  tmpChrName <- sample(chrName, 1)
  tmpChrLength <- chrLength[tmpChrName]
  tmpGRanges <- GRanges( seqnames = tmpChrName, strand = "+",
                         IRanges(start = sample(1:tmpChrLength,1), width = seqLength) )
  if( length(findOverlaps(ctcfPeaks, tmpGRanges)) == 0 ){
    randomGr <- c( randomGr, tmpGRanges)
    print(length(randomGr))
  }else{
    print(paste(length(randomGr), "overlap"))
  }
}
randomGr

# Overlap checking
findOverlaps(ctcfPeaks, randomGr)

# Fasta file generation
getFasta(randomGr, Mmusculus, width = seqLength, 
         filename = paste0(workingPath, "unbound_500.fa"))

# Extract sample sequences for BOUND data--> used this command instead because it worked better for me
fn <- file.path("/Users", "ceciliavillanueva", "Downloads", "BISC481-master", "CTCF", "bound_500.fa")

# Predict DNA shapes
pred <- getShape(fn)

# Generate ensemble plots
plotShape(pred$MGW)
plotShape(pred$ProT)
plotShape(pred$Roll)
plotShape(pred$Helt)

# Extract sample sequences for UNBOUND data --> used this command instead because it worked better for me
fn2 <- file.path("/Users", "ceciliavillanueva", "Downloads", "BISC481-master", "CTCF", "unbound_500.fa")

# Predict DNA shapes
pred2 <- getShape(fn2)

# Generate ensemble plots
plotShape(pred2$MGW)
plotShape(pred2$ProT)
plotShape(pred2$Roll)
plotShape(pred2$Helt)




## Question 6- Logistic regression on ChIP-seq data
##Install packages
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
# DNAshapeR
BiocManager::install("DNAshapeR")
# Biostrings
BiocManager::install("Biostrings")
# Others
install.packages("caret")
install.packages("e1071")
install.packages("ROCR")

## Initialization
library(DNAshapeR)
library(caret)
library(ROCR)
library(Biostrings)
workingPath <- "/Users/ceciliavillanueva/Desktop/BISC481/CTCF/"

## Generate data for the classifcation (assign Y to bound and N to non-bound)
# bound
boundFasta <- readDNAStringSet(paste0(workingPath, "bound_30.fa"))
sequences <- paste(boundFasta)
boundTxt <- data.frame(seq=sequences, isBound="Y")

# non-bound
nonboundFasta <- readDNAStringSet(paste0(workingPath, "unbound_30.fa"))
sequences <- paste(nonboundFasta)
nonboundTxt <- data.frame(seq=sequences, isBound="N")

# merge two datasets
writeXStringSet( c(boundFasta, nonboundFasta), paste0(workingPath, "ctcf.fa"))
exp_data <- rbind(boundTxt, nonboundTxt)


## DNAshapeR prediction
pred <- getShape(paste0(workingPath, "ctcf.fa"))


## Encode feature vectors
featureType <- c("1-mer", "1-shape")
featureVector <- encodeSeqShape(paste0(workingPath, "ctcf.fa"), pred, featureType)
df <- data.frame(isBound = exp_data$isBound, featureVector)


## Logistic regression
# Set parameters for Caret
trainControl <- trainControl(method = "cv", number = 10, 
                             savePredictions = TRUE, classProbs = TRUE)
# Perform prediction
model <- train(isBound~ ., data = df, trControl = trainControl,
               method = "glm", family = binomial, metric ="ROC")
summary(model)

## Plot AUROC
prediction <- prediction( model$pred$Y, model$pred$obs )
performance <- performance( prediction, "tpr", "fpr" )
plot(performance)

## Caluculate AUROC
auc <- performance(prediction, "auc")
auc <- unlist(slot(auc, "y.values"))
auc

####To generate the second plot for only "1-mer", we repeat the same code starting with the new encoding of feature vectors
## Encode feature vectors - "1-mer" only
featureType <- c("1-mer")
featureVector <- encodeSeqShape(paste0(workingPath, "ctcf.fa"), pred, featureType)
df <- data.frame(isBound = exp_data$isBound, featureVector)


## Logistic regression
# Set parameters for Caret
trainControl <- trainControl(method = "cv", number = 10, 
                             savePredictions = TRUE, classProbs = TRUE)
# Perform prediction
model <- train(isBound~ ., data = df, trControl = trainControl,
               method = "glm", family = binomial, metric ="ROC")
summary(model)

## Plot AUROC
prediction <- prediction( model$pred$Y, model$pred$obs )
performance <- performance( prediction, "tpr", "fpr" )
plot(performance)

## Caluculate AUROC
auc <- performance(prediction, "auc")
auc <- unlist(slot(auc, "y.values"))
auc
