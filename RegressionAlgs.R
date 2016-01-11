# clean workspace
rm(list = ls())

# names of datafiles
dataFiles = c("~/Dropbox/AMARETTOdownload/TCGA_BLCA_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData",
              "~/Dropbox/AMARETTOdownload/TCGA_BRCA_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData",
              "~/Dropbox/AMARETTOdownload/TCGA_COADREAD_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData",
              "~/Dropbox/AMARETTOdownload/TCGA_GBM_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData",
              "~/Dropbox/AMARETTOdownload/TCGA_HNSC_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData",
              "~/Dropbox/AMARETTOdownload/TCGA_KIRC_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData",
              "~/Dropbox/AMARETTOdownload/TCGA_LAML_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData",
              "~/Dropbox/AMARETTOdownload/TCGA_LUAD_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData",
              "~/Dropbox/AMARETTOdownload/TCGA_LUSC_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData",
              "~/Dropbox/AMARETTOdownload/TCGA_OV_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData",
              "~/Dropbox/AMARETTOdownload/TCGA_UCEC_ProcessedData_MA20141206_CNV20141017_METMethylMix2015.RData")

# obtain Survival data
load("~/Dropbox/AMARETTOdownload/Survival.RData")
namesSurvival <- data.frame(Survival[,3])
row.names(namesSurvival) <- row.names(Survival)

# choose omic data type
# 1: Methylation data
# 2: Copy number variations data
# 3: Single nucleotide polymorphism (SNP) mutation data
OMIC_DATA_TYPE <- 1

# choose cancer type
# 1.  bladder cancer
# 2.  breast cancer
# 3.  Glioblastoma
# 4.  Ovarian cancer
# 5.  Colorectal cancer
# 6.  Acute myeloid leukemia
# 7.  Head and Neck squamous carcinoma
# 8.  Kidney clear cell carcinoma
# 9.  Lung adencarcinoma 
# 10. Lung squamous carcinoma
# 11. Endometrical cancer
cancerType <- 2

# PARSE DATA SET
# ----------
load(dataFiles[cancerType])
genomicData <- data.frame(ProcessedData[OMIC_DATA_TYPE])
cancerData <- 0
event <- 0
time <- 0
flag <- TRUE
patientId <- list()
j <- 1
for (i in 1:((dim(genomicData))[2])) {
  id <- (colnames(genomicData))[i]
  id <- strsplit(id, "\\.")
  id = paste(id[[1]][2:(length(id[[1]])-1)], collapse = '-')
  
  if(id %in% row.names(namesSurvival)) {
    patientId[length(patientId) + 1] = id
    if(!flag) {
      cancerData <- cbind(cancerData, genomicData[,i]) 
      event[j] <- (Survival[id,])[[3]]
      time[j] <- (Survival[id,])[[2]]
      j <- j + 1
    } else {
      cancerData <- genomicData[,i]
      event <- (Survival[id,])[[3]]
      time <- (Survival[id,])[[2]]
      flag <- FALSE
      j <- 2
    }
  }
}

# store parsed data
NUM_PATIENTS <- length(event)
NUM_GENES <- dim(genomicData)[1]
time <- strtoi(as.vector(time))
event <- strtoi(as.vector(event))
cancerData <- t(cancerData)
cancerData <- round(cancerData, digits = 1)
cancerData <- data.frame(cancerData)
colnames(cancerData) <- row.names(genomicData)
rownames(cancerData) <- patientId[1:NUM_PATIENTS]

# parse gene names
for(j in 1:(dim(cancerData)[2])) {
  if(grepl("`", colnames(cancerData)[j])) {
    colnames(cancerData)[j] <- substr(colnames(cancerData)[j], 2, nchar(colnames(cancerData)[j])-1)
  }
}

# use this for plotting Kaplan-Meier Curves
KM <- cancerData
meanExp <- mean(cancerData[,'ANKRD52'])
for(i in 1:length(cancerData[,'ANKRD52'])){
  if(cancerData[i,'ANKRD52'] >= meanExp){
    KM[i, 'ANKRD52'] <- 0 
  } else{
    KM[i, 'ANKRD52'] <- 1
  }
}

# kaplan meier curve
plot(survfit(Surv(time, event)~ANKRD52,data=KM), conf.int=FALSE, mark.time=TRUE, ylim=c(0.5, 1), xlim=c(0, 3000), xlab='Time (days)', ylab='Cumulative Survival Percentage (%)', main="Kaplan-Meier Survival Curve for ANKRD52 Gene", col=c(2,4, 10))
legend("topright", inset=0.01, legend=c("high ANKRD52 methylation rate", "low ANKRD52 methylation rate"), col=c("red", "blue"), lty=1:2, cex=0.8)

library(survival)    # survival curves library
library(glmnet)      # lasso regression library

TEST_DATA_SIZE <- floor(NUM_PATIENTS * 0.3)

for (q in 1:10) { 
  t <- 1
  OFFSET <- (t-1) * TEST_DATA_SIZE
  cancerDataTest <- cancerData[(1+OFFSET):(TEST_DATA_SIZE + OFFSET), ]
  timeTest <- time[(1+OFFSET):(TEST_DATA_SIZE + OFFSET)]
  eventTest <- event[(1+OFFSET):(TEST_DATA_SIZE + OFFSET)]
  
  # take out data 
  cancerDataTrain <- cancerData[-((1+OFFSET):(TEST_DATA_SIZE + OFFSET)), ]
  timeTrain <- time[-((1+OFFSET):(TEST_DATA_SIZE + OFFSET))]
  eventTrain <- event[-((1+OFFSET):(TEST_DATA_SIZE + OFFSET))]
  
  cancerGenes1Cox <- data.frame(matrix(0,1,1))
  cancerGenes1Lasso <- data.frame(matrix(0,1,1))
  cancerGenes1Ridge <- data.frame(matrix(0,1,1))
  cancerGenes1EN <- data.frame(matrix(0,1,1))
  
  goodCancerGenesCox <- data.frame(matrix(0,1,1))
  goodCancerGenesLasso <- data.frame(matrix(0,1,1))
  goodCancerGenesRidge <- data.frame(matrix(0,1,1))
  goodCancerGenesEN <- data.frame(matrix(0,1,1))
  
  badCancerGenesCox <- data.frame(matrix(0,1,1))
  badCancerGenesLasso <- data.frame(matrix(0,1,1))
  badCancerGenesRidge <- data.frame(matrix(0,1,1))
  badCancerGenesEN <- data.frame(matrix(0,1,1))
  
  genesPerBatch <- 30
  initializationFlag <- TRUE
  for(i in 1:floor(NUM_GENES/(genesPerBatch))) { 
    OFFSET <- (i-1) * genesPerBatch
    batch <- cancerDataTrain[, (1+OFFSET):(OFFSET+genesPerBatch)]
    timeDF <- data.frame(timeTrain)
    colnames(timeDF)[1] <- "time"
    eventDF <- data.frame(eventTrain)
    colnames(eventDF)[1] <- "status"
    
    # Cox regression
    modCoxRegression <- coxph(formula = Surv(timeTrain, eventTrain) ~ ., data = batch)
    # summary(modCoxRegression)
    coxCoeff <- data.frame(sort(coef(modCoxRegression), decreasing = TRUE))
    colnames(coxCoeff) <- "cox regression coefficients"
    
    # lasso regression
    cv.fit <- cv.glmnet(as.matrix(batch), Surv(as.matrix(timeDF), as.matrix(eventDF)), family = "cox", alpha = 1, maxit = 1000)
    lassoCoeff <- data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
    colnames(lassoCoeff) <- "lassoRegressionCoefficients"
    lassoCoeff <- lassoCoeff[order(-lassoCoeff$lassoRegressionCoefficients), , drop = FALSE]
    
    # ridge regression
    cv.fit <- cv.glmnet(as.matrix(batch), Surv(as.matrix(timeDF), as.matrix(eventDF)), family = "cox", alpha = 0, maxit = 1000)
    ridgeCoeff <- data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
    colnames(ridgeCoeff) <- "ridgeRegressionCoefficients"
    ridgeCoeff <- ridgeCoeff[order(-ridgeCoeff$ridgeRegressionCoefficients), , drop = FALSE]
    
    # elastic net 
    cv.fit <- cv.glmnet(as.matrix(batch), Surv(as.matrix(timeDF), as.matrix(eventDF)), family = "cox", alpha = 0.5, maxit = 1000)
    ENCoeff <- data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
    colnames(ENCoeff) <- "ENCoefficients"
    ENCoeff <- ENCoeff[order(-ENCoeff$ENCoefficients), , drop = FALSE]
    
    
    # initialize or update output data based on the coefficients
    if(initializationFlag) {
      # Cox
      badCancerGenesCox <- data.frame(coxCoeff[1:3, ])
      colnames(badCancerGenesCox) <- "Coefficients"
      row.names(badCancerGenesCox) <- row.names(coxCoeff)[1:3]
      goodCancerGenesCox <- data.frame(coxCoeff[(genesPerBatch-2):genesPerBatch, ])
      colnames(goodCancerGenesCox) <- "Coefficients"
      row.names(goodCancerGenesCox) <- row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes1Cox <- data.frame(cbind(batch[, row.names(coxCoeff)[1:3]], batch[, row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch]]))
      
      # Lasso
      badCancerGenesLasso <- data.frame(lassoCoeff[1:3, ])
      colnames(badCancerGenesLasso) <- "Coefficients"
      row.names(badCancerGenesLasso) <- row.names(lassoCoeff)[1:3]
      goodCancerGenesLasso <- data.frame(lassoCoeff[(genesPerBatch-2):genesPerBatch, ])
      colnames(goodCancerGenesLasso) <- "Coefficients"
      row.names(goodCancerGenesLasso) <- row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes1Lasso <- data.frame(cbind(batch[, row.names(lassoCoeff)[1:3]], batch[, row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch]]))
      
      # Ridge
      badCancerGenesRidge <- data.frame(ridgeCoeff[1:3, ])
      colnames(badCancerGenesRidge) <- "Coefficients"
      row.names(badCancerGenesRidge) <- row.names(ridgeCoeff)[1:3]
      goodCancerGenesRidge <- data.frame(ridgeCoeff[(genesPerBatch-2):genesPerBatch, ])
      colnames(goodCancerGenesRidge) <- "Coefficients"
      row.names(goodCancerGenesRidge) <- row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes1Ridge <- data.frame(cbind(batch[, row.names(ridgeCoeff)[1:3]], batch[, row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch]]))
      
      # Elastic net
      badCancerGenesEN <- data.frame(ENCoeff[1:3, ])
      colnames(badCancerGenesEN) <- "Coefficients"
      row.names(badCancerGenesEN) <- row.names(ENCoeff)[1:3]
      goodCancerGenesEN <- data.frame(ENCoeff[(genesPerBatch-2):genesPerBatch, ])
      colnames(goodCancerGenesEN) <- "Coefficients"
      row.names(goodCancerGenesEN) <- row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes1EN <- data.frame(cbind(batch[, row.names(ENCoeff)[1:3]], batch[, row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch]]))
      
      initializationFlag <- FALSE
    } else {
      # Cox
      badCancerGenesCox <- rbind(coxCoeff[1,], coxCoeff[2,], coxCoeff[3,], badCancerGenesCox)
      row.names(badCancerGenesCox)[1:3] <- row.names(coxCoeff)[1:3]
      goodCancerGenesCox <- rbind(coxCoeff[(genesPerBatch-2),], coxCoeff[(genesPerBatch-1),], coxCoeff[genesPerBatch,], goodCancerGenesCox)
      row.names(goodCancerGenesCox)[1:3] <- row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes1Cox <- cbind(batch[row.names(coxCoeff)[1:3]], batch[, row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch]], cancerGenes1Cox)
      colnames(cancerGenes1Cox)[1:6] <- c(row.names(coxCoeff)[1:3], row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch])
      
      # Lasso
      badCancerGenesLasso <- rbind(lassoCoeff[1,], lassoCoeff[2,], lassoCoeff[3,], badCancerGenesLasso)
      row.names(badCancerGenesLasso)[1:3] <- row.names(lassoCoeff)[1:3]
      goodCancerGenesLasso <- rbind(lassoCoeff[(genesPerBatch-2),], lassoCoeff[(genesPerBatch-1),], lassoCoeff[genesPerBatch,], goodCancerGenesLasso)
      row.names(goodCancerGenesLasso)[1:3] <- row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes1Lasso <- cbind(batch[row.names(lassoCoeff)[1:3]], batch[, row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch]], cancerGenes1Lasso)
      colnames(cancerGenes1Lasso)[1:6] <- c(row.names(lassoCoeff)[1:3], row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch])
      
      # Ridge
      badCancerGenesRidge <- rbind(ridgeCoeff[1,], ridgeCoeff[2,], ridgeCoeff[3,], badCancerGenesRidge)
      row.names(badCancerGenesRidge)[1:3] <- row.names(ridgeCoeff)[1:3]
      goodCancerGenesRidge <- rbind(ridgeCoeff[(genesPerBatch-2),], ridgeCoeff[(genesPerBatch-1),], ridgeCoeff[genesPerBatch,], goodCancerGenesRidge)
      row.names(goodCancerGenesRidge)[1:3] <- row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes1Ridge <- cbind(batch[row.names(ridgeCoeff)[1:3]], batch[, row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch]], cancerGenes1Ridge)
      colnames(cancerGenes1Ridge)[1:6] <- c(row.names(ridgeCoeff)[1:3], row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch])
      
      # Elastic net
      badCancerGenesEN <- rbind(ENCoeff[1,], ENCoeff[2,], ENCoeff[3,], badCancerGenesEN)
      row.names(badCancerGenesEN)[1:3] <- row.names(ENCoeff)[1:3]
      goodCancerGenesEN <- rbind(ENCoeff[(genesPerBatch-2),], ENCoeff[(genesPerBatch-1),], ENCoeff[genesPerBatch,], goodCancerGenesEN)
      row.names(goodCancerGenesEN)[1:3] <- row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes1EN <- cbind(batch[row.names(ENCoeff)[1:3]], batch[, row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch]], cancerGenes1EN)
      colnames(cancerGenes1EN)[1:6] <- c(row.names(ENCoeff)[1:3], row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch])
    }
  }
  goodCancerGenesCox <- goodCancerGenesCox[order(-goodCancerGenesCox$Coefficients), , drop = FALSE]
  goodCancerGenesLasso <- goodCancerGenesLasso[order(-goodCancerGenesLasso$Coefficients), , drop = FALSE]
  goodCancerGenesRidge <- goodCancerGenesRidge[order(-goodCancerGenesRidge$Coefficients), , drop = FALSE]
  goodCancerGenesEN <- goodCancerGenesEN[order(-goodCancerGenesEN$Coefficients), , drop = FALSE]
  badCancerGenesCox <- badCancerGenesCox[order(-badCancerGenesCox$Coefficients), , drop = FALSE]
  badCancerGenesLasso <- badCancerGenesLasso[order(-badCancerGenesLasso$Coefficients), , drop = FALSE]
  badCancerGenesRidge <- badCancerGenesRidge[order(-badCancerGenesRidge$Coefficients), , drop = FALSE]
  badCancerGenesEN <- badCancerGenesEN[order(-badCancerGenesEN$Coefficients), , drop = FALSE]
  cancerGenes2Cox <- data.frame(matrix(0,1,1))
  cancerGenes2Lasso <- data.frame(matrix(0,1,1))
  cancerGenes2Ridge <- data.frame(matrix(0,1,1))
  cancerGenes2EN <- data.frame(matrix(0,1,1))
  betterCancerGenesCox <- data.frame(matrix(0,1,1))
  betterCancerGenesLasso <- data.frame(matrix(0,1,1))
  betterCancerGenesRidge <- data.frame(matrix(0,1,1))
  betterCancerGenesEN <- data.frame(matrix(0,1,1))
  worseCancerGenesCox <- data.frame(matrix(0,1,1))
  worseCancerGenesLasso <- data.frame(matrix(0,1,1))
  worseCancerGenesRidge <- data.frame(matrix(0,1,1))
  worseCancerGenesEN <- data.frame(matrix(0,1,1))
  initializationFlag <- TRUE
  genesPerBatch <- 30
  
  for (i in 1:floor((dim(cancerGenes1Cox)[2])/(genesPerBatch))) { 
    OFFSET <- (i-1) * genesPerBatch
    timeDF <- data.frame(timeTrain)
    colnames(timeDF)[1] <- "time"
    eventDF <- data.frame(eventTrain)
    colnames(eventDF)[1] <- "status"
    # Cox regression
    batchCox <- cancerGenes1Cox[, (1+OFFSET):(OFFSET+genesPerBatch)]
    modCoxRegression <- coxph(formula = Surv(timeTrain, eventTrain) ~ ., data = batchCox)
    # summary(modCoxRegression)
    coxCoeff <- data.frame(sort(coef(modCoxRegression), decreasing = TRUE))
    colnames(coxCoeff) <- "cox regression coefficients"
    # lasso regression
    batchLasso <- cancerGenes1Lasso[, (1+OFFSET):(OFFSET+genesPerBatch)]
    cv.fit <- cv.glmnet(as.matrix(batchLasso), Surv(as.matrix(timeDF), as.matrix(eventDF)), family = "cox", alpha = 1, maxit = 1000)
    lassoCoeff <- data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
    colnames(lassoCoeff) <- "lassoRegressionCoefficients"
    lassoCoeff <- lassoCoeff[order(-lassoCoeff$lassoRegressionCoefficients), , drop = FALSE]
    # ridge regression
    batchRidge <- cancerGenes1Ridge[, (1+OFFSET):(OFFSET+genesPerBatch)]
    cv.fit <- cv.glmnet(as.matrix(batchRidge), Surv(as.matrix(timeDF), as.matrix(eventDF)), family = "cox", alpha = 0, maxit = 1000)
    ridgeCoeff <- data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
    colnames(ridgeCoeff) <- "ridgeRegressionCoefficients"
    ridgeCoeff <- ridgeCoeff[order(-ridgeCoeff$ridgeRegressionCoefficients), , drop = FALSE]
    # elastic net
    batchEN <- cancerGenes1EN[, (1+OFFSET):(OFFSET+genesPerBatch)]
    cv.fit <- cv.glmnet(as.matrix(batchEN), Surv(as.matrix(timeDF), as.matrix(eventDF)), family = "cox", alpha = 0.5, maxit = 1000)
    ENCoeff <- data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
    colnames(ENCoeff) <- "ENCoefficients"
    ENCoeff <- ENCoeff[order(-ENCoeff$ENCoefficients), , drop = FALSE]
    
    if(initializationFlag) {
      # Cox
      worseCancerGenesCox <- data.frame(coxCoeff[1:3, ])
      colnames(worseCancerGenesCox) <- "Coefficients"
      row.names(worseCancerGenesCox) <- row.names(coxCoeff)[1:3]
      betterCancerGenesCox <- data.frame(coxCoeff[(genesPerBatch-2):genesPerBatch, ])
      colnames(betterCancerGenesCox) <- "Coefficients"
      row.names(betterCancerGenesCox) <- row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes2Cox <- data.frame(cbind(batchCox[, row.names(coxCoeff)[1:3]], batchCox[, row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch]]))
      # Lasso
      worseCancerGenesLasso <- data.frame(lassoCoeff[1:3, ])
      colnames(worseCancerGenesLasso) <- "Coefficients"
      row.names(worseCancerGenesLasso) <- row.names(lassoCoeff)[1:3]
      betterCancerGenesLasso <- data.frame(lassoCoeff[(genesPerBatch-2):genesPerBatch, ])
      colnames(betterCancerGenesLasso) <- "Coefficients"
      row.names(betterCancerGenesLasso) <- row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes2Lasso <- data.frame(cbind(batchLasso[, row.names(lassoCoeff)[1:3]], batchLasso[, row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch]]))
      # Ridge
      worseCancerGenesRidge <- data.frame(ridgeCoeff[1:3, ])
      colnames(worseCancerGenesRidge) <- "Coefficients"
      row.names(worseCancerGenesRidge) <- row.names(ridgeCoeff)[1:3]
      betterCancerGenesRidge <- data.frame(ridgeCoeff[(genesPerBatch-2):genesPerBatch, ])
      colnames(betterCancerGenesRidge) <- "Coefficients"
      row.names(betterCancerGenesRidge) <- row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes2Ridge <- data.frame(cbind(batchRidge[, row.names(ridgeCoeff)[1:3]], batchRidge[, row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch]]))
      # Elastic net
      worseCancerGenesEN <- data.frame(ENCoeff[1:3, ])
      colnames(worseCancerGenesEN) <- "Coefficients"
      row.names(worseCancerGenesEN) <- row.names(ENCoeff)[1:3]
      betterCancerGenesEN <- data.frame(ENCoeff[(genesPerBatch-2):genesPerBatch, ])
      colnames(betterCancerGenesEN) <- "Coefficients"
      row.names(betterCancerGenesEN) <- row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes2EN <- data.frame(cbind(batchEN[, row.names(ENCoeff)[1:3]], batchEN[, row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch]]))
      initializationFlag <- FALSE
    } else {
      # Cox
      worseCancerGenesCox <- rbind(coxCoeff[1,], coxCoeff[2,], coxCoeff[3,], worseCancerGenesCox)
      row.names(worseCancerGenesCox)[1:3] <- row.names(coxCoeff)[1:3]
      betterCancerGenesCox <- rbind(coxCoeff[(genesPerBatch-2),], coxCoeff[(genesPerBatch-1),], coxCoeff[genesPerBatch,], betterCancerGenesCox)
      row.names(betterCancerGenesCox)[1:3] <- row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes2Cox <- cbind(batchCox[row.names(coxCoeff)[1:3]], batchCox[, row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch]], cancerGenes2Cox)
      colnames(cancerGenes2Cox)[1:6] <- c(row.names(coxCoeff)[1:3], row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch])
      # Lasso
      worseCancerGenesLasso <- rbind(lassoCoeff[1,], lassoCoeff[2,], lassoCoeff[3,], worseCancerGenesLasso)
      row.names(worseCancerGenesLasso)[1:3] <- row.names(lassoCoeff)[1:3]
      betterCancerGenesLasso <- rbind(lassoCoeff[(genesPerBatch-2),], lassoCoeff[(genesPerBatch-1),], lassoCoeff[genesPerBatch,], betterCancerGenesLasso)
      row.names(betterCancerGenesLasso)[1:3] <- row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes2Lasso <- cbind(batchLasso[row.names(lassoCoeff)[1:3]], batchLasso[, row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch]], cancerGenes2Lasso)
      colnames(cancerGenes2Lasso)[1:6] <- c(row.names(lassoCoeff)[1:3], row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch])
      # Ridge
      worseCancerGenesRidge <- rbind(ridgeCoeff[1,], ridgeCoeff[2,], ridgeCoeff[3,], worseCancerGenesRidge)
      row.names(worseCancerGenesRidge)[1:3] <- row.names(ridgeCoeff)[1:3]
      betterCancerGenesRidge <- rbind(ridgeCoeff[(genesPerBatch-2),], ridgeCoeff[(genesPerBatch-1),], ridgeCoeff[genesPerBatch,], betterCancerGenesRidge)
      row.names(betterCancerGenesRidge)[1:3] <- row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes2Ridge <- cbind(batchRidge[row.names(ridgeCoeff)[1:3]], batchRidge[, row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch]], cancerGenes2Ridge)
      colnames(cancerGenes2Ridge)[1:6] <- c(row.names(ridgeCoeff)[1:3], row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch])
      # Elastic net
      worseCancerGenesEN <- rbind(ENCoeff[1,], ENCoeff[2,], ENCoeff[3,], worseCancerGenesEN)
      row.names(worseCancerGenesEN)[1:3] <- row.names(ENCoeff)[1:3]
      betterCancerGenesEN <- rbind(ENCoeff[(genesPerBatch-2),], ENCoeff[(genesPerBatch-1),], ENCoeff[genesPerBatch,], betterCancerGenesEN)
      row.names(betterCancerGenesEN)[1:3] <- row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes2EN <- cbind(batchEN[row.names(ENCoeff)[1:3]], batchEN[, row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch]], cancerGenes2EN)
      colnames(cancerGenes2EN)[1:6] <- c(row.names(ENCoeff)[1:3], row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch])
    }
  }
  
  betterCancerGenesCox <- betterCancerGenesCox[order(-betterCancerGenesCox$Coefficients), , drop = FALSE]
  betterCancerGenesLasso <- betterCancerGenesLasso[order(-betterCancerGenesLasso$Coefficients), , drop = FALSE]
  betterCancerGenesRidge <- betterCancerGenesRidge[order(-betterCancerGenesRidge$Coefficients), , drop = FALSE]
  betterCancerGenesEN <- betterCancerGenesEN[order(-betterCancerGenesEN$Coefficients), , drop = FALSE]
  worseCancerGenesCox <- worseCancerGenesCox[order(-worseCancerGenesCox$Coefficients), , drop = FALSE]
  worseCancerGenesLasso <- worseCancerGenesLasso[order(-worseCancerGenesLasso$Coefficients), , drop = FALSE]
  worseCancerGenesRidge <- worseCancerGenesRidge[order(-worseCancerGenesRidge$Coefficients), , drop = FALSE]
  worseCancerGenesEN <- worseCancerGenesEN[order(-worseCancerGenesEN$Coefficients), , drop = FALSE]
  cancerGenes3Cox <- data.frame(matrix(0,1,1))
  cancerGenes3Lasso <- data.frame(matrix(0,1,1))
  cancerGenes3Ridge <- data.frame(matrix(0,1,1))
  cancerGenes3EN <- data.frame(matrix(0,1,1))
  bestCancerGenesCox <- data.frame(matrix(0,1,1))
  bestCancerGenesLasso <- data.frame(matrix(0,1,1))
  bestCancerGenesRidge <- data.frame(matrix(0,1,1))
  bestCancerGenesEN <- data.frame(matrix(0,1,1))
  worstCancerGenesCox <- data.frame(matrix(0,1,1))
  worstCancerGenesLasso <- data.frame(matrix(0,1,1))
  worstCancerGenesRidge <- data.frame(matrix(0,1,1))
  worstCancerGenesEN <- data.frame(matrix(0,1,1))
  initializationFlag <- TRUE
  genesPerBatch <- 30
  
  for (i in 1:floor((dim(cancerGenes2Cox)[2])/(genesPerBatch))) { 
    OFFSET <- (i-1) * genesPerBatch
    timeDF <- data.frame(timeTrain)
    colnames(timeDF)[1] <- "time"
    eventDF <- data.frame(eventTrain)
    colnames(eventDF)[1] <- "status"
    # Cox regression
    batchCox <- cancerGenes2Cox[, (1+OFFSET):(OFFSET+genesPerBatch)]
    modCoxRegression <- coxph(formula = Surv(timeTrain, eventTrain) ~ ., data = batchCox)
    # summary(modCoxRegression)
    coxCoeff <- data.frame(sort(coef(modCoxRegression), decreasing = TRUE))
    colnames(coxCoeff) <- "cox regression coefficients"
    # lasso regression
    batchLasso <- cancerGenes2Lasso[, (1+OFFSET):(OFFSET+genesPerBatch)]
    cv.fit <- cv.glmnet(as.matrix(batchLasso), Surv(as.matrix(timeDF), as.matrix(eventDF)), family = "cox", alpha = 1, maxit = 1000)
    lassoCoeff <- data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
    colnames(lassoCoeff) <- "lassoRegressionCoefficients"
    lassoCoeff <- lassoCoeff[order(-lassoCoeff$lassoRegressionCoefficients), , drop = FALSE]
    # ridge regression
    batchRidge <- cancerGenes2Ridge[, (1+OFFSET):(OFFSET+genesPerBatch)]
    cv.fit <- cv.glmnet(as.matrix(batchRidge), Surv(as.matrix(timeDF), as.matrix(eventDF)), family = "cox", alpha = 0, maxit = 1000)
    ridgeCoeff <- data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
    colnames(ridgeCoeff) <- "ridgeRegressionCoefficients"
    ridgeCoeff <- ridgeCoeff[order(-ridgeCoeff$ridgeRegressionCoefficients), , drop = FALSE]
    # elastic net
    batchEN <- cancerGenes2EN[, (1+OFFSET):(OFFSET+genesPerBatch)]
    cv.fit <- cv.glmnet(as.matrix(batchEN), Surv(as.matrix(timeDF), as.matrix(eventDF)), family = "cox", alpha = 0.5, maxit = 1000)
    ENCoeff <- data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
    colnames(ENCoeff) <- "ENCoefficients"
    ENCoeff <- ENCoeff[order(-ENCoeff$ENCoefficients), , drop = FALSE]
    if(initializationFlag) {
      # Cox
      worstCancerGenesCox <- data.frame(coxCoeff[1:3, ])
      colnames(worstCancerGenesCox) <- "Coefficients"
      row.names(worstCancerGenesCox) <- row.names(coxCoeff)[1:3]
      bestCancerGenesCox <- data.frame(coxCoeff[(genesPerBatch-2):genesPerBatch, ])
      colnames(bestCancerGenesCox) <- "Coefficients"
      row.names(bestCancerGenesCox) <- row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes3Cox <- data.frame(cbind(batchCox[, row.names(coxCoeff)[1:3]], batchCox[, row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch]]))
      # Lasso
      worstCancerGenesLasso <- data.frame(lassoCoeff[1:3, ])
      colnames(worstCancerGenesLasso) <- "Coefficients"
      row.names(worstCancerGenesLasso) <- row.names(lassoCoeff)[1:3]
      bestCancerGenesLasso <- data.frame(lassoCoeff[(genesPerBatch-2):genesPerBatch, ])
      colnames(bestCancerGenesLasso) <- "Coefficients"
      row.names(bestCancerGenesLasso) <- row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes3Lasso <- data.frame(cbind(batchLasso[, row.names(lassoCoeff)[1:3]], batchLasso[, row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch]]))
      # Ridge
      worstCancerGenesRidge <- data.frame(ridgeCoeff[1:3, ])
      colnames(worstCancerGenesRidge) <- "Coefficients"
      row.names(worstCancerGenesRidge) <- row.names(ridgeCoeff)[1:3]
      bestCancerGenesRidge <- data.frame(ridgeCoeff[(genesPerBatch-2):genesPerBatch, ])
      colnames(bestCancerGenesRidge) <- "Coefficients"
      row.names(bestCancerGenesRidge) <- row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes3Ridge <- data.frame(cbind(batchRidge[, row.names(ridgeCoeff)[1:3]], batchRidge[, row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch]]))
      # Elastic net
      worstCancerGenesEN <- data.frame(ENCoeff[1:3, ])
      colnames(worstCancerGenesEN) <- "Coefficients"
      row.names(worstCancerGenesEN) <- row.names(ENCoeff)[1:3]
      bestCancerGenesEN <- data.frame(ENCoeff[(genesPerBatch-2):genesPerBatch, ])
      colnames(bestCancerGenesEN) <- "Coefficients"
      row.names(bestCancerGenesEN) <- row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes3EN <- data.frame(cbind(batchEN[, row.names(ENCoeff)[1:3]], batchEN[, row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch]]))
      initializationFlag <- FALSE
    }
    else {
      # Cox
      worstCancerGenesCox <- rbind(coxCoeff[1,], coxCoeff[2,], coxCoeff[3,], worstCancerGenesCox)
      row.names(worstCancerGenesCox)[1:3] <- row.names(coxCoeff)[1:3]
      bestCancerGenesCox <- rbind(coxCoeff[(genesPerBatch-2),], coxCoeff[(genesPerBatch-1),], coxCoeff[genesPerBatch,], bestCancerGenesCox)
      row.names(bestCancerGenesCox)[1:3] <- row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes3Cox <- cbind(batchCox[row.names(coxCoeff)[1:3]], batchCox[, row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch]], cancerGenes3Cox)
      colnames(cancerGenes3Cox)[1:6] <- c(row.names(coxCoeff)[1:3], row.names(coxCoeff)[(genesPerBatch-2):genesPerBatch])
      # Lasso
      worstCancerGenesLasso <- rbind(lassoCoeff[1,], lassoCoeff[2,], lassoCoeff[3,], worstCancerGenesLasso)
      row.names(worstCancerGenesLasso)[1:3] <- row.names(lassoCoeff)[1:3]
      bestCancerGenesLasso <- rbind(lassoCoeff[(genesPerBatch-2),], lassoCoeff[(genesPerBatch-1),], lassoCoeff[genesPerBatch,], bestCancerGenesLasso)
      row.names(bestCancerGenesLasso)[1:3] <- row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes3Lasso <- cbind(batchLasso[row.names(lassoCoeff)[1:3]], batchLasso[, row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch]], cancerGenes3Lasso)
      colnames(cancerGenes3Lasso)[1:6] <- c(row.names(lassoCoeff)[1:3], row.names(lassoCoeff)[(genesPerBatch-2):genesPerBatch])
      # Ridge
      worstCancerGenesRidge <- rbind(ridgeCoeff[1,], ridgeCoeff[2,], ridgeCoeff[3,], worstCancerGenesRidge)
      row.names(worstCancerGenesRidge)[1:3] <- row.names(ridgeCoeff)[1:3]
      bestCancerGenesRidge <- rbind(ridgeCoeff[(genesPerBatch-2),], ridgeCoeff[(genesPerBatch-1),], ridgeCoeff[genesPerBatch,], bestCancerGenesRidge)
      row.names(bestCancerGenesRidge)[1:3] <- row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes3Ridge <- cbind(batchRidge[row.names(ridgeCoeff)[1:3]], batchRidge[, row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch]], cancerGenes3Ridge)
      colnames(cancerGenes3Ridge)[1:6] <- c(row.names(ridgeCoeff)[1:3], row.names(ridgeCoeff)[(genesPerBatch-2):genesPerBatch])
      # Elastic net
      worstCancerGenesEN <- rbind(ENCoeff[1,], ENCoeff[2,], ENCoeff[3,], worstCancerGenesEN)
      row.names(worstCancerGenesEN)[1:3] <- row.names(ENCoeff)[1:3]
      bestCancerGenesEN <- rbind(ENCoeff[(genesPerBatch-2),], ENCoeff[(genesPerBatch-1),], ENCoeff[genesPerBatch,], bestCancerGenesEN)
      row.names(bestCancerGenesEN)[1:3] <- row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch]
      cancerGenes3EN <- cbind(batchEN[row.names(ENCoeff)[1:3]], batchEN[, row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch]], cancerGenes3EN)
      colnames(cancerGenes3EN)[1:6] <- c(row.names(ENCoeff)[1:3], row.names(ENCoeff)[(genesPerBatch-2):genesPerBatch])
    }
  }
  bestCancerGenesCox <- bestCancerGenesCox[order(-bestCancerGenesCox$Coefficients), , drop = FALSE]
  bestCancerGenesLasso <- bestCancerGenesLasso[order(-bestCancerGenesLasso$Coefficients), , drop = FALSE]
  bestCancerGenesRidge <- bestCancerGenesRidge[order(-bestCancerGenesRidge$Coefficients), , drop = FALSE]
  bestCancerGenesEN <- bestCancerGenesEN[order(-bestCancerGenesEN$Coefficients), , drop = FALSE]
  worstCancerGenesCox <- worstCancerGenesCox[order(-worstCancerGenesCox$Coefficients), , drop = FALSE]
  worstCancerGenesLasso <- worstCancerGenesLasso[order(-worstCancerGenesLasso$Coefficients), , drop = FALSE]
  worstCancerGenesRidge <- worstCancerGenesRidge[order(-worstCancerGenesRidge$Coefficients), , drop = FALSE]
  worstCancerGenesEN <- worstCancerGenesEN[order(-worstCancerGenesEN$Coefficients), , drop = FALSE]
  
  # END OF DATA PROCESSING
  
  # IDENTIFY GENES
  timeDF <- data.frame(timeTrain)
  colnames(timeDF)[1] <- "time"
  eventDF <- data.frame(eventTrain)
  colnames(eventDF)[1] <- "status"
  
  # Cox regression
  modCoxRegression <- coxph(formula = Surv(timeTrain, eventTrain) ~ ., data = cancerGenes3Cox)
  # summary(modCoxRegression)
  coxCoeffFinal <- data.frame(sort(coef(modCoxRegression), decreasing = TRUE))
  colnames(coxCoeffFinal) <- "cox regression coefficients"
  
  # lasso regression
  cv.fit <- cv.glmnet(as.matrix(cancerGenes3Lasso), Surv(as.matrix(timeDF), as.matrix(eventDF)), family = "cox", alpha = 1, maxit = 1000)
  lassoCoeffFinal <- data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
  colnames(lassoCoeffFinal) <- "lassoRegressionCoefficients"
  lassoCoeffFinal <- lassoCoeff[order(-lassoCoeffFinal$lassoRegressionCoefficients), , drop = FALSE]
  
  # ridge regression
  cv.fit <- cv.glmnet(as.matrix(cancerGenes3Ridge), Surv(as.matrix(timeDF), as.matrix(eventDF)), family = "cox", alpha = 0, maxit = 1000)
  ridgeCoeffFinal <- data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
  colnames(ridgeCoeffFinal) <- "ridgeRegressionCoefficients"
  ridgeCoeffFinal <- ridgeCoeffFinal[order(-ridgeCoeffFinal$ridgeRegressionCoefficients), , drop = FALSE]
  
  # elastic net
  cv.fit <- cv.glmnet(as.matrix(cancerGenes3EN[,1:24]), Surv(as.matrix(timeDF), as.matrix(eventDF)), family = "cox", alpha = 0.5, maxit = 1000)
  ENCoeff1 <- data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
  colnames(ENCoeff1) <- "ENCoefficients"
  ENCoeff1 <- ENCoeff1[order(-ENCoeff1$ENCoefficients), , drop = FALSE]
  cv.fit <- cv.glmnet(as.matrix(cancerGenes3EN[,25:48]), Surv(as.matrix(timeDF), as.matrix(eventDF)), family = "cox", alpha = 0.5, maxit = 1000)
  ENCoeff2 <- data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
  colnames(ENCoeff2) <- "ENCoefficients"
  ENCoeff2 <- ENCoeff2[order(-ENCoeff2$ENCoefficients), , drop = FALSE]
  ENCoeffFinal <- rbind(ENCoeff1,ENCoeff2)
  
  
  ## validation 
  modCox <- coxph(formula = Surv(timeTrain, eventTrain) ~ ., data = cancerGenes3Cox)
  relRiskCox <- sort(predict(modCoxRegression, cbind(cbind(timeTest, eventTest), cancerDataTest), type="risk"))
  modLasso <- coxph(formula = Surv(timeTrain, eventTrain) ~ ., data = cancerGenes3Lasso)
  relRiskLasso <- sort(predict(modLasso, cbind(cbind(timeTest, eventTest), cancerDataTest), type="risk"))
  modRidge <- coxph(formula = Surv(timeTrain, eventTrain) ~ ., data = cancerGenes3Ridge)
  relRiskRidge <- sort(predict(modRidge, cbind(cbind(timeTest, eventTest), cancerDataTest), type="risk"))
  modEN <- coxph(formula = Surv(timeTrain, eventTrain) ~ ., data = cancerGenes3EN)
  relRiskEN <- sort(predict(modEN, cbind(cbind(timeTest, eventTest), cancerDataTest), type="risk"))
  
  riskDeathCox <- 0
  riskSurvivorsCox <- 0
  riskDeathLasso <- 0
  riskSurvivorsLasso <- 0
  riskDeathRidge <- 0
  riskSurvivorsRidge <- 0
  riskDeathEN <- 0
  riskSurvivorsEN <- 0
  
  cntDeaths <- 0
  cntSurvivors <- 0
  
  for (i in 1:TEST_DATA_SIZE){
    if (eventTest[i] == 1){ 
      riskDeathCox <- riskDeathCox + as.numeric(relRiskCox[i])
      riskDeathLasso <- riskDeathLasso + as.numeric(relRiskLasso[i])
      riskDeathRidge <- riskDeathRidge + as.numeric(relRiskRidge[i])
      riskDeathEN <- riskDeathEN + as.numeric(relRiskEN[i])
      cntDeaths <- cntDeaths + 1
    } else{
      riskSurvivorsCox <- riskSurvivorsCox + as.numeric(relRiskCox[i])
      riskSurvivorsLasso <- riskSurvivorsLasso + as.numeric(relRiskLasso[i])
      riskSurvivorsRidge <- riskSurvivorsRidge + as.numeric(relRiskRidge[i])
      riskSurvivorsEN <- riskSurvivorsEN + as.numeric(relRiskEN[i])
      cntSurvivors <- cntSurvivors + 1
    }
  }
  riskDeathCox <- riskDeathCox / cntDeaths
  riskSurvivorsCox <- riskSurvivorsCox / cntSurvivors
  riskDeathLasso <- riskDeathLasso / cntDeaths
  riskSurvivorsLasso <- riskSurvivorsLasso / cntSurvivors
  riskDeathRidge <- riskDeathRidge / cntDeaths
  riskSurvivorsRidge <- riskSurvivorsRidge / cntSurvivors
  riskDeathEN <- riskDeathEN / cntDeaths
  riskSurvivorsEN <- riskSurvivorsEN / cntSurvivors
  
  print(riskDeathCox)
  print(riskSurvivorsCox)
  print(riskDeathLasso)
  print(riskSurvivorsLasso)
  print(riskDeathRidge)
  print(riskSurvivorsRidge)
  print(riskDeathEN)
  print(riskSurvivorsEN)

}
