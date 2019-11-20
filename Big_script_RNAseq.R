setwd('H:/Desktop/PH data analysis')
pheno.all <- read.csv(file = "H:/Desktop/PH data analysis/180614 - pheno all with metabolites and diagnostic genes.csv", row.names = 1)

# Read in with tximport ---------------------------------------------------

source('H:/Desktop/PH data analysis/tximport/tximport.R')
source('H:/Desktop/PH data analysis/tximport/helper.R')
source('H:/Desktop/PH data analysis/tximport/infReps.R')
source('H:/Desktop/PH data analysis/tximport/summarizeToGene.R')

samples=read.table(file = "rnaseq_samples.txt")[,1]
samples_trial=read.table(file = "rnaseq_trial_samples_short.txt")[,1]

files=paste0(samples,"/quant.sf")
files_trial=paste0(samples_trial,"/quant.sf")

quant = read.table(file = files_trial[1], header = T
                 #,nrows = 100
                 )
quant$Name=as.character(quant$Name)
transcripts=do.call(rbind,strsplit(quant$Name,split = "[|]"))
transcripts=as.data.frame(transcripts)
transcripts$TXNAME=quant$Name

colnames(transcripts)=c("TranscriptID","GeneID","OTTGeneID","OTTTranscriptID","TranscriptName","GeneName","Length","Transcript type","TXNAME")
tx2gene=transcripts[,c("TXNAME","GeneName")]
tx2gene$GeneName=as.character(tx2gene$GeneName)
tx2gene$GeneName[is.na(tx2gene$GeneName)]<-"Other"

transcripts[which(transcripts$GeneName=="HLA-DPB1"),]
length(unique(transcripts$GeneName))

transcripts$GeneID2=do.call(rbind,strsplit(as.character(transcripts$GeneID),split="[.]"))[,1]

library(readr)

all_files=c(files,files_trial)
samples_short=do.call(rbind,strsplit(c(files,files_trial), split = "_S"))[,1]
names(all_files)<-samples_short
txi=tximport(all_files, type = "salmon", tx2gene = tx2gene)

# Sample details ----------------------------------------------------------
salmon_samples=colnames(all.tpm)
write.csv(salmon_samples,file="salmon_samples.csv")



# PCA ---------------------------------------------------------------------

all.reads     <- txi$counts
all.tpm       <- txi$abundance
log.reads     <- log10(all.reads)
log.reads[which(is.infinite(log.reads))]<-0

log.tpm       <- log10(all.tpm)
log.tpm[which(is.infinite(log.tpm))]<-0


PCAreads      <- prcomp(x = t(log.reads))
PCAtpm        <- prcomp(x = t(log.tpm))
#PCAqn=prcomp(x = quant.norm.repmin)
#length(which(is.na(ZscoreAL12.nonxeno.repmin)))

#loadings
PCAreads.loadings=as.data.frame(PCAreads$rotation)
PCAtpm.loadings=as.data.frame(PCAtpm$rotation)
#scores
PCAreads.scores=as.data.frame(PCAreads$x)
PCAtpm.scores=as.data.frame(PCAtpm$x)

# PCA plot 1 --------------------------------------------------------------


#plot PCA scores by lane
plot(x = PCAtpm.scores$PC1, y = PCAtpm.scores$PC2,
     col = pheno.all$lane,#Exclude, #pch = as.numeric(AL12pheno$GROUPING)+14,
     #xlim = c(-25,35), ylim = c(-15,15),
     bty="L")
# add legend
legend(x = "right",
       #inset=c(-0.2,0),
       legend = levels(as.factor(pheno.all$lane)),
       col = levels(as.factor(pheno.all$lane)),
       pch = "o"

       #c(rep(15:17,times=2),17,rep(15:17,times=3))
)



# Save basic files --------------------------------------------------------
save(pheno.all,all.reads,all.tpm,tx2gene,txi,
     file = "180423 - salmon pheno and data for 508 samples trial2 and Jan18.RData")

# Detection of transcripts ------------------------------------------------------
transcripts.salmon2=data.frame(gene=rownames(all.reads))
#transcripts.salmon[which(transcripts.salmon$rownames %in% rownames(tpm.filtered.ac1)),]
transcripts.salmon2$detected.controls=apply(X = all.reads, MARGIN = 1, FUN = function(x){
  length(which(x[pheno.all$PAH == 0] > 1))/length(which(pheno.all$PAH == 0))
})
transcripts.salmon2$detected.PAH=apply(X = all.reads, MARGIN = 1, FUN = function(x){
  length(which(x[pheno.all$PAH == 1] > 1))/length(which(pheno.all$PAH == 1))
})
transcripts.salmon2$max.detected.controls.PAH=apply(X = transcripts.salmon2, MARGIN = 1, FUN = function(x){
  max(x['detected.controls'],x['detected.PAH'])
})

length(which(transcripts.salmon2$max.detected.controls.PAH>=0.95))
# [1] 25966
hist(as.numeric(transcripts.salmon2$max.detected.controls.PAH),
     xlab = "Proportion of controls/patients
     gene transcripts detected in", main = "")

# Save image of workspace up to this point --------------------------------


save.image("H:/Desktop/PH data analysis/190608 - salmon RNAseq files up to edgeR.RData")
# edgeR -------------------------------------------------------------------

library(edgeR)

load("C:/Users/po3217/Desktop/Pablo/PhD/Data/RNAseq/Deconvoluted analysis/19 06 10 - New analysis/190608 - salmon RNAseq files up to edgeR.RData")

pheno.all <- read.csv("C:/Users/po3217/Desktop/Pablo/PhD/Data/RNAseq/Deconvoluted analysis/19 06 10 - New analysis/2019-09-04 NEW pheno all with deconvoluted lasso model.csv")


groupA    <- which(pheno.all$Random.group2 == "A" & is.na(pheno.all$Exclude))
groupB    <- which(pheno.all$Random.group2 == "B" & is.na(pheno.all$Exclude))
groupAB   <- which(pheno.all$Random.group2 %in% c("A","B") & is.na(pheno.all$Exclude))
groupC    <- which(pheno.all$Random.group2 == "C" & is.na(pheno.all$Exclude))
groupABC  <- which(pheno.all$Random.group2 %in% c("A","B","C") & is.na(pheno.all$Exclude))


all.reads <- txi$counts
all.tpm <- txi$abundance
cts <- txi$counts
normMat <- txi$length
normMat <- normMat/exp(rowMeans(log(normMat)))

o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))

y=DGEList(counts = all.reads[,groupA], group = pheno.all$PAH[groupA])
y$offset <- t(t(log(normMat[,groupA])) + o[groupA])

cpm=cpm(y$counts)
#rpkm=rpkm(y$counts, gene.length = transcripts.salmon2$Length)
plotMDS(y)
keep <- rowSums(cpm(y)>2) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)


design <- model.matrix(~pheno.all$PAH[groupA]
                       +pheno.all$PC1[groupA]
                       +pheno.all$PC2[groupA]
                       +pheno.all$PC3[groupA]
                       +pheno.all$ciber.T.cells.CD4.naive[groupA]
                       +pheno.all$ciber.B.cells.memory[groupA]
                       +pheno.all$ciber.Mast.cells.resting[groupA]
                       +pheno.all$ciber.Dendritic.cells.resting[groupA]
                       +pheno.all$quant.Tregs[groupA]
                       +pheno.all$quant.T.cells.CD4[groupA]
)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
plotQLDisp(fit)
plotBCV(y)
qlf <- glmQLFTest(fit,coef=2)
tophits=topTags(qlf, n = 25966)
tophits$table->tophits.table
tophits.table$gene=rownames(tophits.table)
tophits.tableA=merge(tophits.table,transcripts.salmon2,by="gene")
write.csv(tophits.tableA,file = paste0(Sys.Date()," - top hits from RNAseq analysisA.csv"))
summary(dt <- decideTestsDGE(qlf))
isDE <- as.logical(dt)
DEnames <- rownames(y)[isDE]
plotSmear(qlf, de.tags=DEnames)
abline(h=c(-1,1), col="blue")

yB=DGEList(counts = all.reads[,groupB], group = pheno.all$PAH[groupB])
yB$offset <- t(t(log(normMat[,groupB])) + o[groupB])

cpmB=cpm(yB$counts)
rpkmB=rpkm(yB$counts, gene.length = transcripts.salmon2$Length)
plotMDS(yB)
keepB <- rowSums(cpm(yB)>2) >= 3
yB <- yB[keepB, , keep.lib.sizes=FALSE]
yB <- calcNormFactors(yB)


designB <- model.matrix(~pheno.all$PAH[groupB]
                        +pheno.all$PC1[groupB]
                        +pheno.all$PC2[groupB]
                        +pheno.all$PC3[groupB]
                        +pheno.all$ciber.T.cells.CD4.naive[groupB]
                        +pheno.all$ciber.B.cells.memory[groupB]
                        +pheno.all$ciber.Mast.cells.resting[groupB]
                        +pheno.all$ciber.Dendritic.cells.resting[groupB]
                        +pheno.all$quant.Tregs[groupB]
                        +pheno.all$quant.T.cells.CD4[groupB]
)
yB <- estimateDisp(yB,designB)
fitB <- glmQLFit(yB,designB)
plotQLDisp(fitB)
plotBCV(yB)
qlfB <- glmQLFTest(fitB,coef=2)
tophitsB=topTags(qlfB, n = 25966)
tophitsB$table->tophits.tableB
tophits.tableB$gene=rownames(tophits.tableB)
tophits.tableB=merge(tophits.tableB,transcripts.salmon2,by="gene")
write.csv(tophits.tableB,file = paste0(Sys.Date()," - top hits from RNAseq analysisB.csv"))
summary(dtB <- decideTestsDGE(qlfB))
isDEB <- as.logical(dtB)
DEnamesB <- rownames(yB)[isDEB]
plotSmear(qlfB, de.tags=DEnamesB)
abline(h=c(-1,1), col="blue")


# In new AB -------------------------------------------------------------------

yAB=DGEList(counts = all.reads[,groupAB], group = pheno.all$PAH[groupAB])
yAB$offset <- t(t(log(normMat[,groupAB])) + o[groupAB])

cpmAB=cpm(yAB$counts)
rpkmAB=rpkm(yAB$counts, gene.length = transcripts.salmon2$Length)
plotMDS(yAB)
keepAB <- rowSums(cpm(yAB)>2) >= 3
yAB <- yAB[keepAB, , keep.lib.sizes=FALSE]
yAB <- calcNormFactors(yAB)


designAB <- model.matrix(~pheno.all$PAH[groupAB]
                         +pheno.all$PC1[groupAB]
                         +pheno.all$PC2[groupAB]
                         +pheno.all$PC3[groupAB]
                         +pheno.all$Age_controls_sample_patients_diagnosis[groupAB]
                         +as.character(pheno.all$sex[groupAB])

                                       +pheno.all$ciber.T.cells.CD4.naive[groupAB]
                                       +pheno.all$ciber.B.cells.memory[groupAB]
                                       +pheno.all$ciber.Mast.cells.resting[groupAB]
                                       +pheno.all$ciber.Dendritic.cells.resting[groupAB]
                                       +pheno.all$quant.Tregs[groupAB]
                                       +pheno.all$quant.T.cells.CD4[groupAB]
                                       )
yAB <- estimateDisp(yAB,designAB)
fitAB <- glmQLFit(yAB,designAB)
plotQLDisp(fitAB)
plotBCV(yAB)
qlfAB <- glmQLFTest(fitAB,coef=2)
tophitsAB=topTags(qlfAB, n = 25966)
tophitsAB$table->tophits.tableAB
tophits.tableAB$gene=rownames(tophits.tableAB)
tophits.tableAB=merge(tophits.tableAB,transcripts.salmon2,by="gene")
write.csv(tophits.tableAB,file = paste0(Sys.Date()," - top hits from RNAseq analysisAB.csv"))
summary(dtAB <- decideTestsDGE(qlfAB))
isDEAB <- as.logical(dtAB)
DEnamesAB <- rownames(yAB)[isDEAB]
plotSmear(qlfAB, de.tags=DEnamesAB)
abline(h=c(-1,1), col="blue")

colnames(tophits.tableA)

#[1] "Geneid"           "logFC"            "logCPM"           "F"                "PValue"           "FDR"
#[7] "Length"           "Gene.description" "Gene.name"        "total.counts"

colnames(tophits.tableA)[2:6]<-c("logFC_A",
                                 "logCPM_A",
                                 "F_A",
                                 "PValue_A",
                                 "FDR_A")
colnames(tophits.tableB)[2:6]<-c("logFC_B",
                                 "logCPM_B",
                                 "F_B",
                                 "PValue_B",
                                 "FDR_B")
colnames(tophits.tableAB)[2:6]<-c("logFC_AB",
                                  "logCPM_AB",
                                  "F_AB",
                                  "PValue_AB",
                                  "FDR_AB")
detectedAB=tophits.tableAB[tophits.tableAB$max.detected.controls.PAH>0.95,]
detectedAB$p.FDR=p.adjust(p = detectedAB$PValue_AB, method = "fdr")
tophits.tableABAB=merge(tophits.tableA,tophits.tableB)
tophits.tableABAB=merge(tophits.tableABAB,tophits.tableAB)
write.csv(tophits.tableABAB,file = paste0(Sys.Date()," - top hits from RNAseq analysisABAB.csv"))


sig.tophits.tableABAB=tophits.tableABAB[which(tophits.tableABAB$PValue_A < 0.05 &
                                                tophits.tableABAB$PValue_B < 0.05),]

sig.tophits.tableABAB$same.direction=sig.tophits.tableABAB$logFC_A*sig.tophits.tableABAB$logFC_B>0



# ROC analysis for each gene ----------------------------------------------
library(pROC)
sig.tophits.tableABAB$ROC.auc=do.call(rbind,lapply(X = sig.tophits.tableABAB$rownames, FUN = function(x){
  roc(pheno.all$PAH[groupAB],as.numeric(all.reads[as.character(x),groupAB]))$auc
}))


# Select well detected and same direction ---------------------------------

well.detected.tophits=sig.tophits.tableABAB[which(sig.tophits.tableABAB$max.detected.controls.PAH>=0.95 &
                                                    sig.tophits.tableABAB$same.direction == TRUE),]
welldetect.tophits.data    <- all.reads[as.character(well.detected.tophits$gene),]
welldetect.tophits.tpm=all.tpm[well.detected.tophits$gene,]
write.csv(welldetect.tophits.data, file=paste0(Sys.Date()," - data for well detected top hits Controls vs PAH in RNAseq.csv"))

welldet.tophits.cor=cor(welldetect.tophits.data)

# Hierarchical clustering of tophits --------------------------------------

clusters <- hclust(dist(welldetect.tophits.data))
plot(clusters)
clusterCut <- cutree(clusters, 15)
clusterCut

well.detected.tophits$cluster=clusterCut
well.detected.tophits.sort=well.detected.tophits[order(well.detected.tophits$PValue_AB),]
well.detected.tophits.sort$p.FDR.AB=p.adjust(well.detected.tophits.sort$PValue_AB, method = "fdr", n=1230)
well.detected.tophits.representative=well.detected.tophits.sort[which(!(duplicated(well.detected.tophits.sort$cluster))),]


write.csv(well.detected.tophits,file = paste0(Sys.Date()," - well detected hits from RNAseq analysisABAB.csv"))

# Residuals from tpm using PC1-3 ------------------------------------------

res.tpm.cohort <- t(apply(X = log.tpm, MARGIN = 1, FUN = function(x){
  model=lm(formula = x ~ PC1 + PC2 + PC3, data = PCAtpm.scores)
  model$residuals+model$coefficients[1][[1]]
}))

save(res.tpm,file="res.tpm.cohort.RData")


# Z-score to PAH patients -------------------------------------------------

PAH.mean  <- apply(res.tpm.cohort[,which(pheno.all$PAH == 1 & is.na(pheno.all$Exclude))], 1, mean)
PAH.sd    <- apply(res.tpm.cohort[,which(pheno.all$PAH == 1 & is.na(pheno.all$Exclude))], 1, sd)
z.res.tpm <- apply(res.tpm.cohort,2,function(x){
  (x-PAH.mean)/PAH.sd
})

write.csv(x = z.res.tpm, file= "z-scored residualised tpm cohort.csv")

save(z.res.tpm, file="res.tpm.cohort.RData")

# LASSO to select independent markers -------------------------------------

#install.packages("glmnet")
library(glmnet)
lasso.phenotypes <- simplepheno2[groupAB,c(3:5, 7:8)]



#Genes from deconvoluted analysis - overlap between new analyses with and without WBC fractions

deconvoluted.genes       <- read.csv("C:/Users/po3217/Desktop/Pablo/PhD/Data/RNAseq/Deconvoluted analysis/19 06 10 - New analysis/2019-08-28 - Shared genes original and WBC-deconvoluted analyses.csv", header = FALSE)
tophits.tableABAB        <- read.csv("C:/Users/po3217/Desktop/Pablo/PhD/Data/RNAseq/Deconvoluted analysis/19 06 10 - New analysis/2019-08-28 - top hits from RNAseq analysis - NoWBC.csv")


sig.tophits.tableABAB                 <- tophits.tableABAB[which(tophits.tableABAB$PValue_A < 0.05 &
                                                                 tophits.tableABAB$PValue_B < 0.05),]

sig.tophits.tableABAB$same.direction  <- sig.tophits.tableABAB$logFC_A*sig.tophits.tableABAB$logFC_B > 0

well.detected.tophits                 <- sig.tophits.tableABAB[which(sig.tophits.tableABAB$max.detected.controls.PAH >= 0.95 &
                                                                     sig.tophits.tableABAB$same.direction == TRUE),]

welldetect.tophits.data               <- all.reads[as.character(well.detected.tophits$gene),]
welldetect.tophits.tpm                <- all.tpm[well.detected.tophits$gene,]


lasso.data                    <- cbind(lasso.phenotypes,t(welldetect.tophits.data[as.character(deconvoluted.genes[,1]),groupAB]))
lasso.data$sex                <- as.numeric(lasso.data$sex)-1
lasso.data2                   <- as.matrix(lasso.data)
lasso.y                       <- as.numeric(simplepheno2[groupAB,1])
lasso.fit                     <- glmnet(x = lasso.data2, y=lasso.y)
plot(lasso.fit, label = T)


##cross-validation fit
cvfit <- glmnet::cv.glmnet(x = lasso.data2, y = lasso.y)

lasso.coef.lambda.1se         <- coef(cvfit, s = "lambda.1se")
lasso.coef.lambda.min         <- coef(cvfit, s = "lambda.min")

lasso.coef.1se.min            <-as.matrix(cbind(lasso.coef.lambda.1se,lasso.coef.lambda.min))

write.csv(x = lasso.coef.1se.min, file="lasso coefficients 1se and min double check.csv")

plot(cvfit)

lasso.tpm=as.matrix(cbind(lasso.phenotypes,t(welldetect.tophits.tpm[,groupAB])))
lasso.tpm[,'sex']=as.numeric(as.factor(lasso.tpm[,'sex']))-1
lasso.tpm2=as.matrix(lasso.tpm)
lasso.tpm.fit=glmnet(x = lasso.tpm2, y=lasso.y)
plot(lasso.tpm.fit, label = T)
##cross-validation fit
cvfit.tpm <- glmnet::cv.glmnet(x = lasso.data2, y = lasso.y)


##predict in C

lasso.phenotypesC    <- simplepheno2[groupC,c(3:5,7:8)]
lasso.dataC          <- cbind(lasso.phenotypesC,t(welldetect.tophits.data[as.character(deconvoluted.genes[,1]),groupC]))
lasso.dataC$sex      <- as.numeric(lasso.dataC$sex)-1
lasso.dataC2         <- as.matrix(lasso.dataC)


lasso.predictionC    <- predict(cvfit, newx = lasso.dataC2, s = "lambda.1se")
write.csv(x = lasso.predictionC, file = "deconvoluted lasso scores group C.csv")

library(pROC)
lasso.predictionC    <- pheno.all$lasso.decon[which(pheno.all$Random.group2 == "C" & is.na(pheno.all$Exclude))]
lassoROC             <- roc(simplepheno3$PAH[groupC],lasso.predictionC)
table(pheno.all$sex[groupC])
plot(lassoROC)
boxplot(lasso.predictionC~pheno.all$PAH[groupC])
wilcox.test(lasso.predictionC~simplepheno3$PAH[groupC])
lassoROC
ci.auc(lassoROC)

#Area under the curve: 0.8655
#95% CI: 0.7911-0.94 (DeLong)


groupABC               <- which(pheno.all$Random.group2 %in% c("A","B","C") & is.na(pheno.all$Exclude))
pheno.all$lassoModel   <- NA

lasso.phenotypesABC    <- simplepheno2[groupABC,c(3:5,7:8)]
lasso.dataABC          <- cbind(lasso.phenotypesABC,t(welldetect.tophits.data[as.character(deconvoluted.genes[,1]),groupABC]))
lasso.dataABC$sex      <- as.numeric(lasso.dataABC$sex)-1
lasso.dataABC2         <- as.matrix(lasso.dataABC)

lasso.predictionABC    <- predict(cvfit, newx = lasso.dataABC2, s = "lambda.1se")
write.csv(x = lasso.predictionABC, file = "deconvoluted lasso scores group ABC.csv")


lassoROC.ABC           <- roc(simplepheno3$PAH[groupABC],pheno.all$lasso.decon[groupABC])


table(pheno.all$sex[groupABC])
plot(lassoROC.ABC)
plot(lasso.predictionABC~simplepheno3$PAH[groupABC])
wilcox.test(lasso.predictionABC~simplepheno3$PAH[groupABC])
lassoROC.ABC
ci.auc(lassoROC.ABC)



lassoROC.ABCdf         <- data.frame(lassoROC.ABC$sensitivities,lassoROC.ABC$specificities,lassoROC.ABC$thresholds)
lassoROC.ABCdf$yi      <- lassoROC.ABCdf$lassoROC.ABC.specificities+lassoROC.ABCdf$lassoROC.ABC.sensitivities


lassoROC.ABCdf[93,]

#lassoROC.ABC.sensitivities lassoROC.ABC.specificities lassoROC.ABC.thresholds       yi
#93               0.8885794                  0.7222222                1.767813 1.610802


pheno.all$lassoModel[groupABC]=predict(cvfit,newx=lasso.dataABC2, s = "lambda.1se")

lasso.phenotypesALL=simplepheno2[,c(3:5,7:8)]
lasso.dataALL=cbind(lasso.phenotypesALL,t(welldetect.tophits.data[,]))
lasso.dataALL$sex=as.numeric(lasso.dataALL$sex)-1
lasso.dataALL2=as.matrix(lasso.dataALL)
pheno.all$lassoModel=predict(cvfit,newx=lasso.dataALL2, s = "lambda.1se")



boxplot(pheno.all$lasso.decon[which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))] ~
          paste0(pheno.all$Random.group2,pheno.all$PAH, pheno.all$sex)[which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))],
        col=2:5, ylab = "lasso Model", lty=1, ylim=c(1,2.4))
abline(h=median(pheno.all$lasso.decon[which(pheno.all$PAH %in% 0)], na.rm = T), lty=2)
abline(v=c(4.5,8.5,12.5), lty=3)
legend("bottomright", legend=c("Control females","Control males","PAH females","PAH males"), col=2:5, pch=15)

# Lasso with residuals ----------------------------------------------------
load("D:/RNAseq/PAXgene samples/lasso genes.RData")
lasso.res.tpm=as.matrix(cbind(simplepheno2[groupAB,c('age','sex')],t(res.tpm[lasso.genes,groupAB])))
lasso.res.tpm[,'sex']=as.numeric(as.factor(lasso.res.tpm[,'sex']))-1
lasso.res.tpm2=apply(lasso.res.tpm,1:2,as.numeric)

lasso.res.tpm.fit=glmnet(x = lasso.res.tpm2, y=lasso.y)

plot(lasso.res.tpm.fit, label = T)
##cross-validation fit
cvfit.res.tpm <- glmnet::cv.glmnet(x = lasso.res.tpm2, y = lasso.y)

save(cvfit.res.tpm,file="cvfit.res.tpm.RData")

## predict (residuals) in C

lasso.res.tpmC=as.matrix(cbind(simplepheno2[groupC,c('age','sex')],t(res.tpm[lasso.genes,groupC])))
lasso.res.tpmC[,'sex']=as.numeric(as.factor(lasso.res.tpmC[,'sex']))-1

lasso.res.tpmC2=apply(lasso.res.tpmC,1:2,as.numeric)
lasso.res.tpm.predictionC=predict(cvfit.res.tpm, newx = lasso.res.tpmC2, s = "lambda.1se")

boxplot(lasso.res.tpm.predictionC~simplepheno2[groupC,1])

library(pROC)
lassoROC.res.tpm=roc(simplepheno3$PAH[groupC],lasso.res.tpm.predictionC[,1])
plot(lassoROC.res.tpm)
plot(lasso.res.tpm.predictionC~simplepheno3$PAH[groupC])
wilcox.test(lasso.res.tpm.predictionC~simplepheno3$PAH[groupC])
lassoROC.res.tpm
ci.auc(lassoROC.res.tpm)



# Lasso with z-scored to PAH residuals ----------------------------------------------------
load("D:/RNAseq/PAXgene samples/lasso genes.RData")
lasso.z.res.tpm=as.matrix(cbind(simplepheno2[groupAB,c('age','sex')],t(z.res.tpm[lasso.genes,groupAB])))
lasso.z.res.tpm[,'sex']=as.numeric(as.factor(lasso.z.res.tpm[,'sex']))-1
lasso.z.res.tpm2=apply(lasso.z.res.tpm,1:2,as.numeric)

lasso.z.res.tpm.fit=glmnet(x = lasso.z.res.tpm2, y=lasso.y)

plot(lasso.z.res.tpm.fit, label = T)
##cross-validation fit
cvfit.z.res.tpm <- glmnet::cv.glmnet(x = lasso.z.res.tpm2, y = lasso.y)

save(cvfit.z.res.tpm,file="cvfit.z.res.tpm.RData")

## predict in ABC

paxgene.zres.lasso.genes.table <- read.csv("C:/Users/po3217/Desktop/Pablo/PhD/Data/RNAseq/PAXgene samples/2019-09-16 - genes z res lasso model PAH vs disease controls in Paxgene RNAseq.csv", header = FALSE)
paxgene.zres.lasso.genes <- as.character(paxgene.zres.lasso.genes.table[, 1])

lasso.z.res.tpmABC            <- as.matrix(cbind(simplepheno2[groupABC, c("PC1", "PC2", "PC3", 'age','sex')], t(z.res.tpm[well.detected.tophits$gene,groupABC])))
lasso.z.res.tpmABC[,'sex']    <- as.numeric(as.factor(lasso.z.res.tpmABC[,'sex']))-1

lasso.z.res.tpmABC2           <- apply(lasso.z.res.tpmABC, 1:2, as.numeric)

colnames(lasso.z.res.tpmABC2)[1:2] <- c("Age.at.sample", "Sex")

lasso.z.res.tpm.predictionABC <- predict(cvfit, newx = data.matrix(lasso.z.res.tpmABC2), s = "lambda.1se")

boxplot(lasso.z.res.tpm.predictionABC ~ simplepheno2[groupABC,1])

library(pROC)
lassoROC.z.res.tpm=roc(simplepheno3$PAH[groupC],lasso.z.res.tpm.predictionC[,1])
plot(lassoROC.z.res.tpm)
plot(lasso.z.res.tpm.predictionC~simplepheno3$PAH[groupC])
wilcox.test(lasso.z.res.tpm.predictionC~simplepheno3$PAH[groupC])
lassoROC.z.res.tpm
ci.auc(lassoROC.z.res.tpm)



# Mann-Whitney test tpm ~ sex ---------------------------------------------

load("C:/Users/po3217/Desktop/Pablo/PhD/Data/RNAseq/Deconvoluted analysis/19 06 10 - New analysis/190608 - salmon RNAseq files up to edgeR.RData")

pheno.all      <- read.csv("C:/Users/po3217/Desktop/Pablo/PhD/Data/RNAseq/Deconvoluted analysis/19 06 10 - New analysis/2019-09-04 NEW pheno all with deconvoluted lasso model.csv")


all.tpm            <- txi$abundance
all.tpm            <- t(all.tpm)
all.tpm            <- as.data.frame(all.tpm)
all.tpm$SampleID   <- rownames(all.tpm)
all.tpm$SampleID   <- substring(all.tpm$SampleID, 2)

pheno.all.mwt      <- merge(pheno.all[, c("SampleID", "sex", "Random.group2")], all.tpm, by = "SampleID", all.x = TRUE)

pheno.female       <- pheno.all.mwt[which(pheno.all.mwt$sex == "female" & pheno.all.mwt$Random.group2 %in% c("A", "B", "C")),]
pheno.male         <- pheno.all.mwt[which(pheno.all.mwt$sex == "male"   & pheno.all.mwt$Random.group2 %in% c("A", "B", "C")),]

mw.genes <- colnames(pheno.all.mwt[, 4:56520])


mw.test <- data.frame(lapply(X = 4:56520, FUN = function(column){
  wilcox.test(x = pheno.female[, column], y = pheno.male[, column])$p.value
} ))

colnames(mw.test) <- mw.genes

mw.test           <- as.data.frame(t(mw.test))



# Summary table of top hit genes ------------------------------------------

MW.lassogenes.PAH=do.call(rbind,apply(all.tpm[lasso.genes,], 1, function(gene){
  data.frame(
    "Controls median A"=median(gene[groupA][which(pheno.all$PAH[groupA] == 0)]),
    "Controls IQR A"=IQR(gene[groupA][which(pheno.all$PAH[groupA] == 0)]),
    "ihPAH median A"=median(gene[groupA][which(pheno.all$PAH[groupA] == 1)]),
    "ihPAH IQR A"=IQR(gene[groupA][which(pheno.all$PAH[groupA] == 1)]),
    "sig. A"=wilcox.test(gene[groupA]~pheno.all$PAH[groupA])$p.value,

    "Controls median B"=median(gene[groupB][which(pheno.all$PAH[groupB] == 0)]),
    "Controls IQR B"=IQR(gene[groupB][which(pheno.all$PAH[groupB] == 0)]),
    "ihPAH median B"=median(gene[groupB][which(pheno.all$PAH[groupB] == 1)]),
    "ihPAH IQR B"=IQR(gene[groupB][which(pheno.all$PAH[groupB] == 1)]),
    "sig. B"=wilcox.test(gene[groupB]~pheno.all$PAH[groupB])$p.value,

    "Controls median C"=median(gene[groupC][which(pheno.all$PAH[groupC] == 0)]),
    "Controls IQR C"=IQR(gene[groupC][which(pheno.all$PAH[groupC] == 0)]),
    "ihPAH median C"=median(gene[groupC][which(pheno.all$PAH[groupC] == 1)]),
    "ihPAH IQR C"=IQR(gene[groupC][which(pheno.all$PAH[groupC] == 1)]),
    "sig. C"=wilcox.test(gene[groupC]~pheno.all$PAH[groupC])$p.value,

    "Controls median"=median(gene[which(pheno.all$PAH == 0)]),
    "Controls IQR"=IQR(gene[which(pheno.all$PAH == 0)]),
    "ihPAH median"=median(gene[which(pheno.all$PAH == 1)]),
    "ihPAH IQR"=IQR(gene[which(pheno.all$PAH == 1)]),
    "sig."=wilcox.test(gene~pheno.all$PAH)$p.value
  )
})
)
MW.lassogenes.PAH$p.FDR=p.adjust(MW.lassogenes.PAH$sig., method = "fdr")
write.csv(MW.lassogenes.PAH,file = "190509 - Mann whitney results lasso genes in Controls v ihPAH.csv")


# Logistic regression model using representative RNAs ---------------------
simplepheno       <- pheno.all[,c("PAH","lane","PC1","PC2","PC3","PC4")] # , "ciber.T.cells.CD4.naive", "ciber.B.cells.memory", "ciber.Mast.cells.resting", "ciber.Dendritic.cells.resting", "quant.Tregs", "quant.T.cells.CD4")]
#
simplepheno$PAH   <- as.factor(simplepheno$PAH)
simplepheno$lane  <- as.factor(simplepheno$lane)
#
simplepheno2      <- simplepheno
simplepheno2$age  <- pheno.all$Age_controls_sample_patients_diagnosis
simplepheno2$sex  <- pheno.all$sex
#
 simplepheno3     <- simplepheno2[,c(1,3:5,7:8)]
 simplepheno3     <- cbind(simplepheno3,t(all.reads[as.character(deconvoluted.genes[,1]),]))
#
# modelAB15RNAs <- glm(PAH ~ .,family=binomial(link='logit'),
#                      data=simplepheno3[groupAB,]
# )
# summary(modelAB15RNAs)
#
# fitted.resultsAB15RNAs <- predict(modelAB15RNAs,newdata=simplepheno3[groupC,])
# library(pROC)
#
# AB15RNAsROC=roc(simplepheno3$PAH[groupC],fitted.resultsAB15RNAs)
# table(pheno.all$sex[groupC])
# plot(AB15RNAsROC)
# plot(fitted.resultsAB15RNAs~simplepheno3$PAH[groupC])
# wilcox.test(fitted.resultsAB15RNAs~simplepheno3$PAH[groupC])
# summary(AB15RNAsROC)
# ci.auc(AB15RNAsROC)
# groupABC=which(pheno.all$Random.group2 %in% c("A","B","C") & is.na(pheno.all$Exclude))
# pheno.all$log.res.model=NA
# pheno.all$log.res.model[groupABC]=predict(modelAB15RNAs,newdata=simplepheno3[groupABC,])

# Survival ----------------------------------------------------------------

library(survival)

survobj   <- Surv(time = pheno.all$Years.since.diagnosis.tested,
                time2 = (pheno.all$Years.since.diagnosis.tested+pheno.all$Years.survived.since.testing),
                event =  pheno.all$Died, type = "counting")

followup  <- pheno.all$Years.since.diagnosis.tested+pheno.all$Years.survived.since.testing


survival.analysis <- function(conditions,xmax){
  survobj  <- Surv(time = pheno.all$Years.since.diagnosis.tested,
                  time2 = (pheno.all$Years.since.diagnosis.tested+pheno.all$Years.survived.since.testing),
                  event =  pheno.all$Died, type = "counting")
  followup <- pheno.all$Years.since.diagnosis.tested+pheno.all$Years.survived.since.testing

  model.limit <- coxph(survobj[which(conditions)]~
                      pheno.all$rs2856830[which(conditions)])

  model.limit.summary <- summary(model.limit)
  print(model.limit.summary)


  plot(survfit(survobj[which(conditions)]~
                 pheno.all$rs2856830[which(conditions)]),
       xmax = xmax,
       #ymin = 0.3,
       col = 1:3, xlab = "Years survived since diagnosis", ylab = "Cumulative survival", axes = F, lty=3:1, lwd=2, cex=2)
  axis(1, cex=2)
  axis(2)

  legend(legend = levels(as.factor(pheno.all$rs2856830))[3:1], bty="n", x = 0.5, y = 0.3, col = 3:1,
         lty = 1:3, title = "rs2856830 C alleles", lwd=2)
  legend(legend = paste0("p = ",round(model.limit.summary$sctest[3][[1]], digits = 3)), bty="n", x = 7, y = 0.1)

  print(survfit(survobj[which(conditions)]~
                  pheno.all$rs2856830[which(conditions)]))
  print(summary(survfit(survobj[which(conditions)]~
                          pheno.all$rs2856830[which(conditions)]), times = c(0,1,2,3,4,5,6,8,10,12,15)
  ))


  model.limit.agesex <- coxph(survobj[which(conditions)]~
                             pheno.all$rs2856830[which(conditions)] + pheno.all$Age_controls_sample_patients_diagnosis[which(conditions)] +
                             pheno.all$sex[which(conditions)])

  model.limit.agesex.summary=summary(model.limit.agesex)
  print(model.limit.agesex.summary)
}
survival.analysis(conditions = followup < 15)
conditions <- followup < 15

transcripts.salmon2[transcripts.salmon2$Gene.stable.ID == "ENSG00000223865",]

model.hladpb1=coxph(survobj[which(conditions)]~
                      all.reads["HLADPB1",which(conditions)])

model.hladpb1.summary=summary(model.hladpb1)
model.hladpb1.summary

model.hladpb1.snp=coxph(survobj[which(conditions)]~
                          all.reads["ENST00000418931.2",which(conditions)] + pheno.all$rs2856830[which(conditions)])

model.hladpb1.snp.summary=summary(model.hladpb1.snp)
model.hladpb1.snp.summary

library(pROC)
hladpb1.roc=roc(pheno.all$Died,all.reads["ENST00000418931.2",])
plot(hladpb1.roc)
df.hladpb1.roc=data.frame(hladpb1.roc$sensitivities,hladpb1.roc$specificities,hladpb1.roc$thresholds)
df.hladpb1.roc$yi=df.hladpb1.roc$hladpb1.roc.specificities+df.hladpb1.roc$hladpb1.roc.sensitivities
summary(roc(pheno.all$Died,all.reads["ENST00000418931.2",]))
model.hladpb1.agesex=coxph(survobj[which(conditions)]~
                             all.reads["ENST00000418931.2",which(conditions)] + pheno.all$Age_controls_sample_patients_diagnosis[which(conditions)] +
                             as.character(pheno.all$sex[which(conditions)]))
class(all.reads["ENST00000418931.2",1:10])
model.hladpb1.agesex.summary=summary(model.hladpb1.agesex)
model.hladpb1.agesex.summary
model.hladpb1cutoff=coxph(survobj[which(conditions)]~
                            all.reads["ENST00000418931.2",which(conditions)]>hladpb1.cutoff)

model.hladpb1cutoff.summary=summary(model.hladpb1cutoff)
model.hladpb1cutoff.summary

model.hladpb1cutoff.summary
hladpb1.cutoff=12.574347

plot(survfit(survobj[which(conditions)]~
               (all.reads["ENST00000418931.2",]>hladpb1.cutoff)[which(conditions)]),
     xmax = 12,
     #ymin = 0.3,
     col = 1:2, xlab = "Years survived since diagnosis", ylab = "Cumulative survival",
     axes = F, lty=2:1, lwd=2, cex=2)
axis(1, cex=2)
axis(2)

legend(legend = levels(as.factor((all.reads["ENST00000418931.2",]>hladpb1.cutoff)[which(conditions)]))[2:1], bty="n", x = 0.5, y = 0.3, col = 2:1,
       lty = 1:2, title = "hladpb1 above cutoff", lwd=2)
legend(legend = paste0("p = ",round(model.hladpb1cutoff.summary$sctest[3][[1]], digits = 3)), bty="n", x = 7, y = 0.1)


# Diagnostic model vs survival --------------------------------------------

library(survival)

survobj       <- Surv(time = pheno.all$Years.since.diagnosis.tested,
                  time2 = (pheno.all$Years.since.diagnosis.tested+pheno.all$Years.survived.since.testing),
                  event =  pheno.all$Died, type = "counting")

followup      <- pheno.all$Years.since.diagnosis.tested+pheno.all$Years.survived.since.testing


model.survroc <- roc(pheno.all$Died,pheno.all$lasso.decon)
plot(model.survroc)


df.model.survroc     <- data.frame(model.survroc$sensitivities,model.survroc$specificities,model.survroc$thresholds)
df.model.survroc$yi  <- df.model.survroc$model.survroc.specificities+df.model.survroc$model.survroc.sensitivities


df.model.survroc[189,]
#    model.survroc.sensitivities model.survroc.specificities model.survroc.thresholds       yi
#189                   0.6458333                   0.7307692                 1.909701 1.376603

ci.auc(model.survroc)
#95% CI: 0.6168-0.7904 (DeLong)

legend("bottomright",legend="95% CI: 0.6168-0.7904 (DeLong)")

model.survroc
#Area under the curve: 0.7036



wilcox.test(pheno.all$Died,pheno.all$lasso.decon)
# p-value < 2.2e-16

plot(pheno.all$lasso.decon[which(pheno.all$Died2 %in% c(0,1)& pheno.all$PAH %in% 0:1)]
     ~ as.factor(paste0(pheno.all$PAH[which(pheno.all$Died2 %in% c(0,1)& pheno.all$PAH %in% 0:1)],
                        pheno.all$Died2[which(pheno.all$Died2 %in% c(0,1)& pheno.all$PAH %in% 0:1)
     ])
      ), ylab = "Lasso diagnostic model", xlab = "          Controls                PAH survivors      PAH non-survivors",
     col=c("darkgreen","orange","red"), lty=1, bty="n")

abline(h=1.909701, lty=2)

pheno.all$lasso1.767813     <- pheno.all$lasso.decon > 1.767813
pheno.all$lasso1.909701     <- pheno.all$lasso.decon > 1.909701


plot(survfit(survobj[which(conditions)]~

               pheno.all$lasso.decon[which(conditions)] > median(pheno.all$lasso.decon[which(conditions)], na.rm = T)),

               #pheno.all$lasso1.909701[which(conditions)]),
     #xmax = 10,
     #ymin = 0.3,
     col = 1:4, xlab = "Years survived since diagnosis", ylab = "Cumulative survival",
     axes = F, lty=2:1, lwd=2, cex=2, xlim=c(0,12)#xmax = 12, xmin=0
     )
axis(1, cex=2)
axis(2)

table(pheno.all$lasso1.909701[which(conditions)])
# FALSE  TRUE
# 157    84

summary(survfit(survobj[which(conditions)]~
          pheno.all$lasso1.909701[which(conditions)]),
        #time=c(5,10,15,20),
        time=c(0,2,4,6,8,10,12)
        )

# #pheno.all$lasso1.909701[which(conditions)]=FALSE
# time n.risk n.event survival std.err lower 95% CI upper 95% CI
# 0     15       0    1.000  0.0000        1.000        1.000
# 2     36       3    0.925  0.0416        0.847        1.000
# 4     33       1    0.896  0.0494        0.804        0.998
# 6     25       2    0.837  0.0615        0.724        0.966
# 8     20       1    0.797  0.0703        0.670        0.947
# 10     17       0    0.797  0.0703        0.670        0.947
# 12     17       0    0.797  0.0703        0.670        0.947
#
# pheno.all$lasso1.909701[which(conditions)]=TRUE
# time n.risk n.event survival std.err lower 95% CI upper 95% CI
# 0     26       0    1.000  0.0000       1.0000        1.000
# 2     29       1    0.963  0.0363       0.8943        1.000
# 4     18      10    0.607  0.0920       0.4508        0.817
# 6     19       7    0.408  0.0873       0.2679        0.620
# 8      8       5    0.267  0.0777       0.1505        0.472
# 10      6       2    0.167  0.0748       0.0691        0.402
# 12      7       1    0.125  0.0667       0.0439        0.356


legend(legend = levels(as.factor(pheno.all$lasso1.909701[which(conditions)]))[2:1],
       bty="n", x = 0.5, y = 0.4, col = 2:1,
       lty = 1:2, title = "Diagnostic lasso model > cutoff", lwd=2)

model.lasso.cox          <- coxph(survobj[which(conditions)]~
                            pheno.all$lasso1.909701[which(conditions)])

model.lasso.cox.summary  <- summary(model.lasso.cox)
model.lasso.cox.summary

model.lasso.cox.summary$sctest[3][[1]]
#[1] 1.036858e-05

legend(legend = "p = 1.03 x 10-5", bty="n", x = 7, y = 0.1)

#legend(legend = paste0("p = ",round(model.lasso.cox.summary$sctest[3][[1]],
#                                    digits = 5)), bty="n", x = 7, y = 0.1)

survobj.fromsample<- Surv(time = pheno.all$Years.survived.since.testing,
                                   event =  pheno.all$Died)

plot(survfit(survobj.fromsample[which(conditions)]~
               pheno.all$lasso1.909701[which(conditions)]),
     #xmax = 10,
     #ymin = 0.3,
     col = 1:4, xlab = "Years survived since sample", ylab = "Cumulative survival",
     axes = F, lty=2:1, lwd=2, cex=2, xmax = 4)
axis(1, cex=2)
axis(2)

table(pheno.all$lasso1.909701[which(conditions)])
# FALSE  TRUE
# 18   223

summary(survfit(survobj.fromsample[which(conditions)]~
                  pheno.all$lasso1.909701[which(conditions)]),
        #time=c(5,10,15,20),
        time=c(0,2,4,6,8,10,12)
)

# pheno.all$lasso1.767813[which(conditions)]=FALSE
# time n.risk n.event survival std.err lower 95% CI upper 95% CI
# 0     18       0    1.000  0.0000        1.000            1
# 2      7       2    0.872  0.0858        0.719            1
# 4      4       0    0.872  0.0858        0.719            1
# 6      1       0    0.872  0.0858        0.719            1
# 8      1       0    0.872  0.0858        0.719            1
# 10      1       0    0.872  0.0858        0.719            1
#
# pheno.all$lasso1.767813[which(conditions)]=TRUE
# time n.risk n.event survival std.err lower 95% CI upper 95% CI
# 0    222       0    1.000  0.0000        1.000        1.000
# 2     98      15    0.900  0.0256        0.851        0.951
# 4     31      16    0.683  0.0549        0.584        0.800
# 6      1       4    0.406  0.1710        0.178        0.927
# 8      1       0    0.406  0.1710        0.178        0.927



legend(legend = levels(as.factor(pheno.all$lasso1.909701[which(conditions)]))[2:1],
       bty="n", x = 0.1, y = 0.4, col = 2:1,
       lty = 1:2, title = "Diagnostic lasso model > cutoff", lwd=2)


model.lasso.cox          <- coxph(survobj.fromsample[which(conditions)]~
                           pheno.all$lasso1.909701[which(conditions)])

model.lasso.cox.summary  <- summary(model.lasso.cox)
model.lasso.cox.summary


model.lasso.cox.summary$sctest[3][[1]]
#[1] 0.4662129

legend(legend = "p = 0.466", bty="n", x = 3, y = 0.1)



#legend(legend = paste0("p = ",round(model.lasso.cox.summary$sctest[3][[1]],
#                                    digits = 3)), bty="n", x = 7, y = 0.1)


pheno.all$lasso1.909701=pheno.all$lasso.decon>1.909701

plot(survfit(survobj.fromsample[which(conditions)]~
               pheno.all$lasso1.909701[which(conditions)]),
     #xmax = 10,
     #ymin = 0.3,
     col = 1:4, xlab = "Years survived since sample", ylab = "Cumulative survival",
     axes = F, lty=2:1, lwd=2, cex=2, xmax = 4)
axis(1, cex=2)
axis(2)

table(pheno.all$lasso1.909701[which(conditions)])

summary(survfit(survobj.fromsample[which(conditions)]~
                  pheno.all$lasso1.909701[which(conditions)]),
        #time=c(5,10,15,20),
        time=c(0,1,2,3,4,6,8,10,12)
)

legend(legend = levels(as.factor(pheno.all$lasso1.889752[which(conditions)]))[2:1],
       bty="n", x = 0.5, y = 0.5, col = 2:1,
       lty = 1:2, title = "Diagnostic lasso model > cutoff", lwd=2)



model.lasso.cox.fromsample=coxph(survobj.fromsample[which(conditions)]~
                        pheno.all$lasso1.909701[which(conditions)])

model.lasso.cox.fromsample.summary=summary(model.lasso.cox.fromsample)
model.lasso.cox.fromsample.summary


model.lasso.cox.fromsample.summary$sctest[3][[1]]
#[1] 4.666215e-06

legend(legend = "p = 4.67 x 10-6", bty="n", x = 7, y = 0.1)

#legend(legend = paste0("p = ",round(model.lasso.cox.fromsample.summary$sctest[3][[1]],
#                                    digits = 3)), bty="n", x = 3, y = 0.1)



# calculate quartiles of score --------------------------------------------
quantile(pheno.all$lassoModel[pheno.all$PAH==1],c(0.25,0.5,0.75), na.rm = T)
pheno.all$lassoModel.quartiles=NA


pheno.all$lassoModel.quartiles[pheno.all$lassoModel<
                                    quantile(pheno.all$lassoModel[pheno.all$PAH==1],
                                             0.25, na.rm = T)]<-1
pheno.all$lassoModel.quartiles[pheno.all$lassoModel>
                                    quantile(pheno.all$lassoModel[pheno.all$PAH==1],
                                             0.25, na.rm = T)]<-2
pheno.all$lassoModel.quartiles[pheno.all$lassoModel>
                                    quantile(pheno.all$lassoModel[pheno.all$PAH==1],
                                             0.5, na.rm = T)]<-3
pheno.all$lassoModel.quartiles[pheno.all$lassoModel>
                                    quantile(pheno.all$lassoModel[pheno.all$PAH==1],
                                             0.75, na.rm = T)]<-4
library(survival)
plot(survfit(survobj[which(conditions)]~
               pheno.all$lassoModel.quartiles[which(conditions)]),
     #xmax = 10,
     #ymin = 0.3,
     col = 1:4, xlab = "Years survived since diagnosis", ylab = "Cumulative survival",
     axes = F, lty=4:1, lwd=2, cex=2)
axis(1, cex=2)
axis(2)

legend(legend = levels(as.factor(pheno.all$lassoModel.quartiles[which(conditions)]))[4:1],
       bty="n", x = 0.5, y = 0.5, col = 4:1,
       lty = 1:4, title = "Quartiles of model", lwd=2)
model.quartiles.cox=coxph(survobj[which(conditions)]~
                            pheno.all$lassoModel.quartiles[which(conditions)])

model.quartiles.cox.summary=summary(model.quartiles.cox)
model.quartiles.cox.summary

legend(legend = paste0("p = ",round(model.quartiles.cox.summary$sctest[3][[1]],
                                    digits = 3)), bty="n", x = 7, y = 0.1)



# Normal range in controls ------------------------------------------------
pheno.all$log.res.model.normal=NA
pheno.all$log.res.model.normal[pheno.all$log.res.model<
                                 quantile(pheno.all$log.res.model[pheno.all$PAH==0],
                                          0.95, na.rm = T)]<-0
pheno.all$log.res.model.normal[pheno.all$log.res.model>
                                 quantile(pheno.all$log.res.model[pheno.all$PAH==0],
                                          0.95, na.rm = T)]<-1

table(pheno.all$log.res.model.normal,
      pheno.all$PAH
)
# survival analysis -------------------------------------------------------


library(survival)

survobj <- Surv(time = pheno.all$Years.since.diagnosis.tested,
                time2 = (pheno.all$Years.since.diagnosis.tested+pheno.all$Years.survived.since.testing),
                event =  pheno.all$Died, type = "counting")
followup=pheno.all$Years.since.diagnosis.tested+pheno.all$Years.survived.since.testing

plot(survfit(survobj[which(conditions)]~
               pheno.all$log.res.model.normal[which(conditions)]),
     xmax = 10,
     #ymin = 0.3,
     col = 1:2, xlab = "Years survived since diagnosis", ylab = "Cumulative survival",
     axes = F, lty=2:1, lwd=2, cex=2)
axis(1, cex=2)
axis(2)

legend(legend = levels(as.factor(pheno.all$log.res.model.normal[which(conditions)]))[2:1],
       bty="n", x = 0.5, y = 0.5, col = 2:1,
       lty = 1:2, title = "Outside of normal range", lwd=2)

legend(legend = paste0("p = ",round(model.cox.summary$sctest[3][[1]], digits = 3)), bty="n", x = 7, y = 0.1)


# cox analysis ------------------------------------------------------------

model.cox=coxph(survobj[which(conditions)]~
                  pheno.all$log.res.model.normal[which(conditions)])

model.cox.summary=summary(model.cox)
model.cox.summary
legend(legend = paste0("p = ",round(model.cox.summary$sctest[3][[1]], digits = 3)), bty="n", x = 7, y = 0.1)

model.cox.agesex=coxph(survobj[which(conditions)]~
                         pheno.all$log.res.model.normal[which(conditions)] +
                         pheno.all$Age_controls_sample_patients_diagnosis[which(conditions)] +
                         as.character(pheno.all$sex[which(conditions)]))

model.cox.agesex.summary=summary(model.cox.agesex)
model.cox.agesex.summary


# WGCNA analysis ----------------------------------------------------------

#install.packages("WGCNA")
#source("https://bioconductor.org/biocLite.R")
#biocLite("impute")
library(WGCNA)
allowWGCNAThreads()

options(stringsAsFactors = FALSE);

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(t(all.reads[which(transcripts.salmon2$max.detected.controls.PAH > 0.95),pheno.all$PAH %in% 0:1]), powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# One-step network construction and module detection
net = blockwiseModules(t(all.reads)[pheno.all$PAH %in% 0:1,
                                which(transcripts.salmon2$max.detected.controls.PAH > 0.95)],
                       power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "PAHrnaseq",
                       verbose = 3, maxBlockSize = 29000)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
transcripts.salmon2$modulesWDsft6=NA
transcripts.salmon2$modulesWDsft6[which(transcripts.salmon2$max.detected.controls.PAH > 0.95)]=moduleLabels

datExpr=t(all.reads)[pheno.all$PAH %in% 0:1,
                 which(transcripts.salmon2$max.detected.controls.PAH > 0.95)]
colnames(datExpr)<-as.character(transcripts.salmon2$Transcript.name[which(transcripts.salmon2$max.detected.controls.PAH > 0.95)])
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

datTraits=simplepheno3[pheno.all$PAH %in% 0:1,c(1,5,6)]
datTraits$sex=as.numeric(datTraits$sex)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "/n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# 3.b Gene relationship to trait and important modules: Gene Significance and Module
# Membership

# Define variable PAH containing the PAH column of datTrait
PAH = as.data.frame(datTraits$PAH);
names(PAH) = "PAH"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, PAH, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(PAH), sep="");
names(GSPvalue) = paste("p.GS.", names(PAH), sep="");

# 3.c Intramodular analysis: identifying genes with high GS and MM
module = "yellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for PAH",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
names(datExpr)[moduleColors=="yellow"]


# GO analysis -------------------------------------------------------------

##need to make allLLIDs a list of entrez IDs for genes LocusLinkIDs

library(biomaRt)
ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
map = getBM(c("hgnc_symbol","entrezgene"),
            filters=c("hgnc_symbol","with_entrezgene"),
            values=list(transcripts.salmon2$Gene.name, TRUE), mart=ensembl)
colnames(map)[1]<-"Gene.name"
transcripts.salmon22=merge(transcripts.salmon2,map,by="Gene.name", all.x=T)

modules.df=data.frame(numbers=moduleLabels,colours=moduleColors)
modules.df=modules.df[which(!(duplicated(modules.df$numbers))),]
write.csv(modules.df,file="180302 - modules and colours.csv")

GOmodules=transcripts.salmon22$modulesWDsft6[which(!(is.na(transcripts.salmon22$modulesWDsft6) |
                                                is.na(transcripts.salmon22$entrezgene)))]
allLLIDs=transcripts.salmon22$entrezgene[which(!(is.na(transcripts.salmon22$modulesWDsft6) |
                                            is.na(transcripts.salmon22$entrezgene)))]

# Run GO enrichment -------------------------------------------------------


GOenr = GOenrichmentAnalysis(GOmodules,
                             allLLIDs,
                             organism = "human", nBestP = 10)
tab = GOenr$bestPTerms[[4]]$enrichment
names(tab)
write.table(tab, file = paste0(Sys.Date(),
                               " - RNAseq WGCNA analysis GOEnrichmentTable.csv"),
            sep = ",", quote = TRUE, row.names = FALSE)

keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Finally, display the enrichment table:
View(screenTab)


# 5. Network visualization using WGCNA functions -------------------------

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
# load(net$TOMFiles)


dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# plotTOM = TOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, PAH))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


# Exporting networks ------------------------------------------------------

#6.a Exporting to VisANT
#The package provides a convenient function for exporting the network to VisANT [1]. We illustrate a simple export
#of the full weighted network of a single module.
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Read in the annotation file
annot = read.csv(file = "GeneAnnotation.csv");
# Select module
module = "black";
# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            #probeToGene = data.frame(gene., annot$gene_symbol)
)
#Because the brown module is rather large, we can restrict the genes in the output to say the 30 top hub genes in the
#module:
nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            #probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol)
)


# Cytoscape ---------------------------------------------------------------

#6.b Exporting to Cytoscape
#Cytoscape [2] allows the user to input an edge file and a node file, allowing the user to specify for example the link
#weights and the node colors. Here we demonstrate the output of two modules, the red and brown ones, to Cytoscape.
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Read in the annotation file
annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
modules = c("yellow");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = transcripts.salmon2$rownames[match(modProbes, transcripts.salmon2$Gene.name)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];

dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);


# Survival analysis in A and B from all RNA species-------------------------------------------------------
write.csv(pheno.all,file="180605 - pheno all.csv")
pheno.all2=read.csv(file = "180605 - pheno all census updated.csv", row.names = 1)
library(survival)
library(pROC)

table(pheno.all$Random.group2, pheno.all$Died)
table(pheno.all2$Random.group2, pheno.all2$Died)
pheno.all=pheno.all2

ROCsurv.AB=apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                      pheno.all$Random.group2 %in% c("A","B")),],
                 MARGIN = 2, FUN = function(x){
                   roc=roc(pheno.all$Died[which(pheno.all$PAH == 1 &
                                                  pheno.all$Random.group2 %in% c("A","B"))]~x)
                   auc=roc$auc
                   sens.spec=roc$sensitivities+roc$specificities
                   if(auc>0.5){
                     best.cutoff=roc$thresholds[which(sens.spec == max(sens.spec))[1]]
                     best.cutoff.sens=roc$sensitivities[which(sens.spec == max(sens.spec))[1]]
                     best.cutoff.spec=roc$specificities[which(sens.spec == max(sens.spec))[1]]}else{
                       best.cutoff=roc$thresholds[which(sens.spec == min(sens.spec))[1]]
                       best.cutoff.sens=roc$sensitivities[which(sens.spec == min(sens.spec))[1]]
                       best.cutoff.spec=roc$specificities[which(sens.spec == min(sens.spec))[1]]
                     }
                   #rocci=ci(roc, of="auc")
                   c(auc,
                     #rocci[1:2],
                     best.cutoff,best.cutoff.sens,best.cutoff.spec)
                 })

ROCsurv.AB.df=as.data.frame(t(ROCsurv.AB))
colnames(ROCsurv.AB.df)<-c("AUC",
                           #"95% lower","95% upper",
                           "Best cut-off","Sensitivity","Specificity")
ROCsurv.AB.df$negAUC=1-ROCsurv.AB.df$AUC
ROCsurv.AB.df$maxAUC=apply(ROCsurv.AB.df[,c(1,#7
                                            5)], MARGIN = 1, FUN = max)
ROCsurv.AB.df$Sig.=apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                              pheno.all$Random.group2 %in% c("A","B")),],
                         MARGIN = 2, FUN = function(gene){
                           test=wilcox.test(gene ~ pheno.all$Died[which(pheno.all$PAH == 1 &
                                                                       pheno.all$Random.group2 %in% c("A","B"))])
                           test$p.value
                           })

ROCsurv.AB.df$p.FDR=p.adjust(p = ROCsurv.AB.df$Sig., method = "fdr")
ROCsurv.AB.df.sig=ROCsurv.AB.df[which(ROCsurv.AB.df$p.FDR<0.05),]

length(which(transcripts.salmon2$max.detected.controls.PAH>=0.95))

ROCsurv.AB.df$logres=apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                                pheno.all$Random.group2 %in% c("A","B")),],
                           MARGIN = 2, FUN = function(x){
                             test=glm(pheno.all$Died[which(pheno.all$PAH == 1 &
                                                             pheno.all$Random.group2 %in% c("A","B"))]~x  +
                                        pheno.all$PC1[which(pheno.all$PAH == 1 &
                                                              pheno.all$Random.group2 %in% c("A","B"))] +
                                        pheno.all$PC2[which(pheno.all$PAH == 1 &
                                                              pheno.all$Random.group2 %in% c("A","B"))],
                                      family=binomial(link='logit'))
                             summary(test)$coefficients[2,4]})
ROCsurv.AB.df$logresp.FDR=p.adjust(p = ROCsurv.AB.df$logres, method = "fdr")


# In A --------------------------------------------------------------------
ROCsurv.A=apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                          pheno.all$Random.group2 %in% c("A")),],
                 MARGIN = 2, FUN = function(x){
                   roc=roc(pheno.all$Died[which(pheno.all$PAH == 1 &
                                                  pheno.all$Random.group2 %in% c("A"))]~x)
                   auc=roc$auc
                   sens.spec=roc$sensitivities+roc$specificities
                   if(auc>0.5){
                     best.cutoff=roc$thresholds[which(sens.spec == max(sens.spec))[1]]
                     best.cutoff.sens=roc$sensitivities[which(sens.spec == max(sens.spec))[1]]
                     best.cutoff.spec=roc$specificities[which(sens.spec == max(sens.spec))[1]]}else{
                       best.cutoff=roc$thresholds[which(sens.spec == min(sens.spec))[1]]
                       best.cutoff.sens=roc$sensitivities[which(sens.spec == min(sens.spec))[1]]
                       best.cutoff.spec=roc$specificities[which(sens.spec == min(sens.spec))[1]]
                     }
                   #rocci=ci(roc, of="auc")
                   c(auc,
                     #rocci[1:2],
                     best.cutoff,best.cutoff.sens,best.cutoff.spec)
                 })

ROCsurv.A.df=as.data.frame(t(ROCsurv.A))
colnames(ROCsurv.A.df)<-c("AUC",
                           #"95% lower","95% upper",
                           "Best cut-off","Sensitivity","Specificity")
ROCsurv.A.df$negAUC=1-ROCsurv.A.df$AUC
ROCsurv.A.df$maxAUC=apply(ROCsurv.A.df[,c(1,#7
                                            5)], MARGIN = 1, FUN = max)
ROCsurv.A.df$Sig.=apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                                  pheno.all$Random.group2 %in% c("A")),],
                         MARGIN = 2, FUN = function(gene){
                           test=wilcox.test(gene ~ pheno.all$Died[which(pheno.all$PAH == 1 &
                                                                          pheno.all$Random.group2 %in% c("A"))])
                           test$p.value
                         })

ROCsurv.A.df$p.FDR=p.adjust(p = ROCsurv.A.df$Sig., method = "fdr")
ROCsurv.A.df.sig=ROCsurv.A.df[which(ROCsurv.A.df$p.FDR<0.05),]

length(which(transcripts.salmon2$max.detected.controls.PAH>=0.95))

ROCsurv.A.df$logres=apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                                    pheno.all$Random.group2 %in% c("A")),],
                           MARGIN = 2, FUN = function(x){
                             test=glm(pheno.all$Died[which(pheno.all$PAH == 1 &
                                                             pheno.all$Random.group2 %in% c("A"))]~x  +
                                        pheno.all$PC1[which(pheno.all$PAH == 1 &
                                                              pheno.all$Random.group2 %in% c("A"))] +
                                        pheno.all$PC2[which(pheno.all$PAH == 1 &
                                                              pheno.all$Random.group2 %in% c("A"))]+
                                        pheno.all$Age_controls_sample_patients_diagnosis[which(pheno.all$PAH == 1 &
                                                              pheno.all$Random.group2 %in% c("A"))]+
                                        pheno.all$sex[which(pheno.all$PAH == 1 &
                                                              pheno.all$Random.group2 %in% c("A"))],
                                      family=binomial(link='logit'))
                             summary(test)$coefficients[2,4]})
ROCsurv.A.df$logresp.FDR=p.adjust(p = ROCsurv.A.df$logres, method = "fdr")


# In B --------------------------------------------------------------------


ROCsurv.B=apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                          pheno.all$Random.group2 %in% c("B")),],
                 MARGIN = 2, FUN = function(x){
                   roc=roc(pheno.all$Died[which(pheno.all$PAH == 1 &
                                                  pheno.all$Random.group2 %in% c("B"))]~x)
                   auc=roc$auc
                   sens.spec=roc$sensitivities+roc$specificities
                   if(auc>0.5){
                     best.cutoff=roc$thresholds[which(sens.spec == max(sens.spec))[1]]
                     best.cutoff.sens=roc$sensitivities[which(sens.spec == max(sens.spec))[1]]
                     best.cutoff.spec=roc$specificities[which(sens.spec == max(sens.spec))[1]]}else{
                       best.cutoff=roc$thresholds[which(sens.spec == min(sens.spec))[1]]
                       best.cutoff.sens=roc$sensitivities[which(sens.spec == min(sens.spec))[1]]
                       best.cutoff.spec=roc$specificities[which(sens.spec == min(sens.spec))[1]]
                     }
                   #rocci=ci(roc, of="auc")
                   c(auc,
                     #rocci[1:2],
                     best.cutoff,best.cutoff.sens,best.cutoff.spec)
                 })

ROCsurv.B.df=as.data.frame(t(ROCsurv.B))
colnames(ROCsurv.B.df)<-c("AUC",
                           #"95% lower","95% upper",
                           "Best cut-off","Sensitivity","Specificity")
ROCsurv.B.df$negAUC=1-ROCsurv.B.df$AUC
ROCsurv.B.df$maxAUC=apply(ROCsurv.B.df[,c(1,#7
                                            5)], MARGIN = 1, FUN = max)
ROCsurv.B.df$Sig.=apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                                  pheno.all$Random.group2 %in% c("B")),],
                         MARGIN = 2, FUN = function(gene){
                           test=wilcox.test(gene ~ pheno.all$Died[which(pheno.all$PAH == 1 &
                                                                          pheno.all$Random.group2 %in% c("B"))])
                           test$p.value
                         })

ROCsurv.B.df$p.FDR=p.adjust(p = ROCsurv.B.df$Sig., method = "fdr")
ROCsurv.B.df.sig=ROCsurv.B.df[which(ROCsurv.B.df$p.FDR<0.05),]

length(which(transcripts.salmon2$max.detected.controls.PAH>=0.95))

ROCsurv.B.df$logres=apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                                    pheno.all$Random.group2 %in% c("B")),],
                           MARGIN = 2, FUN = function(x){
                             test=glm(pheno.all$Died[which(pheno.all$PAH == 1 &
                                                             pheno.all$Random.group2 %in% c("B"))]~x  +
                                        pheno.all$PC1[which(pheno.all$PAH == 1 &
                                                              pheno.all$Random.group2 %in% c("B"))] +
                                        pheno.all$PC2[which(pheno.all$PAH == 1 &
                                                              pheno.all$Random.group2 %in% c("B"))]+
                                        pheno.all$Age_controls_sample_patients_diagnosis[which(pheno.all$PAH == 1 &
                                                                                                 pheno.all$Random.group2 %in% c("B"))]+
                                        pheno.all$sex[which(pheno.all$PAH == 1 &
                                                              pheno.all$Random.group2 %in% c("B"))],
                                      family=binomial(link='logit'))
                             summary(test)$coefficients[2,4]})
ROCsurv.B.df$logresp.FDR=p.adjust(p = ROCsurv.B.df$logres, method = "fdr")

# In C --------------------------------------------------------------------

ROCsurv.C=apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                     pheno.all$Random.group2 %in% c("C")),],
                MARGIN = 2, FUN = function(x){
                  roc=roc(pheno.all$Died[which(pheno.all$PAH == 1 &
                                                 pheno.all$Random.group2 %in% c("C"))]~x)
                  auc=roc$auc
                  sens.spec=roc$sensitivities+roc$specificities
                  if(auc>0.5){
                    best.cutoff=roc$thresholds[which(sens.spec == max(sens.spec))[1]]
                    best.cutoff.sens=roc$sensitivities[which(sens.spec == max(sens.spec))[1]]
                    best.cutoff.spec=roc$specificities[which(sens.spec == max(sens.spec))[1]]}else{
                      best.cutoff=roc$thresholds[which(sens.spec == min(sens.spec))[1]]
                      best.cutoff.sens=roc$sensitivities[which(sens.spec == min(sens.spec))[1]]
                      best.cutoff.spec=roc$specificities[which(sens.spec == min(sens.spec))[1]]
                    }
                  #rocci=ci(roc, of="auc")
                  c(auc,
                    #rocci[1:2],
                    best.cutoff,best.cutoff.sens,best.cutoff.spec)
                })

ROCsurv.C.df=as.data.frame(t(ROCsurv.C))
colnames(ROCsurv.C.df)<-c("AUC",
                          #"95% lower","95% upper",
                          "Best cut-off","Sensitivity","Specificity")
ROCsurv.C.df$negAUC=1-ROCsurv.C.df$AUC
ROCsurv.C.df$maxAUC=apply(ROCsurv.C.df[,c(1,#7
                                          5)], MARGIN = 1, FUN = max)
ROCsurv.C.df$Sig.=apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                             pheno.all$Random.group2 %in% c("C")),],
                        MARGIN = 2, FUN = function(gene){
                          test=wilcox.test(gene ~ pheno.all$Died[which(pheno.all$PAH == 1 &
                                                                      pheno.all$Random.group2 %in% c("C"))])
                          test$p.value})




ROCsurv.C.df$logres=apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                               pheno.all$Random.group2 %in% c("C")),],
                          MARGIN = 2, FUN = function(x){
                            test=glm(pheno.all$Died[which(pheno.all$PAH == 1 &
                                                            pheno.all$Random.group2 %in% c("C"))]~x  +
                                       pheno.all$PC1[which(pheno.all$PAH == 1 &
                                                             pheno.all$Random.group2 %in% c("C"))] +
                                       pheno.all$PC2[which(pheno.all$PAH == 1 &
                                                             pheno.all$Random.group2 %in% c("C"))]+
                                       pheno.all$Age_controls_sample_patients_diagnosis[which(pheno.all$PAH == 1 &
                                                                                                pheno.all$Random.group2 %in% c("C"))]+
                                       pheno.all$sex[which(pheno.all$PAH == 1 &
                                                             pheno.all$Random.group2 %in% c("C"))],
                                     family=binomial(link='logit'))
                            summary(test)$coefficients[2,4]})
ROCsurv.C.df$logresp.FDR=p.adjust(p = ROCsurv.C.df$logres, method = "fdr")



# In  ABC -----------------------------------------------------------------

ROCsurv.ABC  <- apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                         pheno.all$Random.group2 %in% c("C","A","B")),],
                MARGIN = 2, FUN = function(x){
                  roc=roc(pheno.all$Died[which(pheno.all$PAH == 1 &
                                                 pheno.all$Random.group2 %in% c("C","A","B"))]~x)
                  auc=roc$auc
                  sens.spec=roc$sensitivities+roc$specificities
                  if(auc>0.5){
                    best.cutoff=roc$thresholds[which(sens.spec == max(sens.spec))[1]]
                    best.cutoff.sens=roc$sensitivities[which(sens.spec == max(sens.spec))[1]]
                    best.cutoff.spec=roc$specificities[which(sens.spec == max(sens.spec))[1]]}else{
                      best.cutoff=roc$thresholds[which(sens.spec == min(sens.spec))[1]]
                      best.cutoff.sens=roc$sensitivities[which(sens.spec == min(sens.spec))[1]]
                      best.cutoff.spec=roc$specificities[which(sens.spec == min(sens.spec))[1]]
                    }
                  #rocci=ci(roc, of="auc")
                  c(auc,
                    #rocci[1:2],
                    best.cutoff,best.cutoff.sens,best.cutoff.spec)
                })

ROCsurv.ABC.df              <- as.data.frame(t(ROCsurv.ABC))
colnames(ROCsurv.ABC.df)    <- c("AUC",
                                 #"95% lower","95% upper",
                                 "Best cut-off","Sensitivity","Specificity")
ROCsurv.ABC.df$negAUC       <- 1-ROCsurv.ABC.df$AUC
ROCsurv.ABC.df$maxAUC       <- apply(ROCsurv.ABC.df[,c(1,#7
                                          5)], MARGIN = 1, FUN = max)
ROCsurv.ABC.df$Sig.         <- apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                                 pheno.all$Random.group2 %in% c("C","A","B")),],
                        MARGIN = 2, FUN = function(gene){
                          test=wilcox.test(gene ~ pheno.all$Died[which(pheno.all$PAH == 1 &
                                                                         pheno.all$Random.group2 %in% c("C","A","B"))])
                          test$p.value})




ROCsurv.ABC.df$logres       <- apply(X = t(all.reads)[which(pheno.all$PAH == 1 &
                                                   pheno.all$Random.group2 %in% c("C","A","B")),],
                          MARGIN = 2, FUN = function(x){

                            test=  glm(pheno.all$Died[which(pheno.all$PAH == 1 &
                                                            pheno.all$Random.group2 %in% c("C","A","B"))]~x  +
                                       pheno.all$PC1[which(pheno.all$PAH == 1 &
                                                             pheno.all$Random.group2 %in% c("C","A","B"))] +
                                       pheno.all$PC2[which(pheno.all$PAH == 1 &
                                                             pheno.all$Random.group2 %in% c("C","A","B"))]+
                                       pheno.all$Age_controls_sample_patients_diagnosis[which(pheno.all$PAH == 1 &
                                                                                                pheno.all$Random.group2 %in% c("C","A","B"))]+
                                       pheno.all$sex[which(pheno.all$PAH == 1 &
                                                             pheno.all$Random.group2 %in% c("C","A","B"))],
                                     family=binomial(link='logit'))
                            summary(test)$coefficients[2,4]})

ROCsurv.ABC.df$logresp.FDR=p.adjust(p = ROCsurv.ABC.df$logres, method = "fdr")


# Diagnostic genes in survival --------------------------------------------

lasso.coef.1se.min              <- read.csv(file = "C:/Users/po3217/Desktop/Pablo/PhD/Data/RNAseq/Deconvoluted analysis/19 06 10 - New analysis/2019-08-29 - lasso coefficients 1se and min.csv")
lasso.model.1se                 <- lasso.coef.1se.min[which(!(lasso.coef.1se.min[, 2] == 0)), 2]
lasso.model.1se                 <- as.data.frame(lasso.model.1se)
rownames(lasso.model.1se)       <- as.character(lasso.coef.1se.min[which(!(lasso.coef.1se.min[, 2] == 0)), 1])


ROCsurv.ABC.diagnostic          <- ROCsurv.ABC.df[which(rownames(ROCsurv.ABC.df) %in% rownames(lasso.model.1se)),]
ROCsurv.ABC.diagnostic$p.FDR    <- p.adjust(ROCsurv.ABC.diagnostic$Sig., method = "fdr")

write.csv(ROCsurv.ABC.diagnostic, "survival analysis of diagnostic genes.csv")

# Bind results and pick top hits------------------------------------------------------------



ROCsurv.df.ABC=cbind(ROCsurv.AB.df,ROCsurv.C.df)

save(ROCsurv.df.ABC,file = "ROCsurv.df.ABC.RData")

ROCsurv.df.A.B.C=cbind(ROCsurv.A.df,ROCsurv.B.df,ROCsurv.C.df)
ROCsurv.df.A.B.C$gene=rownames(ROCsurv.df.A.B.C)
ROCsurv.df.A.B.C=merge(ROCsurv.df.A.B.C,transcripts.salmon2,by="gene")

ROCsurv.df.A.B.C.sigA=ROCsurv.df.A.B.C[which(ROCsurv.df.A.B.C[,which(colnames(ROCsurv.df.A.B.C)=="logres")[1]] < 0.05 &
                                               ROCsurv.df.A.B.C$max.detected.controls.PAH>0.95),]
ROCsurv.df.A.B.C.sigA$p.FDR.B=p.adjust(ROCsurv.df.A.B.C.sigA[,which(colnames(ROCsurv.df.A.B.C.sigA)=="logres")[2]], method = "fdr")
ROCsurv.df.A.B.C.sigA$p.bon.B=p.adjust(ROCsurv.df.A.B.C.sigA[,which(colnames(ROCsurv.df.A.B.C.sigA)=="logres")[2]], method = "bon")

ROCsurv.df.A.B.C.sigA.Bfdr=ROCsurv.df.A.B.C.sigA[which(ROCsurv.df.A.B.C.sigA$p.FDR.B < 0.05),]
ROCsurv.df.A.B.C.sigA.Bbon=ROCsurv.df.A.B.C.sigA[which(ROCsurv.df.A.B.C.sigA$p.bon.B < 0.05),]


plot(x = -log10(ROCsurv.df.A.B.C$logres), y = -log10(ROCsurv.df.A.B.C[,which(colnames(ROCsurv.df.A.B.C) == "logres")[2]]))
abline(v=-log10(0.05))


# lasso to select best from top hits in A and B ---------------------------
library(glmnet)

welldetect.topSurvhits.data=all.reads[ROCsurv.df.A.B.C.sigA.Bfdr$gene,]

#lassoSurv.phenotypes=simplepheno2[groupAB,c(3:4)]
lassoSurv.phenotypes=simplepheno2[groupAB,c(3:4,7:8)]
lassoSurv.data=cbind(lassoSurv.phenotypes,t(welldetect.topSurvhits.data[,groupAB]))
lassoSurv.data$sex=as.numeric(lassoSurv.data$sex)-1
lassoSurv.data2=as.matrix(lassoSurv.data)
lassoSurv.y=as.numeric(simplepheno2[groupAB,1])
lassoSurv.fit=glmnet(x = lassoSurv.data2, y=lassoSurv.y)
plot(lassoSurv.fit, label = T)
##cross-validation fit
cvfit.Surv <- glmnet::cv.glmnet(x = lassoSurv.data2, y = lassoSurv.y)

lassoSurv.coef.lambda.1se=coef(cvfit.Surv, s = "lambda.1se")
lassoSurv.coef.lambda.min=coef(cvfit.Surv, s = "lambda.min")

lassoSurv.coef.1se.min=as.matrix(cbind(lassoSurv.coef.lambda.1se,lassoSurv.coef.lambda.min))
write.csv(x = lassoSurv.coef.1se.min, file="lassoSurv coefficients 1se and min.csv")
plot(cvfit.Surv,)

##predict in C
#lassoSurv.phenotypesC=simplepheno2[groupC,c(3:4)]
lassoSurv.phenotypesC=simplepheno2[groupC,c(3:4,7:8)]
lassoSurv.dataC=cbind(lassoSurv.phenotypesC,t(welldetect.topSurvhits.data[,groupC]))
lassoSurv.dataC$sex=as.numeric(lassoSurv.dataC$sex)-1
lassoSurv.dataC2=as.matrix(lassoSurv.dataC)

lassoSurv.predictionC=predict(cvfit.Surv, newx = lassoSurv.dataC2, s = "lambda.1se")
library(pROC)
lassoSurvROC=roc(pheno.all$Died[groupC],lassoSurv.predictionC[,1])
table(pheno.all$sex[groupC])
plot(lassoSurvROC)
boxplot(lassoSurv.predictionC~pheno.all$Died[groupC],col=c("red","black"), xlab="Alive or died during follow-up",
        ylab= "Lasso survival score in C")
wilcox.test(lassoSurv.predictionC~pheno.all$Died[groupC])
lassoSurvROC
ci.auc(lassoSurvROC)
groupABC=which(pheno.all$Random.group2 %in% c("A","B","C") & is.na(pheno.all$Exclude))
pheno.all$lassoSurvModel=NA

lassoSurv.phenotypesABC=simplepheno2[groupABC,c(3:4,7:8)]
#lassoSurv.phenotypesABC=simplepheno2[groupABC,c(3:4)]
lassoSurv.dataABC=cbind(lassoSurv.phenotypesABC,t(welldetect.topSurvhits.data[,groupABC]))
lassoSurv.dataABC$sex=as.numeric(lassoSurv.dataABC$sex)-1
lassoSurv.dataABC2=as.matrix(lassoSurv.dataABC)

pheno.all$lassoSurvModel[groupABC]=predict(cvfit.Surv,newx=lassoSurv.dataABC2, s = "lambda.1se")
pheno.all$Died2=pheno.all$Died
pheno.all$Died2[pheno.all$PAH == 0]<-0



boxplot(pheno.all$lassoSurvModel[which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))] ~
          paste0(pheno.all$Random.group2,pheno.all$PAH, pheno.all$Died2)[which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))],
        col=2:5, ylab = "lassoSurv Model", lty=1
        #, ylim=c(1.65,1.9)
        )
abline(h=median(pheno.all$lassoSurvModel[which(pheno.all$PAH %in% 0)], na.rm = T), lty=2)
abline(v=c(4.5,8.5,12.5), lty=3)
legend("bottomright", legend=c("Controls","PAH alive","PAH died","PAH no follow-up"), col=2:5, pch=15)

boxplot(all.reads["TRPC1",which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))] ~
          paste0(pheno.all$Random.group2,pheno.all$PAH, pheno.all$Died2)[which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))],
        col=2:5, ylab = "TRPC1", lty=1
        #, ylim=c(1.65,1.9)
)
abline(h=median(all.reads["TRPC1",which(pheno.all$PAH %in% 0)], na.rm = T), lty=2)
abline(v=c(4.5,8.5,12.5), lty=3)
legend("bottomright", legend=c("Controls","PAH alive","PAH died","PAH no follow-up"), col=2:5, pch=15)



# results in AB and C -----------------------------------------------------




ROCsurv.df.ABC=ROCsurv.df.ABC[order(ROCsurv.df.ABC$Sig.),]

ROCsurv.df.ABC[1:10,]
ROCsurv.df.ABC[c("TIMP1","TIMP2","IL1RL1","APOE","IGFBP1","EPO","PLG","CFD","CFH"),]
ROCsurv.df.ABC[c("GDF2","BMPR2","ACVRL1","BMPR1A","HLA-DPB1","SOX17","CACNA2D2"),]

ROCsurv.df.ABC$p.Bon=ROCsurv.df.ABC$Sig.* 25966
ROCsurv.df.ABC$gene=rownames(ROCsurv.df.ABC)
ROCsurv.df.ABC=merge(ROCsurv.df.ABC,transcripts.salmon2,by="gene")
ROCsurv.df.ABC.sigAB=ROCsurv.df.ABC[which(ROCsurv.df.ABC$p.Bon<0.05),]
ROCsurv.df.ABC.sigAB=ROCsurv.df.ABC.sigAB[order(ROCsurv.df.ABC.sigAB[,which(colnames(ROCsurv.df.ABC.sigAB)=="Sig.")[2]]),]
ROCsurv.df.ABC.sigAB=ROCsurv.df.ABC.sigAB[which(ROCsurv.df.ABC.sigAB$max.detected.controls.PAH>0.95),]
ROCsurv.df.ABC.sigAB[1:10,]

ROCsurv.df.ABC.sigAB$p.FDR.C=p.adjust(p = ROCsurv.df.ABC.sigAB[,which(colnames(ROCsurv.df.ABC.sigAB)=="Sig.")[2]], method = "fdr")

boxplot(all.reads["C6orf52",which(pheno.all$PAH == 1 & pheno.all$Random.group2 %in% c("A","B"))]~
       pheno.all$Died[which(pheno.all$PAH == 1 & pheno.all$Random.group2 %in% c("A","B"))])
test=wilcox.test(t(all.reads)[which(pheno.all$PAH == 1 & pheno.all$Random.group2 %in% c("A","B")),"A3GALT2"]~
              pheno.all$Died[which(pheno.all$PAH == 1 & pheno.all$Random.group2 %in% c("A","B"))])
test$p.value


# KM plot -----------------------------------------------------------------

gene="TRPC1"
cutoff=ROCsurv.ABC.diagnostic$`Best cut-off`[which(rownames(ROCsurv.ABC.diagnostic)==gene)]

table(all.reads[gene,]>cutoff, pheno.all$Died)
conditions=pheno.all$Random.group2 %in% c("A","B","C")
plot(survfit(survobj[which(conditions)]~
               (all.reads[gene,]>cutoff)[which(conditions)]),
     xmax = 12,
     #ymin = 0.3,
     col = 1:2, xlab = "Years survived since diagnosis", ylab = "Cumulative survival",
     axes = F, lty=2:1, lwd=2, cex=2)
axis(1, cex=2)
axis(2)

legend(legend = levels(as.factor((all.reads[gene,]>cutoff)[which(conditions)]))[2:1], bty="n", x = 0.2, y = 0.3, col = 2:1,
       lty = 1:2, title = paste0(gene," above ",cutoff), lwd=2)

legend(legend = paste0("p = ",round(model.hladpb1cutoff.summary$sctest[3][[1]], digits = 3)), bty="n", x = 7, y = 0.1)
legend(legend= "n=165", y=0.4, x=9, bty="n")
legend(legend= "n=202", y=0.95, x=9, bty="n", col="red")

boxplot(all.reads['CACNA2D2',]~pheno.all$PAH, ylab= "CACNA2D2", xlab = "Controls and PAH")

conditions=pheno.all$Random.group2 %in% c("A","B","C") & pheno.all$PAH == 1

for (gene in rownames(ROCsurv.ABC.diagnostic)[which(ROCsurv.ABC.diagnostic$p.FDR < 0.05)]){

    cutoff=round(x = ROCsurv.ABC.diagnostic$`Best cut-off`[which(rownames(ROCsurv.ABC.diagnostic) == gene)], digits = 1)

    gene.model=summary(coxph(survobj.fromsample[which(conditions)]~
                            (all.reads[gene,]>cutoff)[which(conditions)]))
    pdf(paste0(gene, cutoff, " KM.pdf"), width = 5, height = 6)


    plot(survfit(survobj.fromsample[which(conditions)]~
                (all.reads[gene,] > cutoff)[which(conditions)]),
         xmax = 4,
         #ymin = 0.3,
         col = 1:2, xlab = "Years survived since sample", ylab = "Cumulative survival",
         axes = F, lty=2:1, lwd=2, cex=2)
    axis(1, cex=2)
    axis(2)

    legend(legend = levels(as.factor((all.reads[gene, ] > cutoff)[which(conditions)]))[2:1], bty="n", x = 0.2, y = 0.3, col = 2:1,
           lty = 1:2, title = paste0(gene," above ",cutoff), lwd=2)

    legend(legend = paste0("p = ",round(gene.model$sctest[3][[1]], digits = 3)), bty="n", x = 2.5, y = 0.15)

    nabove=length(which(all.reads[gene,conditions]>cutoff))
    nbelow=length(which(all.reads[gene,conditions]<=cutoff))

    legend(legend=paste0("above ",cutoff,", n=",nabove), bty="n", x= 0.2, y=0.1)
    legend(legend=paste0("below ",cutoff,", n=",nbelow), bty="n", x= 0.2, y=0.06)

    dev.off()
  }

  # Plot some data ----------------------------------------------------------

gene="TLR5"
gene="ZNF723"
gene="MAP3K7CL"
gene="ACTB"
gene="TSIX"
gene="KLF10"
gene="XKRX"
gene="DAP"
gene="XIST"
gene="ADGRG7"
gene="OTOF"
gene="GPR15"
gene="FAM132B"
gene="IL6ST"
gene="DNAJB4"
gene="HIF1A"
gene="DMD"
gene="SPSB1"
gene="SMAD5"
gene="TRPC1"
gene="TGFBI"
gene="RBMXL1"
gene="GPR55"
gene="ID1"
gene="ZEB2"
gene="FHIT"
gene="LCK"
gene="H19"
gene="TSPO"
gene="DNMT1"
gene="DNMT3A"
gene="TET2"
gene="SEMA3G"
gene="SNX29"

boxplot(log10(all.reads[gene,which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))]) ~
          paste0(pheno.all$Random.group2,pheno.all$PAH)[which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))],
        col=c("darkgreen","red"), ylab = gene, lty=1)
abline(h=median(log10(all.reads[gene,which(pheno.all$PAH %in% 0)])), lty=2)
abline(v=c(2.5,4.5,12.5), lty=3)
legend("bottomright", legend=c("Controls","PAH"), col=c("darkgreen","red"), pch=15)


boxplot(
  log10(
    all.reads[gene,which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))]
   )
  ~
          paste0(pheno.all$Random.group2,pheno.all$PAH, pheno.all$sex)[which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))],
        col=2:5, ylab = gene, lty=1)
abline(h=median(log10(all.reads[gene,which(pheno.all$PAH %in% 0)])), lty=2)
abline(v=c(4.5,8.5,12.5), lty=3)
legend("bottomright", legend=c("Control females","Control males","PAH females","PAH males"), col=2:5, pch=15)

boxplot(log10(all.reads[gene,which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))]) ~
          paste0(pheno.all$PAH, pheno.all$sex)[which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))],
        col=2:5, ylab = gene, lty=1)
abline(h=median(log10(all.reads[gene,which(pheno.all$PAH %in% 0)])), lty=2)
abline(v=c(4.5,8.5,12.5), lty=3)
legend("bottomright", legend=c("Control females","Control males","PAH females","PAH males"), col=2:5, pch=15)

boxplot(log10(all.reads[gene,which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))]) ~
          paste0(pheno.all$PAH, pheno.all$vasorespond.bin)[which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))],
        col=2:5, ylab = gene, lty=1)
abline(h=median(log10(all.reads[gene,which(pheno.all$PAH %in% 0)])), lty=2)
abline(v=c(4.5,8.5,12.5), lty=3)
legend("bottomright", legend=c("Controls","PAH nonresponders","PAH vasoresponders"), col=2:4, pch=15)



hemnes.genes=c("EPDR1", "DSG2", "SCD5", "LPAR6", "MGAT5", "RHOQ", "UCHL1", "ZNF652", "RALGPS2", "TPD52", "MKLN1", "RAPGEF2", "PIAS1")
length(which(hemnes.genes %in% transcripts.salmon2$gene))
length(which(hemnes.genes %in% well.detected.tophits$gene)) #0
hemnes.genes[which(!(hemnes.genes %in% transcripts.salmon2$gene))]
"MKLN1" %in% transcripts.salmon2$gene
# "P2RY5" "MKNL1"
# P2RY5 == LPAR6
# MKNL1 typo of MKLN1?
for(gene in hemnes.genes){
  boxplot(log10(all.reads[gene,which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))]) ~
            paste0(pheno.all$vasorespond.bin,pheno.all$PAH, pheno.all$sex)[which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))],
          col=2:7, ylab = gene, lty=1)
  legend("topleft", legend=c("Control females","Control males","PAH females","PAH males","Vasoresponder females","Vasoresponder males"), col=2:7, pch=15)
  abline(h=median(log10(all.reads[gene,which(pheno.all$PAH %in% 0)])), lty=2)

}


boxplot(log10(all.reads[gene,which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))]) ~
          pheno.all$PAH[which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))], col=2:5)

plot(log10(all.reads["ZEB2",which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))]) ~
       log10(all.reads["SOX17",which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))]),
       col=pheno.all$PAH[which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))]+1)

zeb2soxcor=cor.test(y = as.numeric(log10(all.reads["ZEB2",which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))]))[which(is.finite(log10(all.reads["SOX17",which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))])))], x=
          as.numeric(log10(all.reads["SOX17",which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))]))[which(is.finite(log10(all.reads["SOX17",which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))])))], use = "c", )

zeb2soxcor
# Files for genetic analysis ----------------------------------------------

wgsid.tpm=t(all.tpm[,which(!(is.na(pheno.all$WGS.ID)))])
rownames(wgsid.tpm)<-pheno.all$WGS.ID[which(!(is.na(pheno.all$WGS.ID)))]

wgsid.tpm=wgsid.tpm[which(!(rownames(wgsid.tpm) == "")),]
wgsid.tpm=wgsid.tpm[which(!(duplicated(rownames(wgsid.tpm)))),]
write.csv(wgsid.tpm,file="RNAseq tpm with wgsid.csv")


# Metabolites -------------------------------------------------------------

pheno.plusmetabs=read.csv(file = "180614 - pheno all with metabolites and diagnostic genes.csv", row.names = 1)
table(pheno.plusmetabs$SampleID == pheno.all$SampleID)

# Metabolites in HV IPAH model in Circ paper
# dehydroisoandrosterone sulfate (DHEA-S)
# methionine sulfone
# N1-methylinosine
# oleoylcarnitine
# palmitoylcholine
# sphingomyelin (d18:1/20:0, d16:1/22:0)*
#   X - 24513

table(pheno.plusmetabs$PAH )
model.metabolites <- glm(pheno.plusmetabs$PAH ~ pheno.plusmetabs$dehydroisoandrosterone.sulfate..DHEA.S. +
                           pheno.plusmetabs$methionine.sulfone +
                           pheno.plusmetabs$N1.methylinosine +
                           pheno.plusmetabs$oleoylcarnitine +
                           pheno.plusmetabs$palmitoylcholine +
                           pheno.plusmetabs$sphingomyelin..d18.1.20.0..d16.1.22.0.. +
                           pheno.plusmetabs$X...24513,
                         family=binomial(link='logit'),
                     )
model.metabolites.summ=summary(model.metabolites)
length(model.metabolites$fitted.values)
plot(model.metabolites$fitted.values ~ pheno.plusmetabs$PAH[which(pheno.plusmetabs$PAH %in% 0:1 &
                                                                    !(is.na(pheno.plusmetabs$methionine.sulfone)))])
pheno.plusmetabs$metabolite.model=NA
pheno.plusmetabs$metabolite.model[which(pheno.plusmetabs$PAH %in% 0:1 &
                                          !(is.na(pheno.plusmetabs$methionine.sulfone)))]<-model.metabolites$fitted.values


model.metab.RNA <- glm(pheno.plusmetabs$PAH ~ pheno.plusmetabs$metabolite.model +
                           pheno.plusmetabs$lassoModel,
                         family=binomial(link='logit'),
)
pheno.plusmetabs$metab.RNA.model=NA
pheno.plusmetabs$metab.RNA.model[which(pheno.plusmetabs$PAH %in% 0:1 &
                                          !(is.na(pheno.plusmetabs$methionine.sulfone)) &
                                         !(is.na(pheno.plusmetabs$lassoModel)))]<-model.metab.RNA$fitted.values




model.metab.RNA2              <- glm(pheno.all$PAH ~ pheno.metabolite.score$IPAH.HV.DA.score.7.metabs +
                                                     pheno.all$lasso.decon,
                                     family = binomial(link='logit'),
)

pheno.all$metab.RNA.model2    <- NA

pheno.all$metab.RNA.model2[which(pheno.all$PAH %in% 0:1 &
                                 !(is.na(pheno.metabolite.score$IPAH.HV.DA.score.7.metabs)) &
                                 !(is.na(pheno.all$lasso.decon)))]                              <-model.metab.RNA2$fitted.values


library(pROC)

metabROC=roc(pheno.plusmetabs$PAH,pheno.plusmetabs$metabolite.model)
plot(metabROC)
wilcox.test(pheno.plusmetabs$metabolite.model~pheno.plusmetabs$PAH)
metabROC
ci.auc(metabROC)

metab.RNAROC=roc(pheno.plusmetabs$PAH,pheno.plusmetabs$metab.RNA.model)
plot(metab.RNAROC)
wilcox.test(pheno.plusmetabs$metab.RNA.model~pheno.plusmetabs$PAH)
metab.RNAROC
ci.auc(metab.RNAROC)



metab7ROC        <- roc(pheno.all$PAH, pheno.metabolite.score$IPAH.HV.DA.score.7.metabs)

wilcox.test(pheno.metabolite.score$IPAH.HV.DA.score.7.metabs ~ pheno.all$PAH)
metab7ROC
ci.auc(metab7ROC)



metab7.RNA2ROC   <- roc(pheno.all$PAH, pheno.all$metab.RNA.model2)
plot(metab7ROC)
plot(metab7.RNA2ROC, add=T, col="red", lty=2)
ci.auc(metab7.RNA2ROC)
legend(x = 0.8, y = 0.25, legend = c("RNA plus metabolites, 95% CI: 0.9702-1",
                                     "Metabolites, 95% CI: 0.8602-0.9622"), col=c("red","black"), lty=2:1, bty="n")
legend("bottomright",legend="p=0.001753",bty="n")

roc.test(metab7ROC,metab7.RNA2ROC,method="delong")


#PRELIMINARY MODEL - CHRIS
# DeLong's test for two correlated ROC curves
#
# data:  metab7ROC and metab7.RNA2ROC
# Z = -3.3972, p-value = 0.0006807
# alternative hypothesis: true difference in AUC is not equal to 0
# sample estimates:
# AUC of roc1 AUC of roc2
# 0.9076875   0.9933837
#

#DECONVOLUTED MODEL - PABLO
# DeLong's test for two correlated ROC curves
#
# data:  metab7ROC and metab7.RNA2ROC
# Z = -3.1291, p-value = 0.001753
# alternative hypothesis: true difference in AUC is not equal to 0
# sample estimates:
# AUC of roc1 AUC of roc2
#   0.9076875   0.9867675



lasso.deconvoluted             <- read.csv("C:/Users/po3217/Desktop/Pablo/PhD/Data/RNAseq/Deconvoluted analysis/19 06 10 - New analysis/deconvoluted lasso scores group ABC.csv")
pheno.metabolite.score         <- read.csv("C:/Users/po3217/Desktop/Pablo/PhD/Data/RNAseq/Deconvoluted analysis/19 06 10 - New analysis/180614 - pheno all with metabolites and diagnostic genes.csv")



pheno.all$lasso.decon          <- lasso.deconvoluted$lasso.deconvoluted
pheno.all$metabolite.score     <- pheno.metabolite.score$IPAH.HV.DA.score.7.metabs

pheno.all$lasso.decon[pheno.all$lasso.decon == "#N/A"] <- NA
pheno.all$lasso.decon <- as.numeric(as.character(pheno.all$lasso.decon))

plot(pheno.all$lasso.decon[groupABC] ~ as.numeric(pheno.all$metabolite.score[groupABC]), col = pheno.all$PAH[groupABC]+1, #pheno.plusmetabs$vasorespond.bin,
     xlab= "7 metabolite discriminant score",
     ylab="RNA lasso model")
abline(h= 1.767813, lty=2)
abline(v=0, lty=2)
legend("bottomright",legend=c("HV","PAH"),col=c("black","red"),pch="o")

write.csv(pheno.plusmetabs, file="180618 - pheno with metabs and scores.csv")
pheno.plusmetabs=read.csv(file="180618 - pheno with metabs and scores.csv")


# SMAD5 snp ---------------------------------------------------------------

smad5snp=read.csv(file = "smad5snp.csv")
colnames(smad5snp)<-c("WGS.ID","rs4146187")

smad5data=data.frame(SMAD5reads=all.reads["SMAD5",], SMAD5tpm=all.tpm["SMAD5",], WGS.ID=pheno.all$WGS.ID, BMPR2=pheno.all$BMPR2, group=pheno.all$Description.1)

smad5data=merge(smad5data,smad5snp)



# plot by SNP -------------------------------------------------------------


#boxplot(smad5data$SMAD5reads ~ smad5data$rs4146187, xlab="rs4146187", ylab="SMAD5 reads", col=c(3,7,2))
boxplot(c(smad5data$SMAD5tpm,all.tpm["SMAD5",which(pheno.all$PAH==0)]) ~
          c(smad5data$rs4146187,rep(-1,times=length(which(pheno.all$PAH==0)))),
        xlab="rs4146187", ylab="SMAD5 tpm", col=c(3,7,2,6), lty=1, frame.plot=F, bty="L", axes=F, ylim=c(0,25))
axis(2)
abline(v=1.5,lty=3)
abline(h=ave(all.tpm["SMAD5",which(pheno.all$PAH==0)], FUN = median), lty=2)
#summary(aov(smad5data$SMAD5reads ~ smad5data$rs4146187))
summary(aov(smad5data$SMAD5tpm ~ smad5data$rs4146187))
legend("topright",legend = "p=1.6x10-8", bty="n")
legend(x = 0.5, y = 2, legend = "Controls", bty="n")
legend(x = 1.5, y = 2, legend = "PAH C/C", bty="n")
legend(x = 2.5, y = 2, legend = "PAH C/A", bty="n")
legend(x = 3.5, y = 2, legend = "PAH A/A", bty="n")


# plot by BMPR2 status ----------------------------------------------------


boxplot(as.numeric(c(smad5data$SMAD5tpm,all.tpm["SMAD5",which(pheno.all$PAH==0)])) ~
          c(as.numeric(as.factor(smad5data$BMPR2)),rep(-1,times=length(which(pheno.all$PAH==0)))),
        xlab="BMPR2", ylab="SMAD5 tpm", col=c(3,7,2), lty=1, frame.plot=F, bty="L", axes=F, ylim=c(0,25))
axis(2)
abline(v=1.5,lty=3)
abline(h=ave(all.tpm["SMAD5",which(pheno.all$PAH==0)], FUN = median), lty=2)
#summary(aov(smad5data$SMAD5reads ~ smad5data$rs4146187))
summary(aov(smad5data$SMAD5tpm ~ smad5data$rs4146187))
legend("topright",legend = "p=1.6x10-8", bty="n")
legend(x = 0.5, y = 2, legend = "Controls", bty="n")
legend(x = 1.5, y = 2, legend = "PAH non-carriers", bty="n")
legend(x = 2.5, y = 2, legend = "PAH BMPR2 carriers", bty="n")



# plot by group -----------------------------------------------------------

boxplot(smad5data$SMAD5tpm ~ smad5data$group,
        xlab="group", ylab="SMAD5 tpm", col=c(3,7,2),
        lty=1, frame.plot=F, bty="L",
        #axes=F,
        ylim=c(0,25)
        )

abline(h=ave(all.tpm["SMAD5",which(pheno.all$PAH==0)], FUN = median), lty=2)

boxplot(pheno.all$lassoModel ~ pheno.all$Description.1,
        xlab="group", ylab="RNA diagnostic score", col=c(3,7,2),
        lty=1, frame.plot=F, bty="L",
        #axes=F,
        #ylim=c(0,25)
)
abline(h=ave(pheno.all$lassoModel[which(pheno.all$PAH==0)], FUN = median), lty=2)

# Vasoresponders ----------------------------------------------------------
# Vasoresponders edgeR -------------------------------------------------------------------

library(edgeR)


PAHgroupA=which(pheno.all$Random.group2 == "A" & is.na(pheno.all$Exclude) & pheno.all$PAH == 1)
PAHgroupB=which(pheno.all$Random.group2 == "B" & is.na(pheno.all$Exclude) & pheno.all$PAH == 1)
PAHgroupAB=which(pheno.all$Random.group2 %in% c("A","B") & is.na(pheno.all$Exclude) & pheno.all$PAH == 1)
PAHgroupC=which(pheno.all$Random.group2 == "C" & is.na(pheno.all$Exclude) & pheno.all$PAH == 1)

#cts <- txi$counts
#normMat <- txi$length
#normMat <- normMat/exp(rowMeans(log(normMat)))
#o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
#y <- DGEList(cts)
#y$offset <- t(t(log(normMat)) + o)


y.vasoA=DGEList(counts = all.reads[,PAHgroupA], group = pheno.plusmetabs$vasorespond.bin[PAHgroupA])
y.vasoA$offset <- t(t(log(normMat[,PAHgroupA])) + o[PAHgroupA])

cpm.vasoA=cpm(y.vasoA$counts)
#rpkm=rpkm(y$counts, gene.length = transcripts.salmon2$Length)
plotMDS(y.vasoA)
keep.vasoA <- rowSums(cpm(y.vasoA)>2) >= 3
y.vasoA <- y.vasoA[keep.vasoA, , keep.lib.sizes=FALSE]
y.vasoA <- calcNormFactors(y.vasoA)


design.vasoA <- model.matrix(~pheno.plusmetabs$vasorespond.bin[PAHgroupA]
                       +pheno.plusmetabs$PC1[PAHgroupA]
                       +pheno.plusmetabs$PC2[PAHgroupA]
                       +pheno.plusmetabs$PC3[PAHgroupA]
)
y.vasoA <- estimateDisp(y.vasoA,design.vasoA)
fit.vasoA <- glmQLFit(y.vasoA,design.vasoA)
plotQLDisp(fit.vasoA)
plotBCV(y.vasoA)
qlf.vasoA <- glmQLFTest(fit.vasoA,coef=2)
tophits.vasoA=topTags(qlf.vasoA, n = 180224)
tophits.vasoA$table->tophits.vasoA.table
tophits.vasoA.table$gene=rownames(tophits.vasoA.table)
tophits.vasoA.tableA=merge(tophits.vasoA.table,transcripts.salmon2,by="gene")
write.csv(tophits.vasoA.tableA,file = paste0(Sys.Date()," - top hits from RNAseq analysis.vasoA.csv"))
summary(dt.vasoA <- decideTestsDGE(qlf.vasoA))
isDE.vasoA <- as.logical(dt.vasoA)
DEnames.vasoA <- rownames(y.vasoA)[isDE.vasoA]
plotSmear(qlf.vasoA, de.tags=DEnames.vasoA)
abline(h=c(-1,1), col="blue")

y.vasoB=DGEList(counts = all.reads[,PAHgroupB], group = pheno.plusmetabs$vasorespond.bin[PAHgroupB])
y.vasoB$offset <- t(t(log(normMat[,PAHgroupB])) + o[PAHgroupB])

cpm.vasoB=cpm(y.vasoB$counts)
rpkm.vasoB=rpkm(y.vasoB$counts, gene.length = transcripts.salmon2$Length)
plotMDS(y.vasoB)
keep.vasoB <- rowSums(cpm(y.vasoB)>2) >= 3
y.vasoB <- y.vasoB[keep.vasoB, , keep.lib.sizes=FALSE]
y.vasoB <- calcNormFactors(y.vasoB)


design.vasoB <- model.matrix(~pheno.plusmetabs$vasorespond.bin[PAHgroupB]
                        +pheno.plusmetabs$PC1[PAHgroupB]
                        +pheno.plusmetabs$PC2[PAHgroupB]
                        +pheno.plusmetabs$PC3[PAHgroupB]
)
y.vasoB <- estimateDisp(y.vasoB,design.vasoB)
fit.vasoB <- glmQLFit(y.vasoB,design.vasoB)
plotQLDisp(fit.vasoB)
plotBCV(y.vasoB)
qlf.vasoB <- glmQLFTest(fit.vasoB,coef=2)
tophits.vasoB=topTags(qlf.vasoB, n = 180224)
tophits.vasoB$table->tophits.table.vasoB
tophits.table.vasoB$gene=rownames(tophits.table.vasoB)
tophits.table.vasoB=merge(tophits.table.vasoB,transcripts.salmon2,by="gene")
write.csv(tophits.table.vasoB,file = paste0(Sys.Date()," - top hits from RNAseq analysis.vasoB.csv"))
summary(dt.vasoB <- decideTestsDGE(qlf.vasoB))
isDE.vasoB <- as.logical(dt.vasoB)
DEnames.vasoB <- rownames(y.vasoB)[isDE.vasoB]
plotSmear(qlf.vasoB, de.tags=DEnames.vasoB)
abline(h=c(-1,1), col="blue")


# Vasoresponders In new AB -------------------------------------------------------------------

y.vasoAB=DGEList(counts = all.reads[,PAHgroupAB], group = pheno.plusmetabs$vasorespond.bin[PAHgroupAB])
y.vasoAB$offset <- t(t(log(normMat[,PAHgroupAB])) + o[PAHgroupAB])

cpm.vasoAB=cpm(y.vasoAB$counts)
rpkm.vasoAB=rpkm(y.vasoAB$counts, gene.length = transcripts.salmon2$Length)
plotMDS(y.vasoAB)
keep.vasoAB <- rowSums(cpm(y.vasoAB)>2) >= 3
y.vasoAB <- y.vasoAB[keep.vasoAB, , keep.lib.sizes=FALSE]
y.vasoAB <- calcNormFactors(y.vasoAB)


design.vasoAB <- model.matrix(~pheno.plusmetabs$vasorespond.bin[PAHgroupAB]
                         +pheno.plusmetabs$PC1[PAHgroupAB]
                         +pheno.plusmetabs$PC2[PAHgroupAB]
                         +pheno.plusmetabs$PC3[PAHgroupAB]
                         +pheno.plusmetabs$Age_controls_sample_patients_diagnosis[PAHgroupAB]
                         +as.character(pheno.plusmetabs$sex[PAHgroupAB])
)
y.vasoAB <- estimateDisp(y.vasoAB,design.vasoAB)
fit.vasoAB <- glmQLFit(y.vasoAB,design.vasoAB)
plotQLDisp(fit.vasoAB)
plotBCV(y.vasoAB)
qlf.vasoAB <- glmQLFTest(fit.vasoAB,coef=2)
tophits.vasoAB=topTags(qlf.vasoAB, n = 180224)
tophits.vasoAB$table->tophits.table.vasoAB
tophits.table.vasoAB$gene=rownames(tophits.table.vasoAB)
tophits.table.vasoAB=merge(tophits.table.vasoAB,transcripts.salmon2,by="gene")
write.csv(tophits.table.vasoAB,file = paste0(Sys.Date()," - top hits from RNAseq analysis.vasoAB.csv"))
summary(dt.vasoAB <- decideTestsDGE(qlf.vasoAB))
isDE.vasoAB <- as.logical(dt.vasoAB)
DEnames.vasoAB <- rownames(y.vasoAB)[isDE.vasoAB]
plotSmear(qlf.vasoAB, de.tags=DEnames.vasoAB)
abline(h=c(-1,1), col="blue")

colnames(tophits.tableA)

#[1] "Geneid"           "logFC"            "logCPM"           "F"                "PValue"           "FDR"
#[7] "Length"           "Gene.description" "Gene.name"        "total.counts"

colnames(tophits.vasoA.tableA)[2:6]<-c("logFC_A",
                                 "logCPM_A",
                                 "F_A",
                                 "PValue_A",
                                 "FDR_A")
colnames(tophits.table.vasoB)[2:6]<-c("logFC_B",
                                 "logCPM_B",
                                 "F_B",
                                 "PValue_B",
                                 "FDR_B")
colnames(tophits.table.vasoAB)[2:6]<-c("logFC_AB",
                                  "logCPM_AB",
                                  "F_AB",
                                  "PValue_AB",
                                  "FDR_AB")
detected.vasoAB=tophits.table.vasoAB[tophits.table.vasoAB$max.detected.controls.PAH>0.95,]
detected.vasoAB$p.FDR=p.adjust(p = detected.vasoAB$PValue_AB, method = "fdr")
sig.detected.vasoAB=detected.vasoAB[detected.vasoAB$p.FDR<0.05,]

tophits.table.vasoABAB=merge(tophits.vasoA.tableA,tophits.table.vasoB)
tophits.table.vasoABAB=merge(tophits.table.vasoABAB,tophits.table.vasoAB)
write.csv(tophits.table.vasoABAB,file = paste0(Sys.Date()," - top hits from RNAseq analysis.vasoABAB.csv"))

sig.tophits.table.vasoABAB=tophits.table.vasoABAB[which(tophits.table.vasoABAB$PValue_A < 0.05 &
                                                tophits.table.vasoABAB$PValue_B < 0.05),]
sig.tophits.table.vasoABAB$same.direction=sig.tophits.table.vasoABAB$logFC_A*sig.tophits.table.vasoABAB$logFC_B>0

hemnes.table.vasoABAB=tophits.table.vasoABAB[which(tophits.table.vasoABAB$gene %in% hemnes.genes),]
write.csv(hemnes.table.vasoABAB,file="RNAseq vasoresponder results for hemnes genes.csv")
# Vasoresponders ROC analysis for each gene ----------------------------------------------
library(pROC)
sig.tophits.table.vasoABAB$ROC.auc=do.call(rbind,lapply(X = sig.tophits.table.vasoABAB$rownames, FUN = function(x){
  roc(pheno.plusmetabs$vasorespond.bin[PAHgroupAB],as.numeric(all.reads[as.character(x),PAHgroupAB]))$auc
}))


# Vasoresponders Select well detected and same direction ---------------------------------

well.detected.vaso.tophits=sig.tophits.table.vasoABAB[which(sig.tophits.table.vasoABAB$max.detected.controls.PAH>=0.95 &
                                                    sig.tophits.table.vasoABAB$same.direction == TRUE),]
welldetect.vaso.tophits.data=all.reads[well.detected.vaso.tophits$gene,]

write.csv(welldetect.vaso.tophits.data, file="data for well detected top hits Controls vs PAH in RNAseq.vaso.csv")

welldet.vaso.tophits.cor=cor(welldetect.vaso.tophits.data)

sig.detected.vasoAB.data=all.reads[sig.detected.vasoAB$gene,]


# Vasoresponders LASSO to select independent markers -------------------------------------

#install.packages("glmnet")
library(glmnet)
lasso.vaso.phenotypes=simplepheno2[PAHgroupAB,c(3:5,7:8)]

lasso.vaso.data=cbind(lasso.vaso.phenotypes,t(sig.detected.vasoAB.data[,PAHgroupAB]))
lasso.vaso.data$sex=as.numeric(lasso.vaso.data$sex)-1
lasso.vaso.data2=as.matrix(lasso.vaso.data)
lasso.vaso.y=as.numeric(pheno.plusmetabs$vasorespond.bin[PAHgroupAB])
lasso.vaso.fit=glmnet(x = lasso.vaso.data2, y=lasso.vaso.y)
plot(lasso.vaso.fit, label = T)
##cross-validation fit
cvfit.vaso <- glmnet::cv.glmnet(x = lasso.vaso.data2, y = lasso.vaso.y)

lasso.vaso.coef.lambda.1se=coef(cvfit.vaso, s = "lambda.1se")
lasso.vaso.coef.lambda.min=coef(cvfit.vaso, s = "lambda.min")

lasso.vaso.coef.1se.min=as.matrix(cbind(lasso.vaso.coef.lambda.1se,lasso.vaso.coef.lambda.min))
write.csv(x = lasso.vaso.coef.1se.min, file="lasso.vaso coefficients 1se and min.csv")
plot(cvfit.vaso,)

##predict in C
lasso.vaso.phenotypesC=simplepheno2[PAHgroupC,c(3:5,7:8)]
lasso.vaso.dataC=cbind(lasso.vaso.phenotypesC,t(sig.detected.vasoAB.data[,PAHgroupC]))
lasso.vaso.dataC$sex=as.numeric(lasso.vaso.dataC$sex)-1
lasso.vaso.dataC2=as.matrix(lasso.vaso.dataC)

lasso.vaso.predictionC=predict(cvfit.vaso, newx = lasso.vaso.dataC2, s = "lambda.min")
library(pROC)
lasso.vasoROC=roc(pheno.plusmetabs$vasorespond.bin[PAHgroupC],lasso.vaso.predictionC[,1])
table(pheno.plusmetabs$sex[PAHgroupC])
plot(lasso.vasoROC)
plot(lasso.vaso.predictionC~pheno.plusmetabs$vasorespond.bin[PAHgroupC])
wilcox.test(lasso.vaso.predictionC~pheno.plusmetabs$vasorespond.bin[PAHgroupC])
lasso.vasoROC
ci.auc(lasso.vasoROC)
PAHgroupABC=which(pheno.plusmetabs$Random.group2 %in% c("A","B","C") & is.na(pheno.plusmetabs$Exclude) & pheno.plusmetabs$PAH == 1)
pheno.plusmetabs$lasso.vasoModel=NA

lasso.vaso.phenotypesABC=simplepheno2[PAHgroupABC,c(3:5,7:8)]
lasso.vaso.dataABC=cbind(lasso.vaso.phenotypesABC,t(welldetect.vaso.tophits.data[,PAHgroupABC]))
lasso.vaso.dataABC$sex=as.numeric(lasso.vaso.dataABC$sex)-1
lasso.vaso.dataABC2=as.matrix(lasso.vaso.dataABC)
pheno.plusmetabs$lasso.vasoModel[PAHgroupABC]=predict(cvfit.vaso,newx=lasso.vaso.dataABC2, s = "lambda.min")
lasso.vasoROC.ABC=roc(pheno.plusmetabs$vasorespond.bin[PAHgroupABC],pheno.plusmetabs$lasso.vasoModel[PAHgroupABC])

lasso.vasoROC.ABCdf=data.frame(lasso.vasoROC.ABC$sensitivities,lasso.vasoROC.ABC$specificities,lasso.vasoROC.ABC$thresholds)
lasso.vasoROC.ABCdf$yi=lasso.vasoROC.ABCdf$lasso.vasoROC.ABC.specificities+lasso.vasoROC.ABCdf$lasso.vasoROC.ABC.sensitivities
lasso.vasoROC.ABCdf[247,]
# lasso.vasoROC.ABC.sensitivities lasso.vasoROC.ABC.specificities lasso.vasoROC.ABC.thresholds       yi
# 89                  0.9108635                  0.7777778                1.740904 1.688641


boxplot(pheno.plusmetabs$lasso.vasoModel[which(pheno.plusmetabs$PAH %in% 0:1 & is.na(pheno.plusmetabs$Exclude))] ~
          paste0(pheno.plusmetabs$Random.group2,pheno.plusmetabs$PAH, pheno.plusmetabs$vasorespond.bin2)[which(pheno.plusmetabs$PAH %in% 0:1 & is.na(pheno.plusmetabs$Exclude))],
        col=2:3, ylab = "lasso.vaso Model", lty=1, ylim=c(-0.25,0.25)
        )
abline(h=median(pheno.plusmetabs$lasso.vasoModel[which(pheno.plusmetabs$vasorespond.bin %in% 0)], na.rm = T), lty=2)
abline(v=c(2.5,4.5), lty=3)
legend("bottomright", legend=c("Non-responders","Vasoresponders"), col=2:3, pch=15)


for(gene in hemnes.genes){
  png(paste0(gene," boxplot.png"), 600, 560)

  boxplot(log10(all.reads[gene,which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))]) ~
            paste0(pheno.plusmetabs$Random.group2,pheno.plusmetabs$PAH, pheno.plusmetabs$vasorespond.bin2)[which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))],
          col=2:4, ylab = gene, lty=1)
  abline(h=median(log10(all.reads[gene,which(pheno.plusmetabs$vasorespond.bin %in% 0)])), lty=2)
  abline(v=c(3.5,6.5,9.5), lty=3)
  legend("topleft", legend=c("Controls","PAH non-responders","PAH vasoresponders"), col=2:4, pch=15)
  dev.off()
}

for(gene in rownames(welldetect.vaso.tophits.data)){
  png(paste0(gene," boxplot.png"), 600, 560)

  boxplot(log10(all.reads[gene,which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))]) ~
            paste0(pheno.plusmetabs$Random.group2,pheno.plusmetabs$PAH, pheno.plusmetabs$vasorespond.bin2)[which(pheno.all$PAH %in% 0:1 & is.na(pheno.all$Exclude))],
          col=2:4, ylab = gene, lty=1)
  abline(h=median(log10(all.reads[gene,which(pheno.plusmetabs$vasorespond.bin %in% 0)])), lty=2)
  abline(v=c(3.5,6.5,9.5), lty=3)
  legend("topleft", legend=c("Controls","PAH non-responders","PAH vasoresponders"), col=2:4, pch=15)
  dev.off()
}
