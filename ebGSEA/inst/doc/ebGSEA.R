## ----vignette-options, echo=FALSE, warning=FALSE------------------------------
require(BiocStyle)

## -----------------------------------------------------------------------------
load("ebgseaDATA.rda");
ls();
dim(dataSMK.m)

## -----------------------------------------------------------------------------
library(ebGSEA)
sgt.m <-doGT(phenoSMK.v,dataSMK.m,array="450k",ncores=4)
head(sgt.m)

## ----warning=FALSE------------------------------------------------------------
data("MSigDB-28Feb14-data")
topGSEA.lm <- doGSEAwt(rankEID.m = sgt.m, ptw.ls = listEZ.lv, ncores = 4, minN = 5,adjPVth = 0.05)

## -----------------------------------------------------------------------------
topGSEA.lm[[1]]

## -----------------------------------------------------------------------------
head(topGSEA.lm$Genestat[[1]])

## -----------------------------------------------------------------------------
plot(x = topGSEA.lm[[1]][,2], y = -log10(topGSEA.lm[[1]][,5]), xlab ="AUC", ylab = "-log10(adjP)", main = "AUC and adjP for each enriched pathway", pch = 21, bg = "red");

## -----------------------------------------------------------------------------
data("SampleCpG")
sigEID.ls <- selEIDfromSelCpG(selCpG.v = sampleCpG.v, allCpG.v = allCpG.v, array = "450k")

## -----------------------------------------------------------------------------
summary(sigEID.ls)

## ----message=FALSE------------------------------------------------------------
topGSEAft.lm <- doGSEAft(selEID.v = sigEID.ls$selEID, ptw.ls = listEZ.lv, allEID.v = names(mapEIDto450k.lv), ncores = 1, adjPVth = 0.05)

## -----------------------------------------------------------------------------
summary(topGSEAft.lm)

## -----------------------------------------------------------------------------
head(topGSEAft.lm$`Rank(P)`)

## -----------------------------------------------------------------------------
plot(x = log2(as.numeric(topGSEAft.lm[[1]][,3])), y = -log10(as.numeric(topGSEAft.lm[[2]][,5])), xlab ="log2(OR)", ylab = "-log10(adjP)", main = "OR and adjP for each enriched pathway", pch = 21, bg = "red")

## ----sessionInfo, echo=T------------------------------------------------------
sessionInfo()

