
library(TCGAbiolinks)
library(maftools)
library(dplyr)
library(DT)

#########################function definition################################
myLegend <- function (leg){
  legend("top", legend=leg, box.lwd = 0, text.font = 2, inset = .04)
  dev.off()
}

kmPlot <- function (genename){
  png(sprintf('%sKM_lihc_muse[%s].png', PlotPath, genename))
  mafSurvival(maf = lihc.muse, genes = genename, time = 'days_to_last_follow_up', Status = 'vital_status', isTCGA = TRUE)
  myLegend(genename)
  png(sprintf('%sKM_lihc_mutect[%s].png', PlotPath, genename))
  mafSurvival(maf = lihc.mutect, genes = genename, time = 'days_to_last_follow_up', Status = 'vital_status', isTCGA = TRUE)
  myLegend(genename)
  png(sprintf('%sKM_lihc_somaticsniper[%s].png', PlotPath, genename))
  mafSurvival(maf = lihc.somaticsniper, genes = genename, time = 'days_to_last_follow_up', Status = 'vital_status', isTCGA = TRUE)
  myLegend(genename)
  png(sprintf('%sKM_lihc_varscan2[%s].png', PlotPath, genename))
  mafSurvival(maf = lihc.varscan2, genes = genename, time = 'days_to_last_follow_up', Status = 'vital_status', isTCGA = TRUE)
  myLegend(genename)
}
##############################################################################

#Output plots path
PlotPath<- "./plots/"

if ( !dir.exists(PlotPath) ) { dir.create(PlotPath) }

if ( !exists("mutatedGenes") ) { mutatedGenes = c() }


if (!exists("lihc.clinical")) {
  lihc.clinical = GDCquery_clinic(project = "TCGA-LIHC", type = "clinical")
}


###### Two Important Modification to Clinical Data in Compliance with MAF data ###########
# Whereas first column has a different name(submitter_id), it is the Tumor_Sample_Barcode. Therefore it should be renamed.
colnames(lihc.clinical)[1] = "Tumor_Sample_Barcode"
# The mafsurvival needs the status field be in boolean type (The event is 'death' in the survival plot)
lihc.clinical$vital_status = lihc.clinical$vital_status=="Dead"


# pipeline options are muse, mutect, somaticsniper, varscan2
# muse
if (!exists("lihc.muse")) {
  lihc.muse     = GDCquery_Maf("LIHC", directory = "GDCdata", pipelines = "muse")
  lihc.muse$Tumor_Sample_Barcode = substr(lihc.muse$Tumor_Sample_Barcode, 1, 12)
  lihc.muse     = read.maf(maf = lihc.muse, clinicalData = lihc.clinical)
  mutatedGenes  = rbind(mutatedGenes, mafSummary(lihc.muse)$gene[1:3,1])
}
# mutect
if (!exists("lihc.mutect")) {
  lihc.mutect   = GDCquery_Maf("LIHC", directory = "GDCdata", pipelines = "mutect")
  lihc.mutect$Tumor_Sample_Barcode = substr(lihc.mutect$Tumor_Sample_Barcode, 1, 12)
  lihc.mutect   = read.maf(maf = lihc.mutect, clinicalData = lihc.clinical)
  mutatedGenes  = rbind(mutatedGenes, mafSummary(lihc.mutect)$gene[1:3,1])
}
# somaticsniper
if (!exists("lihc.somaticsniper")) {
  lihc.somaticsniper  = GDCquery_Maf("LIHC", directory = "GDCdata", pipelines = "somaticsniper")
  lihc.somaticsniper$Tumor_Sample_Barcode = substr(lihc.somaticsniper$Tumor_Sample_Barcode, 1, 12)
  lihc.somaticsniper  = read.maf(maf = lihc.somaticsniper, clinicalData = lihc.clinical)
  mutatedGenes        = rbind(mutatedGenes, mafSummary(lihc.somaticsniper)$gene[1:3,1])
}
# varscan2
if (!exists("lihc.varscan2")) {
  lihc.varscan2       = GDCquery_Maf("LIHC", directory = "GDCdata", pipelines = "varscan2")
  lihc.varscan2$Tumor_Sample_Barcode = substr(lihc.varscan2$Tumor_Sample_Barcode, 1, 12)
  lihc.varscan2       = read.maf(maf = lihc.varscan2, clinicalData = lihc.clinical)
  mutatedGenes        = rbind(mutatedGenes, mafSummary(lihc.varscan2)$gene[1:3,1])
}

mutatedGenes = unique(mutatedGenes)

png(paste0(PlotPath, "lihc_muse.png"))
plotmafSummary(maf = lihc.muse, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
 
png(paste0(PlotPath, "lihc_mutect.png"))
plotmafSummary(maf = lihc.mutect, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

png(paste0(PlotPath, "lihc_somaticsniper.png"))
plotmafSummary(maf = lihc.somaticsniper, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

png(paste0(PlotPath, "lihc_varscan2.png"))
plotmafSummary(maf = lihc.varscan2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()


sapply(unlist(mutatedGenes), kmPlot)