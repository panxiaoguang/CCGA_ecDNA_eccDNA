set.seed(80)
library(maftools)
library(berryFunctions)

pdf("UBC_lollipopPlot_KDM6A.pdf",w=10,h=4.5)
ubc.maf <- "UBC_KDM6A.TSV"

ubc = read.maf(maf = ubc.maf, verbose=FALSE)

lollipopPlot(maf = ubc,
             gene = "KDM6A",
             AACol = "KDM6A:NM_021140",
             refSeqID = "NM_021140",
             labelPos = "all",
             labPosAngle = 90,
             repel = TRUE,
             roundedRect = TRUE,
             defaultYaxis = TRUE,
             domainBorderCol = NA,
             showDomainLabel = FALSE,
             showMutationRate = FALSE)

dev.off()