set.seed(80)
library(maftools)
library(berryFunctions)

pdf("UBC_lollipopPlot_KMT2D.pdf",w=10,h=4.5)
ubc.maf <- "UBC_KMT2D.TSV"

ubc = read.maf(maf = ubc.maf, verbose=FALSE)

lollipopPlot(maf = ubc,
             gene = "KMT2D",
             AACol = "KMT2D:NM_003482",
             refSeqID = "NM_003482",
             labelPos = "all",
             labPosAngle = 90,
             repel = TRUE,
             roundedRect = TRUE,
             defaultYaxis = TRUE,
             domainBorderCol = NA,
             showDomainLabel = FALSE,
             showMutationRate = FALSE)

dev.off()