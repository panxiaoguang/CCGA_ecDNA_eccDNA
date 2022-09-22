set.seed(80)
library(maftools)
library(berryFunctions)

pdf("UBC_lollipopPlot_TP53.pdf",w=10,h=4.5)
ubc.maf <- "UBC_TP53.TSV"

ubc = read.maf(maf = ubc.maf, verbose=FALSE)

lollipopPlot(maf = ubc,
             gene = "TP53", 
             AACol = "TP53:NM_000546", 
             refSeqID = "NM_000546",
             labelPos = "all",
             labPosAngle = 90,
             repel = TRUE,
             roundedRect = TRUE,
             defaultYaxis = TRUE,
             domainBorderCol = NA,
#             showDomainLabel = FALSE,
             showMutationRate = FALSE)

dev.off()