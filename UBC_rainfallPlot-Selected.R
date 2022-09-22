set.seed(80)
library(maftools)
pdf("UBC_rainfallPlot_CCGA-UBC-001T.pdf",w=6,h=3)

ubc.maf <- "UBC_Combine_GW.maf.tsv"
ubc = read.maf(maf = ubc.maf,verbose = FALSE)

rainfallPlot(maf = ubc, detectChangePoints = TRUE, ref.build = "hg38",
             tsb = "CCGA-UBC-001T",
             pointSize = 0.3,
             fontSize = 1,
             savePlot = FALSE)

dev.off()
