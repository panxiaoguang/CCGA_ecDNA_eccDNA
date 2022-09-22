set.seed(80)
library(maftools)
library(berryFunctions)

pdf("UBC_lollipopPlot_TTN_ecDNA.pdf",w=16,h=4.5)

cohort1 <- "UBC_ecDNA_Negative_TTN.tsv"
cohort2 <- "UBC_ecDNA_Positive_TTN.tsv"
C1 = read.maf(maf=cohort1, verbose=FALSE)
C2 = read.maf(maf=cohort2, verbose=FALSE)

lollipopPlot2(m1 = C1, m2 = C2,
              m1_name = "ecDNA-", m2_name = "ecDNA+",
              gene = "TTN", 
              AACol1 = "TTN:NM_001267550",
              AACol2 = "TTN:NM_001267550",
              refSeqID = "NM_001267550",
              m1_label = "all", m2_label = "all",
              labPosAngle = 90,
              roundedRect = TRUE,
              showDomainLabel = FALSE,
              domainBorderCol = NA)

dev.off()
