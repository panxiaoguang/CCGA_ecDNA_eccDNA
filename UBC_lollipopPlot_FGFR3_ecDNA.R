set.seed(80)
library(maftools)
library(berryFunctions)

pdf("UBC_lollipopPlot_FGFR3_ecDNA.pdf",w=10,h=4.5)

cohort1 <- "UBC_ecDNA_Negative_FGFR3.tsv"
cohort2 <- "UBC_ecDNA_Positive_FGFR3.tsv"
C1 = read.maf(maf=cohort1, verbose=FALSE)
C2 = read.maf(maf=cohort2, verbose=FALSE)

lollipopPlot2(m1 = C1, m2 = C2,
              m1_name = "ecDNA-", m2_name = "ecDNA+",
              gene = "FGFR3", 
              AACol1 = "FGFR3:NM_001163213",
              AACol2 = "FGFR3:NM_001163213",
              refSeqID = "NM_001163213",
              m1_label = "all", m2_label = "all",
              labPosAngle = 90,
              roundedRect = TRUE,
              showDomainLabel = FALSE,
              domainBorderCol = NA)

dev.off()
