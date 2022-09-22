set.seed(80)
library(maftools)

ubc.maf <- "UBC_Combine.maf.tsv"
ubc.clin <- "UBC-clinical.tsv"

ubc = read.maf(maf = ubc.maf,clinicalData = ubc.clin,verbose = FALSE)

survGroup(maf = ubc, top = 20, geneSetSize = 1, time = "days_to_last_followup", Status = "Overall_Survival_Status", verbose = FALSE)
survGroup(maf = ubc, top = 20, geneSetSize = 1, time = "Progression_free_survival_days", Status = "Progression_Status", verbose = FALSE)

pdf("UBC_Survival_BIRC6.pdf",w=4.5,h=4)
mafSurvival(maf = ubc, genes = 'BIRC6', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', col = c("#C8493A","#325189"))
dev.off()