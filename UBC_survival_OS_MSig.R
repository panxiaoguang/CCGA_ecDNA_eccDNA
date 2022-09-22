library(survival)
library(survminer)
pdf("UBC_survival_OS_MSig.pdf",w=4.5,h=4,onefile=F)
data <- read.table("UBC-clinical.tsv",sep="\t",header=T)
fit <- survfit(Surv(days_to_last_followup, Overall_Survival_Status) ~ MSig, data = data)
ggsurvplot(fit, data = data, pval=TRUE, pval.method=FALSE,conf.int=FALSE,risk.table=FALSE,palette="npg",tables.height=0.2,xlab="Overall survival (Days)",legend.labs=c("MSig1","MSig2","MSig3"),legend.title="",legend=c(0.9, 0.9))
dev.off()
#time  survival or censoring time
#status  0=censored, 1=death