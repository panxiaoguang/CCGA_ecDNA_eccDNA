set.seed(80)
library(maftools)
library(sigminer)
library(NMF)
library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)

ubc.maf <- "UBC_Combine.maf.tsv"
ubc = read.maf(maf = ubc.maf,verbose = FALSE)

mt_tally_SBS <- sig_tally(ubc, ref_genome="BSgenome.Hsapiens.UCSC.hg38", useSyn=TRUE,mode="SBS")
write.table(mt_tally_SBS,file="UBC_SBS_mt_tally.tsv",sep = "\t", col.names = NA)

# Testing
# mt_est_SBS <- sig_estimate(mt_tally_SBS$nmf_matrix, range = 2:8, nrun =10, use_random = FALSE,cores = 2,verbose = TRUE)
# show_sig_number_survey2(mt_est_SBS$survey)
# show_sig_number_survey(mt_est_SBS$survey)

mt_sig_SBS <- sig_extract(mt_tally_SBS$nmf_matrix, n_sig=3, nrun=100, cores=4)

sim_SBS <- get_sig_similarity(mt_sig_SBS, sig_db="SBS")
sim_SBS <- get_sig_similarity(mt_sig_SBS)
pdf("UBC_SBS_Sig.pdf",w=8,h=4,onefile=F)
show_sig_profile(mt_sig_SBS, mode="SBS", style="cosmic", x_label_angle=90)
dev.off()

pdf("UBC_SBS_exposure.pdf",w=8,h=4,onefile=F)
show_sig_exposure(mt_sig_SBS,hide_samps= FALSE,sig_names = c("APOBEC","Defective DNA mismatch repair","Aristolochic acid"), palette = c("#E64B35FF","#4DBBD5FF","#00A087FF"), rm_panel_border=TRUE)
sig_exposure_SBS <- get_sig_exposure(mt_sig_SBS)
write.table(sig_exposure_SBS,file="UBC_SBS_exposure.tsv",sep = "\t", col.names = NA)
dev.off()
