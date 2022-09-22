# CCGA_circulome
[![DOI](https://zenodo.org/badge/539755936.svg)](https://zenodo.org/badge/latestdoi/539755936)


This repository contains code  from "Lv et al. Chromosomal, Extrachromosomal, and Transcriptomic Alternation Landscape in Urothelial Bladder Carcinoma".

### Code

#### obtain GC contents

- **calGC.py**
- **calGC.sh**

#### obtain repeat ratios

- **calReps.sh**
- **calReps.py**

#### some basic statistics

- **calMappings_nochrM.py** obtain total mapping reads for each sample except chrM
- **get_index_statistic.py** obtain total mapping ratio for each chromosome
- **filter.py** filter raw eccDNA using some metrics 
- **start_anno.jl** annotate eccDNA usings their start junction
- **run_GISTIC2.sh** analyze total CNA characters from all samples
- **getBreaks.jl** obtain all discordant break points from ecDNA graphs
- **PABPC1_CIRCLE_RNA_WGS_alleleCount.sh** example for allelic frequency and coverage 

#### scripts for main figure plots

- **annotate_eccDNA_elements.plotData.R**
- **B-allele_WGS_CircleSeq.plot.R**
- **classify_samples_ecDNA.plotData.R**
- **compare_eccDNA_abundance_heatmap.plot.R**
- **compare_ecDNA_size_copy_BP.plotData.R**
- **compare_EPM.plot.R**
- **detected_chimeric_eccDNA.plotData.R**
- **ecc_abundance_distribution.plot.R**
- **ecc_map_genes.plotData.R**
- **ecc_mRNA_cor.plot.R**
- **ecc_mRNA_cor_rank.plot.R**
- **eccNumbers_amp_noamp.plotData.R**
- **EPM_diff_GSEA.plot.R**
- **epm_mRNA_cor.plot.R**
- **foldchange_amp_noamp.plot.R**
- **genomic_dotplot_allele_frequency.plot.R**
- **length_category_compare.plotData.R**
- **length_distribution_eccDNA.plot.R**
- **ligation_genes_exp_cor.plot.R**
- **mRNA_ecDNA_eccDNA.plot.R**
- **multiomicsViz.plot.R**
- **oncoprint_for_BLCA_UBC.plot.R**
- **PABPC1_circle_coverage.plot.R**
- **statistic_for_BLCA_proportion.plotData.R**

#### SV format and merge

- **format_delly.jl**
- **format_manta.jl**
- **merge.manta.delly.py**

### Contact
If you have any questions concerning code or data, please do not hesitate to contact us at panxiaoguang@genomics.cn.