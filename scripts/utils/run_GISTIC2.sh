basedir=/home/panxiaoguang/Project/WGS/CNV/gistic2_results
echo --- running GISTIC --- 
segfile=/home/panxiaoguang/Project/WGS/CNV/final.seg
refgenefile=/home/panxiaoguang/app/GISTIC_2_0_23/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat
/home/panxiaoguang/app/GISTIC_2_0_23/gistic2 -b $basedir -seg $segfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme