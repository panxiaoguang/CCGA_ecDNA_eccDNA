source /home/panxiaoguang/miniconda3/bin/activate biosoft
### get region-based bi-allele site
bcftools view --with-header -Ov -r chr8:99497918-119476539 -g het -m2 -M2 -v snps -o PABPC1/sample16N.biallele.vcf 16N.vcf.gz

### annotate allele site with exons
bedtools intersect -a PABPC1/sample16N.biallele.vcf -b hg38.protein.exon.bed -wa -wb  > PABPC1/sample16N.biallele.anno.vcf

### process snp.locis
Rscript process_snp.r PABPC1/sample16N.biallele.anno.vcf PABPC1/snpPos.gene.tsv PABPC1/snpPos.bed

### calculate coverage "deeptools"
bamCoverage -b Bca_16_T_subjunc.rmdup.bam --region chr8:99497918:119476539 -of bedgraph -o PABPC1/RNA.cov.bdg
bamCoverage -b 16T.cs.rmdup.sort.bam --region chr8:99497918:119476539 -of bedgraph -o PABPC1/WGS.cov.bdg
bamCoverage -b sorted_cBca_16T_circle.rmdup.bam --region chr8:99497918:119476539 -of bedgraph -o PABPC1/CIRCLE.cov.bdg

source /home/panxiaoguang/miniconda3/bin/activate allelecount
### calculate allele count "allelecounter"
alleleCounter -l PABPC1/snpPos.bed -b Bca_16_T_subjunc.rmdup.bam -o PABPC1/PABPC1.RNA.allele.txt -r /home/panxiaoguang/Project/DataBase/hg38.fa.fai -f 1 -F 1548 -m 0 -q 0
alleleCounter -l PABPC1/snpPos.bed -b 16T.cs.rmdup.sort.bam -o PABPC1/PABPC1.WGS.allele.txt -r /home/panxiaoguang/Project/DataBase/hg38.fa.fai
alleleCounter -l PABPC1/snpPos.bed -b sorted_cBca_16T_circle.rmdup.bam -o PABPC1/PABPC1.CIRCLE.allele.txt -r /home/panxiaoguang/Project/DataBase/hg38.fa.fai -f 1 -F 1548 -m 0 -q 0


