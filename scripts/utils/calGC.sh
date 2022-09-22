#!/usr/bin/bash
chromSize="/home/panxiaoguang/Project/qd-ECC4/S/ECC_report/DataBase/hg38.chromo.size"
Genome="/home/panxiaoguang/Project/DataBase/hg38.fa"
awk -F '\t' '{printf "%s\t%s\t%s\t%s:%s-%s\n",$1,$2,$3,$1,$2,$3}' $1.bed | \
bedtools flank -i - -l 1.0 -r 0 -pct -g $chromSize | \
bedtools nuc -fi $Genome -bed - | \
awk -F '\t' '{print $4,$6} ' OFS='\t' > GCs/$1.upstream.gc.txt

awk -F '\t' '{printf "%s\t%s\t%s\t%s:%s-%s\n",$1,$2,$3,$1,$2,$3}' $1.bed | \
bedtools nuc -fi $Genome -bed - | \
awk -F '\t' '{print $4,$6} ' OFS='\t' > GCs/$1.ecc.gc.txt

awk -F '\t' '{printf "%s\t%s\t%s\t%s:%s-%s\n",$1,$2,$3,$1,$2,$3}' $1.bed | \
bedtools flank -i - -l 0 -r 1.0 -pct -g $chromSize | \
bedtools nuc -fi $Genome -bed - | \
awk -F '\t' '{print $4,$6} ' OFS='\t' > GCs/$1.downstream.gc.txt

python3 calGC.py --sample $1 --out $1.gcContents.txt

rm -rf GCs/$1.upstream.gc.txt
rm -rf GCs/$1.ecc.gc.txt
rm -rf GCs/$1.downstream.gc.txt