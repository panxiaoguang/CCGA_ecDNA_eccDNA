#!/usr/bin/bash
sample=$1

echo '#!/usr/bin/bash'

echo "bedtools multicov -bams alignData/sorted_${sample}_circle.bam -bed dbs/UCSC_repeatMarsker.region.bed > repeats/repeats.in.${sample}.bed"
