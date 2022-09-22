import glob
import pysam as ps
import pandas as pd


chroms=["chr"+str(x) for x in range(1,23)]
chroms.append("chrM")
chroms.append("chrX")
chroms.append("chrY")

def get_mappings(sample,chr):
    samfile=ps.AlignmentFile(f"sorted_{sample}_circle.bam","rb")
    df=pd.DataFrame(samfile.get_index_statistics()).set_index("contig")
    total_reads=sum(df.total)
    fin=df.loc[chr]
    fin["mappedRatio"]=fin["mapped"]/total_reads*100
    fin.to_csv(f"{sample}.mapping.stats.tsv",sep="\t")



files=glob.glob("sorted_*_circle.bam")
for fs in files:
   sp=fs.replace("_circle.bam","").replace("sorted_","")
   get_mappings(sp,chroms)
