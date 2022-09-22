import pysam as ps
import pandas as pd
import fire



def get_repeats(sample="none"):
    df=pd.read_csv(f"repeats/repeats.in.{sample}.bed",header=None,sep="\t")
    df.columns=["chrm","st","ed","name","value"]
    fin=df[df.chrm !="chrM"].groupby("name")["value"].sum().reset_index()
    samfile=ps.AlignmentFile(f"alignData/sorted_{sample}_circle.bam","rb")
    totals=sum([x.mapped for x in samfile.get_index_statistics() if x.contig !="chrM"])
    fin["ratio"]=fin.value/totals
    fin.to_csv(f"repeats/repStats.{sample}.tsv",index=False,sep="\t")


if __name__ == "__main__":
    fire.Fire(get_repeats)
