import pandas as pd
import fire
from functools import reduce

def readIN(sample="BC1C",out="BC1C.gcContents.txt"):
    df1=pd.read_csv(f"GCs/{sample}.upstream.gc.txt",sep="\t")
    df1.columns=["ecc","upstream"]
    df2=pd.read_csv(f"GCs/{sample}.ecc.gc.txt",sep="\t")
    df2.columns=["ecc","self"]
    df3=pd.read_csv(f"GCs/{sample}.downstream.gc.txt",sep="\t")
    df3.columns=["ecc","downstream"]
    fin=reduce(lambda x,y:pd.merge(x,y,on="ecc",how="outer"), [df1,df2,df3])
    fin.to_csv(f"GCs/{out}",sep="\t",index=False)

if __name__ == '__main__':
    fire.Fire(readIN)