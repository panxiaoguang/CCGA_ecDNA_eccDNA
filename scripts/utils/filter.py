import pandas as pd
import numpy as np
import glob
## filter parameters: split >=2 circlescore >=200 sci>=0.33 eci >=0.33 sd<mean coveragecon <=0.1


def filt(sample):
    try:
        df=pd.read_table(f"eccDNAs/{sample}_circle_site.bed")
        df["length"] = df.end- df.start
        df = df.loc[(df["soft-clipped"] >=2) & (df.score >=200) & (df["start_ratio"] >=0.33) & (df["end_ratio"] >=0.33) & (df.continuity <=0.2) &(df["std"] < df["mean"])]
        df['ft']=np.where(df.length<=2000, 100,df.discordants)
        df = df.loc[df["ft"]>0]
        df = df.drop("ft",axis=1)
    except Exception as e:
        print(f"{sample} have error as {e}")
    return df

if __name__ == "__main__":
    for fs in glob.glob("eccDNAs/*_circle_site.bed"):
        sample = fs.replace("_circle_site.bed","").replace("eccDNAs/","")
        rst = filt(sample)
        rst.to_csv(f"filtered_eccs/{sample}_circle_site.filter.tsv",sep="\t",index=None)
