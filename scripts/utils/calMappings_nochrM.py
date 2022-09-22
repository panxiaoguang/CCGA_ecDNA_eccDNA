import glob
import pysam as ps
def get_mappings(sample):
    samfile=ps.AlignmentFile(sample,'rb')
    totals=sum([x.mapped for x in samfile.get_index_statistics() if x.contig !="chrM"])
    return totals

files=glob.glob("/dellfsqd2/ST_LBI/USER/xuzhe/ECC/KidneyCancer/preData/sorted_*_circle.bam")
print("sample\ttotalmappings")
for fs in files:
   sp=fs.replace("_circle.bam","").replace("sorted_","")
   mappings=get_mappings(fs)
   print(sp,"\t",mappings)
