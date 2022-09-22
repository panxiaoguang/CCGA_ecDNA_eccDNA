set.seed(30)
rawdata <- read.table("UBC_SBS_exposure_PCT.tsv",sep="\t",header=T,row.names=1)
data <- rawdata
rownames(data)
colnames(data)
K <- kmeans(data[,c(-1)], centers = 3, nstart = 100)
write.table(K$cluster,file="UBC_SBS_exposure_PCT_kmeans.tsv",sep = "\t", col.names = NA)