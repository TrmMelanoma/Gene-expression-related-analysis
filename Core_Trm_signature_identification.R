[1]Expression merge
quantile_normalization <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

source("../Function/Probe_data_affy.R")


myinf1 <- "../raw/GSE15907_series_matrix.txt"
myinf2 <- "../raw/GSE37448_series_matrix.txt"
myinf3 <- "../raw/GSE47045_series_matrix.txt"

myoutf1 <- "../Meta_data/probe_set_normalization_data.txt"
myoutf2 <- "../Meta_data/Gene_symbol_normalization_data.txt"



data1 <- probe_data_affy(myinf1)
data1 <- log2(data1+1)

data2 <- probe_data_affy(myinf2)
data2 <- log2(data2+1)

data3 <- probe_data_affy(myinf3)

xx1 <- row.names(data1)
xx2 <- row.names(data2)
xx3 <- row.names(data3)

com <- Reduce(intersect, list(xx1,xx2,xx3))

data1 <- data1[com,]
data2 <- data2[com,]
data3 <- data3[com,]

data <- cbind(data1,data2,data3)
metadata <- quantile_normalization(data)

batch <- c(rep(1,ncol(data1)),rep(2,ncol(data2)),rep(3,ncol(data3)))

library(sva)
metadata <- ComBat(dat=as.matrix(data),batch=batch,mod=NULL,par.prior=T)

write.table(metadata,myoutf1,sep="\t",quote=F)

#------------------------------------------------------

myinf1<-"../GSE15907/processed/GSE15907_Sample_info.txt"
myinf2<-"../GSE37448/processed/GSE37448_Sample_info.txt"
myinf3<-"../GSE47045/processed/GSE47045_Sample_info.txt"
info1<-read.csv(myinf1,sep='\t')
info2<-read.csv(myinf2,sep='\t')
info3<-read.csv(myinf3,sep='\t')

info1<-info1[,c("Title","phenotype.markers")]
info2<-info2[,c("Title","phenotype.markers")]
info3<-info3[,c("Title","cell.type")]
colnames(info3)<-c("Title","phenotype.markers")
meta_info<-rbind(info1,info2,info3)
myoutf<-"../Meta_data/Meta_cell_info.txt"
write.table(meta_info,myoutf,sep='\t',quote=F)

myinf<-"../Meta_data/Meta_cell_info.txt"
info<-read.csv(myinf,sep='\t')


###############################
myinf1 <- "../Meta_data/probe_set_normalization_data.txt"
metadata <- read.table(myinf1,sep="\t",quote=NULL)


tmpf = "../GPL/GPL6246.annot"
conIn = file(tmpf, "r")
data = readLines(conIn, -1)
close(conIn)
mysta = grep("platform_table_begin", data)
myend = grep("platform_table_end", data)
data = read.table(tmpf, sep="\t", header=T, row.names=1, skip=mysta, nrows=myend-mysta-2, quote="", comment.char ="")
tmp = row.names(data)
GPL = data[, "Gene.symbol"]
names(GPL) = tmp

platform = 1
data <- metadata
data <- as.data.frame(data)
comGene = intersect(names(GPL), row.names(data))
comGene = comGene[comGene!=""]
GPL = as.character(GPL[comGene])
data = data[comGene,]

data = cbind(GPL, data)
data[,1] = as.character(data[,1])
se = grep("///", data[,1])
if(length(se)>0)
{
  data = data[-grep("///", data[,1]), ]
}
data = data[data[,1]!="", ]

tmp = as.character(data[,1])
tmp = rle(sort(tmp)) 
uni.gen = tmp$values[tmp$lengths==1]
mul.gen = tmp$values[tmp$lengths>1]
dat1 = data[data[,1]%in%uni.gen, ]
tmp = as.character(dat1[,1])
dat1 = dat1[, 2:ncol(dat1)]
row.names(dat1) = tmp
if(length(mul.gen)>0)
{
  dat2 = matrix(0, length(mul.gen), ncol(dat1))
  row.names(dat2) = mul.gen
  colnames(dat2) = colnames(dat1)
  for(k in 1:length(mul.gen))
  {
    xx = data[data[,1]==mul.gen[k], 2:ncol(data)]
    if(platform==1)
    {
      tmp = apply(xx, 1, mean, na.rm=T)
      dat2[k,] = as.numeric(xx[which.max(tmp), ])
    }else
    {
      tmp = apply(xx, 2, mean, na.rm=T)
      dat2[k,] = as.numeric(tmp)  
    }
  }
  data = rbind(dat1, dat2)
}else
{
  data = dat1
}
data = round(data, 6)
write.table(data,myoutf2,sep="\t",quote=F)


###############################

[1]process the data
myinf2 <- "../Meta_data/Meta_cell_info.txt"
myinf1 <- "../Meta_data/Gene_symbol_normalization_data.txt"

data = read.table(myinf1,sep="\t",quote=NULL)
info = read.csv(myinf2,sep="\t",header=T,stringsAsFactors=F)
info$Sample = row.names(info)
info = info[,c("Sample","Title")]

tag1 = grep("T.8",info$Title, fix=T)
tag2 = grep("spleen_naive",info$Title, fix=T)
tag3 = grep("spleen_tcm",info$Title, fix=T)
tag4 = grep("spleen_tem",info$Title, fix=T)
tag5 = grep("skin_trm",info$Title, fix=T)

info1 = info[tag1,]
info2 = info[tag2,]
info3 = info[tag3,]
info4 = info[tag4,]
info5 = info[tag5,]

info = rbind(info1,info2,info3,info4,info5)
info[,"label"] = info$Title

#contruct the group of the cells
info$label = gsub("#1","", info$label,fix=T)
info$label = gsub("#2","", info$label,fix=T)
info$label = gsub("#3","", info$label,fix=T)
info$label = gsub("#4","", info$label,fix=T)
info$label = gsub("#5","", info$label,fix=T)

info$label = gsub("_1","", info$label,fix=T)
info$label = gsub("_2","", info$label,fix=T)
info$label = gsub("_3","", info$label,fix=T)

tar = unique(info$label)

#create the new matrix
res = matrix(0, nrow(data), length(tar))
row.names(res) = row.names(data)
colnames(res) = tar
res = as.data.frame(res)

for(i in 1 : length(tar))
{
  cat("\r",i)
  tag = info$label == tar[i]
  sam = row.names(info)[tag]
  tmp_data = data[,sam]
  exp = apply(tmp_data, 1, mean)
  res[,i] = exp
}

myoutf ="../Meta_data/average_gene_expression_CD8_T.xls"
write.table(res,myoutf,sep="\t",quote=F)

[2]begin to build the signature
myinf1 = "../Meta_data/average_gene_expression_CD8_T.xls"
data = read.table(myinf1,sep="\t",quote=NULL)

med = apply(data, 1, median)
data = data-med

avg = apply(data, 2, mean)
std = apply(data, 2, sd)
for(k in 1:ncol(data))
{
  data[,k] = (data[,k]-avg[k])/std[k]
}


res1 = data
for(k in 1:ncol(res1))
{
  tmp = data[,k]
  tmp[tmp<0]=0
  tmp = -log10(pnorm(-tmp)*2)
  tmp[tmp>10]=10
  res1[,k] = tmp
}
colnames(res1) = paste(colnames(data), "_up", sep="")

res2 = data
for(k in 1:ncol(res2))
{
  tmp = data[,k]
  tmp[tmp>0]=0
  tmp = -log10(pnorm(tmp)*2)
  tmp[tmp>10]=10
  res2[,k] = tmp
}
colnames(res2) = paste(colnames(data), "_dn", sep="")

res = cbind(res1, res2)

minv = min(res)
maxv= max(res)
res = (res-minv)/(maxv-minv)
colnames(res) = gsub(" ", "", colnames(res))

myoutf = "../Meta_data/Weight_for_Trm_across_CD8_T.xls"
write.table(res, myoutf, sep="\t", quote=F)

sum(res$skin_trm_up > 0.3)
sum(res$skin_trm_dn > 0.3)


up_gene<-as.matrix(toupper(row.names(res)[res$skin_trm_up > 0.3]))
colnames(up_gene)<-"up_gene"
dn_gene<-as.matrix(toupper(row.names(res)[res$skin_trm_dn > 0.3]))
colnames(dn_gene)<-"dn_gene"

library(gdata)
gene_list <- cbindX(up_gene, dn_gene)

myoutf<-"../Trm_signature/Mackey_immgen_skin_Trm_signature_UP&DN.txt"
write.table(gene_list,myoutf,sep='\t',quote = F)













