
######Load necessary softwares/packages########

library(Seurat)
library(umap)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(GSVA)
library(snow)
library(fgsea)
library(ggplot2)
library(sctransform)

rm(list=ls()) 

#######Get all cell barcodes in single cell TCR sequencing #############


tar<-c(1,3,5,6)
barcode_info_all<-data.frame()

for(i in 1:length(tar))
{
  myinf1<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/skin/skin_filtered_contig_annotations.csv")
  myinf2<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/blood/blood_filtered_contig_annotations.csv")
  myinf3<-paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scTCRseq/Run_",tar[i],"/TCR/data/tumor/tumor_filtered_contig_annotations.csv")
  data_skin<-read.csv(myinf1)
  data_blood<-read.csv(myinf2)
  data_tumor<-read.csv(myinf3)
  
  data_skin[,"label_PT"]<-rep(paste0("skin_",tar[i]),nrow(data_skin))
  data_blood[,"label_PT"]<-rep(paste0("blood_",tar[i]),nrow(data_blood))
  data_tumor[,"label_PT"]<-rep(paste0("tumor_",tar[i]),nrow(data_tumor))
  
  data_skin[,"label"]<-rep("skin",nrow(data_skin))
  data_blood[,"label"]<-rep("blood",nrow(data_blood))
  data_tumor[,"label"]<-rep("tumor",nrow(data_tumor))
  
  data_skin[,"labeled_clonotype_id"]<-paste0(data_skin$raw_clonotype_id,"_",data_skin$label,tar[i])
  data_blood[,"labeled_clonotype_id"]<-paste0(data_blood$raw_clonotype_id,"_",data_blood$label,tar[i])
  data_tumor[,"labeled_clonotype_id"]<-paste0(data_tumor$raw_clonotype_id,"_",data_tumor$label,tar[i])
  
  barcode_info<-rbind(data_skin,data_blood,data_tumor)
  
  
  tag<-barcode_info$is_cell=="True"
  barcode_info<-barcode_info[tag,]
  tag<-barcode_info$high_confidence=="True"
  barcode_info<-barcode_info[tag,]
  tag<-barcode_info$productive=="True"
  barcode_info<-barcode_info[tag,]
  
  barcode_info$barcode<-gsub("-1","",barcode_info$barcode)
  barcode_info$barcode<-paste0(barcode_info$label,"_",barcode_info$barcode)
  barcode_info$barcode<-gsub("_",paste0("_",tar[i],"_"),barcode_info$barcode)
  barcode_info_all<-rbind(barcode_info,barcode_info_all)
}

tag<-barcode_info_all$v_gene=="TRAV10"&barcode_info_all$j_gene=="TRAJ18"



#############Merge single cell RNA-seq data for 4 patients, 3 tissues/patient##############


myinf1 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/Run_6/skin/"
skin_6.data <- Read10X(data.dir=myinf1)
skin_6 <- CreateSeuratObject(counts = skin_6.data, project = "skin_6",min.cells=3,min.features=200)

myinf2 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/Run_1/skin/"
skin_1.data <- Read10X(data.dir=myinf2)
skin_1 <- CreateSeuratObject(counts = skin_1.data, project = "skin_1",min.cells=3,min.features=200)

myinf3 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/Run_3/skin/"
skin_3.data <- Read10X(data.dir=myinf3)
skin_3 <- CreateSeuratObject(counts = skin_3.data, project = "skin_3",min.cells=3,min.features=200)


myinf5 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/Run_5/skin/"
skin_5.data <- Read10X(data.dir=myinf5)
skin_5 <- CreateSeuratObject(counts = skin_5.data, project = "skin_5",min.cells=3,min.features=200)

skin_1356<- merge(x= skin_6, y = c(skin_1,skin_3,skin_5), 
                  add.cell.ids=c("skin_6","skin_1","skin_3","skin_5")
                  ,project = "skin_1356")


myinf1 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/Run_6/tumor/"
tumor_6.data <- Read10X(data.dir=myinf1)
tumor_6 <- CreateSeuratObject(counts = tumor_6.data, project = "tumor_6",min.cells=3,min.features=200)

myinf2 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/Run_1/tumor/"
tumor_1.data <- Read10X(data.dir=myinf2)
tumor_1 <- CreateSeuratObject(counts = tumor_1.data, project = "tumor_1",min.cells=3,min.features=200)

myinf3 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/Run_3/tumor/"
tumor_3.data <- Read10X(data.dir=myinf3)
tumor_3 <- CreateSeuratObject(counts = tumor_3.data, project = "tumor_3",min.cells=3,min.features=200)


myinf5 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/Run_5/tumor/"
tumor_5.data <- Read10X(data.dir=myinf5)
tumor_5 <- CreateSeuratObject(counts = tumor_5.data, project = "tumor_5",min.cells=3,min.features=200)

tumor_1356<- merge(x= tumor_6, y = c(tumor_1,tumor_3,tumor_5), 
                   add.cell.ids=c("tumor_6","tumor_1","tumor_3","tumor_5")
                   ,project = "tumor_1356")



myinf1 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/Run_6/blood/"
blood_6.data <- Read10X(data.dir=myinf1)
blood_6 <- CreateSeuratObject(counts = blood_6.data, project = "blood_6",min.cells=3,min.features=200)

myinf2 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/Run_1/blood/"
blood_1.data <- Read10X(data.dir=myinf2)
blood_1 <- CreateSeuratObject(counts = blood_1.data, project = "blood_1",min.cells=3,min.features=200)

myinf3 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/Run_3/blood/"
blood_3.data <- Read10X(data.dir=myinf3)
blood_3 <- CreateSeuratObject(counts = blood_3.data, project = "blood_3",min.cells=3,min.features=200)


myinf5 <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/Run_5/blood/"
blood_5.data <- Read10X(data.dir=myinf5)
blood_5 <- CreateSeuratObject(counts = blood_5.data, project = "blood_5",min.cells=3,min.features=200)

blood_1356<- merge(x= blood_6, y = c(blood_1,blood_3,blood_5), 
                   add.cell.ids=c("blood_6","blood_1","blood_3","blood_5")
                   ,project = "blood_1356")






STB_1356<-merge(x=skin_1356,y=c(tumor_1356,blood_1356))

table(STB_1356@meta.data$orig.ident)

#STB_1356
#blood_1 blood_3 blood_5 blood_6  skin_1  skin_3  skin_5  skin_6 tumor_1 tumor_3 
#2769    3563    1129    1661    1388     711    7136    2840    1201    4699 
#tumor_5 tumor_6 
#8312    8914

#get the raw counts
raw.count <- as.matrix(STB_1356@meta.data)

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/before_clean/Seruat_STB_1356_object.RDS"
saveRDS(STB_1356, file = myoutf)



#########Extract CD8+ T cells only and re-clustering all CD8+ T cells#########

myinf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/before_clean/Seruat_STB_1356_object.RDS"
STB_1356 <- readRDS(myinf)

STB_1356[["percent.mt"]] <- PercentageFeatureSet(STB_1356, pattern = "^MT-")


myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/before_clean/QC_control_nGene_UMI_mito.pdf"
pdf(myoutf,width=20,height=5)
VlnPlot(STB_1356, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/before_clean/QC_control_correlation.pdf"
pdf(myoutf,width=10,height=5)
plot1 <- FeatureScatter(STB_1356, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(STB_1356, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

STB_1356 <- subset(STB_1356, subset = nCount_RNA<20000 & nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

table(STB_1356@meta.data$orig.ident)

#blood_1 blood_3 blood_5 blood_6  skin_1  skin_3  skin_5  skin_6 tumor_1 tumor_3 
#2560    2658    1047    1556     565     439    4846    1754     938    3704 
#tumor_5 tumor_6 
#7449    4488

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/before_clean/QC_control_correlation_afterfilter.pdf"
pdf(myoutf,width=10,height=5)
plot1 <- FeatureScatter(STB_1356, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(STB_1356, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/before_clean/nGene_UMI_mito_afterfilter.pdf"
pdf(myoutf,width=15,height=5)
VlnPlot(STB_1356, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


STB_1356 <- SCTransform(STB_1356, vars.to.regress = c("percent.mt",'orig.ident'), verbose = T,return.only.var.genes = F)



#plot for showing the normalization is working
tmp <- GetAssayData(object = STB_1356)
tmp <- as.matrix(tmp)
avg <- apply(tmp, 1, mean)
tag <- avg >=1
tmp <- tmp[tag,]

tag1 <- grep("skin_6",colnames(tmp))
tag2 <- grep("skin_1",colnames(tmp))
tag3 <- grep("skin_3",colnames(tmp))
tag4 <- grep("skin_5",colnames(tmp))
tag5 <- grep("blood_6",colnames(tmp))
tag6 <- grep("blood_1",colnames(tmp))
tag7 <- grep("blood_3",colnames(tmp))
tag8 <- grep("blood_5",colnames(tmp))
tag9 <- grep("tumor_6",colnames(tmp))
tag10 <- grep("tumor_1",colnames(tmp))
tag11 <- grep("tumor_3",colnames(tmp))
tag12 <- grep("tumor_5",colnames(tmp))



res1 <- tmp[,sample(tag1,10)]
res2 <- tmp[,sample(tag2,10)]
res3 <- tmp[,sample(tag3,10)]
res4 <- tmp[,sample(tag4,10)]
res5 <- tmp[,sample(tag5,10)]
res6 <- tmp[,sample(tag6,10)]
res7 <- tmp[,sample(tag7,10)]
res8 <- tmp[,sample(tag8,10)]

res9 <- tmp[,sample(tag9,10)]
res10 <- tmp[,sample(tag10,10)]
res11 <- tmp[,sample(tag11,10)]
res12 <- tmp[,sample(tag12,10)]



res <- cbind(res1,res2,res3,res4,res9,res10,res11,res12,res5,res6,res7,res8)

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/Validation_normalization.pdf"
pdf(myoutf,width=80,height=8)
boxplot(res[,1:ncol(res)])
dev.off()


myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/sctransformed_CD8_STB1356.RDS"
saveRDS(CD8_STB1356, file = myoutf)





myoutf <-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/densityplot_CD8A_STB1356.pdf"
pdf(myoutf,width=15,height=15)
hist(GetAssayData(object = STB_1356)["CD8A",])
dev.off()





CD8_STB1356<-subset(STB_1356,subset = CD8A>0&CD4<0.01&CD79A<0.01&CD19<0.01&FOXP3<0.01)

table(CD8_STB1356@meta.data$orig.ident)

#blood_1 blood_3 blood_5 blood_6  skin_1  skin_3  skin_5  skin_6 tumor_1 tumor_3 
#1087    1167     387    1113     236     196    2461    1268     704    1618 
#tumor_5 tumor_6 
#1330    1210

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/CD8_STB1356_clean_CD4.RDS"
saveRDS(CD8_STB1356, file = myoutf)

myinf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/CD8_STB1356_clean_CD4.RDS"
CD8_STB1356<-readRDS(myinf)







#Extract cells also got paired TCR sequenced#

myinf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/CD8_STB1356_clean_CD4.RDS"
CD8_STB1356<-readRDS(myinf)

tag<-row.names(CD8_STB1356@meta.data)%in%barcode_TCRAB

CD8_STB1356@meta.data[,"TCR"]<-rep("No",nrow(CD8_STB1356@meta.data))

CD8_STB1356@meta.data$TCR[tag]<-"Yes"

#nrow(CD8_STB1356@meta.data) 12777
#sum(CD8_STB1356@meta.data$TCR=="Yes") 10658



Idents(CD8_STB1356)<-"TCR"
CD8_STB1356<-subset(CD8_STB1356,idents="Yes")

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/CD8_STB1356_clean_CD4.RDS"
saveRDS(CD8_STB1356, file = myoutf)

myinf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/CD8_STB1356_clean_CD4.RDS"
CD8_STB1356<-readRDS(myinf)


CD8_STB1356 <- SCTransform(CD8_STB1356, vars.to.regress = c("percent.mt",'orig.ident'), verbose = T,return.only.var.genes = F)


VariableFeatures(object = CD8_STB1356)<-VariableFeatures(object = CD8_STB1356)[!grepl("TRA",VariableFeatures(object = CD8_STB1356))]
VariableFeatures(object = CD8_STB1356)<-VariableFeatures(object = CD8_STB1356)[!grepl("TRB",VariableFeatures(object = CD8_STB1356))]
VariableFeatures(object = CD8_STB1356)<-VariableFeatures(object = CD8_STB1356)[!grepl("IGH",VariableFeatures(object = CD8_STB1356))]
VariableFeatures(object = CD8_STB1356)<-VariableFeatures(object = CD8_STB1356)[!grepl("IGK",VariableFeatures(object = CD8_STB1356))]
VariableFeatures(object = CD8_STB1356)<-VariableFeatures(object = CD8_STB1356)[!grepl("IGL",VariableFeatures(object = CD8_STB1356))]
VariableFeatures(object = CD8_STB1356)<-VariableFeatures(object = CD8_STB1356)[!grepl("TRG",VariableFeatures(object = CD8_STB1356))]



myoutf <-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/densityplot_CD8A_CD8_STB1356.pdf"
pdf(myoutf,width=15,height=15)
hist(GetAssayData(object = CD8_STB1356)["CD8A",])
dev.off()



range(GetAssayData(object = CD8_STB1356))
mean(GetAssayData(object = CD8_STB1356, slot = "scale.data"))




CD8_STB1356 <- RunPCA(CD8_STB1356, features = VariableFeatures(object = CD8_STB1356))


myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/CD8_STB1356_PCA_genes_fig1_.pdf"
pdf(myoutf,width=25,height=25)
VizDimLoadings(CD8_STB1356, dims = 1:30, reduction = "pca")
DimHeatmap(CD8_STB1356, dims = 1:30, cells = 1000, balanced = TRUE)
dev.off()

CD8_STB1356 <- JackStraw(CD8_STB1356, num.replicate = 100,dims=30)
CD8_STB1356 <- ScoreJackStraw(CD8_STB1356, dims = 1:30)

myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/CD8_STB1356_Sig_PCA_components_noregreetissue.pdf"
pdf(myoutf,width=25,height=20)
JackStrawPlot(CD8_STB1356, dims = 1:30)
ElbowPlot(CD8_STB1356,ndims = 30)
dev.off()

print(CD8_STB1356[["pca"]], dims = 1:30)




x<-c(1:15)


eigs <- CD8_STB1356[["pca"]]@stdev^2 #PCA的STBandard dev的平方#
print(paste0("Variance captured by 35 PCs: ",sum(eigs[x] / sum(eigs)))) #0.831#




#Do clustering
CD8_STB1356<-FindNeighbors(CD8_STB1356, dims = x)
CD8_STB1356 <- FindClusters(CD8_STB1356,resolution = seq(0.1,2.0,0.1))

PrintFindClustersParams(object = CD8_STB1356)

clusters_resolution<-sapply(grep("^SCT_snn_res",colnames(CD8_STB1356@meta.data),value = TRUE),
                            function(x) length(unique(CD8_STB1356@meta.data[,x])))
clusters_resolution


#PC=x CD8_STB1356
#SCT_snn_res.0.1 SCT_snn_res.0.2 SCT_snn_res.0.3 SCT_snn_res.0.3 SCT_snn_res.0.5 
#5              10              10              11              13 
#SCT_snn_res.0.6 SCT_snn_res.0.7 SCT_snn_res.0.8 SCT_snn_res.0.9   SCT_snn_res.1 
#16              17              19              20              20 
#SCT_snn_res.1.1 SCT_snn_res.1.2 SCT_snn_res.1.3 SCT_snn_res.1.4 SCT_snn_res.1.5 
#22              23              23              24              24 
#SCT_snn_res.1.6 SCT_snn_res.1.7 SCT_snn_res.1.8 SCT_snn_res.1.9   SCT_snn_res.2 
#25              25              26              28              28 

#Set up a loop to look for markers

tar <- c("SCT_snn_res.0.1","SCT_snn_res.0.2","SCT_snn_res.0.3","SCT_snn_res.0.4",
         "SCT_snn_res.0.5","SCT_snn_res.0.6","SCT_snn_res.0.7","SCT_snn_res.0.8",
         "SCT_snn_res.0.9","SCT_snn_res.1","SCT_snn_res.1.1",'SCT_snn_res.1.2',
         'SCT_snn_res.1.3','SCT_snn_res.1.4',"SCT_snn_res.1.5","SCT_snn_res.1.8",
         "SCT_snn_res.1.9",'SCT_snn_res.2')


tar1 <- c("SCT_snn_res.0.3","SCT_snn_res.0.4","SCT_snn_res.0.2","SCT_snn_res.0.6",
          "SCT_snn_res.0.5","SCT_snn_res.0.1")

for(i in 1 : length(tar1))
{
  cat("\r",i)
  Idents(object = CD8_STB1356) <- tar1[i]
  tmp <- CD8_STB1356
  tmp.markers <- FindAllMarkers(object = tmp, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1,assay="SCT",return.thresh = 0.05)
  myoutf <- paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/cluster_markers_",tar1[i],"_ScRNA.xls")
  write.table(tmp.markers,myoutf,sep="\t",quote=F)
  myoutf1 <- paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/heatmap_markers_",tar1[i],".pdf")
  pdf(myoutf1,30,35)
  top20 <- tmp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
  print(DoHeatmap(tmp, features = top20$gene) + NoLegend())
  dev.off()
}


CD8_STB1356 <- RunUMAP(CD8_STB1356, dims= x)



for(i in 1:length(tar))
{
  Idents(object = CD8_STB1356) <- tar[i]
  myoutf1 <- paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/UMAP_SB6_",tar[i],"_undefined.pdf")
  pdf(myoutf1,7,5)
  print(DimPlot(CD8_STB1356, reduction = "umap"))
  dev.off()
}



tag<-grep("skin",CD8_STB1356@meta.data$orig.ident)

CD8_STB1356@meta.data[,'label']<-rep('tumor',nrow(CD8_STB1356@meta.data))
CD8_STB1356@meta.data$label[tag]<-"skin"
tag1<-grep("blood",CD8_STB1356@meta.data$orig.ident)
CD8_STB1356@meta.data$label[tag1]<-"blood"


tag1<-grep("_1",CD8_STB1356@meta.data$orig.ident)
tag2<-grep("_3",CD8_STB1356@meta.data$orig.ident)
tag3<-grep("_4",CD8_STB1356@meta.data$orig.ident)
tag4<-grep("_5",CD8_STB1356@meta.data$orig.ident)
tag5<-grep("_6",CD8_STB1356@meta.data$orig.ident)


CD8_STB1356@meta.data[,'PT_number']<-rep('PT',nrow(CD8_STB1356@meta.data))
CD8_STB1356@meta.data$PT_number[tag1]<-"PT_1"
CD8_STB1356@meta.data$PT_number[tag2]<-"PT_3"
CD8_STB1356@meta.data$PT_number[tag3]<-"PT_4"
CD8_STB1356@meta.data$PT_number[tag4]<-"PT_5"
CD8_STB1356@meta.data$PT_number[tag5]<-"PT_6"



myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/TSNE_CD8_STB1356_groupby_tissue_noregreetissue.pdf"
pdf(myoutf,width=8,height=5)
DimPlot(object = CD8_STB1356,group.by='label',reduction="umap")
dev.off()


myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/TSNE_CD8_STB1356_groupby_PT_number.pdf"
pdf(myoutf,width=8,height=5)
DimPlot(object = CD8_STB1356,group.by='PT_number',reduction="umap",pt.size=0.6)
dev.off()



myoutf <- "/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/Seruat_CD8_STB1356__finished.RDS"
saveRDS(CD8_STB1356, file = myoutf)
myinf<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/Seruat_CD8_STB1356__finished.RDS"
CD8_STB1356<-readRDS(myinf)


########Feature plot figure 1c########
residency_blood_skin_difference<-c("S1PR1","KLF2","SELL","TCF7","CD69",'RGS1',"ITGAE","PDCD1")


p<-FeaturePlot(object = CD8_STB1356, features = Trm_for_paper,cols= c('#deebf7', "#bd0026"),
               reduction = "umap",slot = "data",ncol=2,max.cutoff=4,pt.size=0.5,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(axis.title = element_text(family="Helvetica",face = "bold", color = "black",size=20), 
                           axis.text=element_text(family="Helvetica",face = "bold", color = "black",size=16),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.ticks = element_line(colour = "black", size = 1),
                           axis.ticks.length=unit(.3, "cm"))
}

myoutf <-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/Trm_for_paper.pdf"
pdf(myoutf,width=11.7,height=10)
print(cowplot::plot_grid(plotlist = p))
dev.off()


##########Heatmap fig 1f############

tar<-c("2",'5','4','8','9','3','6','0','1','7')
all_cell<-character()
for(i in 1:length(tar))
{
  cat("\r",i)
  len<-sum(CD8_STB1356@meta.data$res.0.4_combined==tar[i])
  tag<-sample(1:len,242)
  info<-CD8_STB1356@meta.data[CD8_STB1356@meta.data$res.0.4_combined==tar[i],]
  barcode_use<-row.names(info)[tag]
  all_cell<-combine(all_cell,barcode_use)
}

myoutf1 <- paste0("/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/cluster_identity/heatmap_marker_genes_ordered_242cells_each_cluster.pdf")


Idents(CD8_STB1356)<-"res.0.4_combined"

tmp.markers <- FindAllMarkers(object = CD8_STB1356, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.15,return.thresh=0.05)

ts<-tmp.markers[order(factor(tmp.markers$cluster, levels = tar)),]

top30 <- ts %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)

CD8_STB1356@meta.data[,"ordered_ident"]<-factor(CD8_STB1356@meta.data$res.0.4_combined,levels=tar)

pdf(myoutf1,12,10)
print(DoHeatmap(CD8_STB1356, cells=all_cell,features = top30$gene,group.by="ordered_ident") + 
        scale_fill_gradient2( low = (c('#0702a6',"#4575b4",'#67a9cf')), 
                              mid = "white", high = (c('#fdae61','#f46d43','#ff001f')), 
                              midpoint = 0, guide = "colourbar", aesthetics = "fill"))
dev.off()


############Cluster distribution to different tissues, Ro/e###############

myinf<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/Seruat_CD8_STB1356__finished.RDS"
CD8_STB1356<-readRDS(myinf)

tar<-unique(CD8_STB1356@meta.data$res.0.4_combined)
tar<-tar[order(tar)]
chi<-matrix(0,12,length(tar))
number<-sort(rep((unique(CD8_STB1356@meta.data$PT_number)),3))
tissue<-rep(c('skin','tumor','blood'),4)
row.names(chi)<-paste0(number,tissue)
colnames(chi)<-tar

tar2<-unique(CD8_STB1356@meta.data$PT_number)
tar1<-unique(CD8_STB1356@meta.data$label)

for(j in 1:length(tar2))
{
  yy<-CD8_STB1356@meta.data[CD8_STB1356@meta.data$PT_number==tar2[j],]
  
  for (i in 1:length(tar1))
  {
    tag<-yy$label==tar1[i]
    xx<-yy[tag,]
    for (k in 1:length(tar))
    {
      chi[paste0(tar2[j],tar1[i]),tar[k]]<-sum(xx$res.0.4_combined==tar[k])
    }
    
  }
  
  
}


chi_square<-chisq.test(chi)

chi_observed<-chi_square$observed
chi_expected<-chi_square$expected
ROE<-chi_observed/chi_expected

myoutf<-"/dartfs-hpc/rc/lab/T/TurkLab/Jichang/scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/cluster_tissue_distribution/ROE.xls"
write.table(ROE,myoutf,sep="\t",quote=F)





