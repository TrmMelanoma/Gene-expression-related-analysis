##########Generate marker genes list of each cluster for survival analysis#############
myinf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/Seruat_CD8_STB1356__finished.RDS"
CD8_STB1356<-readRDS(myinf)

Idents(CD8_STB1356)<-'res.0.4_combined'

tmp.markers<-FindMarkers(CD8_STB1356,
                         ident.1=c('0'),logfc.threshold = 0.01)
myoutf <- "../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/non_tumor_targeting_Trm_list.xls"
write.table(tmp.markers,myoutf,sep="\t",quote=F)

tmp.markers<-FindMarkers(CD8_STB1356,
                         ident.1=c('7'),logfc.threshold = 0.01)
myoutf <- "../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/Tex_list.xls"
write.table(tmp.markers,myoutf,sep="\t",quote=F)

tmp.markers<-FindMarkers(CD8_STB1356,
                         ident.1=c('1'),logfc.threshold = 0.01)
myoutf <- "../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/tumor_targeting_Trm_list.xls"
write.table(tmp.markers,myoutf,sep="\t",quote=F)

#########K-M plot and Forest plot using the TCGA data###############
#K-M plot
rm(list=ls())
myinf1 = "../SKCM_RNAseqv2_ALL_Symbol.rda" #gene expression
myinf2 = "../SKCM_Clincial_info.txt" #clinical information
myinf4 <- "../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/survival/signatures/" #gene signature

load(myinf1)
mydata = log10(mydata+1)

tag = combine(grep("06A", colnames(mydata)))
mydata = mydata[,tag]


files = list.files(myinf4)

res = matrix(0, ncol(mydata), 8)
row.names(res) = colnames(mydata)
colnames(res) = gsub("_list.xls", "",files)
res = as.data.frame(res)

for(i in 1 : length(files))
{
  cat("\r",i)
  tmpinf = paste0(myinf4, files[i])
  tmp = read.table(tmpinf,sep="\t",quote=NULL)
  
  tag = tmp$p_val_adj < 0.05
  tmp = tmp[tag,]
  tag1 = tmp$avg_logFC >=0.25
  tmp_UP = tmp[tag1,]
  
  
  if(nrow(tmp_UP) < 10)
  {
    print("Warning Signature less than 100 gene")
  }
  
  
  tmp_UP=tmp_UP[order(-tmp_UP$avg_logFC),]
  
  com_UP = intersect(combine(c("CD8A","CD8B"),row.names(tmp_UP)), row.names(mydata))
  tmp_data_UP = mydata[com_UP,]
  tmp_data_UP = apply(tmp_data_UP,1, function(x) (x-mean(x))/sd(x))
  score_UP = apply(t(tmp_data_UP), 2, mean)
  
  score<-score_UP
  
  res[,i] = score
}

res<-as.matrix(res[, colSums(is.na(res)) != nrow(res)])

#process the clinical information for K-M plot
info = read.table(myinf2,sep="\t",quote=NULL)
tar <- c("days_to_death","days_to_last_followup","person_neoplasm_cancer_status","vital_status")
info <- info[,tar]

tag <- !is.na(info$vital_status)
info <- info[tag,]

info[,"Time"] <- rep("NA",nrow(info))

tag1 <- info$vital_status == "dead"
tag2 <- info$vital_status == "alive"

info[tag1,"Time"] <- info[tag1, "days_to_death"]
info[tag2,"Time"] <- info[tag2, "days_to_last_followup"]

tag <- !is.na(info$Time)
info <- info[tag,]
info$Time = as.numeric(info$Time)
info$Time = info$Time/365

info[,"status"] <- rep(0,nrow(info))
tag <- info$vital_status == "dead"
info[tag,"status"] = "1"
info$status = as.numeric(info$status)

library(survminer)

for(i in 1 : ncol(res))
{
  cat("\r",i)
  
  tmp_res = res[,i]
  tmp_res = ifelse(tmp_res > median(tmp_res), 1, 0)
  sam = names(tmp_res)
  
  for(k in 1 : length(sam))
  {
    tmp_sam = strsplit(sam[k],"-")[[1]]
    names(tmp_res)[k] = paste(tmp_sam[1],tmp_sam[2],tmp_sam[3],sep="-")
  }
  
  tag1 = tmp_res == 1 
  tag2 = tmp_res == 0 
  
  sam1 = names(tmp_res)[tag1]
  sam2 = names(tmp_res)[tag2]
  
  com = intersect(row.names(info), names(tmp_res))
  
  xx = info[com,]
  xx[sam1,"group"] = "1"
  xx[sam2,"group"] = "2"
  
  fit2 <- survfit( Surv(Time, status) ~ group,
                   data = xx )
  
  labs = c(paste0(colnames(res)[i],"_Hi"),paste0(colnames(res)[i],"_Lo"))
  
  myoutf = paste0("../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/survival/TCGA/Mets/cebersort/0.25/TCGA_",colnames(res)[i],".pdf")
  pdf(myoutf, width=15, height=30)
  p= ggsurvplot(fit2, pval = TRUE, 
                break.time.by = 10,
                risk.table = TRUE,
                risk.table.height = 0.5,
                palette = c("royalblue4", "tomato3"),
                legend.labs = labs)
  p = p + guides(colour = guide_legend(nrow = 2))
  print(p)
  dev.off()
  
}	

#Forest Plot
rm(list=ls())
myinf1 = "../SKCM_RNAseqv2_ALL_Symbol.rda" #gene expression
myinf2 = "../SKCM_Clincial_info.txt" #clinical information
myinf4 <- "../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/survival/signatures/" #gene signature

load(myinf1)
mydata = log10(mydata+1)

tag = grep("06A", colnames(mydata))
mydata = mydata[,tag]

files = list.files(myinf4)

res = matrix(0, ncol(mydata), 8)
row.names(res) = colnames(mydata)
colnames(res) = gsub("_list.xls", "",files)
res = as.data.frame(res)

for(i in 1 : length(files))
{
  cat("\r",i)
  tmpinf = paste0(myinf4, files[i])
  tmp = read.table(tmpinf,sep="\t",quote=NULL)
  tag = tmp$p_val_adj < 0.05
  tmp = tmp[tag,]
  tag = tmp$avg_logFC >= 0.25
  tmp = tmp[tag,]
  if(nrow(tmp) < 10)
  {
    print("Warning Signature less than 100 gene")
  }
  com = intersect(combine(row.names(tmp),c("CD8A","CD8B")), row.names(mydata))
  
  tmp_data = mydata[com,]
  
  tmp_data = apply(tmp_data,1, function(x) (x-mean(x))/sd(x))
  
  score = apply(t(tmp_data), 2, mean)
  res[,i] = score
}


com1<-combine("CD6", "CD3D", "CD3E", "SH2D1A", "TRAT1", "CD3G")
tmp_data1 = mydata[com1,]
tmp_data1 = apply(tmp_data1,1, function(x) (x-mean(x))/sd(x))
score1 = apply(tmp_data1, 1, mean)
res[,9] = score1

colnames(res)[9]<-c("CD3")
data<-as.matrix(res)
tag = grep("06A", row.names(data))
data = data[tag,]
sam = row.names(data)

for(i in 1 : length(sam))
{
  tmp_sam = strsplit(sam[i],"-")[[1]]
  sam[i] = paste(tmp_sam[1],tmp_sam[2],tmp_sam[3],sep="-")
}

row.names(data) = sam
myinf2 <- "../SKCM_Clincial_info.txt" 
info <- read.table(myinf2,sep="\t",quote=NULL,row.names=1,header=T,stringsAsFactors=F)
com = intersect(row.names(info), row.names(data))
info = info[com,]
data = data[com,]
tar <- c("days_to_death","days_to_last_followup","person_neoplasm_cancer_status","vital_status","age_at_initial_pathologic_diagnosis","gender","stage_event.pathologic_stage")
info <- info[,tar]
tag <- !is.na(info$vital_status)
info <- info[tag,]
xx<-ifelse(is.na(info$days_to_death),0,info$days_to_death)
yy<-ifelse(is.na(info$days_to_last_followup),0,info$days_to_last_followup)
info[,"Time"] <- xx+yy
tag <- info$Time!=0
info <- info[tag,]
info[,"status"] <- rep(0,nrow(info))
tag <- info$vital_status == "dead"
info[tag,"status"] = "1"
comSam = intersect(row.names(data), row.names(info))
data = data[comSam,]
info = info[comSam,]
mytf1 <- as.numeric(data[,"tumor_targeting_Trm"])
mytf2 <- as.numeric(data[,"Tex"])
mytf5 <- as.numeric(data[,"CD3"])


xx = cbind(info, mytf1,mytf2,mytf5)
info$stage_event.pathologic_stage<-ifelse(is.na(info$stage_event.pathologic_stage),"Unknown",info$stage_event.pathologic_stage)
tag1 <- info$stage_event.pathologic_stage == "stage ia" | info$stage_event.pathologic_stage == "stage i" | info$stage_event.pathologic_stage == "stage ib"
tag2 <- info$stage_event.pathologic_stage == "stage ii" | info$stage_event.pathologic_stage == "stage iia" | info$stage_event.pathologic_stage == "stage iib" | info$stage_event.pathologic_stage == "stage iic"
tag3 <- info$stage_event.pathologic_stage == "stage iii" | info$stage_event.pathologic_stage == "stage iiia" | info$stage_event.pathologic_stage == "stage iiib" | info$stage_event.pathologic_stage == "stage iiic"
tag4 <- info$stage_event.pathologic_stage == "stage iv" | info$stage_event.pathologic_stage == "stage iva" | info$stage_event.pathologic_stage == "stage ivb" | info$stage_event.pathologic_stage == "stage ivc"

xx1 <-xx[tag1|tag2,]
xx1[,"stage"] <- "I/II"
xx2 <- xx[tag3|tag4,]
xx2[,"stage"] <- "III/IV"
xx <- rbind(xx1,xx2)
time <- as.numeric(xx$Time)
xx = cbind(time,xx)
xx$status <- as.numeric(xx$status)
xx$Time <- as.numeric(xx$Time)
xx$Time <- xx$Time/365
xx$mytf5 <- ifelse(xx$mytf5>median(xx$mytf5),1,0)
tag <- !is.na(xx$Time)
xx <- xx[tag,]
xx[,"age"] <- ifelse(xx$age_at_initial_pathologic_diagnosis>=65,1,0)
mycox <- coxph(formula = Surv(Time, status) ~ mytf1+mytf2+mytf5+stage+gender, data = xx)
cox <- mycox
beta <- coef(cox)
beta <- exp(beta)
beta <- signif(beta,2)
se   <- sqrt(diag(cox$var))
p    <- 1 - pchisq((coef(cox)/se)^2, 1)
p    <- signif(p,1)
CI   <- round(confint(cox), 4)
CI   <- exp(CI)
CI   <- signif(CI,2)
res <- cbind(beta, se = exp(beta), CI, p)
res <- as.data.frame(res)
res[,"label"] <- row.names(res)
colnames(res)[1:4] <- c("mean","se","lower","upper")
res[,"label"] <- factor(res[,"label"],level=c("gendermale","stageIII/IV","mytf5","mytf2","mytf1"))

library(ggplot2)
myoutf <-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/survival/TCGA/Mets/nocebersort/0.25/forest_TrmTOX_Tcell.pdf"
pdf(myoutf, width= 4, height= 3)
p <- ggplot(data=res, aes(x=label, y=mean,ymin=lower,ymax=upper))  +
  geom_pointrange(shape=22,fatten=6) + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)+
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Hazard ratio") + ggtitle("TCGA")+
  theme_bw() + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line.x = element_line(colour = "black"),axis.title=element_text(size=10),plot.title=element_text(size=10,hjust = 0.5,vjust=0.5),axis.text.x=element_text(size= 10,hjust=0.5))# use a white background
p
dev.off()







