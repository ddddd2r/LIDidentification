#####################################################################################################################
## Code: Identification of LID targets and dependency validation

## Contents:
## 0. Preparation and Load data 
## 1. Data filtering and processing
## 2. Analysis for identification of LID targets
## 3. Broad DepMap validation
## 4. LID target selection and data output

#####################################################################################################################
## 0. Preparation and Load data  #####
#####################################################################################################################

rm(list = ls())
options(stringsAsFactors = F)


library(stringr)
library(dplyr)
library(tidyr)


load("isozyme_data.Rdata")
load("FUSCC_RNA_data.Rdata")
load("Breast_CRISPR_depscore.Rdata")


#####################################################################################################################
## 1. Data filtering and processing  #####
#####################################################################################################################

data_cancer_iso<-FUSCC_RNA_matched_T[isozyme_info$gene,]
data_normal_iso<-FUSCC_RNA_matched_N[isozyme_info$gene,]

enz_num <- length(enzyme_info[[1]])
iso_num <- length(isozyme_info[,1])

iso_comp<-matrix(nrow=iso_num,ncol=5)
rownames(iso_comp)<-rownames(isozyme_info)
colnames(iso_comp)<-c("mean_cancer","mean_normal","l2FC","p_wilcox","fdr_wilcox")
iso_comp<-as.data.frame(iso_comp)

for (i in 1:iso_num){
  iso_comp$mean_cancer[i]<-mean(as.numeric(data_cancer_iso[i,]))
  iso_comp$mean_normal[i]<-mean(as.numeric(data_normal_iso[i,]))
  iso_comp$l2FC[i] <- mean(c(t(data_cancer_iso[i,]))) - mean(c(t(data_normal_iso[i,])))
  
  data_table<-matrix(nrow=ncol(data_cancer_iso),ncol=2)
  colnames(data_table)<-c("cancer","normal")
  rownames(data_table)<-colnames(data_cancer_iso)
  data_table<-as.data.frame(data_table)
  data_table$normal <- c(t(data_normal_iso[i,]))
  data_table$cancer <- c(t(data_cancer_iso[i,]))
  
  res.wil <- wilcox.test(data_table$cancer, data_table$normal, paired=TRUE)
  iso_comp$p_wilcox[i]<-res.wil$p.value
}
iso_comp$fdr_wilcox<-p.adjust(iso_comp$p_wilcox, method ="fdr")


#####################################################################################################################
## 2. Analysis for identification of LID targets  #####
#####################################################################################################################

### 2.1 find enzymes with only one increased/not change isozyme in cancer #####

sig_lev <- 0.05
iso_decr<- rownames(iso_comp[which(iso_comp$fdr_wilcox < sig_lev & iso_comp$l2FC < 0),])

enz_tmp<-c()
for (i in iso_decr){
  enz_iso<-unlist(strsplit(isozyme_info[i,"EC_number"],", "))
  enz_tmp<-c(enz_tmp,enz_iso)
}
enz_tmp<-unique(enz_tmp)

enz_LID<-c()
for(i in enz_tmp){
  iso_tmp<-unlist(enzyme_info$isozyme_list[which(enzyme_info$EC_number == i)])
  if(length(intersect(iso_tmp,iso_decr)) == length(iso_tmp)-1){
    enz_LID<-c(enz_LID,i)
  }
}


add_LID<-c()
add_tar<-c()
for(i in enzyme_info$EC_number){
  iso_tmp<-unlist(enzyme_info$isozyme_list[which(enzyme_info$EC_number == i)])
  comp_tmp<-iso_comp[iso_tmp,1:3]
  num_tar<-length(comp_tmp[comp_tmp$l2FC > -0.5,1])
  num_oth<-length(comp_tmp[comp_tmp$l2FC < -0.5,1])
  if (num_tar == 1 & num_oth == (length(iso_tmp)-1)){
    if (max(comp_tmp[,"mean_cancer"])== comp_tmp[comp_tmp$l2FC > -0.5,1]){
      add_LID<-c(add_LID,i)
      add_tar<-c(add_tar,rownames(comp_tmp[comp_tmp$l2FC < -0.5,]))
    }
  }
}

enz_LID<-c(enz_LID,add_LID)

enz_LID_info<-list(enzyme_info$EC_number[which(enzyme_info$EC_number %in% enz_LID)],
                   enzyme_info$Name[which(enzyme_info$EC_number %in% enz_LID)],
                   enzyme_info$iso_num[which(enzyme_info$EC_number %in% enz_LID)],
                   enzyme_info$isozyme_list[which(enzyme_info$EC_number %in% enz_LID)]
)
names(enz_LID_info)<-c("EC_number","Name","iso_num","isozyme_list")

iso_LID<-c()
for (i in 1:length(enz_LID)){
  iso_LID<-c(iso_LID,unlist(enz_LID_info$isozyme_list[i]))
}
iso_LID<-unique(iso_LID)

iso_LID_info<-isozyme_info[iso_LID,]
for (i in rownames(iso_LID_info)){
  iso_LID_info[i,"EC_number"]<-paste(enz_LID_info$EC_number[grep(i,enz_LID_info$isozyme_list)],collapse=", ")
  iso_LID_info[i,"enz_number"]<-length(enz_LID_info$EC_number[grep(i,enz_LID_info$isozyme_list)])
}


###  2.2 prioritize LID targets #####
LID_breast <- find_tagets(enz_LID_info,iso_LID_info,iso_decr,add_tar,data_cancer_iso,data_normal_iso)
LID_breast<-LID_breast %>% distinct(EC_number,.keep_all = T)


#####################################################################################################################
## 3. Broad DepMap validation  #####
#####################################################################################################################

Cus_taget_iso<-intersect(LID_breast$target_iso,rownames(Cus_CRISPR_dep_score))
LID_CRISPR<-cbind(Cus_taget_iso,Cus_CRISPR_dep_score[Cus_taget_iso,])

for(i in 1:nrow(LID_CRISPR)){
  a<-wilcox.test(c(t(LID_CRISPR[i,2:ncol(LID_CRISPR)])), mu = 0, alternative = "less")
  LID_CRISPR$p_0[i]<-a$p.value
}
LID_CRISPR$FDR_0<-p.adjust(LID_CRISPR$p_0,method="fdr")

LID_breast$FDR_0_CRISPR<-LID_CRISPR[LID_breast$target_iso,"FDR_0"]

LID_breast$dep_score_tar_CRISPR<-NA
LID_breast$dep_score_oth_CRISPR<-NA
LID_breast$p_tar_other_CRISPR<-NA
LID_breast$FDR_tar_other_CRISPR<-NA
for (i in 1:nrow(LID_breast)){
  if (LID_breast$target_iso[i] %in% rownames(Cus_CRISPR_dep_score)){
    other_iso<-unlist(strsplit(LID_breast[i,"other_iso"],","))
    if (any(other_iso %in% rownames(Cus_CRISPR_dep_score))){
      other_dep_score<-c(t(rowMeans(Cus_CRISPR_dep_score[other_iso,],na.rm=T)))
      min_other_iso<-other_iso[which.min(other_dep_score)]
      
      LID_breast$dep_score_tar_CRISPR[i]<-rowMeans(Cus_CRISPR_dep_score[LID_breast$target_iso[i],],na.rm = T)
      LID_breast$dep_score_oth_CRISPR[i]<-other_dep_score[which.min(other_dep_score)]
      
      a<-wilcox.test(c(t(Cus_CRISPR_dep_score[min_other_iso,])), c(t(Cus_CRISPR_dep_score[LID_breast$target_iso[i],])))
      LID_breast$p_tar_other_CRISPR[i]<-a$p.value
    }
  }
}
LID_breast$FDR_tar_other_CRISPR<-p.adjust(LID_breast$p_tar_other_CRISPR,method="fdr")

write.csv(LID_breast,"LID_breast.csv")

#####################################################################################################################
## 4. LID target selection and data output  #####
#####################################################################################################################

TableS3_LID_breast<-LID_breast[,1:6]
TableS3_LID_breast$S1<-LID_breast$S1
TableS3_LID_breast$S2<-LID_breast$S2
TableS3_LID_breast$S3<-LID_breast$S3
TableS3_LID_breast$S4<-LID_breast$S4
TableS3_LID_breast$step1<-(TableS3_LID_breast$S1 > 0.5 & TableS3_LID_breast$S2 > 0.5
                           & TableS3_LID_breast$S3 > -0.5 & TableS3_LID_breast$S4 > 0.7)
TableS3_LID_breast$TCGA_LID<-TableS3_LID_breast$EC_number %in% LID_TCGA$EC_number
TableS3_LID_breast$TCGA_step1<-TableS3_LID_breast$EC_number %in% LID_TCGA$EC_number[LID_TCGA$S1 > 0.5 & LID_TCGA$S2 > 0.5  
                                                                                    & LID_TCGA$S3 > -0.5 & LID_TCGA$S4 > 0.7]

TableS3_LID_breast$dep_target<-LID_breast$dep_score_tar_CRISPR
TableS3_LID_breast$dep_other<-LID_breast$dep_score_oth_CRISPR
TableS3_LID_breast$dep_target_fdr<-LID_breast$FDR_0_CRISPR
TableS3_LID_breast$dep_target_other_fdr<-LID_breast$FDR_tar_other_CRISPR

write.csv(TableS3_LID_breast,"TableS3_LID_breast.csv")
