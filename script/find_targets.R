
find_tagets<-function(enz_LID_info,iso_LID_info,iso_decr,add_tar,data_cancer_iso,data_normal_iso){
  
  
  target_expr_descr<-matrix(nrow=length(enz_LID_info$EC_number),ncol=6)
  rownames(target_expr_descr)<-enz_LID_info$EC_number
  colnames(target_expr_descr)<-c("EC_number","name","iso_num","target_iso","uniprot","other_iso")
  target_expr_descr<-as.data.frame(target_expr_descr)
  
  target_expr_descr$EC_number<-enz_LID_info$EC_number
  target_expr_descr$name<-enz_LID_info$Name
  target_expr_descr$iso_num<-enz_LID_info$iso_num
  for (i in 1:nrow(target_expr_descr)){
    if(length(setdiff(unlist(enz_LID_info$isozyme_list[i]),iso_decr)) >0){
      target_expr_descr$target_iso[i]<-setdiff(unlist(enz_LID_info$isozyme_list[i]),iso_decr)
      target_expr_descr$uniprot[i]<-iso_LID_info[target_expr_descr$target_iso[i],"swissprot"]
      target_expr_descr$other_iso[i]<-paste(intersect(unlist(enz_LID_info$isozyme_list[i]),iso_decr),collapse=",")
    }
    else {
      target_expr_descr$target_iso[i]<-intersect(unlist(enz_LID_info$isozyme_list[i]),add_tar)
      target_expr_descr$uniprot[i]<-iso_LID_info[target_expr_descr$target_iso[i],"swissprot"]
      target_expr_descr$other_iso[i]<-paste(setdiff(unlist(enz_LID_info$isozyme_list[i]),add_tar),collapse=",")
    }
  }

  
  target_expr_mean<-matrix(nrow=length(target_expr_descr$target_iso),ncol=2)
  colnames(target_expr_mean)<-c("target_normal","target_cancer")
  target_expr_mean<-as.data.frame(target_expr_mean)
  for (i in 1:nrow(target_expr_mean)){
    target_expr_mean[i,"target_normal"]<-mean(c(t(data_normal_iso[target_expr_descr$target_iso[i],])))
    target_expr_mean[i,"target_cancer"]<-mean(c(t(data_cancer_iso[target_expr_descr$target_iso[i],])))
  }
  
  
  other_expr_max<-matrix(nrow=length(target_expr_descr$other_iso),ncol=2)
  colnames(other_expr_max)<-c("other_normal","other_cancer")
  other_expr_max<-as.data.frame(other_expr_max)
  for (i in 1:nrow(other_expr_max)){
    iso_list<-unlist(strsplit(target_expr_descr[i,"other_iso"],","))
    other_expr_max[i,"other_normal"]<-max(apply(data_normal_iso[iso_list,],1,mean))
    other_expr_max[i,"other_cancer"]<-max(apply(data_cancer_iso[iso_list,],1,mean))
  }
  
  
  l2FC_iso<-matrix(nrow=length(target_expr_descr$EC_number),ncol=3)
  colnames(l2FC_iso)<-c("target_other_normal","target_other_cancer","cancer_other")
  l2FC_iso<-as.data.frame(l2FC_iso)
  l2FC_iso$target_other_normal<-target_expr_mean$target_normal - other_expr_max$other_normal
  l2FC_iso$target_other_cancer<-target_expr_mean$target_cancer - other_expr_max$other_cancer
  l2FC_iso$cancer_other<-l2FC_iso$target_other_cancer - l2FC_iso$target_other_normal
  
  
  score_all<-matrix(nrow=length(target_expr_descr$EC_number),ncol=4)
  rownames(score_all)<-target_expr_descr$EC_number
  colnames(score_all)<-c("S1","S2","S3","S4")
  score_all<-as.data.frame(score_all)
  score_all$S1<-l2FC_iso$target_other_cancer
  score_all$S2<-other_expr_max$other_normal - other_expr_max$other_cancer
  score_all$S3<-target_expr_mean$target_cancer - target_expr_mean$target_normal
  score_all$S4<-l2FC_iso$cancer_other
  
  score<-score_all
  score$S1<-rank(score$S1,ties.method = "average")
  score$S2<-rank(score$S2,ties.method = "average")
  score$S3<-rank(score$S3,ties.method = "average")
  score$S4<-rank(score$S4,ties.method = "average")
  score<-score%>%
    mutate(score_sum = rowSums(.[1:4]))
  score$score_sum <- 100*(score$score_sum - min(score$score_sum))/(max(score$score_sum)-min(range(score$score_sum)))
  
  
  LID_res<-cbind(target_expr_descr,score$score_sum,score_all)
  colnames(LID_res)[7]<-"score"
  LID_res<-LID_res[order(LID_res$score,decreasing = T),]
  LID_res<-LID_res[,-7]
  
  return(LID_res)
}
