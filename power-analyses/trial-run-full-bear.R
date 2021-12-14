# set library path
.libPaths("/home/adrianf/project-folder/nobackup/ADRIAN/metagenomics-power-analyses-tutorial/all-bears/r-env/lib/R/library")

## ---- install, eval=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## install.packages('devtools')
## devtools::install_github("micropower")
## install.packages('knitr')
## install.packages('kableExtra')
## install.packages('dplyr')
## install.packages('ggplot2')
## install.packages('parallel')

## ---- setup------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(micropower)
library(tidyverse)
library(parallel)
library(phyloseq)
cores = 20
set.seed(515087345)

## ----load-otu----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  otu<-readRDS("/home/adrianf/project-folder/nobackup/ADRIAN/metagenomics-power-analyses-tutorial/all-bears/RDS/amr-taxa-bear-0.01-af-bracken-max-highoral-0.01-otu-env-filtered.rds")
  samples<-data.frame(sample_data(otu)) %>% 
    filter(Spec.host=="Bear"&!is.na(Spec.HalfCentury)) %>%
    mutate(antib.era=ifelse(Spec.HalfCentury%in%c("c20.2","c21.1"),"post","pre"))

# iterate through simulations of the full sample list, or random sample of individuals from each group
for (i in c(1:10)) {
  for (j in c(5,10,15,20)){
  ss<-samples %>% 
    select(SampleID,antib.era) %>% 
    group_by(antib.era) %>% 
    slice_sample(n = j) %>% 
    pull(SampleID)
  
  df<-otu %>% 
    otu_table() %>% as.data.frame() %>%
    select(ss) %>% 
    filter(rowSums(.)!=0)
     
  ## ----within_dist,warnings=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------
  jaccard<-calcWJstudy(df)
  m<-mean(jaccard)
  s<-sd(jaccard)
  
  ## ----hashMean----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  weight_hm<-hashMean(rare_levels=runif(1000,0,1),rep_per_level=10,otu_number=1000,sequence_depth=1000)
  means<-mclapply(weight_hm, function(x) (mean(lowerTriDM(calcWJstudy(x)))),mc.cores=cores)
  names(means) <- sapply(strsplit(names(means),"_"),FUN=function(x) {x[[1]]})
  save(means,file=paste0("/home/adrianf/project-folder/nobackup/ADRIAN/metagenomics-power-analyses-tutorial/all-bears/data/p",j,".",i,"means.Rdata"))
  load(paste0("/home/adrianf/project-folder/nobackup/ADRIAN/metagenomics-power-analyses-tutorial/all-bears/data/p",j,".",i,"means.Rdata"))
  mean_df <- data.frame(subsampling=as.numeric(names(means)),means=as.numeric(means),target=abs(as.numeric(means)-m))
  subsample <- mean_df[which(mean_df$target==min(mean_df$target)),]$subsampling
  p1<-qplot(data=mean_df,x=subsampling,y=means)
  ggsave(plot=p1,file=paste0("/home/adrianf/project-folder/nobackup/ADRIAN/metagenomics-power-analyses-tutorial/all-bears/figures/p",j,".",i,"-means.png"),dpi=300)
  
  ## ---- hashSD, warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------
  weight_hsd<-hashSD(rare_depth=subsample,otu_number_range=10^runif(n = 100,min = 1,max = 5),sim_number=length(ss), sequence_depth=1000)
  sds<-mclapply(weight_hsd, function(x) (sd(lowerTriDM(calcWJstudy(x)))), mc.cores=cores)
  names(sds) <- sapply(names(sds), function(x) substring(x,4))
  save(sds,file=paste0("/home/adrianf/project-folder/nobackup/ADRIAN/metagenomics-power-analyses-tutorial/all-bears/data/p",j,".",i,"sds.Rdata"))
  load(paste0("/home/adrianf/project-folder/nobackup/ADRIAN/metagenomics-power-analyses-tutorial/all-bears/data/p",j,".",i,"sds.Rdata"))
  sds_df <- data.frame(otunum=as.numeric(names(sds)),sd=as.numeric(sds),target=abs(as.numeric(sds)-s))
  otunum <- sds_df[which(sds_df$target==min(sds_df$target)),]$otunum
  p2<-qplot(data=sds_df,x=otunum,y=sd) + geom_smooth(method="lm") + scale_x_log10() + scale_y_log10() + xlab("log10(otunum)") + ylab("log10(sd)")
  ggsave(plot=p2,file=paste0("/home/adrianf/project-folder/nobackup/ADRIAN/metagenomics-power-analyses-tutorial/all-bears/figures/p",j,".",i,"-otu-sd.png"),dpi=300)
  
  ## ---- simPower---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  sp<-simPower(group_size_vector=c(length(ss),length(ss)), otu_number=otunum, rare_depth=subsample,sequence_depth=1000,effect_range=seq(0,0.3,length.out=100))
  wj<-mclapply(sp,function(x) calcWJstudy(x),mc.cores=cores)
  # check if mean & sd is approximately OK
  mean(as.numeric(mclapply(wj, mean, mc.cores=cores)))
  mean(as.numeric(mclapply(wj, sd, mc.cores=cores)))
  
  ## ---- bootPower--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  bp80<-bootPower(wj, boot_number=10, subject_group_vector=c(length(ss),length(ss)),alpha=0.05)
  bp80_model <- subset(bp80, power < 0.95 & power > 0.2)
  bp80_model <- data.frame(log_omega2=log10(bp80_model$simulated_omega2),log_power=log10(bp80_model$power))
  bp80_model <- subset(bp80_model, log_omega2>-Inf)
  bp80_lm <- lm(log_omega2 ~ log_power, data=bp80_model)
  power80 <- 10^predict(bp80_lm, newdata=data.frame(log_power=log10(0.8)))
  power90 <- 10^predict(bp80_lm, newdata=data.frame(log_power = log10(0.9)))
  write.table(cbind(rbind(80,90),rbind(power80,power90)),file=paste0("/home/adrianf/project-folder/nobackup/ADRIAN/metagenomics-power-analyses-tutorial/all-bears/results/p",j,".",i,"-bootstrap-power.txt"),row.names = F,quote = F)
  p3<-ggplot2::qplot(data=bp80_model,x=log_power,y=log_omega2) + 
  geom_smooth(method="lm") + 
  ggtitle(paste0("log10(Omega2) to log10(Power) with",length(ss), "Subjects per Group")) + xlab("log10(Power)") + ylab("log10(Omega2)")
  ggsave(plot=p3,file=paste0("/home/adrianf/project-folder/nobackup/ADRIAN/metagenomics-power-analyses-tutorial/all-bears/figures/p",j,".",i,"-omega-power.png"),dpi=300)
}
}