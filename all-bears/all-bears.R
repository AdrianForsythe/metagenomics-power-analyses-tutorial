## ---- setup------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
.libPaths("/home/adrianf/project-folder/nobackup/ADRIAN/R_library")
# devtools::install_github("brendankelly/micropower")
library(micropower)
library(tidyverse)
library(parallel)
library(phyloseq)
set.seed(515087345)

# The necessary input for micropower is within-group mean and standard deviations.
# These can be computed using the chosen distance metric from OTU tables.
# Thus any dataset that is used should either provide within-group mean and standard deviations already, 
# provide a within-group distance matrix, or provide data that can be transformed into an OTU table.

# 3 different OTU table simulations to determine different parameters
  # 1) simulate OTU tables with different levels of subsampling (rarefaction)
    # in order to find the subsampling level that corresponds to the within-group mean
  # 2) simulate OTU tables with different numbers of OTUs and the given rarefaction 
    # level in order to find the number of OTUs corresponding to the within-group standard deviation
  # 3) simulate OTU tables with different effect sizes and the given rarefaction level and number of OTUs.
    # This finally allows us to answer questions relating power, sample size, and effect size
      # how much power we will have given an effect size and sample size?
      # what is the minimum detectible effect size for a given power level and sample size?
      # how many samples we need to achieve at least x power and be able to detect an effect size of y?

# in this script, we are going to skip over the first two steps
# we have within group mean and standard deviation
# we need to figure out the level of subsampling correspoding to this sd

## ----load-otu----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
otu<-readRDS("RDS/amr-taxa-bear-0.01-af-bracken-max-highoral-0.01-otu-env-filtered.rds")

samples<-data.frame(sample_data(otu)) %>% 
  filter(Spec.host=="Bear"&!is.na(Spec.HalfCentury) & Spec.HalfCentury != "c21.1") %>%
  mutate(antib.era=ifelse(Spec.HalfCentury == "c20.2","post","pre"))

# how many individuals in each group?
nsamples<-samples %>% 
  group_by(antib.era) %>% 
  summarise(n=n_distinct(SampleID))
nsamples %>% knitr::kable()

pre<-samples %>% 
  filter(antib.era=="pre") %>% 
  distinct(SampleID) %>% 
  pull()

post<-samples %>% 
  filter(antib.era=="post") %>% 
  distinct(SampleID) %>% 
  pull()

load("all-bears/RDS/p0means.Rdata")
means %>% as.data.frame() %>% filter()

# make sure the otu table is clean
df<-otu %>% 
  otu_table() %>% as.data.frame() %>%
  select(all_of(pre),all_of(post)) %>% 
  filter(rowSums(.)!=0) %>% as.matrix()

# Weighted Jaccard Distances
# total mean and sd
m.t<-df %>% as.data.frame() %>% calcWJstudy(.) %>% mean() # 0.845
sd.t<-df %>% as.data.frame() %>% calcWJstudy(.) %>% sd() # 0.17

# 'pre' mean and sd
m.pre<-df %>% as.data.frame() %>% select(all_of(pre)) %>% calcWJstudy(.) %>% mean()
sd.pre<-df %>% as.data.frame() %>% select(all_of(pre)) %>% calcWJstudy(.) %>% sd()

# 'post' mean and sd
m.post<-df %>% as.data.frame() %>% select(all_of(post)) %>% calcWJstudy(.) %>% mean()
sd.post<-df %>% as.data.frame() %>% select(all_of(post)) %>% calcWJstudy(.) %>% sd()

load("currbio-bears/RDS/means.Rdata")
mean_df <- data.frame(subsampling=as.numeric(names(means)),means=as.numeric(means),target=abs(as.numeric(means)-m.t))
subsample <- mean_df[which(mean_df$target==min(mean_df$target)),]$subsampling
# qplot(data=mean_df,x=subsampling,y=means)

# construct a distance matrix on all samples
named_df<-calcWJstudy(df)

# make col and row names reflect group and individual identity
colnames(named_df)<-ifelse(colnames(named_df) %in% pre, 
                           paste0("g1s",which(colnames(named_df) %in% pre)),
                           paste0("g2s",which(colnames(named_df) %in% pre)))
rownames(named_df)<-colnames(named_df)

# we can estimate the effect size from this data using the calcOmega2 function
# i.e. permanova power!

# what size of an effect ca we detect a difference between groups?
# with no subsampling
calcOmega2(named_df)
calcPERMANOVAp(perm = PERMANOVA(named_df))
calcR2(perm = PERMANOVA(named_df))

obs.results<-NULL
for (s in c(5,10,15,20)) {
  # for (i in 1:100) {
  print(s)
  ndfb<-bootDM(named_df,subject_group_vector = c(s,s))
  om<-calcOmega2(ndfb)
  perm.p<-calcPERMANOVAp(perm = PERMANOVA(ndfb))
  obs.results<-rbind(obs.results,cbind(s,om,perm.p))
  # }
}

subsample=0.
results<-NULL
for (s in 5) {
  sp<-simPower(group_size_vector=c(s,s),
               otu_number=100,
               rare_depth=subsample,
               sequence_depth=1000,
               effect_range=seq(0,0.1,
                                length.out=100))
  wj<-mclapply(sp,function(x) calcWJstudy(x),mc.cores=cores)
  
  bp<-bootPower(wj, boot_number=100, subject_group_vector=c(s,s),alpha=0.05)
  bp_model <- subset(bp, power < 0.95 & power > 0.2)
  bp_model <- data.frame(log_omega2=log10(bp_model$simulated_omega2),log_power=log10(bp_model$power))
  bp_model <- subset(bp_model, log_omega2>-Inf)
  bp_lm <- lm(log_omega2 ~ log_power, data=bp_model)
  power80 <- 10^predict(bp_lm, newdata=data.frame(log_power=log10(0.8)))
  power90 <- 10^predict(bp_lm, newdata=data.frame(log_power = log10(0.9)))
  p<-ggplot2::qplot(data=bp_model,x=log_power,y=log_omega2) + geom_smooth(method="lm") + ggtitle(paste0("log10(Omega2) to log10(Power) with", s ,"Subjects/Group")) + xlab("log10(Power)") + ylab("log10(Omega2)")
  print(p)
  results<-rbind(results,cbind(s,bp,power80,power90))
}

results %>% distinct(s,power80,power90)

results %>% ggplot(aes(x=simulated_omega2,y=power,color=as.factor(s)))+
  geom_point()+
  geom_smooth()+
  coord_cartesian(ylim=c(0,1),xlim = c(0,0.08))
  # geom_vline(xintercept = 0.039)+theme_classic()
ggsave(plot = last_plot(),filename = "currbio-bears/sim-effect-sizes.png",dpi = 300)