## ---- setup------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(micropower)
library(tidyverse)
library(parallel)
library(phyloseq)
set.seed(515087345)

## ----load-otu----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# load in the AMR taxa from bears
amr.otu<-readRDS("RDS/amr-taxa-bear-0.01-af-bracken-max-highoral-0.01-otu-env-filtered.rds")
# to subset the full otu table
otu<-readRDS("~/Lab-Notes/calculus/DC2/all-samples/RDS/full-0.01-af-bracken-max-highoral-0.01-otu-env-filtered.rds") %>% prune_taxa(taxa = taxa_names(amr.otu),.)

# make a new sample category for diet
samples<-data.frame(sample_data(otu)) %>% 
  filter(Spec.host2 %in% c("Bear","Reindeer","Gorilla")) %>% 
  mutate(diet=ifelse(Spec.host2=="Bear","carnivore","herbivore")) %>% 
  select(SampleID,Spec.host2,diet)

# count how many individuals are in each group
nsamples<-samples %>% group_by(diet,Spec.host2) %>% summarise(n=n_distinct(SampleID))
nsamples %>% knitr::kable()

# make list for each diet/species
# herbivores
herb<-samples %>% 
  filter(diet=="herbivore") %>% 
  distinct(SampleID) %>% 
  pull()

# carinvores
carn<-samples %>% filter(diet=="carnivore") %>% 
  distinct(SampleID) %>% 
  pull()

# reindeer
reindeer<-samples %>% filter(Spec.host2=="Reindeer") %>% 
  distinct(SampleID) %>% 
  pull()

# bear
bear<-samples %>% filter(Spec.host2=="Bear") %>% 
  distinct(SampleID) %>% 
  pull()

# gorilla
gorilla<-samples %>% filter(Spec.host2=="Gorilla") %>% 
  distinct(SampleID) %>% 
  pull()

# make sure the otu table is clean
df<-otu %>% 
  otu_table() %>% as.data.frame() %>%
  select(all_of(herb),all_of(carn)) %>% 
  select(which(colSums(.)>0)) %>% 
  filter(rowSums(.)!=0) %>% as.matrix()

# calculate total mean and sd
df %>% as.data.frame() %>% calcWJstudy(.) %>% mean() # 0.899
df %>% as.data.frame() %>% calcWJstudy(.) %>% sd() # 0.124

# create a distance matrix on all samples
named_df<-calcWJstudy(df)

# add unique group and individual labels, by changing col and row names
colnames(named_df)<-ifelse(colnames(named_df) %in% reindeer, 
                           paste0("g1s",which(colnames(named_df) %in% reindeer)),
                           ifelse(colnames(named_df) %in% gorilla,paste0("g2s",which(colnames(named_df) %in% gorilla)),
                                  paste0("g3s",which(colnames(named_df) %in% bear))))
rownames(named_df)<-colnames(named_df)

# we can estimate the effect size from this data using the calcOmega2 function
# i.e. permanova power!

# create a dummy variable to collect results
obs.results<-NULL
# iterate through all combinations of species (and dietary groups)
for (s in list(c(length(gorilla),length(reindeer),length(bear)),
            c(length(gorilla),length(reindeer),0),
            c(length(gorilla),0,length(bear)),
            c(0,length(reindeer),length(bear)))) {
  print(s)
  ndfb<-bootDM(named_df,subject_group_vector = s)
  om<-calcOmega2(ndfb)
  perm.p<-calcPERMANOVAp(perm = PERMANOVA(ndfb))
  obs.results<-rbind(obs.results,cbind(s,om,perm.p))
  
}

results<-NULL
for (s in list(c(length(gorilla),length(reindeer),length(bear)),
               c(length(gorilla),length(reindeer),0),
               c(length(gorilla),0,length(bear)),
               c(0,length(reindeer),length(bear)))) {
  sp<-simPower(group_size_vector=c(s,s),
               otu_number=100,
               rare_depth=subsample,
               sequence_depth=1000,
               effect_range=seq(0,0.1,
                                length.out=100))
  wj<-mclapply(sp,function(x) calcWJstudy(x),mc.cores=cores)
  
  bp<-bootPower(wj, boot_number=10, subject_group_vector=c(s,s),alpha=0.05)
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

