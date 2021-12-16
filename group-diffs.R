library(tidyverse)
library(phyloseq)
library(microbiome)
library(microbiomeutilities)

otu<-readRDS("~/Downloads/amr-taxa-bear-0.01-af-bracken-max-highoral-0.01-otu-env-filtered.rds")
otu@sam_data$amr.era<-ifelse(otu@sam_data$Spec.HalfCentury %in% c("c20.2","c21.1"),"post","pre")
ord<-ordinate(transform(otu,"compositional"),method = "NMDS",distance = "bray")
plot_ordination(otu,ord,color = "amr.era")+stat_ellipse(type="norm",aes(color=amr.era))
pairwise.adonis(distance(transform(otu,transform = "compositional"),method = "jaccard"),factors = sample_data(otu)$amr.era) %>% knitr::kable() 

otu<-readRDS("calculus/DC2/all-samples/RDS/full-0.01-af-bracken-max-highoral-0.01-otu-env-contaminants.rds")
otu@sam_data$diet<-ifelse(otu@sam_data$Spec.host %in% c("Gorilla","Reindeer"),"herb",
                          ifelse(otu@sam_data$Spec.host=="Bear","omni",NA))

is.na(otu@otu_table)<-0
otu.sub<-prune_samples(rownames(otu@sam_data[!is.na(otu@sam_data$diet)]),otu) %>% 
  filter_taxa(.,function(x)mean(x)!=0,TRUE)
ord<-ordinate(otu.sub,method = "NMDS",distance = "bray")
plot_ordination(otu.sub,ord,shape = "diet",color="Spec.host")+stat_ellipse(aes(color=diet))

library(pairwiseAdonis)
pairwise.adonis(distance(transform(otu.sub,transform = "compositional"),method = "jaccard"),factors = sample_data(otu.sub)$diet) %>% knitr::kable() 
pairwise.adonis(distance(transform(otu.sub,transform = "compositional"),method = "jaccard"),factors = sample_data(otu.sub)$Spec.host) %>% knitr::kable()
