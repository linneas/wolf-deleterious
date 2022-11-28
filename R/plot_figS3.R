#################
### Plot Figure S3
### PLOTTING GERP SCORE DISTRIBUTIONS FOR ALL CLASSES OF VEP, AND GENOME WIDE


#Clear the workspace
rm(list=ls())

###############################################################################
# Make DENSITY plot of synonymous and deleterious in the same figure
# (on top of eachother)
require(tidyverse)

rm(list=ls())
plotdir="plots/"
filt="mac2"
anc=paste("Pol.2out.",filt, sep="")

#Input files
all_gerp_file=paste("gerp/canFam3/",anc,"/100S95F14R.chr1-38.allSNPs.rates.unique.bed", sep="")
vep_gerp_file=paste("gerp/canFam3/",anc,"/100S95F14R.chr1-38.vepfinal.rates.unique.bed", sep="")
mod_gerp_file=paste("gerp/canFam3/",anc,"/100S95F14R.chr1-38.modifier.rates.unique.bed", sep="")
vep_file=paste("vep/",anc,"/veptypes.chr1-X.txt", sep="")
sift_file=paste("vep/",anc,"/sifttypes.chr1-X.txt", sep="")

allgerp_tib <- all_gerp_file %>% read.table(header=FALSE) %>% as_tibble() %>%
                rename(CHROM=V1, POS=V3, GERP=V5) %>%
                mutate(CHROM = as.character(CHROM))
gerp_tib<-vep_gerp_file %>% read.table(header=FALSE) %>% as_tibble() %>% rename(CHROM=V1, POS=V3, GERP=V5) %>%
                mutate(CHROM = as.character(CHROM))
mod_tib<-mod_gerp_file %>% read.table(header=FALSE) %>% as_tibble() %>% rename(CHROM=V1, POS=V3, GERP=V5) %>%
                  mutate(CHROM = as.character(CHROM))
vep_tib <- vep_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))
sift_tib <- sift_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))

comb_tib<-gerp_tib %>% add_row(mod_tib)

tidy_tab <- comb_tib %>% inner_join(vep_tib) %>% left_join(sift_tib) %>%
      rename(Type=VEP_TYPE) %>%
      mutate(Type=case_when((SIFT_TYPE=="deleterious" | SIFT_TYPE=="tolerated") ~ SIFT_TYPE,
                      (is.na(Type) ~ 'modifier'),
                      TRUE ~ Type)) %>%
      mutate(Region="Genic") %>%
      select(CHROM, POS, GERP, Type, Region)

tidy_full <- allgerp_tib  %>% mutate(Type="all", Region="WholeGenome") %>%
            select(CHROM, POS, GERP, Type, Region)

comb_tib=bind_rows(tidy_tab, tidy_full)

#Add factors
comb_tib$Type <- factor(comb_tib$Type, levels = c("modifier", "synonymous","tolerated", "deleterious", "nonsense","all"))
comb_tib$Region <- factor(comb_tib$Region, levels = c("Genic", "WholeGenome"))

################################################################################
# PLOT
outfile=paste(plotdir, "FigureS4.revised.pdf", sep="")

p<-ggplot(data=comb_tib, aes(x=GERP, fill=Type)) +
    geom_density(alpha = 0.5, color=NA)+
    facet_wrap(~Region, scales="free") +
    scale_fill_manual(values=c("darkgrey", "lightblue", "seagreen4", "dodgerblue4","chocolate", "gray50")) +
    theme(panel.grid.major = element_line(colour = 'white'),
          panel.background = element_rect(fill = '#f5f4e6', colour = '#FFFDF8'),
          legend.box = "vertical",
          panel.spacing = unit(2, "lines"),
          legend.position="bottom",
          legend.title=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    labs(title="", x="GERP score", y="Density") +
    geom_vline(xintercept=4, colour="darkred",
             linetype="longdash", size=1)

ggsave(
  outfile,
  plot = p,
  scale = 1,
  dpi = 300,
  width=7.5,
  height=4,
  limitsize = TRUE)
