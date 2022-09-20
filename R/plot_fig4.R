#################
### Plotting GERP scores distribution for all and for deleterious and synonymous

# Load libraries
require(tidyverse)

rm(list=ls())
plotdir="plots/"
filt="mac2"
anc=paste("Pol.2out.",filt, sep="")

#Input files
all_gerp_file=paste("gerp/canFam3/",anc,"/100S95F14R.chr1-38.allSNPs.rates.unique.bed", sep="")
vep_gerp_file=paste("gerp/canFam3/",anc,"/100S95F14R.chr1-38.vepfinal.rates.unique.bed", sep="")
vep_file=paste("vep/",anc,"/veptypes.chr1-X.txt", sep="")
sift_file=paste("vep/",anc,"/sifttypes.chr1-X.txt", sep="")
# Read in data
all_gerp_tib <- all_gerp_file %>% read.table(header=FALSE) %>% as_tibble() %>%
                rename(CHROM=V1, POS=V3, GERP=V5) %>%
                mutate(CHROM = as.character(CHROM))
gerp_tib<-vep_gerp_file %>% read.table(header=FALSE) %>% as_tibble() %>% rename(CHROM=V1, POS=V3, GERP=V5) %>%
                mutate(CHROM = as.character(CHROM))
vep_tib <- vep_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))
sift_tib <- sift_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))

# Coding variants
coding_tib <- gerp_tib %>% inner_join(vep_tib) %>% left_join(sift_tib) %>%
      rename(Type=VEP_TYPE) %>% filter(Type=='synonymous' | SIFT_TYPE == 'deleterious') %>%
      mutate(Type=ifelse(Type=='missense', 'deleterious',Type)) %>%
      mutate(Region="Coding") %>%
      select(CHROM, POS, GERP, Type, Region)

# Whole genome variants
wg_tib <- all_gerp_tib  %>% mutate(Type="all", Region="WholeGenome") %>%
            select(CHROM, POS, GERP, Type, Region)

# Combine
comb_tib=bind_rows(coding_tib, wg_tib)

#Add factors
comb_tib$Type <- factor(comb_tib$Type, levels = c("synonymous","deleterious", "all"))
comb_tib$Region <- factor(comb_tib$Region, levels = c("Coding", "WholeGenome"))

################################################################################
# PLOT
# Output
outfile=paste(plotdir, "Figure4.pdf", sep="")

p<-ggplot(data=comb_tib, aes(x=GERP, fill=Type)) +
    geom_density(alpha = 0.5, color=NA)+
    facet_wrap(~Region, scales="free") +
    scale_fill_manual(values=c("lightblue", "dodgerblue4", "darkgrey")) +
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
