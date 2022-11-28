#################
### Plotting violin plots of different GERP loads

# Load libraries
require(tidyverse)

rm(list=ls())
plotdir="plots/"
filt="mac2"
anc=paste("Pol.2out.",filt, sep="")
gthres=4


# Check Gerp distr for thresholds
all_gerp=paste("gerp/canFam3/",anc,"/100S95F14R.chr1-38.allSNPs.rates.unique.bed", sep="")
all_gerp_tib <- all_gerp %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Gerp=V5)
# Check some different quantiles
all_gerp_tib %>% summarize(mean=mean(Gerp), max=max(Gerp), min=min(Gerp),
                qs=quantile(Gerp, c(0.1,0.25,0.5,0.75,0.9, 0.95, 0.99)),
                q=c(0.1,0.25,0.5,0.75,0.9, 0.95, 0.99))

# Input files
gerp_file=paste("gerp/canFam3/",anc,"/LoadPerInd.thres",gthres,".100S95F14R.chr1-38.allSNPs.wfm.txt", sep="")
meta_file="help_files/metadata.txt"

# Reading in the data
gerp_tib <- gerp_file %>% read.table(header=TRUE) %>% as_tibble() %>% rename(Indiv=IND, Realized=REALIZED_LOAD, Masked=MASKED_LOAD)
meta_tib<-meta_file %>% read.table(header=TRUE) %>% as_tibble() %>% rename(Indiv=UU_ID, Cat=Category)


# Adding all relevant columns
full_tib <- gerp_tib %>%
  mutate(SCORE=TOT_SUM/TOT_CALLS*mean(gerp_tib$TOT_CALLS)) %>%
  left_join(meta_tib) %>%
  mutate(Grp=case_when(GenClass=='Finland' | GenClass=='Russia' | GenClass=='NR_immigrants' ~ 'Reference', TRUE ~ 'Scandinavia'))  %>%
  mutate(GenClass=case_when(GenClass=="R_immigrants" ~ "Reproducing immigrants",
                            GenClass=="NR_immigrants" ~ "Non-reproducing immigrants",
                            TRUE ~ GenClass))

# Extracting only important columns and weight load
extr_tib <- full_tib %>% pivot_longer(cols=c(Realized, Masked),
                                      names_to="TYPE", values_to="LOAD") %>%
                      mutate(WLOAD=LOAD*mean(TOT_CALLS)) %>%
                      select(Indiv, GenClass, Grp, TYPE, LOAD, WLOAD)

# reorder boxes and make factor
extr_tib$TYPE <- factor(extr_tib$TYPE, levels = c("Masked", "Realized"))
extr_tib$GenClass <- factor(extr_tib$GenClass, levels = c("Finland", "Russia", "Non-reproducing immigrants", "Founders", "F1", "F2", "F3", "F4", "F5", "F6", "Reproducing immigrants", "L1", "L2", "L3"))


################################################################################
# MASKED (OR POTENTIAL) LOAD (assuming all deleterious are fully recessive)

#Output files
out=paste(plotdir,"Figure5.pdf", sep="")

# Plot
p<- ggplot(extr_tib, aes(x=GenClass, y=LOAD))+
  geom_violin(fill="white", alpha=0.5)+
  facet_grid(TYPE~factor(Grp, level=c("Reference","Scandinavia")), scales="free", space="free_x") +
  labs(x=NULL,y="GERP load per call")+
  theme(panel.grid.major = element_line(colour = 'white'),
        panel.background = element_rect(fill = '#f5f4e6', colour = '#FFFDF8'),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        legend.box = "vertical") +
  theme(legend.key=element_rect(fill='white', colour='NA')) +
  geom_point(colour="grey", size=1, alpha=0.5)

ggsave(
  out,
  plot = p,
  scale = 1,
  width = 7.5,
  height = 6,
  unit = "in",
  dpi = 300,
  limitsize = TRUE,
)
