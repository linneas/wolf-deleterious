#################
### Allele frequency changes in Scandinavian wolves after 5-6 generations
# Only sites with data for at least 8 ind(>=16 alleles) are used, and no
# singletons (MAC=2)

rm(list=ls())

# Setting up, loading R libraries and set working directory
require(vcfR)
require(tidyverse)
plotdir="plots/"

filt="mac2"
anc=paste("Pol.2out.",filt, sep="")

# Reading in the data
vcf_file=paste("vcf/",anc,"/100S95F14R.chr1-38.vepfinal.wfm.vcf", sep="")
meta_file="help_files/metadata.txt"
vep_file=paste("vep/",anc,"/veptypes.chr1-X.txt", sep="")
sift_file=paste("vep/",anc,"/sifttypes.chr1-X.txt", sep="")
vcf <- read.vcfR(vcf_file)
meta_tib<-meta_file %>% read.table(header=TRUE) %>% as_tibble() %>% rename(Indiv=UU_ID, Cat=Category)
vep_tib <- vep_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))
sift_tib <- sift_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))

# Convert to tidy format
tidy_vcf <- vcfR2tidy(vcf,
                    info_fields=c("AA"),
                    format_fields=c("GT"),
                    dot_is_NA=TRUE)


# Tibble with allele counts for populations of interest
gt_tib <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
	inner_join(tidy_vcf$fix) %>% inner_join(meta_tib) %>%
	inner_join(vep_tib) %>% left_join(sift_tib) %>%
	rename(Type=VEP_TYPE) %>% filter(Type=='synonymous' | SIFT_TYPE == 'deleterious') %>%
	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1', (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
	mutate(Type=ifelse(Type=='missense', 'deleterious','synonymous')) %>%
	select(CHROM, POS, Indiv, new_gt, Type, GenClass, Cat) %>%
  filter(GenClass=="Finland" | GenClass=="F5" | GenClass=="F6" | GenClass=="Founders" |Cat=="2007-2014I") %>%
  mutate(Group=case_when((GenClass=="F5" | GenClass=="F6") ~ "LastScand",TRUE ~ Cat)) %>%
	group_by(CHROM, POS, new_gt, Type, Group) %>% summarize(count=n()) %>%
	mutate(ancestral=case_when((new_gt=='0/0') ~ (count*2.0), (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0), derived=case_when((new_gt=='1/1') ~ (count*2.0), (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0)) %>%
	group_by(CHROM, POS, Type, Group) %>% summarize(totanc=sum(ancestral), totder=sum(derived))


# Calculate derived allele frequencies
allfreq <- gt_tib %>%
  mutate_if(is.character, str_replace_all, pattern = "-", replacement = "to") %>%
  pivot_wider(names_from=Group, values_from=c(totanc,totder)) %>%
  filter(totanc_Founders+totder_Founders==6 & totder_Founders>0 &
          totanc_LastScand+totder_LastScand>=16) %>%
  mutate(Sc_derFreq=totder_LastScand/(totder_LastScand+totanc_LastScand),
        Fin_derFreq=totder_Finland/(totder_Finland+totanc_Finland),
        Fou_derFreq=totder_Founders/(totder_Founders+totanc_Founders),
        Id_derFreq=totder_2007to2014I/(totder_2007to2014I+totanc_2007to2014I)) %>%
  select(CHROM, POS, Type, Sc_derFreq, Id_derFreq, Fin_derFreq, Fou_derFreq, totder_Founders)
allfreq$totder_Founders<-as.factor(allfreq$totder_Founders)
allfreq$Type <- factor(allfreq$Type, levels = c("synonymous", "deleterious"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DOT PLOT WITH FOUNDER VS LAST SCAND TIME CLASS

# plot both next to eachother
p<-ggplot(allfreq, aes(x=totder_Founders, y=Sc_derFreq)) +
  geom_boxplot(aes(colour=Type, fill=Type), alpha=0.5) +
  scale_color_manual(values=c("lightblue", "dodgerblue4")) +
  scale_fill_manual(values=c("lightblue", "dodgerblue4")) +
  theme(panel.grid.major = element_line(colour = 'white'),
        panel.background = element_rect(fill = '#f5f4e6', colour = '#FFFDF8'),
        legend.position="bottom",
        legend.title=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key = element_rect(fill = '#f5f4e6', color = NA)) +
  labs(x="Derived alleles in original founders", y="Derived allele frequency after 5-6 generations")

outfile=paste(plotdir,"Figure2.revision.pdf", sep="")
show(outfile)

ggsave(
  outfile,
  plot = p,
  scale = 1,
  dpi = 300,
  limitsize = TRUE,
  width=4,
  height=4,
)

###########################################################
# Some numbers for paper
# Sites variable in founder
allfreq %>% filter(totder_Founders!=6)
# All sites variable in founder that has fixed in F5-F6
allfreq %>% filter(totder_Founders!=6 & (Sc_derFreq==1 |Sc_derFreq==0))
# Sites fixed for the least common variant in founders
allfreq %>% filter((totder_Founders==4 | totder_Founders==5) & (Sc_derFreq==0))
allfreq %>% filter((totder_Founders==1 | totder_Founders==2) & (Sc_derFreq==1))
# Startfreq of 0.5, now fixed for ancestral
allfreq %>% filter((totder_Founders==3) & (Sc_derFreq==0))

# Start freq of 0.5, no fixed for Derived
allfreq %>% filter((totder_Founders==3) & (Sc_derFreq==1))

# Sites fixed in F5-F6 that are not fixed in L1-L3
allfreq %>% filter(totder_Founders!=6 & (Sc_derFreq==1 |Sc_derFreq==0) &(Id_derFreq!=1 &Id_derFreq!=0) )
