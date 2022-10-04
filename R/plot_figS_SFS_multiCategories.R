#################
### R code for plotting Figure S X
### Make SFS for Russian and Founders (compare all classes of SNPs)

rm(list=ls())

# Setting up, loading R libraries and set working directory
require(vcfR)
require(tidyverse)
plotdir="plots/"

# For SFS, we use also low frequency alleles
filt="mac1"
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

# Only need founders and russians
# (only sites where we have no missing data (3 ind = 6 alleles for founders
# and 28 alleles for Russians) AND at least 1 derived allele!)
founder_tib <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
	inner_join(tidy_vcf$fix) %>% inner_join(meta_tib) %>% filter(Cat=="Founders") %>%
	inner_join(vep_tib) %>% left_join(sift_tib) %>%
	rename(Type=VEP_TYPE) %>% filter(Type=='synonymous' | SIFT_TYPE == 'deleterious') %>%
	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1', (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
	mutate(Type=ifelse(Type=='missense', 'deleterious','synonymous')) %>%
	select(CHROM, POS, Indiv, new_gt, Type) %>%
	group_by(CHROM, POS, new_gt, Type) %>% summarize(count=n()) %>%
	mutate(ancestral=case_when((new_gt=='0/0') ~ (count*2.0), (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0), derived=case_when((new_gt=='1/1') ~ (count*2.0), (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0)) %>%
	group_by(CHROM, POS, Type) %>%  summarize(totanc=sum(ancestral), totder=sum(derived)) %>%
	filter(totanc+totder==6 && totder>0) %>%  group_by(Type, totder) %>% summarize(count=n()) %>%
  mutate(frac=ifelse(Type=='deleterious',
                    count/sum(count[Type=='deleterious']),
                    count/sum(count[Type=='synonymous']))) %>%
  ungroup() %>% mutate(Pop="Founders")

russian_tib <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
  	inner_join(tidy_vcf$fix) %>% inner_join(meta_tib) %>%
  	filter(Cat=='Russia')  %>% inner_join(vep_tib) %>% left_join(sift_tib) %>%
  	rename(Type=VEP_TYPE) %>% filter(Type=='synonymous' | SIFT_TYPE == 'deleterious') %>%
  	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1', (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
  	mutate(Type=ifelse(Type=='missense', 'deleterious','synonymous')) %>%
  	select(CHROM, POS, Indiv, new_gt, Type) %>%
  	group_by(CHROM, POS, new_gt, Type) %>% summarize(count=n()) %>%
  	mutate(ancestral=case_when((new_gt=='0/0') ~ (count*2.0),
                                (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0),
          derived=case_when((new_gt=='1/1') ~ (count*2.0),
                            (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0)) %>%
  	group_by(CHROM, POS, Type) %>%
    summarize(totanc=sum(ancestral), totder=sum(derived)) %>%
  	filter(totanc+totder==28 && totder>0) %>%  group_by(Type, totder) %>%
    summarize(count=n()) %>%
    mutate(frac=ifelse(Type=='deleterious',
          count/sum(count[Type=='deleterious']),
          count/sum(count[Type=='synonymous']))) %>%
    ungroup() %>% mutate(Cat="Russia")

# Combine both and factorize
comb_tib<-bind_rows(founder_tib, russian_tib)
comb_tib$Cat<-factor(comb_tib$Cat, levels = c("Russia", "Founders"))
comb_tib$Type <- factor(comb_tib$Type, levels = c("synonymous", "deleterious"))



################################################################################
# PLOT (combine with facet)
# Don't plot the bars for "zero derived"

outfile=paste(plotdir,"Figure1.pdf", sep="")

p<-ggplot(comb_tib, aes(x=totder, y=frac, fill=Type)) +
    geom_bar(position="dodge", stat="identity", alpha = 0.5) +
    scale_fill_manual(values=c("lightblue", "dodgerblue4")) +
    facet_wrap(~Pop, scales="free") +
    labs(x="Derived alleles",y="Fraction of sites")+
    theme(panel.grid.major = element_line(colour = 'white'),
      panel.background = element_rect(fill = '#f5f4e6', colour = '#FFFDF8'),
      panel.spacing = unit(2, "lines"),
      legend.position="bottom",
      legend.title=element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_blank())

ggsave(outfile,
	plot = p,
	scale = 1,
	dpi = 300,
	limitsize = TRUE,
  width=7.5,
  height=4)
