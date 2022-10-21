#################
### R code for plotting Figure S5, same as figure 1 but comparing low and high gerp
### SFS for Russian and Founders (compare deleterious and synonymous)

# Loading R libraries
require(vcfR)
require(tidyverse)

m(list=ls())
plotdir<-"plots/"
filt="mac1"
anc=paste("Pol.2out.",filt, sep="")


# Reading in the data
vcf_file=paste("vcf/",anc,"/100S95F14R.chr1-38.vepfinal.wfm.vcf", sep="")
meta_file="help_files/metadata.txt"
vep_file=paste("vep/",anc,"/veptypes.chr1-X.txt", sep="")
sift_file=paste("vep/",anc,"/sifttypes.chr1-X.txt", sep="")
gerp_file=paste("gerp/canFam3/",anc,"/100S95F14R.chr1-38.vepfinal.rates.unique.bed", sep="")
vcf <- read.vcfR(vcf_file)
meta_tib<-meta_file %>% read.table(header=TRUE) %>% as_tibble() %>% rename(Indiv=UU_ID, Cat=Category)
vep_tib <- vep_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))
sift_tib <- sift_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))
gerp_tib <- gerp_file %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Gerp=V5, CHROM=V1, POS=V3) %>% mutate(CHROM = as.character(CHROM))

# Convert to tidy format
tidy_vcf <- vcfR2tidy(vcf,
  info_fields=c("AA"),
  format_fields=c("GT"),
  dot_is_NA=TRUE)

# Only need founders and russians
# (only sites where we have no missing data (3 ind = 6 alleles for founders
# and 28 alleles for Russians) AND at least 1 derived allele!)
founder_tib <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
	inner_join(tidy_vcf$fix) %>% filter(Indiv=="D-85-01" | Indiv=="FM1" | Indiv == "FM2") %>%
	inner_join(vep_tib) %>% left_join(sift_tib) %>%
	rename(Type=VEP_TYPE) %>% filter(Type=='synonymous' | SIFT_TYPE == 'deleterious') %>%
	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1', (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
  inner_join(gerp_tib) %>%
	mutate(Type=case_when((Type=='missense' & Gerp>4) ~ 'deleterious high', (Type=='missense' & Gerp<=4) ~ 'deleterious low', TRUE ~'synonymous')) %>%
	select(CHROM, POS, Indiv, new_gt, Type) %>%
	group_by(CHROM, POS, new_gt, Type) %>% summarize(count=n()) %>%
	mutate(ancestral=case_when((new_gt=='0/0') ~ (count*2.0), (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0), derived=case_when((new_gt=='1/1') ~ (count*2.0), (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0)) %>%
	group_by(CHROM, POS, Type) %>%  summarize(totanc=sum(ancestral), totder=sum(derived)) %>%
	filter(totanc+totder==6 && totder>0) %>%  group_by(Type, totder) %>% summarize(count=n()) %>%
  mutate(frac=case_when(Type=='deleterious high' ~ count/sum(count[Type=='deleterious high']),
                        Type=='deleterious low' ~ count/sum(count[Type=='deleterious low']),
                        Type=='synonymous' ~ count/sum(count[Type=='synonymous']))) %>%
  ungroup() %>% mutate(Cat="Founders")

russian_tib <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
  	inner_join(tidy_vcf$fix) %>% inner_join(meta_tib) %>%
  	filter(Cat=='Russia')  %>% inner_join(vep_tib) %>% inner_join(sift_tib) %>%
  	rename(Type=VEP_TYPE) %>% filter(Type=='synonymous' | SIFT_TYPE == 'deleterious') %>%
  	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1', (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
    inner_join(gerp_tib) %>%
  	mutate(Type=case_when((Type=='missense' & Gerp>4) ~ 'deleterious high', (Type=='missense' & Gerp<=4) ~ 'deleterious low', TRUE ~'synonymous')) %>%
	  select(CHROM, POS, Indiv, new_gt, Type) %>%
  	group_by(CHROM, POS, new_gt, Type) %>% summarize(count=n()) %>%
  	mutate(ancestral=case_when((new_gt=='0/0') ~ (count*2.0), (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0), derived=case_when((new_gt=='1/1') ~ (count*2.0), (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0)) %>%
  	group_by(CHROM, POS, Type) %>%  summarize(totanc=sum(ancestral), totder=sum(derived)) %>%
  	filter(totanc+totder==28 && totder>0) %>%  group_by(Type, totder) %>% summarize(count=n()) %>%
    mutate(frac=case_when(Type=='deleterious high' ~ count/sum(count[Type=='deleterious high']),
                          Type=='deleterious low' ~ count/sum(count[Type=='deleterious low']),
                          Type=='synonymous' ~ count/sum(count[Type=='synonymous']))) %>%
    ungroup() %>% mutate(Cat="Russia")

# Combine both and factorize
comb_tib<-bind_rows(founder_tib, russian_tib)
comb_tib$Cat<-factor(comb_tib$Cat, levels = c("Russia", "Founders"))
comb_tib$Type <- factor(comb_tib$Type, levels = c("synonymous", "deleterious low", "deleterious high"))



################################################################################
# PLOT (combine with facet)
# Don't plot the bars for "zero derived"

outfile=paste(plotdir,"FigureS5.pdf", sep="")

p<-ggplot(comb_tib, aes(x=totder, y=frac, fill=Type)) +
    geom_bar(position="dodge", stat="identity", alpha = 0.5) +
    scale_fill_manual(values=c("lightblue", "dodgerblue4", "cadetblue4")) +
    facet_wrap(~Cat, scales="free") +
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
