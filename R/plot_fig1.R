#################
### R code for plotting Figure 1
### Make SFS for Russian and Founders (all VEP categories incl modifier)

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
vcf_file2=paste("vcf/",anc,"/100S95F14R.chr1-38.modifier.wfm.vcf", sep="")
meta_file="help_files/metadata.txt"
vep_file=paste("vep/",anc,"/veptypes.chr1-X.txt", sep="")
sift_file=paste("vep/",anc,"/sifttypes.chr1-X.txt", sep="")
vcf <- read.vcfR(vcf_file)
vcf2 <- read.vcfR(vcf_file2)
meta_tib<-meta_file %>% read.table(header=TRUE) %>% as_tibble() %>% rename(Indiv=UU_ID, Cat=Category)
vep_tib <- vep_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))
sift_tib <- sift_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))

# Convert to tidy format
tidy_vcf <- vcfR2tidy(vcf,
  info_fields=c("AA"),
  format_fields=c("GT"),
  dot_is_NA=TRUE)
tidy_vcf2 <- vcfR2tidy(vcf2,
  info_fields=c("AA"),
  format_fields=c("GT"),
  dot_is_NA=TRUE)

# Combine coding and modifier
  combined_gt<-bind_rows(tidy_vcf$gt, tidy_vcf2$gt)
  combined_fix<-bind_rows(tidy_vcf$fix,tidy_vcf2$fix)

# Only need founders and russians
# (only sites where we have no missing data (3 ind = 6 alleles for founders
# and 28 alleles for Russians) AND at least 1 derived allele!)
founder_tib <- combined_gt %>% filter(!is.na(gt_GT)) %>%
	inner_join(combined_fix) %>% inner_join(meta_tib) %>%
	filter(Cat=="Founders") %>% left_join(vep_tib) %>% left_join(sift_tib) %>%
	rename(Type=VEP_TYPE)  %>%
	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
                          (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
	mutate(Type=case_when((SIFT_TYPE=="deleterious" | SIFT_TYPE=="tolerated") ~ SIFT_TYPE,
                        is.na(Type) ~ 'modifier', TRUE ~ Type)) %>%
	select(CHROM, POS, Indiv, new_gt, Type) %>%
	group_by(CHROM, POS, new_gt, Type) %>% summarize(count=n()) %>%
	mutate(ancestral=case_when((new_gt=='0/0') ~ (count*2.0), (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0), derived=case_when((new_gt=='1/1') ~ (count*2.0), (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0)) %>%
	group_by(CHROM, POS, Type) %>%  summarize(totanc=sum(ancestral), totder=sum(derived)) %>%
	filter(totanc+totder==6 && totder>0 && Type!='missense') %>%  group_by(Type, totder) %>% summarize(count=n()) %>%
  mutate(frac=case_when((Type=='deleterious') ~ (count/sum(count[Type=='deleterious'])),
                        (Type=='synonymous') ~ (count/sum(count[Type=='synonymous'])),
                        (Type=='tolerated') ~ (count/sum(count[Type=='tolerated'])),
                        (Type=='nonsense') ~ (count/sum(count[Type=='nonsense'])),
                        (Type=='modifier') ~ (count/sum(count[Type=='modifier'])),
                        TRUE ~ 0)) %>%
  ungroup() %>% mutate(Cat="Founders")


russian_tib <- combined_gt %>% filter(!is.na(gt_GT)) %>%
	inner_join(combined_fix) %>% inner_join(meta_tib) %>%
  	filter(Cat=='Russia')  %>% left_join(vep_tib) %>% left_join(sift_tib) %>%
  	rename(Type=VEP_TYPE) %>%
  	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
                            (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
    mutate(Type=case_when((SIFT_TYPE=="deleterious" | SIFT_TYPE=="tolerated") ~ SIFT_TYPE,
                          is.na(Type) ~ 'modifier', TRUE ~ Type)) %>%
  	select(CHROM, POS, Indiv, new_gt, Type) %>%
  	group_by(CHROM, POS, new_gt, Type) %>% summarize(count=n()) %>%
  	mutate(ancestral=case_when((new_gt=='0/0') ~ (count*2.0),
                                (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0),
          derived=case_when((new_gt=='1/1') ~ (count*2.0),
                            (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0)) %>%
  	group_by(CHROM, POS, Type) %>%
    summarize(totanc=sum(ancestral), totder=sum(derived)) %>%
  	filter(totanc+totder==28 && totder>0 && Type!='missense') %>%  group_by(Type, totder) %>%
    summarize(count=n()) %>%
    mutate(frac=case_when((Type=='deleterious') ~ (count/sum(count[Type=='deleterious'])),
                          (Type=='synonymous') ~ (count/sum(count[Type=='synonymous'])),
                          (Type=='tolerated') ~ (count/sum(count[Type=='tolerated'])),
                          (Type=='nonsense') ~ (count/sum(count[Type=='nonsense'])),
                          (Type=='modifier') ~ (count/sum(count[Type=='modifier'])),
                          TRUE ~ 0)) %>%
    ungroup() %>% mutate(Cat="Russia")

# Combine both and factorize
comb_tib<-bind_rows(founder_tib, russian_tib)
comb_tib$Cat<-factor(comb_tib$Cat, levels = c("Russia", "Founders"))
comb_tib$Type <- factor(comb_tib$Type, levels = c("modifier", "synonymous", "tolerated", "deleterious", "nonsense"))



################################################################################
# PLOT (combine with facet)
# Don't plot the bars for "zero derived"

outfile=paste(plotdir,"Figure1.revision.pdf", sep="")

p<-ggplot(comb_tib, aes(x=totder, y=frac, fill=Type)) +
    geom_bar(position="dodge", stat="identity", alpha = 0.5) +
    scale_fill_manual(values=c("darkgray","lightblue", "seagreen4", "dodgerblue4", "chocolate")) +
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

################################################################################
# STATS - COMPARE DISTRIBUTIONS WITH CHI-SQUARE GOODNESS OF FIT TEST
# Compare distribution of deleterious and synonymous
founder_pivot <- founder_tib %>% select(Type, totder, count, frac) %>%
              pivot_wider(names_from=Type, values_from=c(count, frac)) %>%
              mutate(Cat="Founders", max=max(totder))
russian_pivot <- russian_tib %>% select(Type, totder, count, frac) %>%
              pivot_wider(names_from=Type, values_from=c(count, frac)) %>%
              mutate(Cat="Russia", max=max(totder))
comb_pivot<-bind_rows(founder_pivot, russian_pivot)
comb_pivot$Cat<-factor(comb_pivot$Cat, levels = c("Russia", "Founders"))


# Test if the distribution of non-synonymous are deviating significantly from
# the synonymous distribution:
chisq.test(founder_pivot$count_deleterious, p=founder_pivot$frac_synonymous)
#data:  founder_pivot$count_deleterious
#X-squared = 128.29, df = 5, p-value < 2.2e-16
# Yes!
chisq.test(russian_pivot$count_deleterious, p=russian_pivot$frac_synonymous)
#data:  russian_pivot$count_deleterious
#X-squared = 640.55, df = 27, p-value < 2.2e-16
# Yes!

# Modifier vs synonymous
chisq.test(founder_pivot$count_modifier, p=founder_pivot$frac_synonymous)
  # X-squared = 28.331, df = 5, p-value = 3.136e-05
chisq.test(russian_pivot$count_modifier, p=russian_pivot$frac_synonymous)
  #X-squared = 234.45, df = 27, p-value < 2.2e-16

#Nonsense vs synonymous
chisq.test(founder_pivot$count_nonsense, p=founder_pivot$frac_synonymous)
  # X-squared = 0.98077, df = 5, p-value = 0.9641
chisq.test(russian_pivot$count_nonsense, p=russian_pivot$frac_synonymous)
  # X-squared = 38.673, df = 27, p-value = 0.06782

# Tolerated vs synonymous
chisq.test(founder_pivot$count_tolerated, p=founder_pivot$frac_synonymous)
  #X-squared = 9.3002, df = 5, p-value = 0.09767
chisq.test(russian_pivot$count_tolerated, p=russian_pivot$frac_synonymous)
  #X-squared = 70.661, df = 27, p-value = 9.024e-06
