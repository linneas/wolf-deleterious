#################
### Calculating GERP stats for tables and text


rm(list=ls())

# Setting up, loading R libraries and set working directory
require(vcfR)
require(tidyverse)

###############################################################################
# CALCULATING NUMBERS FOR GERP SECTION IN PAPER
# PLUS NUMBERS FOR TABLE S3
filt="mac2"
anc=paste("Pol.2out.", filt, sep="")
gthres=4

# Check Gerp distr for thresholds
all_gerp=paste("gerp/canFam3/",anc,"/100S95F14R.chr1-38.rates.unique.bed", sep="")
allgerp_tib <- all_gerp %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Gerp=V5)
# Check some different quantiles
allgerp_tib %>% summarize(mean=mean(Gerp), max=max(Gerp), min=min(Gerp),
                qs=quantile(Gerp, c(0.1,0.25,0.5,0.75,0.9, 0.95, 0.99)),
                q=c(0.1,0.25,0.5,0.75,0.9, 0.95, 0.99))

#The number of sites with GERP values are given from the length of the tibble
allgerp_tib

# To get the numbers of deleterious sites per individual, we need to read in
# the full table
# FILE NAMES
gerp_file=paste("gerp/canFam3/",anc,"/LoadPerInd.thres",gthres,".allSNPs.withMaleFounders.chr1-38.txt", sep="")
meta_file="help_files/metadata.txt"
vep_file=paste("vep/",anc,"/veptypes.chr1-X.txt", sep="")
sift_file=paste("vep/",anc,"/sifttypes.chr1-X.txt", sep="")

# READ IN FILES
gerp_tib <- gerp_file %>% read.table(header=TRUE) %>% as_tibble() %>% rename(Indiv=IND, Realized=REALIZED, Masked=MASKED)
meta_tib<-meta_file %>% read.table(header=TRUE) %>% as_tibble() %>% rename(Indiv=UU_ID, Cat=Category)
vep_tib <- vep_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))
sift_tib <- sift_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))

# Adding all relevant columns and remove the nonbreeding immigrants
data_tib <- gerp_tib %>%
  mutate(SCORE=TOT_SUM/TOT_CALLS*mean(gerp_tib$TOT_CALLS)) %>%
  left_join(meta_tib) %>%
  mutate(Grp=case_when(Cat=='Finland' | Cat=='Russia' ~ 'Source', TRUE ~ 'Scandinavia'))  %>%
  mutate(m_sites=TOT_DEL-HOM_DEL, r_sites=HOM_DEL/2)


# Extract only relevant and summarize per Population
extr_tib <- data_tib %>% select(Indiv, GenClass, m_sites, r_sites, TOT_CALLS) %>%
        group_by(GenClass) %>% pivot_longer(cols=c(m_sites, r_sites),
                                          names_to="TYPE", values_to="sites") %>%
        mutate(cov=TOT_CALLS/4995746, norm_sites=sites/cov)

# Mean and SD for table S3
res_masked<- extr_tib %>% group_by(GenClass, TYPE) %>%
      summarize(n=n(), avg=mean(norm_sites), sd=sd(norm_sites), se=sd/sqrt(n),
                min=min(norm_sites), max=max(norm_sites)) %>% filter(TYPE=="m_sites") %>%
                print(n=50)

res_realized<- extr_tib %>% group_by(GenClass, TYPE) %>%
      summarize(n=n(), avg=mean(norm_sites), sd=sd(norm_sites), se=sd/sqrt(n),
                min=min(norm_sites), max=max(norm_sites)) %>% filter(TYPE=="r_sites") %>%
                print(n=50)

## COMBINE VEP AND GERP, BUT USE COMB (NOT SICOMB AS FOR THE SFS IN SUPPL FIG)
gt_tib <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
	inner_join(tidy_vcf$fix) %>%
	inner_join(vep_tib) %>% left_join(sift_tib) %>%
	rename(Type=VEP_TYPE) %>% filter(Type=='synonymous' | SIFT_TYPE == 'deleterious') %>%
	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1', (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
  inner_join(allgerp_tib) %>%
	mutate(Type=case_when((Type=='missense' & Gerp>4) ~ 'deleterious_high', (Type=='missense' & Gerp<=4) ~ 'deleterious_low', TRUE ~'synonymous')) %>%
  group_by(Indiv, new_gt, Category) %>% summarize(count=n()) %>% ungroup() %>%
  mutate_if(is.character, str_replace_all, pattern = "/", replacement = "") %>%
   pivot_wider(names_from=new_gt, names_prefix="gt", values_from=count) %>%
  mutate(sum=gt00+gt01+gt11) %>%
  mutate(hetfrq=gt01/sum, homderfrq=gt11/sum, deralfreq=(gt11*2+gt01)/(sum*2))  %>%

# How many in total?
tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>% inner_join(tidy_vcf$fix) %>%
      inner_join(vep_tib) %>% left_join(sift_tib) %>%rename(Type=VEP_TYPE) %>%
      filter(Type=='synonymous' | SIFT_TYPE == 'deleterious') %>%
      mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
                              (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
      inner_join(gerp_tib) %>%
      mutate(Type=case_when((Type=='missense' & Gerp>4) ~ 'deleterious_high',
                                (Type=='missense' & Gerp<=4) ~ 'deleterious_low',
                                 TRUE ~'synonymous')) %>%
      group_by(Indiv, new_gt, Category) %>% filter(Type=="deleterious_high") %>%
       ungroup() %>% select(CHROM, POS, Indiv) %>% group_by(CHROM, POS) %>% summarize(n=n())


# And summarize highly deleterious per group
gt_tib %>% select(Type, gt11, sum, GenClass) %>%
						filter(Category=="deleterious_high") %>% mutate(fracCov=sum/2002) %>%
						mutate(newGt11=gt11/fracCov) %>% group_by(GenClass) %>%
						summarize(n=n(), mean=mean(newGt11), sd=sd(newGt11), max=max(newGt11), min=min(newGt11))
