#################
### Calculating VEP stats for tables and text


# Setting up, loading R libraries and set working directory
require(vcfR)
require(tidyverse)

# NOTE! For this first part, we use a file WITHOUT the imputed male founders
# NOTE! Also, the female founder is counted with "Original Scandinavia!"

rm(list=ls())
# Set filtering here
filt<-"mac2"
anc<-paste("Pol.2out.",filt, sep="")

# Reading in the data
vcf_file=paste("vcf/",anc,"/100S95F14R.chr1-38.vepfinal.vcf", sep="")
vcf_file2=paste("vcf/",anc,"/100S95F14R.chr1-38.modifier.wfm.vcf", sep="")
Xvcf_file=paste("vcf/",anc,"/74Females.chrX.vepfinal.vcf", sep="")
meta_file="help_files/metadata.txt"
vep_file=paste("vep/",anc,"/veptypes.chr1-X.txt", sep="")
sift_file=paste("vep/",anc,"/sifttypes.chr1-X.txt", sep="")

vcf <- read.vcfR(vcf_file)
Xvcf <- read.vcfR(Xvcf_file)
vcf2 <- read.vcfR(vcf_file2)
meta_tib<-meta_file %>% read.table(header=TRUE) %>% as_tibble() %>% rename(Indiv=UU_ID, Cat=Category)
vep_tib <- vep_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))
sift_tib <- sift_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))


# Convert to tidy format
tidy_vcf <- vcfR2tidy(vcf,
  info_fields=c("AA"),
  format_fields=c("GT"),
  dot_is_NA=TRUE)
tidy_Xvcf <- vcfR2tidy(Xvcf,
	info_fields=c("AA"),
	format_fields=c("GT"),
	dot_is_NA=TRUE)
tidy_vcf2 <- vcfR2tidy(vcf2,
  info_fields=c("AA"),
  format_fields=c("GT"),
  dot_is_NA=TRUE)



# Extracting relevant data, numbers per individual
# VEP categories
gt_vep <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
              inner_join(tidy_vcf$fix)  %>% inner_join(vep_tib) %>%
              mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
                        (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
              group_by(Indiv, new_gt, VEP_TYPE) %>% summarize(count=n()) %>% ungroup() %>%
              mutate_if(is.character, str_replace_all, pattern = "/", replacement = "") %>%
              pivot_wider(names_from=new_gt, names_prefix="gt", values_from=count) %>%
              inner_join(meta_tib) %>% mutate(Grp=ifelse((Cat=="1983-1990" | Cat=="1991-1998" | Cat=="1999-2006" | Cat=="2007-2014S" | Indiv=="D-85-01"), "Original", Cat))

# SIFT categories
gt_sift <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
              inner_join(tidy_vcf$fix)  %>% inner_join(vep_tib) %>% inner_join(sift_tib) %>%
              mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
                        (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
              group_by(Indiv, new_gt, VEP_TYPE, SIFT_TYPE) %>% summarize(count=n()) %>% ungroup() %>%
              mutate_if(is.character, str_replace_all, pattern = "/", replacement = "") %>%
              pivot_wider(names_from=new_gt, names_prefix="gt", values_from=count) %>%
              inner_join(meta_tib) %>% mutate(Grp=ifelse((Cat=="1983-1990" | Cat=="1991-1998" | Cat=="1999-2006" | Cat=="2007-2014S" | Indiv=="D-85-01"), "Original", Cat))


# Numbers per population instead
# (to save space, only save lines with the derived variant (new_gt!="0/0"))
gt_der <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
              inner_join(tidy_vcf$fix)  %>% inner_join(vep_tib) %>% left_join(sift_tib) %>%
              mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
                                    (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
              filter(new_gt!="0/0") %>% select(CHROM, POS, Indiv, VEP_TYPE, SIFT_TYPE, new_gt) %>% inner_join(meta_tib)

# And modifier category
gt_mod_der <- tidy_vcf2$gt %>% filter(!is.na(gt_GT)) %>%
              inner_join(tidy_vcf2$fix)  %>% mutate(VEP_TYPE="modifier") %>%
              mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
                                    (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
              filter(new_gt!="0/0") %>% select(CHROM, POS, Indiv, VEP_TYPE, new_gt) %>% inner_join(meta_tib)

################################################################################
# FOR TABLE 1:
# Table with counts per vep type and group
gt_der %>% mutate(Grp=ifelse((Cat=="1983-1990" | Cat=="1991-1998" | Cat=="1999-2006" | Cat=="2007-2014S" | Indiv=="D-85-01"), "Original", Cat)) %>%
          select(CHROM, POS, Cat, VEP_TYPE, Grp) %>%
          group_by(CHROM, POS, VEP_TYPE, Grp) %>%
          summarize(indCount=n()) %>% group_by(VEP_TYPE, Grp) %>%
          summarize(count=n()) %>% print(n=50)

# Table with counts for modifier
gt_mod_der %>% mutate(Grp=ifelse((Cat=="1983-1990" | Cat=="1991-1998" | Cat=="1999-2006" | Cat=="2007-2014S" | Indiv=="D-85-01"), "Original", Cat)) %>%
          select(CHROM, POS, Cat, VEP_TYPE, Grp) %>%
          group_by(CHROM, POS, VEP_TYPE, Grp) %>%
          summarize(indCount=n()) %>% group_by(VEP_TYPE, Grp) %>%
          summarize(count=n()) %>% print(n=50)

# Table with counts per sift type and group
gt_der %>% mutate(Grp=ifelse((Cat=="1983-1990" | Cat=="1991-1998" | Cat=="1999-2006" | Cat=="2007-2014S" | Indiv=="D-85-01"), "Original", Cat)) %>%
          select(CHROM, POS, Cat, VEP_TYPE, SIFT_TYPE, Grp) %>%
          group_by(CHROM, POS,VEP_TYPE,SIFT_TYPE, Grp) %>% summarize(indCount=n()) %>% group_by(VEP_TYPE, SIFT_TYPE, Grp) %>% summarize(count=n()) %>% print(n=50)

# Table with total counts per group
gt_der %>% mutate(Grp=ifelse((Cat=="1983-1990" | Cat=="1991-1998" | Cat=="1999-2006" | Cat=="2007-2014S" | Indiv=="D-85-01"), "Original", Cat)) %>%
          select(CHROM, POS, Cat, Grp) %>% group_by(CHROM, POS, Grp) %>%
          summarize(indCount=n()) %>% group_by(Grp) %>%
          summarize(count=n()) %>% print(n=50)

# Per individual:
# VEP types
gt_vep %>% mutate(sum=gt00+gt01+gt11, hetfrq=gt01/sum, homderfrq=gt11/sum, deralfreq=(gt11*2+gt01)/(sum*2), der=gt01+gt11) %>% group_by(Grp, VEP_TYPE) %>% summarize(Mean=mean(der), n=n(), sd=sd(der), se=sd/sqrt(n), min=min(der), max=max(der))
# SIFT types
gt_sift %>% mutate(sum=gt00+gt01+gt11, hetfrq=gt01/sum, homderfrq=gt11/sum, deralfreq=(gt11*2+gt01)/(sum*2), der=gt01+gt11) %>% group_by(Grp, VEP_TYPE, SIFT_TYPE) %>% summarize(Mean=mean(der), n=n(), sd=sd(der), se=sd/sqrt(n), min=min(der), max=max(der))

# For text, downsample the number of Original wolves to 14 (same as the Russian pop)
big_table<- gt_der %>% mutate(Grp=ifelse((Cat=="1983-1990" | Cat=="1991-1998" |
                                        Cat=="1999-2006" | Cat=="2007-2014S" |
                                        Indiv=="D-85-01"), "Original", Cat)) %>%
          select(CHROM, POS, Indiv, VEP_TYPE, SIFT_TYPE, Grp)

# Table of original Scand
orig<- big_table %>% filter(Grp=="Original") %>% select(Indiv) %>%
                    group_by(Indiv) %>% summarize()

del_vector<-c()
for (i in 1:100) {
  subset<- orig %>% slice_sample(n=14)
  temp<-big_table %>% inner_join(subset) %>% group_by(CHROM, POS,VEP_TYPE,SIFT_TYPE, Grp) %>%
            summarize(indCount=n()) %>% group_by(VEP_TYPE, SIFT_TYPE, Grp) %>%
            summarize(count=n())
            del_vector<-c(del_vector, temp$count[1])
}


################################################################################
# FOR TABLE S1:

# COUNT NUMBER OF VARIANTS SEEN PER GENERATION
gt_der %>% filter(VEP_TYPE=="synonymous" | SIFT_TYPE=="deleterious") %>%
          mutate(Type=case_when(VEP_TYPE=="synonymous" ~ VEP_TYPE,
                                            TRUE ~ SIFT_TYPE)) %>%
          select(CHROM, POS, Type, Indiv, GenClass) %>%
          group_by(CHROM, POS, Type, GenClass)%>% summarize(I_cnt=n()) %>%
          group_by(GenClass, Type) %>%
          summarize(count=n()) %>% print(n=30)

# NUMBERS PER INDIVIDUAL
# SYNONYMOUS
syn_all_ind<- gt_der %>% filter(VEP_TYPE=="synonymous") %>%
          mutate(GenClass=ifelse(GenClass=="F5" | GenClass=="F6", "F5-6", GenClass)) %>%
          select(CHROM, POS, Indiv, GenClass) %>%
          group_by(Indiv, GenClass)%>% summarize(count=n())
#DELETERIOUS
del_all_ind<- gt_der %>% filter(SIFT_TYPE=="deleterious") %>%
          mutate(GenClass=ifelse(GenClass=="F5" | GenClass=="F6", "F5-6", GenClass)) %>%
          select(CHROM, POS, Indiv, GenClass, SIFT_TYPE) %>%
          group_by(Indiv, GenClass, SIFT_TYPE)%>% summarize(count=n())


# To get proportions, also ancestral 0/0 are needed - make a gt table per generation
# TABLE PER GENERATION
gt_gen <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
        inner_join(tidy_vcf$fix)  %>% inner_join(vep_tib) %>%
        inner_join(sift_tib) %>% inner_join(meta_tib) %>%
        filter(VEP_TYPE=="synonymous" | SIFT_TYPE=="deleterious") %>%
        mutate(Type=case_when(VEP_TYPE=="synonymous" ~ VEP_TYPE,
                                                TRUE ~ SIFT_TYPE)) %>%
        mutate(GenClass=ifelse(GenClass=="F5" | GenClass=="F6", "F5-6", GenClass)) %>%
        mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
                                    (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
        select(CHROM, POS, Indiv, Type, GenClass, new_gt) %>%
        group_by(Indiv, new_gt, Type, GenClass) %>% summarize(count=n()) %>% ungroup() %>%
        mutate_if(is.character, str_replace_all, pattern = "/", replacement = "") %>%
        pivot_wider(names_from=new_gt, names_prefix="gt", values_from=count)

# GET PROP DELETERIOUS DERIVED HOMOZYGOTES, MEAN PER GEN CLASS
gt_gen %>% mutate(sum=gt00+gt01+gt11, homderfrq=gt11/sum) %>%
        group_by(GenClass, Type) %>%
        summarize(Mean=mean(homderfrq), n=n(), sd=sd(homderfrq), se=sd/sqrt(n), min=min(homderfrq), max=max(homderfrq))

# GET (NORMALIZED) NUMBER OF DELETERIOUS DERIVED HOMOZYGOTES, MEAN PER GEN CLASS
gt_gen %>% mutate(sum=gt00+gt01+gt11, callability=sum/max(sum), norm11=gt11/callability) %>%
      select(-Indiv, -gt00, -gt01) %>% group_by(Type, GenClass) %>%
      summarize(mean=mean(norm11), sd=sd(norm11), min=min(norm11), max=max(norm11), n=n())
# Number for F1: 176, and for F5-6: 255. Increase over these generations: 45%!

################################################################################
# FOR SUPPLEMENTARY TABLE 2:

# COUNT NUMBER OF VARIANTS SEEN PER FRoH class
gt_der %>% filter(VEP_TYPE=="synonymous" | SIFT_TYPE=="deleterious") %>%
                mutate(Type=case_when(VEP_TYPE=="synonymous" ~ VEP_TYPE,
                                                        TRUE ~ SIFT_TYPE)) %>%
                select(CHROM, POS, Indiv, Type, RoHClass) %>%
                group_by(CHROM, POS, Type, RoHClass)%>% summarize(I_cnt=n()) %>%
                group_by(RoHClass, Type) %>%
                summarize(count=n())

# PROP DD
# First need a new table with 0/0, just as above
gt_froh <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
            inner_join(tidy_vcf$fix)  %>% inner_join(vep_tib) %>%
            inner_join(sift_tib) %>% inner_join(meta_tib) %>%
            filter(VEP_TYPE=="synonymous" | SIFT_TYPE=="deleterious") %>%
            mutate(Type=case_when(VEP_TYPE=="synonymous" ~ VEP_TYPE,
                                                TRUE ~ SIFT_TYPE)) %>%
            mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
                                    (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
            select(CHROM, POS, Indiv, Type, RoHClass, new_gt) %>%
            group_by(Indiv, new_gt, Type, RoHClass) %>% summarize(count=n()) %>% ungroup() %>%
            mutate_if(is.character, str_replace_all, pattern = "/", replacement = "") %>%
            pivot_wider(names_from=new_gt, names_prefix="gt", values_from=count)

# Get Proportion of DD
gt_froh %>% mutate(sum=gt00+gt01+gt11, homderfrq=gt11/sum) %>%
                group_by(RoHClass, Type) %>% summarize(Mean=mean(homderfrq), n=n(), sd=sd(homderfrq), se=sd/sqrt(n), min=min(homderfrq), max=max(homderfrq))


################################################################################
# FOR NUMBERS OF VARIANTS IN TEXT, COMPARING 2007-2014I and 2007-2014S

gt_der %>% select(CHROM, POS, Cat, SIFT_TYPE) %>%
                filter(SIFT_TYPE=="deleterious") %>%
                group_by(CHROM, POS, SIFT_TYPE,Cat) %>%
                summarize(indCount=n()) %>% group_by(SIFT_TYPE, Cat) %>%
                summarize(count=n())

# Or maybe all functional variants?
gt_der %>% select(CHROM, POS, Cat) %>%
                group_by(CHROM, POS, Cat) %>%
                summarize(indCount=n()) %>% group_by(Cat) %>%
                summarize(count=n())


# Check how many of the deleterious among 2007-2014I are new
# (not seen in original scandinavia)
gt_der %>% filter(SIFT_TYPE=="deleterious") %>%
            mutate(Grp=ifelse((Cat=="1983-1990" | Cat=="1991-1998" |
                              Cat=="1999-2006" | Cat=="2007-2014S" |
                              Indiv=="D-85-01"), "Original", Cat)) %>%
          select(CHROM, POS, Grp) %>% group_by(CHROM, POS, Grp) %>%
          summarize(n=n()) %>%
          pivot_wider(names_prefix= "n_", names_from=Grp, values_from=n) %>%
          select(CHROM, POS, `n_2007-2014I`, n_Original) %>%
          filter(!is.na(`n_2007-2014I`) & is.na(n_Original))

# Deleterious in reproducing immigrants, not seen in original Scandinavia
gt_der %>% filter(SIFT_TYPE=="deleterious") %>%
            mutate(Grp=ifelse((Cat=="1983-1990" | Cat=="1991-1998" |
                              Cat=="1999-2006" | Cat=="2007-2014S" |
                              Indiv=="D-85-01"), "Original", Cat)) %>%
          select(CHROM, POS, Grp) %>% group_by(CHROM, POS, Grp) %>%
          summarize(n=n()) %>%
          pivot_wider(names_prefix= "n_", names_from=Grp, values_from=n) %>%
          select(CHROM, POS, n_R_immigrants, n_Original) %>%
          filter(!is.na(n_R_immigrants) & is.na(n_Original))

# Deleterious in L3, not previously seen in Scandinavia
gt_der %>% filter(SIFT_TYPE=="deleterious") %>%
            mutate(Grp=ifelse((Cat=="1983-1990" | Cat=="1991-1998" |
                              Cat=="1999-2006" | Cat=="2007-2014S" |
                              Indiv=="D-85-01"), "Original", GenClass)) %>%
          select(CHROM, POS, Grp) %>% group_by(CHROM, POS, Grp) %>%
          summarize(n=n()) %>%
          pivot_wider(names_prefix= "n_", names_from=Grp, values_from=n) %>%
          select(CHROM, POS, n_L3, n_Original) %>%
          filter(!is.na(n_L3) & is.na(n_Original)) %>% print(n=50)

# Homozygous in L3, not present at all before
gt_der %>% filter(!(GenClass=="L3" & new_gt=="0/1")) %>% filter(SIFT_TYPE=="deleterious") %>%
            mutate(Grp=ifelse((Cat=="1983-1990" | Cat=="1991-1998" |
                              Cat=="1999-2006" | Cat=="2007-2014S" |
                              Indiv=="D-85-01"), "Original", GenClass)) %>%
          select(CHROM, POS, Grp) %>% group_by(CHROM, POS, Grp) %>%
          summarize(n=n()) %>%
          pivot_wider(names_prefix= "n_", names_from=Grp, values_from=n) %>%
          select(CHROM, POS, n_L3, n_Original) %>%
          filter(!is.na(n_L3) & is.na(n_Original)) %>% print(n=50)


 ###############################################################################
 # LOOKING AT DRIFT! For the section "The effect of drift"
 # Checking numbers of fixed etc among the sites where all founders have data.

# NOTE! Now using the file WITH male founders, and first female is counted as
# a founder, not Original Scandinavian!
 filt<-"mac2"
 anc<-paste("Pol.2out.",filt, sep="")

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

# A tibble with allele counts for Finland, Immigrant offspring and two last
# generations of original Scandinavia
 gt_tib <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
 	inner_join(tidy_vcf$fix) %>%
 	inner_join(vep_tib) %>% left_join(sift_tib) %>%
 	rename(Type=VEP_TYPE) %>% filter(Type=='synonymous' | SIFT_TYPE == 'deleterious') %>%
 	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1', (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
 	mutate(Type=ifelse(Type=='missense', 'deleterious','synonymous')) %>%
 	select(CHROM, POS, Indiv, new_gt, Type) %>%
  inner_join(meta_tib) %>%
  filter(Cat=="Finland" | GenClass=="F5" | GenClass=="F6" | Cat=="Founder" |Cat=="2007-2014I") %>%
   mutate(Group=case_when((GenClass=="F5" | GenClass=="F6") ~ "LastScand",TRUE ~ Cat)) %>%
 	group_by(CHROM, POS, new_gt, Type, Group) %>% summarize(count=n()) %>%
 	mutate(ancestral=case_when((new_gt=='0/0') ~ (count*2.0), (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0), derived=case_when((new_gt=='1/1') ~ (count*2.0), (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0)) %>%
 	group_by(CHROM, POS, Type, Group) %>% summarize(totanc=sum(ancestral), totder=sum(derived))

# First, look at numbers only for the sites where we have all founders
counts_allF <- gt_tib %>%
    mutate_if(is.character, str_replace_all, pattern = "-", replacement = "to") %>%
    pivot_wider(names_from=Group, values_from=c(totanc,totder)) %>%
    filter(totanc_Founder+totder_Founder==6 & totder_Founder>0)

# How many sites are fixed for deleterious derived in the founders:
counts_allF %>% filter(Type=="deleterious" & totder_Founder==6)
# How many were polymorphic in the founders?
counts_allF %>% filter(Type=="synonymous" & totder_Founder!=6)
counts_allF %>% filter(Type=="deleterious" & totder_Founder!=6)

# For the polymorphic, we don't need to demand that we have data for all founders
gt_tib %>%
    mutate_if(is.character, str_replace_all, pattern = "-", replacement = "to") %>%
    pivot_wider(names_from=Group, values_from=c(totanc,totder)) %>%
    filter(Type=="deleterious" & totder_Founder>0 & totanc_Founder>0)
gt_tib %>%
    mutate_if(is.character, str_replace_all, pattern = "-", replacement = "to") %>%
    pivot_wider(names_from=Group, values_from=c(totanc,totder)) %>%
    filter(Type=="synonymous" & totder_Founder>0 & totanc_Founder>0)

gt_tib %>%
    mutate_if(is.character, str_replace_all, pattern = "-", replacement = "to") %>%
    pivot_wider(names_from=Group, values_from=c(totanc,totder)) %>%
    filter(Type=="deleterious" & totder_LastScand>0 & totanc_LastScand>0)

# And how many of those were fixed in the last Scand gen?
# Synonymous (first check fixed for derived, then fixed for ancestral):
counts_allF %>% filter(Type=="synonymous" & totder_Founder!=6 & totanc_LastScand==0)
counts_allF %>% filter(Type=="synonymous" & totder_Founder!=6 & totder_LastScand==0)
# Deleterious:
counts_allF %>% filter(Type=="deleterious" & totder_Founder!=6 & totanc_LastScand==0)
counts_allF %>% filter(Type=="deleterious" & totder_Founder!=6 & totder_LastScand==0)

# Check how many deleterious derived that were regained with the immigrants:
counts_allF %>% filter(Type=="deleterious" & totder_Founder!=6 & totanc_LastScand==0 & totanc_2007to2014I!=0)
