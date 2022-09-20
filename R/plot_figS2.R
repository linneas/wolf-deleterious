#################
### Make SFS for Russians and Founders, only X chromosome!

#loading R libraries
require(vcfR)
require(tidyverse)

rm(list=ls())
plotdir<-"plots/"
filt="mac1"
anc=paste("Pol.2out.",filt, sep="")

# Reading in the data
Xvcf_file=paste("vcf/",anc,"/74Females.chrX.vepfinal.wfm.vcf", sep="")
meta_file="help_files/metadata.txt"
vep_file=paste("vep/",anc,"/veptypes.chr1-X.txt", sep="")
sift_file=paste("vep/",anc,"/sifttypes.chr1-X.txt", sep="")
Xvcf <- read.vcfR(Xvcf_file)
meta_tib<-meta_file %>% read.table(header=TRUE) %>% as_tibble() %>% rename(Indiv=UU_ID, Cat=Category)
vep_tib <- vep_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))
sift_tib <- sift_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))

# Convert to tidy format
tidy_Xvcf <- vcfR2tidy(Xvcf,
  info_fields=c("AA"),
  format_fields=c("GT"),
  dot_is_NA=TRUE)

# Get numbers of SNPs per group and type
Xgt_der <- tidy_Xvcf$gt %>% filter(!is.na(gt_GT)) %>%
              inner_join(tidy_Xvcf$fix)  %>% inner_join(vep_tib) %>% left_join(sift_tib) %>%
              mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
                                    (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
              filter(new_gt!="0/0") %>% select(CHROM, POS, Indiv, VEP_TYPE, SIFT_TYPE,new_gt) %>%
              inner_join(meta_tib)
# Table with counts per sift type and group
Xgt_der %>% select(CHROM, POS, Cat, VEP_TYPE, SIFT_TYPE) %>%
                mutate(Grp=ifelse((Cat=="1983-1990" | Cat=="1991-1998" | Cat=="1999-2006" | Cat=="2007-2014S"), "Original", Cat)) %>%
                group_by(CHROM, POS,VEP_TYPE,SIFT_TYPE,Grp) %>% summarize(indCount=n()) %>%
                group_by(VEP_TYPE, SIFT_TYPE, Grp) %>% summarize(count=n()) %>% print(n=50)



# Use the Russian females, only using sites where we have no missing data
#(8 females = 16 alleles)), and only look at non-PAR (Pos >7Mb)
russianX_tib <- tidy_Xvcf$gt %>% filter(!is.na(gt_GT)) %>%
	inner_join(tidy_Xvcf$fix) %>% inner_join(meta_tib) %>%
  filter(CHROM=="X" & POS>7000000) %>%
	filter(Cat=='Russia') %>%
  inner_join(vep_tib) %>% left_join(sift_tib) %>%
	rename(Type=VEP_TYPE) %>% filter(Type=='synonymous' | SIFT_TYPE == 'deleterious') %>%
	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1', (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
	mutate(Type=ifelse(Type=='missense', 'deleterious','synonymous')) %>%
	select(CHROM, POS, Indiv, new_gt, Type) %>%
	group_by(CHROM, POS, new_gt, Type) %>% summarize(count=n()) %>%
  mutate(ancestral=case_when((new_gt=='0/0') ~ (count*2.0),
                              (new_gt=='0/1') ~ (count*1.0),TRUE ~ 0),
        derived=case_when((new_gt=='1/1') ~ (count*2.0),
                              (new_gt=='0/1') ~ (count*1.0),TRUE ~ 0)) %>%
	group_by(CHROM, POS, Type) %>%  summarize(totanc=sum(ancestral), totder=sum(derived)) %>%
	filter(totanc+totder==16 & totder>0) %>%
  group_by(Type, totder) %>% summarize(count=n()) %>%
  mutate(frac=ifelse(Type=='deleterious',
    count/sum(count[Type=='deleterious']),
    count/sum(count[Type=='synonymous']))) %>% ungroup() %>%
    mutate(Cat="Russia")

# Add missing lines
for(i in 1:16) {
  v<-(russianX_tib$totder[russianX_tib$Type=="deleterious"]==i)
  if(any(v)==FALSE) {
    show(paste("there is no value for", i))
    russianX_tib <- russianX_tib %>% add_row(Type="deleterious", totder=i, count=0, frac=0.0, Cat="Russia", .before=i)
  }
}

# Also check founders (ONLY non par!!)
founderX_tib <- tidy_Xvcf$gt %>% filter(!is.na(gt_GT)) %>%
  	inner_join(tidy_Xvcf$fix) %>% filter(Indiv=="D-85-01" | Indiv=="FM1" | Indiv == "FM2") %>%
    filter(CHROM=="X" & POS>7000000) %>%
  	inner_join(vep_tib) %>% left_join(sift_tib) %>%
  	rename(Type=VEP_TYPE) %>% filter(Type=='synonymous' | SIFT_TYPE == 'deleterious') %>%
  	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1', (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
  	mutate(Type=ifelse(Type=='missense', 'deleterious','synonymous')) %>%
  	select(CHROM, POS, Indiv, new_gt, Type) %>%
    mutate(new_gt=case_when((Indiv=="FM1" | Indiv == "FM2") & new_gt=='0/0' ~ '0',
                            (Indiv=="FM1" | Indiv == "FM2") & new_gt=='1/1' ~ '1', TRUE ~ new_gt)) %>%
  	group_by(CHROM, POS, new_gt, Type) %>% summarize(count=n()) %>%
  	mutate(ancestral=case_when((new_gt=='0/0') ~ (count*2.0), (new_gt=='0/1') ~ (count*1.0),
                              new_gt=='0' ~ count*1, TRUE ~ 0),
          derived=case_when((new_gt=='1/1') ~ (count*2.0), (new_gt=='0/1') ~ (count*1.0),
                              new_gt=='1' ~ count*1, TRUE ~ 0)) %>%
  	group_by(CHROM, POS, Type) %>%  summarize(totanc=sum(ancestral), totder=sum(derived)) %>%
  	filter(totanc+totder==4 && totder>0) %>%  group_by(Type, totder) %>% summarize(count=n()) %>%
    mutate(frac=ifelse(Type=='deleterious',
                      count/sum(count[Type=='deleterious']),
                      count/sum(count[Type=='synonymous']))) %>%
    ungroup() %>% mutate(Cat="Founders")

    # Add missing rows
for(i in 1:4) {
  v<-(founderX_tib$totder[founderX_tib$Type=="deleterious"]==i)
  if(any(v)==FALSE) {
    show(paste("there is no value for", i))
    founderX_tib <- founderX_tib %>% add_row(Type="deleterious", totder=i, count=0, frac=0.0, Cat="Founders", .before=i)
  }
}

combX_tib<-bind_rows(russianX_tib, founderX_tib)
combX_tib$Cat<-factor(combX_tib$Cat, levels = c("Russia", "Founders"))
combX_tib$Type <- factor(combX_tib$Type, levels = c("synonymous", "deleterious"))


################################################################################
# COMBINE IN SINGLE PLOT USING FACET

outfile=paste(plotdir,"FigureS2.pdf", sep="")

p<-ggplot(combX_tib, aes(x=totder, y=frac, fill=Type)) +
    geom_bar(position="dodge", stat="identity", alpha = 0.5) +
    scale_fill_manual(values=c("lightblue", "dodgerblue4")) +
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
