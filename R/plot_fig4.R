#################
### Plotting violin plots of heterozygous and homozygous derived
### alleles, plus numbers related to this figure in the text

# loading R libraries
require(vcfR)
require(tidyverse)

rm(list=ls())
plotdir<-"plots/"
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

# Extracting relevant data and combine tibbles
gt_tib <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
	inner_join(tidy_vcf$fix) %>% inner_join(vep_tib) %>% left_join(sift_tib) %>%
	filter(VEP_TYPE=='synonymous' | SIFT_TYPE == 'deleterious') %>%
	mutate(Type=ifelse(VEP_TYPE=='missense', 'deleterious','synonymous')) %>%
 	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
		(ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
	group_by(Indiv, new_gt, Type) %>% summarize(count=n()) %>% ungroup() %>%
	mutate_if(is.character, str_replace_all, pattern = "/", replacement = "") %>%
	 pivot_wider(names_from=new_gt, names_prefix="gt", values_from=count) %>%
	mutate(sum=gt00+gt01+gt11) %>%
	mutate(hetfrq=gt01/sum, homderfrq=gt11/sum, deralfreq=(gt11*2+gt01)/(sum*2)) %>%
	inner_join(meta_tib)

gt_tib$Type <- factor(gt_tib$Type, levels = c("synonymous", "deleterious"))


################################ PLOTTING ######################################
# Plot synonymous and deleterious together with generation class on X-axis
extr_gt<-gt_tib %>% select(Indiv,Type,Cat,hetfrq,homderfrq,deralfreq,GenClass) %>%
				mutate(Grp=case_when(GenClass=='Finland' | GenClass=='Russia' | GenClass=='NR_immigrants' ~ 'Reference',
				TRUE ~ 'Scandinavia')) %>%
				mutate(GenClass=case_when(GenClass=="R_immigrants" ~ "Reproducing immigrants",
																	GenClass=="NR_immigrants" ~ "Non-reproducing immigrants",
																	TRUE ~ GenClass))


# Make factor
extr_gt$GenClass <- factor(extr_gt$GenClass, levels = c("Finland", "Russia", "Non-reproducing immigrants", "Founders", "F1", "F2", "F3", "F4", "F5", "F6", "Reproducing immigrants", "L1", "L2", "L3"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# # # # # Heterozygous sites
output <- paste(plotdir,"Figure4A.pdf", sep="")

pA<- ggplot(extr_gt, aes(x=GenClass, y=hetfrq))+
  geom_violin(fill="white", alpha=0.5)+
  facet_grid(Type~factor(Grp, level=c("Reference","Scandinavia")), scales="free", space="free_x") +
  labs(x=NULL,y="proportion of heterozygous genotypes")+
  theme(panel.grid.major = element_line(colour = 'white'),
        panel.background = element_rect(fill = '#f5f4e6', colour = '#FFFDF8'),
        axis.text.x = element_text(angle = 45, hjust=1),
				strip.text.x = element_text(size=12),
				strip.text.y = element_text(size=12),
        legend.box = "vertical") +
  theme(legend.key=element_rect(fill='white', colour='NA')) +
	geom_point(data=extr_gt, colour="grey", size=1, alpha=0.5)

ggsave(
  output,
  plot = pA,
  scale = 1,
  width = 7.5,
  height = 6,
  unit = "in",
  dpi = 300,
  limitsize = TRUE,
)

# Numbers for text:
extr_gt %>% select(Type, hetfrq, GenClass) %>% group_by(Type, GenClass) %>% summarize(mean=mean(hetfrq), sd=sd(hetfrq)) %>% print(n=50)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# # # # # Homozygous derived sites
pB<-ggplot(extr_gt,aes(x=GenClass, y=homderfrq))+
geom_violin(fill="white", alpha=0.5)+
facet_grid(Type~factor(Grp, level=c("Reference","Scandinavia")), scales="free", space="free_x") +
  labs(x=NULL,y="proportion of homozygous derived genotypes")+
  theme(panel.grid.major = element_line(colour = 'white'),
        panel.background = element_rect(fill = '#f5f4e6', colour = '#FFFDF8'),
        axis.text.x = element_text(angle = 45, hjust=1),
				strip.text.x = element_text(size=12),
				strip.text.y = element_text(size=12),
        legend.box = "vertical") +
  theme(legend.key=element_rect(fill='white', colour='NA')) +
	geom_point(colour="grey", size=1, alpha=0.5)

output<-paste(plotdir,"Figure4B.pdf", sep="")
ggsave(
  output,
  plot = pB,
  scale = 1,
  width = 7.5,
  height = 6,
  unit = "in",
  dpi = 300,
  limitsize = TRUE,
)
