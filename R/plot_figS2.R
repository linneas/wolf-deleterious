#################
### SUPPLEMENTARY FIGURE 2
# Comparing Vep/Sift deleteriousness with Sneath and Miyata
# Thresholds 1.85 for Miyata and 25 for Sneath (correspond to ~35% and 30% respectively)

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
aa_file=paste("aa_prop/100S95F14R.mac2.properties.all.txt")
vcf <- read.vcfR(vcf_file)
meta_tib<-meta_file %>% read.table(header=TRUE) %>% as_tibble() %>% rename(Indiv=UU_ID, Cat=Category)
vep_tib <- vep_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))
sift_tib <- sift_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(CHROM = as.character(CHROM))
aa_tib <- aa_file %>% read.table(header=TRUE)%>% as_tibble() %>% mutate(CHROM = as.character(CHROM))



# Convert to tidy format
tidy_vcf <- vcfR2tidy(vcf,
	info_fields=c("AA"),
	format_fields=c("GT"),
	dot_is_NA=TRUE)

# Extracting relevant data and combine tibbles
# First a big general tibble
full_tib <- tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
	inner_join(tidy_vcf$fix) %>% left_join(aa_tib) %>%
	left_join(vep_tib) %>% left_join(sift_tib) %>%
	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1', (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT), Exchgb=if_else(ALT==AA, Exchgb2, Exchgb1))

# Then a table for deleterious SIFT (like before, just for comparison)
sift_tib <- full_tib %>% filter(SIFT_TYPE == 'deleterious') %>%
 	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
		(ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
	group_by(Indiv, new_gt) %>% summarize(count=n()) %>% ungroup() %>%
	mutate_if(is.character, str_replace_all, pattern = "/", replacement = "") %>%
	 pivot_wider(names_from=new_gt, names_prefix="gt", values_from=count) %>%
	mutate(sum=gt00+gt01+gt11) %>% inner_join(meta_tib)  %>%
	mutate(hetfrq=gt01/sum, homderfrq=gt11/sum, deralfreq=(gt11*2+gt01)/(sum*2)) %>%
	select(Indiv,Cat,hetfrq,homderfrq,deralfreq,GenClass) %>% mutate(DEL_TYPE="SIFT")

# And for each different AA repacement score
# Miyata
miyata_tib <- full_tib %>% filter(Miyata>1.85) %>%
 	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
		(ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
	group_by(Indiv, new_gt) %>% summarize(count=n()) %>% ungroup() %>%
	mutate_if(is.character, str_replace_all, pattern = "/", replacement = "") %>%
	 pivot_wider(names_from=new_gt, names_prefix="gt", values_from=count) %>%
	mutate(sum=gt00+gt01+gt11) %>% inner_join(meta_tib)  %>%
	mutate(hetfrq=gt01/sum, homderfrq=gt11/sum, deralfreq=(gt11*2+gt01)/(sum*2))  %>%
	select(Indiv,Cat,hetfrq,homderfrq,deralfreq,GenClass) %>% mutate(DEL_TYPE="Miyata")

#Sneath
sneath_tib <- full_tib %>% filter(Sneath>25) %>%
 	mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
		(ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
	group_by(Indiv, new_gt) %>% summarize(count=n()) %>% ungroup() %>%
	mutate_if(is.character, str_replace_all, pattern = "/", replacement = "") %>%
	 pivot_wider(names_from=new_gt, names_prefix="gt", values_from=count) %>%
	mutate(sum=gt00+gt01+gt11) %>% inner_join(meta_tib)  %>%
	mutate(hetfrq=gt01/sum, homderfrq=gt11/sum, deralfreq=(gt11*2+gt01)/(sum*2))  %>%
	select(Indiv,Cat,hetfrq,homderfrq,deralfreq,GenClass) %>% mutate(DEL_TYPE="Sneath")

# COMBINE
comb_tib <- bind_rows(sift_tib,sneath_tib,miyata_tib) %>%
						mutate(Grp=case_when(GenClass=='Finland' | GenClass=='Russia' | GenClass=='NR_immigrants' ~ 'Reference',
						TRUE ~ 'Scandinavia')) %>%
						mutate(GenClass=case_when(GenClass=="R_immigrants" ~ "Reproducing immigrants",
					                            GenClass=="NR_immigrants" ~ "Non-reproducing immigrants",
					                            TRUE ~ GenClass))
# Make factors
comb_tib$GenClass <- factor(comb_tib$GenClass, levels = c("Finland", "Russia", "Non-reproducing immigrants", "Founders", "F1", "F2", "F3", "F4", "F5", "F6", "Reproducing immigrants", "L1", "L2", "L3"))
comb_tib$DEL_TYPE <- factor(comb_tib$DEL_TYPE, levels = c("SIFT", "Miyata", "Sneath"))

################################# PLOTTING
# # # # # Heterozygous sites

phet<- ggplot(comb_tib, aes(x=GenClass, y=hetfrq))+
  geom_violin(fill="white", alpha=0.5)+
  facet_grid(DEL_TYPE~factor(Grp, level=c("Reference","Scandinavia")), scales="free", space="free_x") +
  labs(x=NULL,y="proportion of heterozygous genotypes")+
  theme(panel.grid.major = element_line(colour = 'white'),
        panel.background = element_rect(fill = '#f5f4e6', colour = '#FFFDF8'),
        axis.text.x = element_text(angle = 45, hjust=1),
				strip.text.x = element_text(size=12),
				strip.text.y = element_text(size=12),
        legend.box = "vertical") +
  theme(legend.key=element_rect(fill='white', colour='NA')) +
	geom_point(colour="grey", size=1, alpha=0.5)

output <- paste(plotdir,"FigureS2A.pdf", sep="")
ggsave(
  output,
  plot = phet,
  scale = 1,
  width = 7,
  height = 9,
  unit = "in",
  dpi = 300,
  limitsize = TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# # # # # Homozygous derived sites
phom<-ggplot(comb_tib,aes(x=GenClass, y=homderfrq))+
geom_violin(fill="white", alpha=0.5)+
facet_grid(DEL_TYPE~factor(Grp, level=c("Reference","Scandinavia")), scales="free", space="free_x") +
  labs(x=NULL,y="proportion of homozygous derived genotypes")+
  theme(panel.grid.major = element_line(colour = 'white'),
        panel.background = element_rect(fill = '#f5f4e6', colour = '#FFFDF8'),
        axis.text.x = element_text(angle = 45, hjust=1),
				strip.text.x = element_text(size=12),
				strip.text.y = element_text(size=12),
        legend.box = "vertical") +
  theme(legend.key=element_rect(fill='white', colour='NA')) +
	geom_point(colour="grey", size=1, alpha=0.5)

output <- paste(plotdir,"FigureS2B.pdf", sep="")
ggsave(
  output,
  plot = phom,
  scale = 1,
  width = 7,
  height = 9,
  unit = "in",
  dpi = 300,
  limitsize = TRUE)
