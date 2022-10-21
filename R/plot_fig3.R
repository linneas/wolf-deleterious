#################
### R code for plotting Figure 3
### R_xy analysis to compare relative numbers of alleles between populations

rm(list=ls())

# Setting up, loading R libraries and set working directory
require(vcfR)
require(tidyverse)
plotdir="plots/"

# For this analysis we use MAC>=2
filt="mac2"
anc=paste("Pol.2out.",filt, sep="")

# Reading in the data
vcf_file=paste("vcf/",anc,"/100S95F14R.chr1-38.vepfinal.wfm.vcf", sep="")
vcf_file2=paste("vcf/",anc,"/100S95F14R.chr1-38.modifier.wfm.vcf", sep="")
vcf_file3=paste("vcf/",anc,"/100S95F14R.chr1-38.intron.wfm.vcf", sep="")
meta_file="help_files/metadata.txt"
vep_file=paste("vep/",anc,"/veptypes.chr1-X.txt", sep="")
sift_file=paste("vep/",anc,"/sifttypes.chr1-X.txt", sep="")
vcf <- read.vcfR(vcf_file)
vcf2 <- read.vcfR(vcf_file2)
vcf3 <- read.vcfR(vcf_file3)
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
tidy_vcf3 <- vcfR2tidy(vcf3,
  info_fields=c("AA"),
  format_fields=c("GT"),
  dot_is_NA=TRUE)


# Combine coding and modifier
combined_gt<-bind_rows(tidy_vcf$gt, tidy_vcf2$gt)
combined_fix<-bind_rows(tidy_vcf$fix, tidy_vcf2$fix)

# Start with intron SNPs to calculate the sum to normalize on
intron_tib <- tidy_vcf3$gt %>% filter(!is.na(gt_GT)) %>%
      inner_join(tidy_vcf3$fix) %>% mutate(Type="intron") %>% inner_join(meta_tib) %>%
      mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
                              (ALT==AA & gt_GT=='1/1') ~'0/0', TRUE ~ gt_GT)) %>%
      mutate(Pop=ifelse((GenClass=="F5" | GenClass=="F6"), "Last", GenClass)) %>%
      select(CHROM, POS, Indiv, new_gt, Pop) %>%
    	group_by(CHROM, POS, new_gt, Pop) %>% summarize(count=n()) %>%
      mutate(ancestral=case_when((new_gt=='0/0') ~ (count*2.0),
                                  (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0),
            derived=case_when((new_gt=='1/1') ~ (count*2.0),
                              (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0)) %>%
    	group_by(CHROM, POS, Pop) %>%
      summarize(totanc=sum(ancestral), totder=sum(derived))

# Calculate derived allele freq and pivot_wider to have one site per line
intr_dfreq_tib <- intron_tib %>% mutate(dfreq=totder/(totanc+totder)) %>% ungroup() %>%
                    select(CHROM, POS, Pop, dfreq) %>% group_by(CHROM,POS) %>%
                    pivot_wider(names_from=Pop, names_prefix="dfreq", values_from=dfreq)


# COMPARE LAST WITH RUSSIA
norm_LR <- (intr_dfreq_tib %>%  mutate(fact=dfreqLast*(1-dfreqRussia)) %>%
            ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
            summarize(s=sum(fact)))[[1]]

norm_RL <- (intr_dfreq_tib %>%  mutate(fact=dfreqRussia*(1-dfreqLast)) %>%
            ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
            summarize(s=sum(fact)))[[1]]

# COMPARE FOUNDER WITH RUSSIA
norm_FR <- (intr_dfreq_tib %>%  mutate(fact=dfreqFounders*(1-dfreqRussia)) %>%
            ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
            summarize(s=sum(fact)))[[1]]

norm_RF <- (intr_dfreq_tib %>%  mutate(fact=dfreqRussia*(1-dfreqFounders)) %>%
            ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
            summarize(s=sum(fact)))[[1]]

# COMPARE FINISH WITH RUSSIA
norm_FiR <- (intr_dfreq_tib %>%  mutate(fact=dfreqFinland*(1-dfreqRussia)) %>%
            ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
            summarize(s=sum(fact)))[[1]]

norm_RFi <- (intr_dfreq_tib %>%  mutate(fact=dfreqRussia*(1-dfreqFinland)) %>%
            ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
            summarize(s=sum(fact)))[[1]]

# COMPARE LAST WITH FOUNDERS
norm_LF <- (intr_dfreq_tib %>%  mutate(fact=dfreqLast*(1-dfreqFounders)) %>%
            ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
            summarize(s=sum(fact)))[[1]]

norm_FL <- (intr_dfreq_tib %>%  mutate(fact=dfreqFounders*(1-dfreqLast)) %>%
            ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
            summarize(s=sum(fact)))[[1]]



##########################################################################
### CALCULATE RXY FOR VEP CATEGORIES, USING BLOCK JACKKNIFE REPLICATES

# First, make combined table to start with
big_tib <- combined_gt %>% filter(!is.na(gt_GT)) %>%
      inner_join(combined_fix) %>% left_join(vep_tib) %>% left_join(sift_tib) %>%
    	rename(Type=VEP_TYPE) %>% inner_join(meta_tib) %>%
      mutate(new_gt=case_when((ALT==AA & gt_GT=='0/0') ~'1/1',
                              (ALT==AA & gt_GT=='1/1') ~'0/0',
                              TRUE ~ gt_GT)) %>%
      mutate(Type=case_when((SIFT_TYPE=="deleterious" | SIFT_TYPE=="tolerated") ~ SIFT_TYPE,
                            (is.na(Type) ~ 'modifier'),
                            TRUE ~ Type)) %>%
      mutate(Pop=ifelse(GenClass=="F5" | GenClass=="F6", "Last", GenClass)) %>%
      select(CHROM, POS, Indiv, new_gt, Pop, Type) %>%
    	group_by(CHROM, POS, new_gt, Pop, Type) %>% summarize(count=n()) %>%
      mutate(ancestral=case_when((new_gt=='0/0') ~ (count*2.0),
                                  (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0),
            derived=case_when((new_gt=='1/1') ~ (count*2.0),
                              (new_gt=='0/1') ~ (count*1.0), TRUE ~ 0))

# Then, set the categories to look at
type<-c("nonsense","deleterious","synonymous","modifier","tolerated")

# Number of Jackknife replicates
num_groups = 100

# Create tibble for saving all the values
all_val_tib <- tibble(i=integer(), type=character(), LR=numeric(),
                      RL=numeric(), FR=numeric(), RF=numeric(),
                      FiR=numeric(), RFi=numeric(), LF=numeric(), FL=numeric())

# Run The replicates for each category
for (t in 1:length(type)) {

  show(paste("subsetting big tib for", type[t]))
  subset_tib <- new_big_tib %>% filter(Type==type[t])

  wanted_sites <- subset_tib %>% ungroup() %>%
                  select(CHROM,POS) %>% group_by(CHROM,POS) %>% summarize()

  gl<- wanted_sites %>%
     group_by((row_number()-1) %/% (n()/num_groups)) %>%
     nest %>% pull(data)

  for (i in 1:num_groups) {
    show(paste("Type:", type[t], "Iteration:", i))
    temp_tib <- subset_tib %>% anti_join(gl[[i]]) %>% group_by(CHROM, POS, Pop) %>%
          summarize(totanc=sum(ancestral), totder=sum(derived)) %>%
          mutate(dfreq=totder/(totanc+totder)) %>% ungroup() %>%
          select(CHROM, POS, Pop, dfreq) %>% group_by(CHROM,POS) %>%
          pivot_wider(names_from=Pop, names_prefix="dfreq", values_from=dfreq)

    # COMPARE LAST WITH RUSSIA
    tmp_LR <- (temp_tib %>%  mutate(fact=dfreqLast*(1-dfreqRussia)) %>%
              ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
              summarize(s=sum(fact)))[[1]]

    tmp_RL <- (temp_tib %>%  mutate(fact=dfreqRussia*(1-dfreqLast)) %>%
              ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
              summarize(s=sum(fact)))[[1]]

    # COMPARE FOUNDER WITH RUSSIA
    tmp_FR <- (temp_tib %>%  mutate(fact=dfreqFounders*(1-dfreqRussia)) %>%
              ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
              summarize(s=sum(fact)))[[1]]

    tmp_RF <- (temp_tib %>%  mutate(fact=dfreqRussia*(1-dfreqFounders)) %>%
              ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
              summarize(s=sum(fact)))[[1]]

    # COMPARE FINNISH WITH RUSSIA
    tmp_FiR <- (temp_tib %>%  mutate(fact=dfreqFinland*(1-dfreqRussia)) %>%
              ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
              summarize(s=sum(fact)))[[1]]

    tmp_RFi <- (temp_tib %>%  mutate(fact=dfreqRussia*(1-dfreqFinland)) %>%
              ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
              summarize(s=sum(fact)))[[1]]

    # COMPARE LAST WITH FOUNDERS
    tmp_LF <- (temp_tib %>%  mutate(fact=dfreqLast*(1-dfreqFounders)) %>%
              ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
              summarize(s=sum(fact)))[[1]]

    tmp_FL <- (temp_tib %>%  mutate(fact=dfreqFounders*(1-dfreqLast)) %>%
              ungroup() %>% select(fact) %>% filter(!is.na(fact)) %>%
              summarize(s=sum(fact)))[[1]]


    #ADD ALL NUMBERS
    all_val_tib<- add_row(all_val_tib, i=i, type=type[t], LR=tmp_LR,
                          RL=tmp_RL, FR=tmp_FR, RF=tmp_RF,
                          FiR=tmp_FiR, RFi=tmp_RFi, LF=tmp_LF, FL=tmp_FL)
  }
}


# Calculate Rxy
rxy_tib<-all_val_tib %>% mutate(Last=(LR/norm_LR)/(RL/norm_RL),
                                Founder=(FR/norm_FR)/(RF/norm_RF),
                                Finland=(FiR/norm_FiR)/(RFi/norm_RFi),
                                Last2Founder=(LF/norm_LF)/(FL/norm_FL)) %>%
        pivot_longer(Last:Last2Founder, names_to="Pop", values_to="Rxy")

# Make factors
rxy_tib$Pop<-as.factor(rxy_tib$Pop)
rxy_tib$type <- factor(rxy_tib$type, levels = c("modifier", "synonymous", "tolerated", "deleterious", "nonsense"))




################################################################################
# Plot Last scandinavians vs Russia and Founders

red_rxy <- rxy_tib %>% filter(Pop=="Last2Founder" | Pop=="Last")
red_rxy$Pop <- factor(red_rxy$Pop, levels = c("Last","Last2Founder"))
red_rxy$type <- factor(red_rxy$type, levels = c("nonsense", "deleterious", "tolerated", "synonymous", "modifier"))

outfile=paste(plotdir,"Figure3.revision.pdf", sep="")

# plot boxes vertically
p<-ggplot(red_rxy, aes(x=type, y=Rxy)) +
  geom_boxplot(aes(colour=type, fill=type), alpha=0.5) +
  coord_flip() +
  facet_wrap(~Pop) +
  scale_color_manual(values=c("chocolate", "dodgerblue4", "seagreen4", "lightblue", "darkgray")) +
  scale_fill_manual(values=c("chocolate", "dodgerblue4", "seagreen4", "lightblue", "darkgray")) +
  theme(panel.grid.major = element_line(colour = 'white'),
        panel.background = element_rect(fill = '#f5f4e6', colour = '#FFFDF8'),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position="none") +
  labs(x="", y=expression(R[XY]))

ggsave(
  outfile,
  plot = p,
  scale = 1,
  dpi = 300,
  limitsize = TRUE,
  width=7,
  height=3.5,
)

# SOME STATS
lastvsrus_tib <- red_rxy %>% filter(Pop=="Last") %>% select(type, Rxy)
lastvsfou_tib <- red_rxy  %>% filter(Pop=="Last2Founder") %>% select(type, Rxy)
res1<-aov(Rxy ~ type, data=lastvsrus_tib)
res2<-aov(Rxy ~ type, data=lastvsfou_tib)