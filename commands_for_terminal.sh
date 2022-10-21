# BASH COMMANDS

# NOTE THAT MOST PROGRAMS ARE RUN AS BATCH JOBS SUBMITTED TO A SLURM CLUSTER
# (settings as number of nodes, requested time and memory are customized and
# might need to be tweaked if run on other systems)

# the commands assumes a working directory with sub directories for all steps
# (vcf, reference, vep, gerp, outgroup etc). The "help_files" should also be
# placed here. The scripts directory can be placed anywhere (specified here):
scrdir="./wolf-deleterious/scripts"
Rdir="./wolf-deleterious/R"
smkdir="./wolf-deleterious/snakemake"

# Most files have the following prefix, standing for "100 Scandinavian,
# 95 Finnish and 14 Russian wolves"
p1="100S95F14R"
# for the X chromosome, we only use the females with prefix
p2="74Females"

# Note that this project stared with gvcf files, produced for another project.
# The code can be found in: https://github.com/linneas/fennoscandian_wolf

########################### MAIN SNAKEMAKE PIPELINE ############################
# Runs joint genotyping, SNP extraction, filtering and VEP
#Make sure the following programs are available on the system:
#snakemake/5.30.1
#GATK/3.8-0
#htslib/1.12
#vcftools/0.1.16
#vep/99
# Plot the DAG (directed acyclic graph)
snakemake --snakefile $smkdir/main.snakefile --dag | dot -Tsvg > dag/dag.$p1.svg
# Run full pipeline
snakemake --snakefile $smkdir/main.snakefile -p -j 64  --cluster "sbatch -p {cluster.partition} -n {cluster.n} -C {cluster.C} -t {cluster.time} -e {cluster.error} -o {cluster.output}" --cluster-config $smkdir/cluster_main.json

# For the X chromosomes we only want females to avoid ploidy problems. Extract
# the X and the 74 females from the original vcf, before any filtering
sbatch -J chrX -t 2-00:00:00 -p core $scrdir/run_vcftools_extractIndListAndChrom.sh vcf/$p1.vcf.gz help_files/females.list "X" "vcf/$p2"
# Then the main snakemake pipeline can be used for filtering and VEP

# Make bed file with autosomes and ChrX
mkdir -p bed
for i in "1" "2"
do
 zcat vcf/$p1.SNPs.HF.mac$i.vcf.gz | awk '(/^[0-9]/){s=$2-1; print $1"\t"s"\t"$2}' >bed/$p1.chr1-38.mac$i.bed
 zcat vcf/$p2.SNPs.HF.mac$i.vcf.gz | awk '(/^X/){s=$2-1; print $1"\t"s"\t"$2}' >bed/74Females.chrX.mac$i.bed
done

# Note: The two filterings, called mac2 and mac1, only differ in the number of
# minor alleles - mac2 demands at least 2 minor alleles (stricter, will remove
# potential calling errors), while mac1 only demands a single allele, which
# includes all singletons but also potentially some call errors. The last is
# only used for site frequency spectrum analysis.

################################# POLARIZING ###################################
# Commands for polarizing SNPs using 2 or 5 outgroups (2 was used in the paper)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DOWNLOAD OUTGROUPS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mkdir fastq/
for ind in $(grep -v "#" help_files/SRA_accession_outgroups.txt |cut -f2)
do
  echo $ind
  sbatch -p core -t 10:00:00  $scrdir/run_downloadFastq.sh $ind help_files/SRA_accession_outgroups.txt
done

# ~~~~~~~~~~~~~~~~~~~~~~~~ OUTGROUP SNAKEMAKE PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~
# Pipeline that maps the outgroups onto the dog reference, makes allsites files
# and extracts the relevant sites (polymorphic in our wolves)
# Software needed:
#snakemake/5.30.1
#bwa/0.7.17
#samtools/1.10
#picard/2.23.4
#GATK/3.8-0
#BEDTools/2.29.2

# Plot DAG (directed acyclic graph)
snakemake --dag --snakefile $smkdir/outgroup.snakefile | dot -Tsvg > dag/dag.outgroups.svg
# Then start the actual pipeline
snakemake --snakefile $smkdir/outgroup.snakefile -p -j 64  --cluster "sbatch -p {cluster.partition} -n {cluster.n} -C {cluster.C} -t {cluster.time} -e {cluster.error} -o {cluster.output} --mail-type {cluster.mail-type} --mail-user {cluster.mail-user}" --cluster-config $smkdir/snake/cluster_outgroup.json


# ~~~~~~~~~~~~~~~~~~~ ASSIGN GENOTYPES (PSEUDO-HAPLOIDIZED) ~~~~~~~~~~~~~~~~~~~~

suffix="extract.noIndel"
d=5
for filt in "mac1" "mac2"
do
 for set in "$p2.chrX" "$p1.chr1-38"
 do
   s="$set.$filt"
  for i in  "AfGWo" "BlBJa" #"Dhole" "SiSJa" "EthWo"
  do
    # Remove indels
    #awk '{if(/^#/){print}else{if(length($4)==1 && length($5)==1 && $5!="*"){print}}}' outgroup/$i.$s.extract.vcf >outgroup/$i.$s.$suffix.vcf
    # Assign genotypes to outgroups (min depth=5 required)
    python3 $scrdir/pseudo_haploidize.py -v outgroup/$i.$s.$suffix.vcf -d $d -o outgroup/$i.$s.weightedAF.DP$d
    # Add Ns (for sites missing, for example if they fall in a deletion)
    intersectBed -wao -a bed/$s.bed -b outgroup/$i.$s.weightedAF.DP$d.bed |awk -v OFS="\t" '{if($8==0){$7="N"}; print}' |cut -f 7 >outgroup/$i.$s.weightedAF.DP$d.addN.txt
  done
 done
done

# ~~~~~~~~~~~~~~~~~~~~~ COMBINE AND ASSIGN ANCESTRAL STATE ~~~~~~~~~~~~~~~~~~~~~

# Using two outgroups (AfGWo and BlBJa, the two with highest coverage)
d=5
for filt in "mac1" "mac2"
do
  for set in "$p2.chrX" "$p1.chr1-38"
  do
    s="$set.$filt"
    # Combine sites with genotypes from the wanted outgroups
    paste bed/$s.bed outgroup/AfGWo.$s.weightedAF.DP$d.addN.txt outgroup/BlBJa.$s.weightedAF.DP$d.addN.txt >outgroup/2outgroups.$s.DP$d.txt
    # Assign ancestral using a python script, this prints the most common allele
    # together with the support (fraction of outgroups supporting this allele) and
    # number of outgroups with data
    python3 $scrdir/assign_ancestral.py -i outgroup/2outgroups.$s.DP$d.txt -c 4 -o outgroup/Ancestral.2outgroups.$s.DP$d 2>stderr.2out.$s.DP$d.txt
    # Extract sites with 100% agreement and no missing outgroups
    awk -v OFS="\t" '($4!="N" && $5=="1.00" && $6==2){print $1,$2,$3,$4}' outgroup/Ancestral.2outgroups.$s.DP$d.ancestral.txt >outgroup/Pol.2out.$s.bed
    # Also output as CHR:POS BASE to be used with old perl script
    awk '($4!="N" && $5=="1.00" && $6==2){print $1":"$3"\t"$4}' outgroup/Ancestral.2outgroups.$s.DP$d.ancestral.txt >outgroup/Pol.2out.$s.txt
  done
done

# Also tried using all 5 outgroups, this resulted in fewer calls due to lower
# coverage in some samples (demanding 100% agreement and data for at least 3 out
# of 5 outgroups)
for filt in "mac1" "mac2"
do
  for set in "$p2.chrX" "$p1.chr1-38"
  do
    s="$set.$filt"
  paste  bed/$s.bed outgroup/AfGWo.$s.weightedAF.DP$n.addN.txt outgroup/EthWo.$s.weightedAF.DP$d.addN.txt outgroup/Dhole.$s.weightedAF.DP$d.addN.txt outgroup/SiSJa.$s.weightedAF.DP$d.addN.txt outgroup/BlBJa.$s.weightedAF.DP$d.addN.txt >outgroup/5outgroups.$s.DP$d.txt
  python3 $scrdir/assign_ancestral.py -i outgroup/5outgroups.$s.DP$d.txt -c 4 -o outgroup/Ancestral.5outgroups.$s.DP$d 2>stderr.5out.$s.DP$d.txt
  awk -v OFS="\t" '($4!="N" && $5=="1.00" && $6>=3){print $1,$2,$3,$4}' outgroup/Ancestral.5outgroups.$s.DP$d.ancestral.txt >outgroup/Pol.5out.$s.bed
done




####################### DELETERIOUSNESS BASED ON VEP/SIFT ######################
# VEP was run as part of the main snakemake pipeline above. Below follow
# extractions of the different categories.

# Note: vep output includes annotation for everything in the input file
# (for example unplaced scaffolds). We filter this out later.
for filt in "mac1" "mac2"
do
  for set in $p1 $p2
  do
    extset="$set.SNPs.HF.$filt"
    s="$set.$filt"
    # Make a list with all different vep types and combinations present in the set
#    grep -v "#" vep/$extset.txt |cut -f7 |sort |uniq -c >vep/$s.vep_types.txt
    # Extract each set separately
    grep "IMPACT=MODIFIER" vep/$extset.txt |grep -v "intergenic_variant" |grep -v "intron" |grep -v "downstream" |grep -v "upstream" | awk '{split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' |uniq >vep/$s.firstExtract.modifier.bed
#    grep "synonymous" vep/$extset.txt |awk '($7=="synonymous_variant"){split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' |uniq >vep/$s.firstExtract.synonymous.bed
#    grep "missense" vep/$extset.txt |awk '($7=="missense_variant"){split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' |uniq >vep/$s.firstExtract.missense.bed
#    grep "missense" vep/$extset.txt |grep "tolerated" |grep -v "low_confidence" | awk '($7=="missense_variant"){split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' |uniq >vep/$s.firstExtract.tolerated.missense.bed
#    grep "missense" vep/$extset.txt |grep "deleterious" |grep -v "low_confidence" | awk '($7=="missense_variant"){split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' |uniq >vep/$s.firstExtract.deleterious.missense.bed
#    grep "IMPACT=HIGH" vep/$extset.txt |awk '{if($7~/stop_gained/){if($11~/*$/){print}}else if($7~/stop_lost/){if($11~/^*/){print}}else if($7~/start_lost/){if($11~/^M/){print}}else{print}}'|awk '{split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' |uniq >vep/$s.firstExtract.nonsense.bed
  done
done

# Note: Some sites have two or even more consequenses listed. For synonymous and
# missense above, I keep only SNPs with these consequences (if they are eg. both
#"synonymous" and "splice_region_variant", they are not included in synonymous).
# For nonsense, I save sites with multiple annotations. However, I only save
# start_lost/stop_lost/stop_gained if they make sense (I check that they lost an
# "M" or "*" or gained a "*")

# However, the files above can still overlap slightly. A snp can be synonymous
# in one gene and nonsynonymous in an overlapping gene (listed multiple times in
# the original file - not with multiple consequenses in one line)

# Keep only the most severe consequence (and combine to single bed)
# In the same time, save only chr1-38 and chrX
for filt in "mac1" "mac2"
do
  for set in "$p1.chr1-38" "$p2.chrX"
  do
    pref=`echo $set | cut -f1 -d"."`
    pf="$pref.$filt"
    s="$set.$filt"
#    intersectBed -v -a vep/$pf.firstExtract.synonymous.bed -b vep/$pf.firstExtract.missense.bed vep/$pf.firstExtract.nonsense.bed | intersectBed -a - -b bed/$s.bed >vep/$s.final.synonymous.bed
#    intersectBed -v -a vep/$pf.firstExtract.missense.bed -b vep/$pf.firstExtract.nonsense.bed | intersectBed -a - -b bed/$s.bed >vep/$s.final.missense.bed
#    intersectBed -v -a vep/$pf.firstExtract.tolerated.missense.bed -b vep/$pf.firstExtract.nonsense.bed vep/$pf.firstExtract.deleterious.missense.bed | intersectBed -a - -b bed/$s.bed >vep/$s.final.tolerated.missense.bed
#    intersectBed -v -a vep/$pf.firstExtract.deleterious.missense.bed -b vep/$pf.firstExtract.nonsense.bed | intersectBed -a - -b bed/$s.bed >vep/$s.final.deleterious.missense.bed
#    intersectBed -a vep/$pf.firstExtract.nonsense.bed -b bed/$s.bed >vep/$s.final.nonsense.bed
#    cat vep/$s.final.missense.bed vep/$s.final.nonsense.bed vep/$s.final.synonymous.bed |sort -k1,1 -k2,2n >bed/$s.vepfinal.bed
    intersectBed -v -a vep/$pf.firstExtract.modifier.bed -b bed/$s.vepfinal.bed >bed/$s.modifier.bed
  done
done

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXTRACT AND POLARIZE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract vep sites from the whole genome vcfs and polarize them

anc="Pol.2out"
for filt in "mac2" #"mac1"
do+
  for set in "$p1.chr1-38" "$p2.chrX"
  do
    pref=`echo $set | cut -f1 -d"."`
    chr=`echo $set | cut -f2 -d"."`
    s="$set.$filt"
#   mkdir -p vcf/$anc.$filt/
#   mkdir -p bed/$anc.$filt/
#    intersectBed -header -a vcf/$pref.SNPs.HF.$filt.vcf.gz -b bed/$s.vepfinal.bed >vcf/$s.vepfinal.vcf
#    perl $scrdir/addAAInfoToVCF_fromList.pl vcf/$s.vepfinal.vcf outgroup/$anc.$s.txt vcf/$anc.$filt/$set.vepfinal.vcf
    # Make bed file:
#    awk -v OFS="\t" '($1!~/^#/){s=$2-1; print $1,s,$2}' vcf/$anc.$filt/$set.vepfinal.vcf >bed/$anc.$filt/$set.vepfinal.bed
    #MODIFIER
    intersectBed -header -a vcf/$pref.SNPs.HF.$filt.vcf.gz -b bed/$s.modifier.bed >vcf/$s.modifier.vcf
    perl $scrdir/addAAInfoToVCF_fromList.pl vcf/$s.modifier.vcf outgroup/$anc.$s.txt vcf/$anc.$filt/$set.modifier.vcf
    awk -v OFS="\t" '($1!~/^#/){s=$2-1; print $1,s,$2}' vcf/$anc.$filt/$set.modifier.vcf >bed/$anc.$filt/$set.modifier.bed
  done
done

# Make tables with vep and sift sites (that could be polarized) to be used in R
anc="Pol.2out"
for filt in "mac1" "mac2"
do
 for set in "$p1.chr1-38" "$p2.chrX"
do
  pref=`echo $set |cut -f1 -d"."`
  s="$set.$filt"
  poldir="$anc.$filt"
  mkdir -p vep/$poldir
  echo "CHROM POS VEP_TYPE" |sed 's/ /\t/g' >vep/$poldir/$set.veptypes.txt
  rm -f $set.veptypes.tmp
  for type in "missense" "nonsense" "synonymous"
  do
    intersectBed -a vep/$s.final.$type.bed -b vcf/$poldir/$set.vepfinal.vcf |awk -v t=$type '{print $1"\t"$3"\t"t}' >>$set.veptypes.tmp
  done

  cat $set.veptypes.tmp | sort -k1,1 -k2,2n  >> vep/$poldir/$set.veptypes.txt
  rm $set.veptypes.tmp
  # SIFT
  echo "CHROM POS SIFT_TYPE" |sed 's/ /\t/g' >vep/$poldir/$set.sifttypes.txt
  rm -f $set.sifttypes.tmp
  intersectBed -v -a bed/$s.vepfinal.bed -b vep/$s.final.deleterious.missense.bed vep/$s.final.tolerated.missense.bed |awk '{print $1"\t"$3"\tNA"}' >$s.sifttypes.tmp
  intersectBed -a vep/$s.final.deleterious.missense.bed -b vcf/$poldir/$set.vepfinal.vcf |awk '{print $1"\t"$3"\tdeleterious"}' >>$set.sifttypes.tmp
  intersectBed -a vep/$s.final.tolerated.missense.bed -b vcf/$poldir/$set.vepfinal.vcf |
  awk '{print $1"\t"$3"\ttolerated"}' >>$set.sifttypes.tmp
  sort -k1,1 -k2,2n $set.sifttypes.tmp >>vep/$poldir/$set.sifttypes.txt
  rm $set.sifttypes.tmp
  done
done

# Merge the autosome and X
anc="Pol.2out"
for filt in "mac1" "mac2"
do
  poldir="$anc.$filt"
#  cat vep/$poldir/$p1.chr1-38.veptypes.txt <(tail -n+2 vep/$poldir/$p2.chrX.veptypes.txt) >vep/$poldir/veptypes.chr1-X.txt
#  cat vep/$poldir/$p1.chr1-38.sifttypes.txt  <(tail -n+2 vep/$poldir/$p2.chrX.sifttypes.txt) >vep/$poldir/sifttypes.chr1-X.txt
  echo "CHROM POS SIFT_TYPE" |sed 's/ /\t/g' >vep/$poldir/modifier.chr1-X.txt
  cat bed/$poldir/100S95F14R.chr1-38.modifier.bed bed/$poldir/74Females.chrX.modifier.bed |awk -v OFS='\t' '{print $1,$3,"modifier"}' >>vep/$poldir/modifier.chr1-X.txt
done

# FOR REVISION: Divide nonsense further into three classes; splice, start and stop
anc="Pol.2out"
for filt in "mac1" "mac2"
do
  poldir="$anc.$filt"
  echo "CHROM POS NONS_TYPE" |sed 's/ /\t/g' >vep/$poldir/nonsensetypes.chr1-38.txt
  grep "nonsense" vep/$poldir/100S95F14R.chr1-38.veptypes.txt |cut -f1,2 | sed 's/\t/:/' |sort |join -1 1 -2 2 - <(grep "HIGH" vep/100S95F14R.SNPs.HF.$filt.txt |sort -k2,2)  | cut -f1,7 -d" " |uniq |awk '{if(/splice_acceptor_variant/){t="splice"}else if(/splice_donor_variant/){t="splice"}else if(/start_lost/){t="start"}else if(/stop/){t="stop"}else{t="other"}; print $1"\t"t}' |sed 's/:/\t/' |sort -k1,1n -k2,2n >>vep/$poldir/nonsensetypes.chr1-38.txt
done

# ~~~~~~~~~~~~~~~~~~~~~~~~~~ INFERRING FOUNDER MALES ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Autosomes
anc="Pol.2out"
set=$p1.chr1-38
for filt in "mac2" "mac1"
do
  poldir=$anc.$filt
#  mkdir -p inferred/$poldir/
  # Connect haplotypes to genotypes for each SNP position
#  python3 $scrdir/assign_genotype_to_haplotype.py -v vcf/$poldir/$set.vepfinal.vcf -p help_files/76Scand_haplotypes.txt -o inferred/$poldir/$set.vepfinal
  # Use this and the infered haplotypes for founder males from Viluma et al.
  # to infer "fake" genotypes.
#  python3 $scrdir/infer_founder_genotypes.py -a inferred/$poldir/$set.vepfinal.assigned.haplotypes.txt -p help_files/3Founders_haplotypes.txt -o inferred/$poldir/3Founders.chr1-38.vepfinal.fake.genotypes.txt
  #Add these to the vcf file
#  awk '(NR>1){s=$2-1; print $1"\t"s"\t"$2"\t"$4"::"$5}' inferred/$poldir/3Founders.chr1-38.vepfinal.fake.genotypes.txt |intersectBed -wo -a vcf/$poldir/$set.vepfinal.vcf -b - |cut -f1-218,222 |sed 's/::/\t/g' |cat <(grep "#" vcf/$poldir/$set.vepfinal.vcf|awk '{if(/#CHROM/){print $0"\tFM1\tFM2"}else{print}}') - >vcf/$poldir/$set.vepfinal.wfm.vcf
  #MODIFIER
  python3 $scrdir/assign_genotype_to_haplotype.py -v vcf/$poldir/$set.modifier.vcf -p help_files/76Scand_haplotypes.txt -o inferred/$poldir/$set.modifier
  python3 $scrdir/infer_founder_genotypes.py -a inferred/$poldir/$set.modifier.assigned.haplotypes.txt -p help_files/3Founders_haplotypes.txt -o inferred/$poldir/3Founders.chr1-38.modifier.fake.genotypes.txt
  awk '(NR>1){s=$2-1; print $1"\t"s"\t"$2"\t"$4"::"$5}' inferred/$poldir/3Founders.chr1-38.modifier.fake.genotypes.txt |intersectBed -wo -a vcf/$poldir/$set.modifier.vcf -b - |cut -f1-218,222 |sed 's/::/\t/g' |cat <(grep "#" vcf/$poldir/$set.modifier.vcf|awk '{if(/#CHROM/){print $0"\tFM1\tFM2"}else{print}}') - >vcf/$poldir/$set.modifier.wfm.vcf
done

# X chromomosome
# The script adds the genotype in diploid form to match with the other genotypes
# (so that if one male has allele 1, it will get the genotype 1/1. This is taken
# care of later so that the males are treated as haploid in R)
anc="Pol.2out"
set=$p2.chrX
for filt in "mac1" "mac2"
do
  poldir=$anc.$filt
  python3 $scrdir/assign_genotype_to_haplotype_haploid.py -v vcf/$poldir/$set.vepfinal.vcf -p help_files/76Scand_haplotypes.txt -o inferred/$poldir/$set.vepfinal
  python3 $scrdir/infer_founder_genotypes_haploid.py -a inferred/$poldir/$set.vepfinal.assigned.haplotypes.txt -p help_files/3Founders_haplotypes.txt -o inferred/$poldir/3Founders.chrX.vepfinal.fake.genotypes.txt
  #Add these to the vcf file
  awk '(NR>1){s=$2-1; print $1"\t"s"\t"$2"\t"$4"::"$5}' inferred/$poldir/3Founders.chrX.vepfinal.fake.genotypes.txt |intersectBed -wo -a vcf/$poldir/$set.vepfinal.vcf -b - |cut -f1-83,87 |sed 's/::/\t/g' |cat <(grep "#" vcf/$poldir/$set.vepfinal.vcf|awk '{if(/#CHROM/){print $0"\tFM1\tFM2"}else{print}}') - >vcf/$poldir/$set.vepfinal.wfm.vcf
done

################ DELETERIOUSNESS BASED ON AMINO-ACID PROPERTIES ################

# The following python script assumes three tables named table.ex.csv,
# table.sneath.csv and table.miyata.csv placed in the subdirectory help_files

mkdir -p aa_prop/
python3 $scrdir/replacement_impact.py -i vep/$p1.SNPs.HF.mac2.txt -o aa_prop/$p1.mac2.properties.all.txt

# This script needs info on amino acid changes and therefore takes the VEP
# output as input (which nicely annotates all changes, including amino acids
# for non-synonymous changes). Except for this it is not related to the VEP
# analysis in any way.

# Just as for VEP, there are some SNPs with multiple annotations (due to over-
# lapping genes). I only save the change with the highest impact (that means,
# highest for Miyata/Sneath, lowest for Exchgb)
cut -f1,2,6- aa_prop/$p1.mac2.properties.all.txt| awk -v OFS="\t" '{if(NR==1){print}else{if(NR==2){lc=$1; lpos=$2; lmi=$3; lexr=$4; lexv=$5; lsn=$6}else{if(lpos==$2){if($3>lmi){lmi=$3}; if($4<lexr){lexr=$4}; if($5<lexv){lexv=$5}; if($6>lsn){lsn=$6}}else{print lc,lpos,lmi,lexr,lexv,lsn; lc=$1; lpos=$2; lmi=$3; lexr=$4; lexv=$5; lsn=$6}}}}'  >aa_prop/$p1.mac2.properties.filt.txt



######################## DELETERIOUSNESS BASED ON GERP #########################
# For the GERP part, we only use the stricter SNP filtering (mac2), and only
# autosomes.
mkdir -p maf/
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/hg38.100way.nh maf/
cat maf/hg38.100way.nh |sed 's/\s//g' |tr "\n" "+" |sed 's/+//g' >maf/hg38.tree

# Run snakemake pipeline that downloads maf alignments from UCSC
# Make sure the following are installed:
#snakemake/5.30.1
#MafFilter/1.1.2
#GERP++/20110522

# Plot the DAG
snakemake --snakefile $smkdir/gerp.snakefile --dag | dot -Tsvg > dag/gerp_pipeline.svg
# Run pipeline
snakemake --snakefile $smkdir/gerp.snakefile -p -j 64  --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time} -e {cluster.error} -o {cluster.output}" --cluster-config $smkdir/cluster_gerp.json

# ~~~~~~~~~~~~~~~~~~~~~~~~~ TRANSLATE TO DOG REFERENCE ~~~~~~~~~~~~~~~~~~~~~~~~~

# Files are sorted after human chromosomes now.
# The following code sorts the files after the dog genome.
# Note! If run several times files chr*.rates.bed.chr* have to be deleted
# before the rerun, since the script only appends, not overwrites.
for hchr in {1..22}
do
  infile="gerp/canFam3/hg38.chr$hchr.rates.bed"
  job=$(sbatch -J h.chr$hchr.divide -t 1:00:00 -p core $scrdir/run_divideLiftover.sh $infile |cut -f4 -d" ")
  for dchr in {1..38}
  do
    sbatch -J h.chr$hchr.c.chr$dchr.sort -d afterok:$job -t 4:00:00 -p core $scrdir/run_sortLiftover.sh $infile.chr$dchr
  done
done

for dchr in {1..38}
do
   outf="gerp/canFam3/chr$dchr.rates.bed"
   echo '#!/bin/bash
   sort -m -k2,2n --buffer-size=6G gerp/canFam3/hg38.chr*.rates.bed.chr'$dchr'.sorted >'$outf'
   '| sbatch -J merge.$dchr -t 1:00:00 -p core
done

# Because of some regions can have multiple hits in the alignment, some positions
# can occur twice (or more) in the output file.
# Remove all multi-occuring positions, to reduce bias caused by alignment issues
for dchr in {1..38}
do
  sbatch -J removedupl.$dchr -t 4:00:00 -p core $scrdir/run_remove_liftover_duplicates.sh gerp/canFam3/chr$dchr.rates.bed gerp/canFam3/chr$dchr.rates.unique.bed
done

 # Make a merged version
 echo '#!/bin/bash
 rm -f gerp/canFam3/Autosomes.rates.unique.bed
for dchr in {1..38}
do
 cat chr$dchr.rates.unique.bed >>gerp/canFam3/Autosomes.rates.unique.bed
done
'| sbatch -J merge -t 5:00:00 -p core


# ~~~~~~~~~~~~~~~~~~~~~~~~~  POLARIZE SNPs GENOME-WIDE ~~~~~~~~~~~~~~~~~~~~~~~~~
anc="Pol.2out"
chr="chr1-38"
for filt in "mac1" "mac2"
do
#  mkdir -p bed/$anc.$filt
#  mkdir -p gerp/canFam3/$anc.$filt
#  perl $scrdir/addAAInfoToVCF_fromList.pl <(zcat vcf/$p1.SNPs.HF.$filt.vcf.gz) outgroup/$anc.$p1.$chr.$filt.txt vcf/$anc.$filt/$p1.$chr.allSNPs.vcf
  # And make a bed with these...
#  awk '(/^[0-9]/){s=$2-1; print $1"\t"s"\t"$2}' vcf/$anc.$filt/$p1.$chr.allSNPs.vcf >bed/$anc.$filt/$p1.$chr.allSNPs.bed
  # .. so we can extract the wanted GERP scores
#  intersectBed -a <(sed "s/chr//" gerp/canFam3/Autosomes.rates.unique.bed) -b bed/$anc.$filt/$p1.$chr.allSNPs.bed >gerp/canFam3/$anc.$filt/$p1.$chr.allSNPs.rates.unique.bed
  # And intersect with VEP (for plotting)
  #intersectBed -a <(sed "s/chr//" gerp/canFam3/Autosomes.rates.unique.bed) -b bed/$anc.$filt/$p1.$chr.vepfinal.bed >gerp/canFam3/$anc.$filt/$p1.$chr.vepfinal.rates.unique.bed
  intersectBed -a <(sed "s/chr//" gerp/canFam3/Autosomes.rates.unique.bed) -b bed/$anc.$filt/$p1.$chr.modifier.bed >gerp/canFam3/$anc.$filt/$p1.$chr.modifier.rates.unique.bed
done

# ~~~~~~~~~~~~~~~~~~~~~~~~~~ INFERRING FOUNDER MALES ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Autosomes (For the GERP part, we only use autosomes)
anc="Pol.2out"
set=$p1.chr1-38
for filt in "mac2" "mac1"
do
  poldir=$anc.$filt
  mkdir -p inferred/$poldir/
  # Connect haplotypes to genotypes for each SNP position
  python3 $scrdir/assign_genotype_to_haplotype.py -v vcf/$poldir/$set.allSNPs.vcf -p help_files/76Scand_haplotypes.txt -o inferred/$poldir/$set.allSNPs
  # Use this and the infered haplotypes for founder males from Viluma et al.
  # to infer the genotypes.
  python3 $scrdir/infer_founder_genotypes.py -a inferred/$poldir/$set.allSNPs.assigned.haplotypes.txt -p help_files/3Founders_haplotypes.txt -o inferred/$poldir/3Founders.allSNPs.fake.genotypes.txt
  #Add these to the vcf file
  awk '(NR>1){s=$2-1; print $1"\t"s"\t"$2"\t"$4"::"$5}' inferred/$poldir/3Founders.allSNPs.fake.genotypes.txt |intersectBed -wo -a vcf/$poldir/$set.allSNPs.vcf -b - |cut -f1-218,222 |sed 's/::/\t/g' |cat <(grep "#" vcf/$poldir/$set.allSNPs.vcf|awk '{if(/#CHROM/){print $0"\tFM1\tFM2"}else{print}}') - >vcf/$poldir/$set.allSNPs.wfm.vcf
done


# ~~~~~~~~~~~~~~~~~~~~  CALCULATE GERP LOAD PER INDIVIDUAL ~~~~~~~~~~~~~~~~~~~~~
anc="Pol.2out"
filt="mac2"
chr="chr1-38"
# Calculating GERP mutation load using different GERP threshold (4 is used for
# the paper)
for gerp in 10 # 2 4 6 8 10
do
  python3 $scrdir/gerp_load_polarized.py -l <(cut -f1,3 help_files/metadata.txt |tail -n+2) -v vcf/$anc.$filt/$p1.$chr.allSNPs.wfm.vcf -b outgroup/$anc.$p1.$chr.$filt.bed -g gerp/canFam3/$anc.$filt/$p1.$chr.rates.unique.bed -t $gerp -o gerp/canFam3/$anc.$filt/LoadPerInd.thres$gerp.$p1.chr1-38.allSNPs.wfm.txt
done


################################ R_xy ANALYSIS #################################
# Added for revision 2022.10.11

 # For the Rxy analysis we want to normalize the statistics based on some non
 # genic control sites. Just like in Grossen et al we use intronic sequence,
 # and further we want to use the same gene set as where we have the coding
 # variants, and only choose among sites with a GERP around 0 (for example >-1
 # and <1.)

 anc="Pol.2out"
 filt="mac2"
 chr="chr1-38"

# List of genes with VEP variants (from above):
intersectBed -a gtf/CanFam3.1.102.CDS.with.gene.names.unmerged.bed -b bed/$anc.$filt/$p1.$chr.vepfinal.bed |cut -f4 |sort |uniq >$anc.$filt.$chr.list_of_genes.txt
# Bed file with all intronic variants in those genes (>3.5M)
grep -f list_of_genes.txt vep/$p1.SNPs.HF.$filt.txt |grep "intron_variant" |awk -v OFS="\t" '{split($s1,t,"_"); s=t[2]-1; print t[1],s,t[2]}' |uniq >$anc.$filt.$chr.list_of_intron_SNPs.bed
# Intersect with gerp scores and only keep sites with GERP >-1 && <1
intersectBed -a gerp/canFam3/$anc.$filt/$p1.$chr.allSNPs.rates.unique.bed -b $anc.$filt.$chr.list_of_intron_SNPs.bed |awk '($5>-1 && $5<1){print}' >$anc.$filt.$chr.introns_with_neutral_gerp.bed
# Extract those variants from VCF file
intersectBed -header -a vcf/$anc.$filt/$p1.$chr.allSNPs.wfm.vcf -b $anc.$filt.$chr.introns_with_neutral_gerp.bed >vcf/$anc.$filt/$p1.$chr.intron.wfm.vcf
########################## STATS AND PLOTTING WITH R ###########################
#R version used: R/4.1.1

# All the stats and numbers in the text and the tables are calculated with R.
# There are two general scripts, and one script per figure (that sometimes also
# contains code for extracting numbers). Everything was run interactively, one
# command at the time.

# R scripts for calculating VEP and GERP stats (for the tables in the paper)
Rscript $Rdir/calculate_Vep_stats_for_tables.R
Rscript $Rdir/calculate_GERP_stats_for_tables.R

# R scripts for all the figures
for i in 1 2 3 4 5 S1 S2 S3 S4
do
Rscript $Rdir/plot_fig$i.R
done
