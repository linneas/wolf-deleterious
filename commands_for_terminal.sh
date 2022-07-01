# BASH COMMANDS

# NOTE THAT MOST PROGRAMS ARE RUN AS BATCH JOBS SUBMITTED TO A SLURM CLUSTER
# (settings as number of nodes, requested time and memory are customized and
# might need to be tweaked if run on other systems)

# the commands assumes a working directory with sub directories for all steps
# (vcf, reference, vep, gerp, outgroup etc.) The script repository can be placed
# within the working directory, or elsewhere:
scrdir="./wolf-deleterious/"

# Most files have the following prefix, standing for "100 Scandinavian,
# 95 Finnish and 14 Russian wolves"
prefix="100S95F14R"

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
# First dry run
snakemake --snakefile $scrdir/main.snakefile -np
# Then plot the DAG
snakemake --snakefile $scrdir/main.snakefile --dag | dot -Tsvg > dag/dag.$prefix.svg
# Run full pipeline
snakemake --snakefile $scrdir/main.snakefile -p -j 64  --cluster "sbatch -p {cluster.partition} -n {cluster.n} -C {cluster.C} -t {cluster.time} -e {cluster.error} -o {cluster.output}" --cluster-config $scrdir/cluster_main.json


# For the X chromosomes we only want females to avoid ploidy problems. Extract
# the X and the 74 females from the original vcf, before any filtering
sbatch -J chrX -t 2-00:00:00 -p core $scrdir/run_vcftools_extractIndListAndChrom.sh vcf/$prefix.vcf.gz help_files/females.list "X" "vcf/74Females"
# Then the main snakemake pipeline can be used for filtering and VEP

# Make bed file with autosomes and ChrX
mkdir -p bed
for i in "1" "2"
do
 zcat vcf/$prefix.SNPs.HF.mac$i.vcf.gz | awk '(/^[0-9]/){s=$2-1; print $1"\t"s"\t"$2}' >bed/$prefix.mac$i.chr1-38.bed
 zcat vcf/74Females.SNPs.HF.mac$i.vcf.gz | awk '(/^[0-9]/){s=$2-1; print $1"\t"s"\t"$2}' >bed/74Females.mac$i.chrX.bed
done



################################# POLARIZING ###################################
# Commands for polarizing SNPs using 2 or 5 outgroups (2 was used in the paper)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DOWNLOAD OUTGROUPS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for ind in $(grep -v "#" help_files/SRA_accession_outgroups.txt |cut -f2)
do
  echo $ind
  sbatch -p core -t 10:00:00  $scrdir/run_downloadFastq.sh $ind help_files/SRA_accession_outgroups.txt
done

# ~~~~~~~~~~~~~~~~~~~~~~~~ OUTGROUP SNAKEMAKE PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~
# pipeine that maps the five outgroups onto the dog reference

# I AM HERE!!!!


# Do this with a Snakemake pipeline! (The outgroups are defined in the file)
screen
module load bioinfo-tools snakemake/5.30.1 bwa/0.7.17 samtools/1.10 picard/2.23.4 GATK/3.8-0 BEDTools/2.29.2
# first try a dryrun
snakemake -np --snakefile $scrdir/outgroup.snakefile
# And try plotting a DAG (directed acyclic graph)
snakemake --dag --snakefile $scrdir/outgroup.snakefile | dot -Tsvg > dag/dag_genotype_outgroups.svg
# Then start the actual pipeline
snakemake --snakefile $scrdir/snake/outgroup.snakefile -p -j 64  --cluster "sbatch -A p2018002 --qos=p2018002_8nodes -p {cluster.partition} -n {cluster.n} -t {cluster.time} -e {cluster.error} -o {cluster.output} --mail-type {cluster.mail-type} --mail-user {cluster.mail-user}" --cluster-config $scrdir/snake/cluster_outgroup.json
# removed this when restarting other steps than the first -C {cluster.C}
# added this for etracting sites step --qos=p2018002_8nodes


################################## GENOTYPES ###################################
# 2021-04-19
# 2022-01-25 Update: rerun for all wolf SNP-sites, not only vep
# 2022-03-23 Rerun on new sets without duplicate ind V138

# I realised I haven't removed indels from the vcf! This can create multiple
# lines for the same site in the script below. Instead of using GATK on the
# original file, I filter on variant length in the extracted vcf files.
for set in "74Females.CombFilter.chrX" "74Females.KeepSingleCombFilter.chrX" #"74Females.CombFilter.vepfinal.chrX" "74Females.KeepSingleCombFilter.vepfinal.chrX" #"100S95F14R.CombFilter.chr1-38" "100S95F14R.KeepSingleCombFilter.chr1-38" # "100S95F14R.CombFilter.vepfinalchr1-38" "100S95F14R.KeepSingleCombFilter.vepfinalchr1-38" #"BiDepMac2.PASS"
do
 for i in   "AfGWo" "BlBJa" #"Dhole" "SiSJa" "EthWo"
 do
  #awk '{if(/^#/){print}else{if(length($4)==1 && length($5)==1 && $5!="*"){print}}}' outgroup/$i.$filt.vep.vcf >outgroup/$i.$filt.vep.SNPs.chr1-X.vcf
  #awk '{if(/^#/){print}else{if(length($4)==1 && length($5)==1 && $5!="*"){print}}}' outgroup/$i.$filt.allSNPs.vcf >outgroup/$i.$filt.allSNPs.SNPs.chr1-38.vcf
  awk '{if(/^#/){print}else{if(length($4)==1 && length($5)==1 && $5!="*"){print}}}' outgroup/$i.$set.extract.vcf >outgroup/$i.$set.extract.noIndel.vcf
 done
done

# Convert the vcf to bed (require minimum depth=5/GQ=15)
# Rerun 2021-05-18, try DP=3 instead
# rerun 2021-05-24, try pseudohaploidization using allele freq as weight
# 2022-03-23 Rerun on new sets without duplicate ind V138
module load python3
for set in "74Females.CombFilter.chrX" "74Females.KeepSingleCombFilter.chrX" "74Females.CombFilter.vepfinal.chrX" "74Females.KeepSingleCombFilter.vepfinal.chrX" #"100S95F14R.CombFilter" "100S95F14R.KeepSingleCombFilter" #"100S95F14R.CombFilter.vepfinal" "100S95F14R.KeepSingleCombFilter.vepfinal" #"CombFilter" #"BiDepMac2.PASS"
do
 for n in 5 #3
 do
  for i in  "AfGWo" "BlBJa" #"Dhole" "SiSJa" "EthWo"
  do
  # python3 $scrdir/convertOutgroupVCFtoBED.py -v outgroup/$i.$filt.vep.SNPs.chr1-X.vcf -d $n -o outgroup/$i.$filt.vep.SNPs.chr1-X.DP$n
#   python3 $scrdir/convertOutgroupVCFtoBED_useAlleleFreq.py -v outgroup/$i.$filt.vep.SNPs.chr1-X.vcf -d $n -o outgroup/$i.$filt.vep.SNPs.weightedAF.chr1-X.DP$n
   #calls=`grep -v "N" outgroup/$i.$filt.vep.SNPs.weightedAF.chr1-X.DP$n.bed |wc -l`
   python3 $scrdir/convertOutgroupVCFtoBED_useAlleleFreq.py -v outgroup/$i.$set.extract.noIndel.vcf -d $n -o outgroup/$i.$set.weightedAF.DP$n
   calls=`grep -v "N" outgroup/$i.$set.weightedAF.DP$n.bed |wc -l`
   echo $i $n $calls
  done
 done
done

# Some sites might be missing due to indels/"<NON_REF>" in vcf file - replace
# missing with Ns!
module load bioinfo-tools  BEDTools/2.29.2
for set in "74Females.CombFilter.vepfinal.chrX" "74Females.KeepSingleCombFilter.vepfinal.chrX" #"100S95F14R.CombFilter.vepfinal" "100S95F14R.KeepSingleCombFilter.vepfinal" #"CombFilter" #"BiDepMac2.PASS"
do
 for n in 5 #3
 do
  for i in "AfGWo" "BlBJa" #"Dhole" "SiSJa" "EthWo"
  do
  # intersectBed -wao -a vep/100S95F15R.$filt.final.all.chr1-X.bed -b outgroup/$i.$filt.vep.SNPs.chr1-X.DP$n.bed |awk -v OFS="\t" '{if($9==0){$8="N"}; print}' |cut -f 8 >outgroup/$i.$filt.DP$n.tmp.txt
   #intersectBed -wao -a vep/100S95F15R.$filt.final.all.chr1-X.bed -b outgroup/$i.$filt.vep.SNPs.weightedAF.chr1-X.DP$n.bed |awk -v OFS="\t" '{if($9==0){$8="N"}; print}' |cut -f 8 >outgroup/$i.$filt.weightedAF.DP$n.tmp.txt
   intersectBed -wao -a bed/$set.bed -b outgroup/$i.$set.weightedAF.DP$n.bed |awk -v OFS="\t" '{if($9==0){$8="N"}; print}' |cut -f 8 >outgroup/$i.$set.weightedAF.DP$n.addN.txt
  done
 done
done
n=5
for set in "74Females.CombFilter.chrX" "74Females.KeepSingleCombFilter.chrX" #"100S95F14R.CombFilter" "100S95F14R.KeepSingleCombFilter"
do
  for i in "AfGWo" "BlBJa" #"Dhole" "SiSJa" "EthWo"
  do
   intersectBed -wao -a bed/$set.bed -b outgroup/$i.$set.weightedAF.DP$n.bed |awk -v OFS="\t" '{if($8==0){$7="N"}; print}' |cut -f 7 >outgroup/$i.$set.weightedAF.DP$n.addN.txt
  done
done


################ ALL 5 outgroups, not used in the end!
# Combine all outgroups and assign ancestral
# + check stats
module load python3
for filt in "CombFilter" "BiDepMac2.PASS"
do
 for n in 5 3
 do
  #paste  vep/100S95F15R.$filt.final.all.chr1-X.bed outgroup/AfGWo.$filt.DP$n.tmp.txt outgroup/EthWo.$filt.DP$n.tmp.txt outgroup/Dhole.$filt.DP$n.tmp.txt outgroup/SiSJa.$filt.DP$n.tmp.txt outgroup/BlBJa.$filt.DP$n.tmp.txt >outgroup/5outgroups.$filt.vep.SNPs.chr1-X.DP$n.txt
  #python3 $scrdir/assignAncestral.py -i outgroup/5outgroups.$filt.vep.SNPs.chr1-X.DP$n.txt -c 5 -o outgroup/Ancestral.$filt.vep.SNPs.chr1-X.DP$n 2>stderr.$filt.vep.SNPs.chr1-X.DP$n.txt
  #echo "From GT: Filt $filt and DP $n:"
  #cut -f5  outgroup/Ancestral.$filt.vep.SNPs.chr1-X.DP$n.ancestral.txt |grep -v "N" |wc -l
  #awk '($6=="1.00"){print}' outgroup/Ancestral.$filt.vep.SNPs.chr1-X.DP$n.ancestral.txt |wc -l
  #awk '($6=="1.00" && $7>=3){print}' outgroup/Ancestral.$filt.vep.SNPs.chr1-X.DP$n.ancestral.txt |wc -l
  #paste vep/100S95F15R.$filt.final.all.chr1-X.bed outgroup/AfGWo.$filt.weightedAF.DP$n.tmp.txt outgroup/EthWo.$filt.weightedAF.DP$n.tmp.txt outgroup/Dhole.$filt.weightedAF.DP$n.tmp.txt outgroup/SiSJa.$filt.weightedAF.DP$n.tmp.txt outgroup/BlBJa.$filt.weightedAF.DP$n.tmp.txt >outgroup/5outgroups.$filt.weightedAF.vep.SNPs.chr1-X.DP$n.txt
  #python3 $scrdir/assignAncestral.py -i outgroup/5outgroups.$filt.weightedAF.vep.SNPs.chr1-X.DP$n.txt -c 5 -o outgroup/Ancestral.$filt.vep.SNPs.weightedAF.chr1-X.DP$n 2>stderr.$filt.vep.SNPs.weightedAF.chr1-X.DP$n.txt
  echo "AF weighted: Filt $filt and DP $n:"
  cut -f5  outgroup/Ancestral.$filt.vep.SNPs.weightedAF.chr1-X.DP$n.ancestral.txt |grep -v "N" |wc -l
  awk '($6=="1.00"){print}' outgroup/Ancestral.$filt.vep.SNPs.weightedAF.chr1-X.DP$n.ancestral.txt |wc -l
  awk '($6=="1.00" && $7>=3){print}' outgroup/Ancestral.$filt.vep.SNPs.weightedAF.chr1-X.DP$n.ancestral.txt |wc -l
 done
done
###########

# Also try using only 2 species (AfGWo and BlBJa, the two best covered outgroups)
for set in "100S95F14R.CombFilter.vepfinal.chr1-38" "100S95F14R.KeepSingleCombFilter.vepfinal.chr1-38" #"74Females.CombFilter.vepfinal.chrX" "74Females.KeepSingleCombFilter.vepfinal.chrX" #" #"CombFilter" #"BiDepMac2.PASS"
do
 for n in 5 #3
 do
  paste  bed/$set.bed outgroup/AfGWo.$set.weightedAF.DP$n.addN.txt outgroup/BlBJa.$set.weightedAF.DP$n.addN.txt >outgroup/2outgroups.$set.DP$n.txt
  echo "start of loop"
  python3 $scrdir/assignAncestral.py -i outgroup/2outgroups.$set.DP$n.txt -c 5 -o outgroup/Ancestral.2outgroups.$set.DP$n 2>stderr.2out.$set.DP$n.txt
  echo "From GT: Filt $set and DP $n:"
  cut -f5  outgroup/Ancestral.2outgroups.$set.DP$n.ancestral.txt |grep -v "N" |wc -l
  awk '($6=="1.00"){print}' outgroup/Ancestral.2outgroups.$set.DP$n.ancestral.txt |wc -l
  awk '($6=="1.00" && $7>=2){print}' outgroup/Ancestral.2outgroups.$set.DP$n.ancestral.txt |wc -l
 done
done
# Need to do full SNP set separately due to less columns in bed file
n=5
for set in "74Females.CombFilter.chrX" "74Females.KeepSingleCombFilter.chrX" #"100S95F14R.CombFilter.chr1-38" "100S95F14R.KeepSingleCombFilter.chr1-38"
do
  paste  bed/$set.bed outgroup/AfGWo.$set.weightedAF.DP$n.addN.txt outgroup/BlBJa.$set.weightedAF.DP$n.addN.txt >outgroup/2outgroups.$set.DP$n.txt
  python3 $scrdir/assignAncestral.py -i outgroup/2outgroups.$set.DP$n.txt -c 4 -o outgroup/Ancestral.2outgroups.$set.DP$n 2>stderr.2out.$set.DP$n.txt
  echo "From GT: Filt $set and DP $n:"
  cut -f5  outgroup/Ancestral.2outgroups.$set.DP$n.ancestral.txt |grep -v "N" |wc -l
  awk '($5=="1.00"){print}' outgroup/Ancestral.2outgroups.$set.DP$n.ancestral.txt |wc -l
  awk '($5=="1.00" && $6>=2){print}' outgroup/Ancestral.2outgroups.$set.DP$n.ancestral.txt |wc -l
done



# FINAL BED FILES WITH ANCESTRAL STATES
# AUTOSOMES
awk -v OFS="\t" '($5!="N" && $6=="1.00" && $7==2){print $1,$2,$3,$5}' outgroup/Ancestral.2outgroups.100S95F14R.CombFilter.vepfinal.chr1-38.DP5.ancestral.txt >outgroup/Anc.2out.comb.vepfinal.chr1-38.bed
awk -v OFS="\t" '($5!="N" && $6=="1.00" && $7>=2){print $1,$2,$3,$5}' outgroup/Ancestral.2outgroups.100S95F14R.KeepSingleCombFilter.vepfinal.chr1-38.DP5.ancestral.txt >outgroup/Anc.2out.sicomb.vepfinal.chr1-38.bed
awk -v OFS="\t" '($4!="N" && $5=="1.00" && $6==2){print $1,$2,$3,$4}' outgroup/Ancestral.2outgroups.100S95F14R.CombFilter.chr1-38.DP5.ancestral.txt >outgroup/Anc.2out.comb.chr1-38.bed
awk -v OFS="\t" '($4!="N" && $5=="1.00" && $6>=2){print $1,$2,$3,$4}' outgroup/Ancestral.2outgroups.100S95F14R.KeepSingleCombFilter.chr1-38.DP5.ancestral.txt >outgroup/Anc.2out.sicomb.chr1-38.bed
# CHR X
awk -v OFS="\t" '($5!="N" && $6=="1.00" && $7==2){print $1,$2,$3,$5}' outgroup/Ancestral.2outgroups.74Females.CombFilter.vepfinal.chrX.DP5.ancestral.txt >outgroup/Anc.2out.comb.vepfinal.chrX.bed
awk -v OFS="\t" '($5!="N" && $6=="1.00" && $7>=2){print $1,$2,$3,$5}' outgroup/Ancestral.2outgroups.74Females.KeepSingleCombFilter.vepfinal.chrX.DP5.ancestral.txt >outgroup/Anc.2out.sicomb.vepfinal.chrX.bed
awk -v OFS="\t" '($4!="N" && $5=="1.00" && $6==2){print $1,$2,$3,$4}' outgroup/Ancestral.2outgroups.74Females.CombFilter.chrX.DP5.ancestral.txt >outgroup/Anc.2out.comb.chrX.bed
awk -v OFS="\t" '($4!="N" && $5=="1.00" && $6>=2){print $1,$2,$3,$4}' outgroup/Ancestral.2outgroups.74Females.KeepSingleCombFilter.chrX.DP5.ancestral.txt >outgroup/Anc.2out.sicomb.chrX.bed








##################################### VEP ######################################
# 2022-03-23
# VEP was run as part of the Snakemake pipeline above. Below follow extractions
# of the different categories.

for prefix in "74Females.chrX" #"100S95F14R"
do
 for filt in "CombFilter" "KeepSingleCombFilter"
 do
  grep -v "#" vep/$prefix.SNPs.HF.$filt.txt |cut -f7 |sort |uniq -c >vep/$prefix.$filt.vep_types.txt
 done
done

# Extract sites
for prefix in "74Females.chrX" #"100S95F14R"
do
 for filt in "CombFilter" "KeepSingleCombFilter"
 do
  grep "synonymous" vep/$prefix.SNPs.HF.$filt.txt |awk '($7=="synonymous_variant"){split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' |uniq >vep/$prefix.$filt.firstExtract.synonymous.bed
  grep "missense" vep/$prefix.SNPs.HF.$filt.txt |awk '($7=="missense_variant"){split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' |uniq >vep/$prefix.$filt.firstExtract.missense.bed
  grep "missense" vep/$prefix.SNPs.HF.$filt.txt |grep "tolerated" |grep -v "low_confidence" | awk '($7=="missense_variant"){split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' |uniq >vep/$prefix.$filt.firstExtract.tolerated.missense.bed
  grep "missense" vep/$prefix.SNPs.HF.$filt.txt |grep "deleterious" |grep -v "low_confidence" | awk '($7=="missense_variant"){split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' |uniq >vep/$prefix.$filt.firstExtract.deleterious.missense.bed
  grep "IMPACT=HIGH" vep/$prefix.SNPs.HF.$filt.txt |awk '{if($7~/stop_gained/){if($11~/*$/){print}}else if($7~/stop_lost/){if($11~/^*/){print}}else if($7~/start_lost/){if($11~/^M/){print}}else{print}}'|awk '{split($2,s,":"); start=s[2]-1; print s[1]"\t"start"\t"s[2]"\t"$1}' |uniq >vep/$prefix.$filt.firstExtract.nonsense.bed
 done
done

# Note: Some sites have two or even more consequenses listed. For synonymous and
# missense above, I keep only SNPs with these consequences (if they are eg. both
#"synonymous" and "splice_region_variant", they are not included in synonymous).
# For nonsense, I save sites with multiple annotations. However, I only save
# start_lost/stop_lost/stop_gained if they make sense (I check that they lost an
# "M" or "*" or gained a "*")

# Warning - noticed that the files above still overlap slightly! A snp can be
# synonymous in one gene and nonsynonymous in an overlapping gene (listed
# multiple times in the original file - not with multiple consequenses in one
# line)

# Keep only the most severe consequence:
module load bioinfo-tools BEDTools/2.29.2
for prefix in "74Females.chrX" #"100S95F14R"
do
 for filt in "CombFilter" "KeepSingleCombFilter"
 do
  intersectBed -v -a vep/$prefix.$filt.firstExtract.synonymous.bed -b vep/$prefix.$filt.firstExtract.missense.bed vep/$prefix.$filt.firstExtract.nonsense.bed >vep/$prefix.$filt.final.synonymous.bed
  intersectBed -v -a vep/$prefix.$filt.firstExtract.missense.bed -b vep/$prefix.$filt.firstExtract.nonsense.bed >vep/$prefix.$filt.final.missense.bed
  intersectBed -v -a vep/$prefix.$filt.firstExtract.deleterious.missense.bed -b vep/$prefix.$filt.firstExtract.nonsense.bed >vep/$prefix.$filt.final.deleterious.missense.bed
  intersectBed -v -a vep/$prefix.$filt.firstExtract.tolerated.missense.bed -b vep/$prefix.$filt.firstExtract.nonsense.bed vep/$prefix.$filt.final.deleterious.missense.bed >vep/$prefix.$filt.final.tolerated.missense.bed
  cp vep/$prefix.$filt.firstExtract.nonsense.bed vep/$prefix.$filt.final.nonsense.bed
 done
done
# Now there are no duplicates! Checked with
for prefix in "74Females.chrX" #"100S95F14R"
do
 for filt in "CombFilter" "KeepSingleCombFilter"
 do
  cat vep/$prefix.$filt.final.missense.bed vep/$prefix.$filt.final.nonsense.bed vep/$prefix.$filt.final.synonymous.bed |wc -l
  cat vep/$prefix.$filt.final.missense.bed vep/$prefix.$filt.final.nonsense.bed vep/$prefix.$filt.final.synonymous.bed |cut -f4 |sort |uniq |wc -l
  cat vep/$prefix.$filt.final.deleterious.missense.bed vep/$prefix.$filt.final.tolerated.missense.bed vep/$prefix.$filt.final.nonsense.bed vep/$prefix.$filt.final.synonymous.bed  |wc -l
  cat vep/$prefix.$filt.final.deleterious.missense.bed vep/$prefix.$filt.final.tolerated.missense.bed vep/$prefix.$filt.final.nonsense.bed vep/$prefix.$filt.final.synonymous.bed |cut -f4 |sort |uniq |wc -l
 done
done

# From the whole genome file, extract only autosomes:
# (Put final file in subdir "bed" to make it easier for outgroup snake pipe!)
prefix="100S95F14R"
for filt in "CombFilter"  "KeepSingleCombFilter"
do
 cat vep/$prefix.$filt.final.missense.bed vep/$prefix.$filt.final.nonsense.bed vep/$prefix.$filt.final.synonymous.bed |awk '(/^[0-9]/){print}' |sort -k1,1 -k2,2n >bed/$prefix.$filt.vepfinal.chr1-38.bed
done
# Make a combined file for Females X chromosome and rename
oldprefix="74Females.chrX"
newprefix="74Females"
for filt in "CombFilter"  "KeepSingleCombFilter"
do
 cat vep/$oldprefix.$filt.final.missense.bed vep/$oldprefix.$filt.final.nonsense.bed vep/$oldprefix.$filt.final.synonymous.bed |grep "^X" |sort -k1,1 -k2,2n  > bed/$newprefix.$filt.vepfinal.chrX.bed
done
# Create new vcf files based on these sites only
# Autosomes
prefix="100S95F14R"
for filt in "CombFilter" "KeepSingleCombFilter" # "BiDepMac2.PASS"
do
 echo '#!/bin/bash -l
 module load bioinfo-tools BEDTools/2.29.2
 intersectBed -header -a vcf/'$prefix'.SNPs.HF.'$filt'.vcf.gz -b bed/'$prefix'.'$filt'.vepfinal.chr1-38.bed >vcf/'$prefix'.'$filt'.vepfinal.chr1-38.vcf
 ' |sbatch -J $set -A p2018002 --qos=p2018002_8nodes -t 5:00:00
done
# ChrX (only females)
prefix="74Females"
for filt in "CombFilter" "KeepSingleCombFilter" # "BiDepMac2.PASS"
do
 echo '#!/bin/bash -l
 module load bioinfo-tools BEDTools/2.29.2
 intersectBed -header -a vcf/'$prefix'.chrX.SNPs.HF.'$filt'.vcf.gz -b bed/'$prefix'.'$filt'.vepfinal.chrX.bed >vcf/'$prefix'.'$filt'.vepfinal.chrX.vcf
 ' |sbatch -J $set -A p2018002 --qos=p2018002_8nodes -t 5:00:00
done
