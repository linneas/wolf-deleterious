
# The following rules will not be executed as jobs
localrules: all

rule all:
    input:
    #    "vep/100S95F14R.SNPs.HF.mac2.txt",
    #    "vep/100S95F14R.SNPs.HF.mac1.txt"
        "vep/74Females.SNPs.HF.mac2.txt",
        "vep/74Females.SNPs.HF.mac1.txt"


rule GenotypeGVCFs:
    input:
        fa="reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa",
        list="{prefix}.gvcf.list"
    output:
        "vcf/{prefix}.vcf"
    params:
        gatk="/sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar",
        mem="1000g"
    log:
        "logs/GATK/GenotypeGVCFs.{prefix}.log"
    threads: 20
    shell:
        "java -Xmx{params.mem} -Djava.io.tmpdir=$SNIC_TMP -jar {params.gatk} \
        -T GenotypeGVCFs -R {input.fa} --variant {input.list} -nt {threads} -o {output}.tmp && mv {output}.tmp {output}"

rule zip:
    input:
        "vcf/{prefix}.vcf"
    output:
        "vcf/{prefix}.vcf.gz"
    log:
        "logs/zip/{prefix}.log"
    threads: 1
    shell:
        "bgzip {input}"

rule tabix:
    input:
        "vcf/{prefix}.vcf.gz"
    output:
        "vcf/{prefix}.vcf.gz.tbi"
    params:
        "-p vcf"
    log:
        "logs/tabix/{prefix}.log"
    threads: 1
    shell:
        "tabix {params} {input}"

rule idepth:
    input:
        index="vcf/{prefix}.vcf.gz.tbi",
        vcf="vcf/{prefix}.vcf.gz"
    output:
        "coverage/{prefix}.vcf.idepth"
    log:
        "logs/coverage/{prefix}.log"
    threads: 1
    shell:
        "vcftools --gzvcf {input.vcf} --stdout --depth >{output}"

rule extractSNPs:
    input:
        fa="reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa",
        index="vcf/{prefix}.vcf.gz.tbi",
        vcf="vcf/{prefix}.vcf.gz"
    output:
        "vcf/{prefix}.SNPs.vcf"
    log:
        "logs/GATK/{prefix}.extractSNPs.log"
    params:
        gatk="/sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar",
        mem="60g"
    threads: 10
    shell:
        "java -Xmx{params.mem} -jar {params.gatk} -T SelectVariants \
        -R {input.fa} -V {input.vcf} -selectType SNP -o {output}"


rule filterSNPs:
    input:
        fa="reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa",
        index="vcf/{prefix}.SNPs.vcf.gz.tbi",
        vcf="vcf/{prefix}.SNPs.vcf.gz"
    output:
        "vcf/{prefix}.SNPs.HF.vcf"
    log:
        "logs/GATK/{prefix}.HardFiltSNPs.log"
    params:
        gatk="/sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar",
        mem="60g",
        filt=r"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
    threads: 10
    shell:
        "java -Xmx{params.mem} -jar {params.gatk} -T VariantFiltration \
        -R {input.fa} -V {input.vcf} --filterExpression '{params.filt}' \
        --filterName 'hard_filt' -o {output}"


rule vcfFilter:
    input:
        index="vcf/{prefix}.vcf.gz.tbi",
        vcf="vcf/{prefix}.vcf.gz",
        dep="coverage/{prefix}.vcf.idepth"
    output:
        "vcf/{prefix}.mac2.vcf"
    log:
        "logs/filter/{prefix}.mac2.log"
    params:
        min=10,
        minal=2,
        maxal=2,
        mac=2,
        maxmiss=0.9,
        GQ=30,
    threads: 1
    shell:
        "m=`awk '(NR>1){{n++; sum+=$3}}END{{mean=sum/n; twice=2*mean; print twice}}' {input.dep}` && "
        "vcftools --gzvcf {input.vcf} --out vcf/{wildcards.prefix} --mac {params.mac}  \
        --max-alleles {params.maxal}  --min-alleles {params.minal}  --min-meanDP {params.min} \
        --max-meanDP $m  --minGQ {params.GQ} --max-missing {params.maxmiss} --remove-filtered-all \
        --recode --recode-INFO-all && mv vcf/{wildcards.prefix}.recode.vcf {output}"

rule keepSingletonFilter:
    input:
        index="vcf/{prefix}.vcf.gz.tbi",
        vcf="vcf/{prefix}.vcf.gz",
        dep="coverage/{prefix}.vcf.idepth"
    output:
        "vcf/{prefix}.mac1.vcf"
    log:
        "logs/filter/{prefix}.mac1.log"
    params:
        min=10,
        minal=2,
        maxal=2,
        maxmiss=0.9,
        GQ=30,
    threads: 1
    shell:
        "m=`awk '(NR>1){{n++; sum+=$3}}END{{mean=sum/n; twice=2*mean; print twice}}' {input.dep}` && "
        "vcftools --gzvcf {input.vcf} --out vcf/{wildcards.prefix}  \
        --max-alleles {params.maxal}  --min-alleles {params.minal}  --min-meanDP {params.min} \
        --max-meanDP $m  --minGQ {params.GQ} --max-missing {params.maxmiss} --remove-filtered-all \
        --recode --recode-INFO-all && mv vcf/{wildcards.prefix}.recode.vcf {output}"

rule vep:
    input:
        vcf="vcf/{prefix}.vcf.gz",
        index="vcf/{prefix}.vcf.gz.tbi",
    output:
        "vep/{prefix}.txt"
    log:
        "logs/vep/{prefix}.log"
    params:
        sp="canis_familiaris"
    threads: 4
    shell:
        "vep --cache --dir $VEP_CACHE -i {input.vcf} -o {output} \
        --species {params.sp} --fork {threads} --force_overwrite --sift b"
