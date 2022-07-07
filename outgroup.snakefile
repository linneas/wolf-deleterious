SAMPLES = ["AfGWo", "EthWo", "Dhole", "BlBJa", "SiSJa"]
FILTER = ["mac2", "mac1"]

# The following rules will not be executed as jobs
localrules: all, makeGVCFList

rule all:
    input:
        expand("outgroup/{prefix}.100S95F14R.chr1-38.{prefix2}.extract.vcf", prefix=SAMPLES, prefix2=FILTER),
        expand("outgroup/{prefix}.74Females.chrX.{prefix2}.extract.vcf", prefix=SAMPLES, prefix2=FILTER)

rule bwa:
    input:
        ref="reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa",
        r1="fastq/{prefix}_1.fastq.gz",
        r2="fastq/{prefix}_2.fastq.gz"
    output:
        "bam/{prefix}.bam"
    log:
        "logs/bwa/{prefix}.log"
    params:
        rg=r"@RG\tID:{prefix}\tSM:{prefix}\tLB:{prefix}\tPL:Illumina"
    threads: 20
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} -M {input.ref} {input.r1} {input.r2} \
         |samtools sort -m 6G -@{threads} -T $SNIC_TMP/{wildcards.prefix} - >{output}.tmp) 2>{log} \
         && mv {output}.tmp {output}"

rule index_bam:
    input:
        "bam/{prefix}.bam"
    output:
        "bam/{prefix}.bam.bai"
    log:
        "logs/samtools/index.{prefix}.log"
    threads: 1
    shell:
        "samtools index {input} {output}"

rule mkdupl:
    input:
        "bam/{prefix}.bam"
    output:
        bam="bam/{prefix}.md.bam",
        metrics="bam/{prefix}.metrics"
    log:
        "logs/picard/{prefix}.mkdupl.log"
    params:
        picard="/sw/apps/bioinfo/picard/2.23.4/rackham/picard.jar",
        mem="60g"
    threads: 10
    shell:
        "(java -Xmx{params.mem} -jar {params.picard} MarkDuplicates INPUT={input} \
         METRICS_FILE={output.metrics} TMP_DIR=$SNIC_TMP ASSUME_SORTED=true \
         VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE OUTPUT={output.bam}.tmp) 2>{log} \
         && mv {output.bam}.tmp {output.bam} && mv {output.bam}.tmp.bai {output.bam}.bai"

rule haplotype_caller:
    input:
        bam="bam/{prefix}.md.bam",
        ref="reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa"
    output:
        "gvcf/{prefix}.g.vcf.gz"
    log:
        "logs/GATK/HaplotypeCaller.{prefix}.log"
    params:
        gatk="/sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar",
        mem="120g"
    threads: 20
    shell:
        "(java -Xmx{params.mem} -Djava.io.tmpdir=$SNIC_TMP -jar {params.gatk} \
        -T HaplotypeCaller -R {input.ref} -I {input.bam} --emitRefConfidence GVCF \
        --variant_index_type LINEAR --variant_index_parameter 128000 -nct {threads} \
        -jdk_deflater -jdk_inflater -o gvcf/{wildcards.prefix}.g.vcf.tmp.gz) 2>{log} \
        && mv gvcf/{wildcards.prefix}.g.vcf.tmp.gz {output} \
        && mv gvcf/{wildcards.prefix}.g.vcf.tmp.gz.tbi {output}.tbi"

rule makeGVCFList:
    input:
        "gvcf/{prefix}.g.vcf.gz"
    output:
        "{prefix}.gvcf.list"
    log:
        "logs/misc/create.GVCFList.{prefix}.log"
    shell:
        "echo {input} >{output}"

rule AllsitesGenotypeGVCFs:
    input:
        fa="reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa",
        list="{prefix}.gvcf.list"
    output:
        "vcf/{prefix}.vcf"
    params:
        gatk="/sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar",
        mem="120g"
    log:
        "logs/GATK/GenotypeGVCFs.{prefix}.log"
    threads: 20
    shell:
        "java -Xmx{params.mem} -Djava.io.tmpdir=$SNIC_TMP -jar {params.gatk} \
        -T GenotypeGVCFs -R {input.fa} --variant {input.list} -nt {threads} \
        -allSites -o {output}.tmp && mv {output}.tmp {output}"

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

rule extractBedSites:
    input:
        vcf="vcf/{species, \w+}.vcf.gz",
        index="vcf/{species, \w+}.vcf.gz.tbi",
        bed="bed/{filter}.bed"
    output:
        "outgroup/{species, \w+}.{filter}.extract.vcf"
    log:
        "logs/vcf/{species, \w+}.{filter}.extract.log"
    threads: 1
    shell:
        "zcat {input.vcf} |grep -v '<NON_REF>' | \
        intersectBed -header -a - -b {input.bed} >{output}"
