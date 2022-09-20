CHROM = [i for i in range(1,22+1)]

# The following rules will not be executed as jobs
localrules: all, maffilter_options

rule all:
    input:
        expand("gerp/canFam3/hg38.chr{prefix}.rates.bed", prefix=CHROM)

# Download data from UCSC
rule download_data:
    output:
        "maf/original/chr{prefix}.maf.gz"
    log:
        "logs/maf/download.chr{prefix}.log"
    threads: 1
    shell:
        "wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/maf/chr{wildcards.prefix}.maf.gz -O {output}"

rule maffilter_options:
    input:
        inmaf="maf/original/chr{prefix}.maf.gz",
        tree="maf/hg38.tree"
    output:
        "maf/optionfiles/chr{prefix}.options",
    log:
        "logs/maf/options.chr{prefix}.log"
    params:
        o="maf/noCanFam/chr{prefix}.maf",
        olog="maf/noCanFam/logs/chr{prefix}.maffilter.log"
    threads: 1
    shell:
        "make_maffilter_options.sh {input.tree} {input.inmaf} {output} {params.olog} {params.o}"

rule maffilter:
    input:
        inmaf="maf/original/chr{prefix}.maf.gz",
        inopt="maf/optionfiles/chr{prefix}.options"
    output:
        "maf/noCanFam/chr{prefix}.maf"
    log:
        "logs/maf/filter.chr{prefix}.log"
    threads: 1
    shell:
        "maffilter param={input.inopt}"

rule gerp:
    input:
        inmaf="maf/noCanFam/chr{prefix}.maf",
        tree="maf/hg38.tree"
    output:
        "maf/noCanFam/chr{prefix}.maf.rates"
    log:
        "logs/gerp/chr{prefix}.log"
    params:
        ref="hg38"
    threads: 1
    shell:
        "/sw/bioinfo/GERP++/20110522/rackham/gerpcol -t {input.tree} -f {input.inmaf} -z -e {params.ref}"

rule rates2bed:
    input:
        "maf/noCanFam/chr{prefix}.maf.rates"
    output:
        "gerp/hg38/chr{prefix}.rates.bed"
    log:
        "logs/gerp/toBed.chr{prefix}.log"
    threads: 1
    shell:
        "awk -v p='chr'{wildcards.prefix} -v OFS='\t' '{{s=NR-1; print p,s,NR,$1,$2}}' {input} >{output}"

rule liftover:
    input:
        rates="gerp/hg38/chr{prefix}.rates.bed",
        lift="liftOver/hg38ToCanFam3.over.chain"
    output:
        o1="gerp/canFam3/hg38.chr{prefix}.rates.bed",
        o2="gerp/hg38/unmapped/chr{prefix}.rates.bed"
    log:
        "logs/gerp/liftover.chr{prefix}.log"
    threads: 1
    shell:
        "~/software/UCSC/liftOver {input.rates} {input.lift} {output.o1} {output.o2}"
