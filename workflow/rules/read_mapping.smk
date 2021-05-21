rule salmon_index:
    input:
        get_reference
    output:
        directory("resources/salmon-index/{sample}")
    log:
        "logs/salmon-index/{ref}.log"
    threads: 2
    params:
        extra="--features -k7"
    wrapper:
        "0.74.0/bio/salmon/index"


rule salmon_alevin:
    input:
        unpack(get_reads),
        index=get_salmon_index,
    output:
        directory("results/salmon-alevin/{sample}")
    params:
        unit=get_unit
    log:
        "logs/salmon-alevin/{sample}-{type}.log"
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon alevin "
        "-l ISR -i {input.index} -1 {input.fq1} -2 {input.fq2} "
        "--read-geometry {params.unit[read_geometry]} --bc-geometry {params.unit[barcode_geometry]} "
        "--umi-geometry {params.unit[umi_geometry]} -o {output} --sketch -p 16 "
        "2> {log}"


rule alevin_fry_preprocess:
    input:
        "results/salmon-alevin/{sample}"
    output:
        directory("results/alevin-fry/rad/{sample}")
    log:
        "logs/alevin-fry/rad/{sample}.log"
    conda:
        "../envs/alevin-fry.yaml"
    threads: 16
    shell:
        "(alevin-fry generate-permit-list -d fw -i {input} -o {output} -k &&"
        " alevin-fry collate -r {input} -i {output} -t {threads})"
        " 2> {log}"


rule alevin_fry_quant:
    input:
        rad="results/alevin-fry/rad/{sample}",
        t2g="resources/reference/{sample}/t2g.tsv"
    output:
        directory("results/alevin-fry/quant/{sample}")
    log:
        "logs/alevin-fry/quant/{sample}.log"
    conda:
        "../envs/alevin-fry.yaml"
    threads: 16
    shell:
        "alevin-fry quant -m {input.t2g} -i {input.rad} "
        "-o {output} -r cr-like -t {threads} --use-mtx 2> {log}"