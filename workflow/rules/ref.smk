rule get_genome:
    output:
        "resources/reference/genome.fasta",
    params:
        species=config["reference"]["species"],
        datatype="dna",
        build=config["reference"]["build"],
        release=config["reference"]["release"],
    log:
        "logs/get-genome.log",
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "0.74.0/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/reference/annotation.gtf",
    params:
        species=config["reference"]["species"],
        build=config["reference"]["build"],
        release=config["reference"]["release"],
        fmt="gtf",
        flavor="",
    log:
        "logs/get-annotation.log",
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "0.74.0/bio/reference/ensembl-annotation"


rule get_geneid2name:
    output:
        "resources/reference/geneid2name.tsv",
    params:
        species=config["reference"]["species"],
    conda:
        "../envs/pybiomart.yaml"
    script:
        "../scripts/get-geneid2name.py"


rule build_splici_transcriptome:
    input:
        fasta="resources/reference/genome.fasta",
        gtf="resources/reference/annotation.gtf",
    output:
        seq="resources/reference/rna/transcriptome.fasta",
        t2g="resources/reference/rna/t2g.tsv",
    log:
        "logs/build-splici-transcriptome.log",
    conda:
        "../envs/r-splici-transcriptome.yaml"
    notebook:
        "../notebooks/build-splici-transcriptome.r.ipynb"


rule spoof_t2g:
    input:
        get_reference,
    output:
        "resources/reference/{type}/t2g.tsv",
    log:
        "logs/spoof-t2g/{type}.log",
    wildcard_constraints:
        type="adt|hto",
    conda:
        "../envs/unix.yaml"
    shell:
        "awk '{{print $1\"\\t\"$1;}}' {input} > {output} 2> {log}"
