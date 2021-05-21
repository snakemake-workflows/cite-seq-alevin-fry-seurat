rule seurat:
    input:
        rna="results/alevin-fry/quant/rna",
        adt="results/alevin-fry/quant/adt",
        hto="results/alevin-fry/quant/hto",
        geneid2name="resources/reference/geneid2name.tsv"
    output:
        "results/seurat/{sample}.all.rds"
    log:
        "logs/seurat/{sample}.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat.r"


rule plot_initial_hto_counts:
    input:
        "results/seurat/all.rds"
    output:
        "results/plots/hto-counts.initial.pdf"
    log:
        "logs/plot-hto-counts-initial.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/plot-hto-counts.r"


rule filter_normalize_demux:
    input:
        "results/seurat/all.rds"
    output:
        "results/seurat/hto-filtered.rds"
    params:
        max_hto_count=config["thresholds"]["max-hto-count"]
    log:
        "logs/filter-normalize-demux.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/filter-normalize-demux.r"


rule plot_counts_hto_filtered:
    input:
        "results/seurat/hto-filtered.rds"
    output:
        hto="results/plots/hto-counts.pdf",
        rna="results/plots/rna-counts.pdf",
    log:
        "logs/plot-counts-hto-filtered.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/plot-counts-hto-filtered.r"


rule filter_negatives:
    input:
        "results/seurat/hto-filtered.rds"
    output:
        "results/seurat/non-negatives.rds"
    log:
        "logs/fitler-negatives.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/filter-negatives.r"


rule plot_umap_singlets_doublets:
    input:
        "results/seurat/non-negatives.rds"
    output:
        by_xlet="results/plots/umap-singlets-doublets.col~xlet.pdf",
        by_hashtag="results/plots/umap-singlets-doublets.col~hashtag.pdf",
    log:
        "logs/plot-umap-singlets-doublets.log"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/plot-umap-singlets-doublets.r"


rule filter_to_singlets:
    input:
        "results/seurat/non-negatives.rds"
    output:
        rds="results/seurat/singlets.rds",
        pdf="results/plots/adt-counts.pdf"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/filter-to-singlets.r"