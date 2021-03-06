def get_reference(wildcards):
    if wildcards.sample == "adt":
        return config["antibodies"]["adt-seqs"]
    elif wildcards.sample == "hto":
        return config["antibodies"]["hto-seqs"]
    elif wildcards.sample == "rna":
        return "resources/reference/rna/transcriptome.fasta"


def get_reads(wildcards):
    sample = config["samples"][wildcards.sample]
    print(sample)
    if "fq1" in sample and "fq2" in sample:
        return {"fq1": sample["fq1"], "fq2": sample["fq2"]}
    else:
        accession = sample["sra"]
        return {
            "fq1": f"resources/sra/{accession}_1.fastq",
            "fq2": f"resources/sra/{accession}_2.fastq",
        }


def get_sample(wildcards):
    return config["samples"][wildcards.sample]


def get_targets():
    if config["thresholds"].get("max-hto-count") is None:
        return "results/plots/hto-counts.initial.pdf"
    if config["thresholds"].get("max-adt-count") is None:
        return "results/plots/adt-counts.pdf"
    # TODO go on with later steps
    return []
