import sys
sys.stderr = open(snakemake.log[0], "w")

from pybiomart import Dataset

prefix, suffix = snakemake.config["reference"]["species"].split("_", 1)
species = prefix[1] + suffix

dataset = Dataset(name=f"{species}_gene_ensembl", host="http://www.ensembl.org")

dataset.query(attributes=["ensembl_gene_id", "external_gene_name"]).to_csv(
    snakemake.output[0], sep="\t"
)
