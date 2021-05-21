rule get_sra:
    output:
        "resources/sra/{accession}_1.fastq",
        "resources/sra/{accession}_2.fastq",
    log:
        "logs/get-sra/{accession}.log",
    wrapper:
        "0.74.0/bio/sra-tools/fasterq-dump"