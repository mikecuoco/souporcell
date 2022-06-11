# Mike's fork of [wheaton5/souporcell](https://github.com/wheaton5/souporcell)

<img src="https://github.com/wheaton5/souporcell/blob/master/souporcell_star.png" width="100">

souporcell is a method for clustering mixed-genotype scRNAseq experiments by individual genotypes.

See [wheaton5/souporcell](https://github.com/wheaton5/souporcell) and his [publication](https://www.nature.com/articles/s41592-020-0820-1)
for complete details.

## Install

```bash
git clone https://github.com/wheaton5/souporcell.git

# create conda environment
conda env create -f /path/to/souporcell/souporcell_env.yaml
conda activate souporcell

# install rust and build
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
cd /path/to/souporcell/souporcell && cargo build --release
cd /path/to/souporcell/troublet && cargo build --release
```

## My Usage

This fork uses `relatedness.sh` to calculate the relatedness between cluster genotypes and known genotypes. I have added this step to `souporcell_pipeline.py`. 

Before running `souporcell`, I process the 10x scRNA-seq and genotyping data in the following way:

1. `cellranger count` on the 10x scRNA-seq data.
2. [CellBender](https://github.com/broadinstitute/CellBender) on the `cellranger count` output
3. `minimac4` imputation on the genotyping data, with the 1000 genomes reference panel, done using the [Michigan imputation server](https://imputationserver.sph.umich.edu/index.html#!)

```bash
conda activate souporcell

mkdir -p souporcell && cd souporcell

# setup variables
PIPELINE="/path/to/souporcell/souporcell_pipeline.py"

REF="/path/to/cellranger/fasta/genome.fa"
BAM="/path/to/cellranger_count_dir/outs/possorted_genome_bam.bam"
BARCODES="/path/to/cellbender_out_dir/output_cell_barcodes.csv" 
VCF="/path/to/snp_imputation_output/variants.vcf"

n_inds=\#
n_threads=\# 

$PIPELINE \
--bam $BAM \
--barcodes $BARCODES \
--fasta $FASTA \
--threads $n_threads \
--min_alt 5 \
--min_ref 5 \
--common_variants $VCF \
--out_dir . \
--clusters $n_lines \
--aligner HISAT2    
```
