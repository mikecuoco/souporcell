#!/usr/bin/env bash
# Author: Mike Cuoco
# Created on: 6/11/22, 1:39 PM
#
# Description: compute relatedness of clusters and genotypes from souporcell run
# Usage: bash annotate_clusters.sh <souporcell output directory> <num_threads>

# exit if any non-zero, exit if undefined var
set -uo pipefail

inDir=$1
cd $inDir
threads=$2

clusters_vcf=cluster_genotypes.vcf
clusters_vcf_fixnames=cluster_genotypes_fixnames.vcf 
genotypes_vcf=common_variants_covered.vcf
samples_file=cluster_names.txt

echo "Comparing $clusters_vcf to $genotypes_vcf"
{
    # fix cluster names
    clusters=$(bcftools query -l $clusters_vcf)
    rm -f $samples_file
    for c in $clusters; do
        echo "cluster$c" >> $samples_file
    done

    bcftools reheader -s $samples_file -o $clusters_vcf_fixnames $clusters_vcf

    # Merge VCF files
    echo "Merging $clusters_vcf with $genotypes_vcf"
    bgzip -f $clusters_vcf_fixnames && tabix -f $clusters_vcf_fixnames.gz
    bgzip -f $genotypes_vcf && tabix -f $genotypes_vcf.gz
    bcftools merge -f PASS --threads $threads -Oz --output $inDir/merged.vcf.gz \
        -m none \
        $clusters_vcf_fixnames.gz $genotypes_vcf.gz

    # compute relatedness
    echo "Computing all pairwise relatedness"
    plink \
        --vcf merged.vcf.gz \
        --genome \
        --threads $threads \
        --silent \
        -allow-extra-chr \
        --out relatedness

    # save cluster matches
    grep cluster relatedness.genome | awk '{ OFS = ","; print $1,$3,$12,$15,$16,$17}' > relatedness.csv
    # print cluster matches > 90% related
    grep cluster relatedness.genome | awk '{if ($12 > 0.90 && $12 != "nan") print $1,$3,$12,$15,$16,$17}'

    # cleanup
    echo "Cleaning up"
    bgzip -d $genotypes_vcf.gz
    rm -f $samples_file $clusters_vcf_fixnames.gz merged.vcf.gz *tbi
} 1> relatedness.out 2> relatedness.err

cd -