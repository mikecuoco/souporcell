#!/usr/bin/env bash
# Author: Mike Cuoco
# Created on: 6/11/22, 1:39 PM
#
# Description: compute relatedness of clusters and genotypes from souporcell run
# Usage: bash annotate_clusters.sh <souporcell output directory>

# exit if any non-zero, exit if undefined var
set -uo pipefail

inDir=$1
threads=$2
clusters_vcf=$inDir/cluster_genotypes.vcf
clusters_vcf_fixnames=$inDir/cluster_genotypes_fixnames.vcf 
genotypes_vcf=$inDir/common_variants_covered.vcf
samples_file=$inDir/cluster_names.txt

echo "Comparing $clusters_vcf to $genotypes_vcf"

# fix cluster names
clusters=$(bcftools query -l $clusters_vcf)
touch $samples_file
for c in $clusters; do
    echo "cluster$c" >> $samples_file
done

bcftools reheader -s $samples_file -o $clusters_vcf_fixnames $clusters_vcf

# Merge VCF files
echo "Merging"
bgzip $clusters_vcf_fixnames && tabix -f $clusters_vcf_fixnames.gz
bgzip $genotypes_vcf && tabix -f $genotypes_vcf.gz
bcftools merge -f PASS --threads $threads -Oz --output $inDir/merged.vcf.gz \
    -m none \
    $clusters_vcf_fixnames.gz $genotypes_vcf.gz

# compute relatedness
echo "Computing relatedness"
plink \
    --vcf $inDir/merged.vcf.gz \
    --genome \
    -allow-extra-chr \
    --out $inDir/relatedness

# take matches > 95% related
grep cluster $inDir/relatedness.genome | awk '{if ($12 > 0.90 && $12 != "nan") print $1,$3,$12,$15,$16,$17}'

# cleanup
echo "Cleaning up"
bgzip -d $genotypes_vcf.gz
rm -f $samples_file $clusters_vcf_fixnames.gz $inDir/merged.vcf.gz $inDir/*tbi