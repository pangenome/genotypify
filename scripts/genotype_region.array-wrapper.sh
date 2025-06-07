#!/bin/bash
#SBATCH -c 48

dir_base=$1

# Read the region info for this array task
region_line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $dir_base/genotyping/regions.txt)
chrom=$(echo "$region_line" | cut -f1)
start=$(echo "$region_line" | cut -f2)
end=$(echo "$region_line" | cut -f3)
name=$(echo "$region_line" | cut -f4)
len=$(echo "$region_line" | cut -f5)

echo "Array task $SLURM_ARRAY_TASK_ID: Processing region: $chrom:$start-$end ($name, len: $len)"

region=${chrom}_${start}_${end}
echo -e "$chrom\t$start\t$end\t$name" > $dir_base/genotyping/$region.bed

hostname
$dir_base/scripts/genotype_region.sh \
    $dir_base/genotyping/$region.bed \
    $dir_base/wfmash \
    $dir_base/reference/GRCh38.fa.gz \
    $dir_base/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    $dir_base/genotyping/samples-to-consider.txt \
    /lizardfs/guarracino/pangenomes/HPRCv2 \
    $dir_base/data/HPRCv2/illumina \
    48 \
    /scratch \
    $dir_base/genotyping/$region