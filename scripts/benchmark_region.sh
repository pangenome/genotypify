#!/bin/bash
# process_region.sh - Script to process a single region

chrom=$1
start=$2
end=$3

echo "Processing region: $chrom:$start-$end"

region=${chrom}_${start}_${end}
mkdir -p benchmark/$chrom/$region

# Step 1: Make TPR table
Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/calc_tpr.r \
    odgi/dissimilarity/$chrom/$region.tsv.gz \
    clusters/$chrom/$region.clusters.json \
    clusters/$chrom/$region.clusters.hapdist.tsv \
    benchmark/$chrom/$region/tpr.tsv \
    cosigt/*/$chrom/$region/sorted_combos.tsv.gz

# # Step 2: Flip PGGB graph
# odgi flip \
#     -i pggb/$chrom/$region/$region.final.og \
#     -o benchmark/$chrom/$region/$region.flip.og -P \
#     --ref-flips <(grep '^P' odgi/view/$chrom/$region.gfa | cut -f 2 | grep "GRCh38#0")

# Step 3: Convert OG to FASTA
odgi paths \
    -i pggb/$chrom/$region/$region.final.og \
    -f | sed 's/_inv$//g' > benchmark/$chrom/$region/$region.fasta

# Step 4: Prepare combinations for QV
mkdir -p benchmark/$chrom/$region/qv_prep
bash /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/prepare_qv.sh \
    benchmark/$chrom/$region/tpr.tsv \
    benchmark/$chrom/$region/$region.fasta \
    benchmark/$chrom/$region/qv_prep

# Step 5: Calculate QV for each sample
for sample_dir in benchmark/$chrom/$region/qv_prep/*/ ; do
    sample=$(basename "$sample_dir")
    echo "Calculating QV for sample $sample"
    bash /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/calculate_qv.sh \
        benchmark/$chrom/$region/qv_prep/$sample \
        benchmark/$chrom/$region/qv_prep/$sample/qv.tsv
done

# Step 6: Combine QV results
cat benchmark/$chrom/$region/qv_prep/*/qv.tsv > benchmark/$chrom/$region/bestqv.tsv

# Step 7: Combine TPR and QV
Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/combine_tpr_qv.r \
    benchmark/$chrom/$region/tpr.tsv \
    benchmark/$chrom/$region/bestqv.tsv \
    $region \
    benchmark/$chrom/$region/tpr_qv.tsv