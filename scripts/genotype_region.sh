#!/bin/bash
#
# Usage: ./genotype_region.sh <bed_file_path> <dir_alignments> <path_reference> <path_reference_cram> <path_samples_txt> <dir_pangenome> <dir_reads> <threads>

set -e
set -o pipefail

# Function to parse command line arguments
parse_arguments() {
    if [ $# -lt 8 ]; then
        echo "Usage: $0 <bed_file_path> <dir_alignments> <path_reference> <path_reference_cram> <path_samples_txt> <dir_pangenome> <dir_reads> <threads>"
        exit 1
    fi
    
    # Return all arguments
    echo "$1 $2 $3 $4 $5 $6 $7 $8"
}

# Function to extract region name from BED file
get_region_name() {
    local bed_file="$1"
    awk '{if(NF>=4) {print $1"_"$2"_"$3"_"$4} else {print $1"_"$2"_"$3}}' "$bed_file" | head -n 1
}

# Function to create and setup working directory
setup_environment() {
    local region="$1"
    
    # Use SLURM JOB ID as part of output directory if available
    if [ -n "$SLURM_JOB_ID" ]; then
        local scratch_dir="/scratch/gt-${SLURM_JOB_ID}"
    else
        local scratch_dir="/scratch/gt-$$"  # Use PID if not running under SLURM
    fi
    
    # Create scratch directory and subdirectories
    mkdir -p "$scratch_dir"
    mkdir -p "$scratch_dir/impg/$region"
    mkdir -p "$scratch_dir/pggb/$region"
    mkdir -p "$scratch_dir/odgi"
    mkdir -p "$scratch_dir/odgi/dissimilarity"
    mkdir -p "$scratch_dir/clusters"
    mkdir -p "$scratch_dir/alignments/$region"
    
    # Setup logging
    local log_file="$scratch_dir/${region}_processing.log"
    exec > >(tee -a "$log_file") 2>&1
    
    # Return directory paths
    echo "$scratch_dir $scratch_dir/impg/$region $scratch_dir/pggb/$region $scratch_dir/odgi $scratch_dir/odgi/dissimilarity $scratch_dir/clusters $scratch_dir/alignments/$region"
}

# Function to process sample PAF files
process_paf_file() {
    local paf_file="$1"
    local bed_file_path="$2"
    local impg_dir="$3"
    
    local sample=$(basename "$paf_file" -vs-grch38.aln.paf)
    echo "  Projecting alignments for sample $sample..."
    
    impg query -p "$paf_file" -b "$bed_file_path" > "$impg_dir/$sample.projected.bedpe"
    bedtools sort -i "$impg_dir/$sample.projected.bedpe" | bedtools merge -d 100000 > "$impg_dir/$sample.merged.bed"
    
    echo "  Completed projection for sample $sample"
}

# Function to process merged BED files
process_merged_bed() {
    local bed_file="$1"
    local dir_pangenome="$2"
    local impg_dir="$3"
    
    local sample=$(basename "$bed_file" .merged.bed)
    echo "  Getting fasta for sample $sample..."
    
    bedtools getfasta -fi "$dir_pangenome/$sample.fa.gz" -bed "$bed_file" > "$impg_dir/$sample.fasta"
}

# Function to process CRAM files
process_cram_file() {
    local cram_file="$1"
    local path_reference_cram="$2"
    local bed_file_path="$3"
    local threads_per_sample="$4"
    local impg_dir="$5"
    local region="$6"
    local odgi_dir="$7"
    local align_dir="$8"
    local scratch_dir="$9"
    local clusters_dir="${10}"
    
    local sample=$(basename "$cram_file" .cram)
    echo "  Processing alignment for sample $sample..."
    
    # Set up working threads
    local bwa_threads=$(( threads_per_sample > 1 ? threads_per_sample - 2 : 1 ))
    local samtools_threads=$(( threads_per_sample > 1 ? 2 : 1 ))
    
    # Create cosigt directory for this sample
    mkdir -p "$scratch_dir/cosigt/$sample"
    
    # Extract reads from region
    samtools view \
        -T "$path_reference_cram" \
        -L "$bed_file_path" \
        -M \
        -b \
        "$cram_file" | \
        samtools sort -n | \
        samtools fasta | \
        bwa-mem2.avx mem -t $bwa_threads "$impg_dir/$region.pangenome.fa.gz" - | \
        samtools view -b -F 4 -@ $samtools_threads - \
        > "$align_dir/$sample.reads-vs-pangenome.bam"
    
    # Inject and genotype
    gfainject \
        --gfa "$odgi_dir/$region.chopped.gfa" \
        --bam "$align_dir/$sample.reads-vs-pangenome.bam" \
        > "$align_dir/$sample.reads-vs-graph.gaf"
    
    gafpack \
        --gfa "$odgi_dir/$region.chopped.gfa" \
        --gaf "$align_dir/$sample.reads-vs-graph.gaf" \
        --len-scale | \
        gzip > "$align_dir/$sample.coverage.gafpack.gz"
    
    cosigt \
        -i "$sample" \
        -p "$odgi_dir/$region.paths_matrix.tsv.gz" \
        -g "$align_dir/$sample.coverage.gafpack.gz" \
        -c "$clusters_dir/$region.clusters.json" \
        -o "$scratch_dir/cosigt/$sample/$region"
        
    echo "  Completed alignment process for sample $sample"
}

# Main function to orchestrate the workflow
main() {
    # Parse arguments
    local args=($(parse_arguments "$@"))
    local bed_file_path="${args[0]}"
    local dir_alignments="${args[1]}"
    local path_reference="${args[2]}"
    local path_reference_cram="${args[3]}"
    local path_samples_txt="${args[4]}"
    local dir_pangenome="${args[5]}"
    local dir_reads="${args[6]}"
    local threads="${args[7]}"
    
    # Get region name
    local region=$(get_region_name "$bed_file_path")
    
    # Setup environment and get directories
    local dirs=($(setup_environment "$region"))
    local scratch_dir="${dirs[0]}"
    local impg_dir="${dirs[1]}"
    local pggb_dir="${dirs[2]}"
    local odgi_dir="${dirs[3]}"
    local diss_dir="${dirs[4]}"
    local clusters_dir="${dirs[5]}"
    local align_dir="${dirs[6]}"
    
    cd "$scratch_dir"
    
    # Print configuration
    echo "=== Starting processing of region $region at $(date) ==="
    echo "Working directory: $scratch_dir"
    echo "BED file: $bed_file_path"
    echo "Alignments directory: $dir_alignments"
    echo "Reference: $path_reference"
    echo "Reference for CRAM: $path_reference_cram"
    echo "Samples list: $path_samples_txt"
    echo "Pangenome directory: $dir_pangenome"
    echo "Reads directory: $dir_reads"
    echo "Threads: $threads"
    
    echo "Processing $region..."
    
    # Export all functions for use with GNU Parallel
    export -f process_paf_file process_merged_bed process_cram_file
    
    # 1. Project alignments (parallel)
    echo "  Projecting alignments in parallel..."
    ls $dir_alignments/*-vs-grch38.aln.paf | grep -Ff $path_samples_txt > $scratch_dir/samples_paf.txt
    
    # Use parallel to process PAF files
    cat $scratch_dir/samples_paf.txt | parallel -j $threads process_paf_file {} "$bed_file_path" "$impg_dir"
    
    # 2. Collect sequences
    echo "  Collecting sequences..."
    # Get reference sequence
    bedtools getfasta -fi $path_reference -bed "$bed_file_path" | \
        sed 's/^>chr/>GRCh38#0#chr/g' > $impg_dir/$region.pangenome.fa
    
    # Process merged bed files in parallel
    ls $impg_dir/*.merged.bed | parallel -j $threads process_merged_bed {} "$dir_pangenome" "$impg_dir"
    
    # Concatenate results
    echo "  Combining all fasta files..."
    cat $impg_dir/*.fasta >> $impg_dir/$region.pangenome.fa
    bgzip -@ $threads $impg_dir/$region.pangenome.fa
    samtools faidx $impg_dir/$region.pangenome.fa.gz
    
    # 3. Build pangenome graph
    echo "  Pangenome graph building..."
    pggb -i $impg_dir/$region.pangenome.fa.gz -o $pggb_dir -t $threads -D $scratch_dir -c 2
    mv $pggb_dir/*smooth.final.og $pggb_dir/$region.final.og
    
    # 4. Process path coverage matrix
    echo "  Getting the path coverage matrix..."
    odgi chop -i $pggb_dir/$region.final.og -c 30 -t $threads -o $odgi_dir/$region.chopped.og
    odgi paths -i $odgi_dir/$region.chopped.og -H -t $threads | \
        cut -f 1,4- | \
        gzip > $odgi_dir/$region.paths_matrix.tsv.gz
    
    # 5. Clustering
    echo "  Pangenome clustering..."
    odgi similarity -i $odgi_dir/$region.chopped.og -d -t $threads > $diss_dir/$region.tsv
    Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/cluster.r \
        $diss_dir/$region.tsv $clusters_dir/$region.clusters.json automatic 0
    
    # 6. Prepare graph for alignment
    echo "  Pangenome indexing..."
    bwa-mem2.avx index $impg_dir/$region.pangenome.fa.gz
    odgi view -i $odgi_dir/$region.chopped.og -g -t $threads > $odgi_dir/$region.chopped.gfa
    
    # 7. Alignment, injection, genotyping (parallel)
    echo "  Alignment, injection, genotyping in parallel..."
    ls $dir_reads/*cram | grep -Ff $path_samples_txt > $scratch_dir/samples_cram.txt
    sample_count=$(cat $scratch_dir/samples_cram.txt | wc -l)
    
    # Calculate per-sample thread allocation
    local threads_per_sample=6  # Minimum threads per sample
    local parallel_samples=$(( threads / threads_per_sample ))
    [[ $parallel_samples -lt 1 ]] && parallel_samples=1
    [[ $parallel_samples -gt $sample_count ]] && parallel_samples=$sample_count
    threads_per_sample=$(( threads / parallel_samples ))
    
    echo "  Processing $sample_count samples with $parallel_samples in parallel ($threads_per_sample threads each)"
    
    # Process CRAM files in parallel
    cat $scratch_dir/samples_cram.txt | parallel -j $parallel_samples \
        process_cram_file {} "$path_reference_cram" "$bed_file_path" "$threads_per_sample" \
        "$impg_dir" "$region" "$odgi_dir" "$align_dir" "$scratch_dir" "$clusters_dir"
    
    echo "=== Completed processing of region $region at $(date) ==="
    echo "Results are in $scratch_dir"
}

# Execute main function with all arguments
main "$@"