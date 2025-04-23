#!/bin/bash
#
# Usage: ./genotype_region.sh <bed_file_path> <dir_alignments> <path_reference> <path_reference_cram> <path_samples_txt> <dir_pangenome> <dir_reads> <threads> <work_dir> <output_dir>

set -e
set -o pipefail

# Function to parse command line arguments
parse_arguments() {
    if [ $# -lt 10 ]; then
        echo "Usage: $0 <bed_file_path> <dir_alignments> <path_reference> <path_reference_cram> <path_samples_txt> <dir_pangenome> <dir_reads> <threads> <work_dir> <output_dir>"
        exit 1
    fi
    
    # Return all arguments
    echo "$1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}"
}

# Function to extract region name and chromosome from BED file
get_region_info() {
    local bed_file="$1"
    local chrom=$(awk '{print $1}' "$bed_file" | head -n 1)
    local region=$(awk '{print $1"_"$2"_"$3}' "$bed_file" | head -n 1)
    echo "$chrom $region"
}

# Function to create and setup working directory
setup_environment() {
    local chrom="$1"
    local region="$2"
    local work_dir="$3"
    
    # Use SLURM JOB ID as part of output directory if available
    if [ -n "$SLURM_JOB_ID" ]; then
        local scratch_dir="${work_dir}/pangt-${SLURM_JOB_ID}"
    else
        local scratch_dir="${work_dir}/pangt-$$"  # Use PID if not running under SLURM
    fi
    
    # Create scratch directory and subdirectories
    mkdir -p "$scratch_dir"
    mkdir -p "$scratch_dir/impg/$chrom/$region"
    mkdir -p "$scratch_dir/pggb/$chrom/$region"
    
    # New directory structure for odgi
    mkdir -p "$scratch_dir/odgi/paths/matrix/$chrom"
    mkdir -p "$scratch_dir/odgi/view/$chrom"
    mkdir -p "$scratch_dir/odgi/dissimilarity/$chrom"
    
    mkdir -p "$scratch_dir/clusters/$chrom"
    mkdir -p "$scratch_dir/cosigt"
    mkdir -p "$scratch_dir/bwa-mem2"
    mkdir -p "$scratch_dir/gfainject"
    mkdir -p "$scratch_dir/gafpack"

    # Return directory paths
    echo "$scratch_dir $scratch_dir/impg/$chrom/$region $scratch_dir/pggb/$chrom/$region $scratch_dir/odgi/paths/matrix/$chrom $scratch_dir/odgi/view/$chrom $scratch_dir/odgi/dissimilarity/$chrom $scratch_dir/clusters/$chrom"
}

# Function to process sample PAF files
process_paf_file() {
    local paf_file="$1"
    local bed_file_path="$2"
    local impg_dir="$3"
    
    local sample=$(basename "$paf_file" -vs-grch38.aln.paf)
    echo "  Projecting alignments for sample $sample..."
    
    impg query -p "$paf_file" -b "$bed_file_path" > "$impg_dir/$sample.projected.bedpe"
    bedtools sort -i "$impg_dir/$sample.projected.bedpe" | bedtools merge -d 200000 > "$impg_dir/$sample.merged.bed"
    
    echo "  Completed projection for sample $sample"
}

# Function to filter sequences
filter_sequences() {
    local impg_dir="$1"
    local threads="$2"
    
    echo "  Filtering sequences..."
    # Combine all merged bed files
    cat $impg_dir/*.merged.bed > $impg_dir/ALL.tmp.bed
    
    # Filter outliers
    Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/outliers.r \
        $impg_dir/ALL.tmp.bed \
        $impg_dir/ALL.merged.filtered.bed
    
    rm $impg_dir/ALL.tmp.bed
    
    # Note: Additional QC could be implemented here as noted in the code comments
    # - Check if sequences fully span the locus (with 5% buffer on both sides)
    # - Filter based on flagger-flagged regions (max 5-10% of bad regions)
    # - Merqury-based filtering
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
    local threads="$4"
    local impg_dir="$5"
    local region="$6"
    local chrom="$7"
    local odgi_view_dir="$8"
    local odgi_paths_matrix_dir="$9"
    local clusters_dir="${10}"
    local scratch_dir="${11}"
    
    local sample=$(basename "$cram_file" .cram)
    echo "  Processing alignment for sample $sample..."
    
    # Set up working threads
    local bwa_threads=$(( threads > 4 ? threads - 6 : 1 ))
    local samtools_threads=$(( threads > 4 ? 3 : 1 ))
    
    # Create directories for this sample
    mkdir -p "$scratch_dir/bwa-mem2/$sample/$chrom"
    mkdir -p "$scratch_dir/gfainject/$sample/$chrom"
    mkdir -p "$scratch_dir/gafpack/$sample/$chrom"
    mkdir -p "$scratch_dir/cosigt/$sample/$chrom"
    
    # Extract reads from region
    samtools view \
        -T "$path_reference_cram" \
        -L "$bed_file_path" \
        -M \
        -b \
        -@ $samtools_threads \
        "$cram_file" | \
        samtools sort -n -T $scratch_dir | \
        samtools fasta | \
        bwa-mem2.avx mem -t $bwa_threads -p -h 10000 "$impg_dir/$region.pangenome.fa.gz" - | \
        samtools view -b -F 4 -@ $samtools_threads - \
        > "$scratch_dir/bwa-mem2/$sample/$chrom/$region.realigned.bam"
    
    # Inject and genotype
    gfainject \
        --gfa "$odgi_view_dir/$region.gfa" \
        --bam "$scratch_dir/bwa-mem2/$sample/$chrom/$region.realigned.bam" \
        --alt-hits 10000 | \
        gzip > "$scratch_dir/gfainject/$sample/$chrom/$region.gaf.gz"
    
    gafpack \
        --gfa "$odgi_view_dir/$region.gfa" \
        --gaf "$scratch_dir/gfainject/$sample/$chrom/$region.gaf.gz" \
        --len-scale \
        --weight-queries | \
        gzip > "$scratch_dir/gafpack/$sample/$chrom/$region.gafpack.gz"
    
    cosigt \
        -i "$sample" \
        -p "$odgi_paths_matrix_dir/$region.paths_matrix.tsv.gz" \
        -g "$scratch_dir/gafpack/$sample/$chrom/$region.gafpack.gz" \
        -c "$clusters_dir/$region.clusters.json" \
        -o "$scratch_dir/cosigt/$sample/$chrom/$region"
        
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
    local work_dir="${args[8]}"
    local output_dir="${args[9]}"
    
    # Get region info
    local region_info=($(get_region_info "$bed_file_path"))
    local chrom="${region_info[0]}"
    local region="${region_info[1]}"
    
    # Setup environment and get directories
    local dirs=($(setup_environment "$chrom" "$region" "$work_dir"))
    local scratch_dir="${dirs[0]}"
    local impg_dir="${dirs[1]}"
    local pggb_dir="${dirs[2]}"
    local odgi_paths_matrix_dir="${dirs[3]}"
    local odgi_view_dir="${dirs[4]}"
    local odgi_diss_dir="${dirs[5]}"
    local clusters_dir="${dirs[6]}"
    
    # Start logging for everything that follows
    log_file="$scratch_dir/${region}_processing.log"
    exec > >(tee -a "$log_file") 2>&1

    cd "$scratch_dir"
    export TMPDIR="$scratch_dir/tmp"
    mkdir -p "$TMPDIR"

    # Print configuration
    echo "=== Starting processing of region $region at $(date) ==="
    echo "BED file: $bed_file_path"
    echo "Alignments directory: $dir_alignments"
    echo "Reference: $path_reference"
    echo "Reference for CRAM: $path_reference_cram"
    echo "Samples list: $path_samples_txt"
    echo "Pangenome directory: $dir_pangenome"
    echo "Reads directory: $dir_reads"
    echo "Threads: $threads"
    echo "Working directory: $scratch_dir"
    echo "Output directory: $output_dir"
    
    echo "Processing $region..."
    
    # Export all functions for use with GNU Parallel
    export -f process_paf_file process_merged_bed process_cram_file
    
    # 1. Project alignments (parallel)
    echo "  Projecting alignments in parallel..."
    ls $dir_alignments/*-vs-grch38.aln.paf | grep -Ff $path_samples_txt > $scratch_dir/samples_paf.txt
    
    # Use parallel to process PAF files
    cat $scratch_dir/samples_paf.txt | parallel --tmpdir $scratch_dir -j $threads process_paf_file {} "$bed_file_path" "$impg_dir"
    
    # 2. Filter sequences
    #filter_sequences "$impg_dir" "$threads"
    
    # 3. Collect sequences
    echo "  Collecting sequences..."
    # Get reference sequence
    bedtools getfasta -fi $path_reference -bed "$bed_file_path" | sed 's/^>chr/>GRCh38#0#chr/g' > $impg_dir/$region.pangenome.fa
    
    # Process merged bed files in parallel
    ls $impg_dir/*.merged.bed | parallel --tmpdir $scratch_dir -j $threads process_merged_bed {} "$dir_pangenome" "$impg_dir"
    
    # Concatenate results
    echo "  Combining all fasta files..."
    cat $impg_dir/*.fasta >> $impg_dir/$region.pangenome.fa
    bgzip -@ $threads $impg_dir/$region.pangenome.fa
    samtools faidx $impg_dir/$region.pangenome.fa.gz
    
    # 4. Build pangenome graph
    echo "  Pangenome graph building..."
    pggb_param="-c 2"
    pggb -i $impg_dir/$region.pangenome.fa.gz -o $pggb_dir -t $threads -D $scratch_dir $pggb_param
    mv $pggb_dir/*smooth.final.og $pggb_dir/$region.final.og
   
    # 5. Process path coverage matrix
    echo "  Getting the path coverage matrix..."
    odgi paths -i $pggb_dir/$region.final.og -H -t $threads | cut -f 1,4- | gzip > $odgi_paths_matrix_dir/$region.paths_matrix.tsv.gz
    odgi view -i $pggb_dir/$region.final.og -g -t $threads > $odgi_view_dir/$region.gfa
    
    # 6. Clustering
    echo "  Pangenome clustering..."
    grep '^S' $odgi_view_dir/$region.gfa | awk '{{print("node."$2,length($3))}}' OFS="\t" > $odgi_view_dir/$region.node.length.tsv
    Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/filter.r \
        $odgi_paths_matrix_dir/$region.paths_matrix.tsv.gz \
        $odgi_view_dir/$region.node.length.tsv \
        no_filter \
        $odgi_paths_matrix_dir/$region.shared.tsv
    region_similarity=$(cut -f 3 $odgi_paths_matrix_dir/$region.shared.tsv | tail -1)
    odgi similarity -i $pggb_dir/$region.final.og --distances --all -t $threads > $odgi_diss_dir/$region.tsv
    
    Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/cluster.r \
        $odgi_diss_dir/$region.tsv \
        $clusters_dir/$region.clusters.json \
        automatic \
        $region_similarity
    
    # 7. Prepare graph for alignment
    echo "  Pangenome indexing..."
    bwa-mem2.avx index $impg_dir/$region.pangenome.fa.gz
    
    # 8. Alignment, injection, genotyping (parallel)
    echo "  Alignment, injection, genotyping in parallel..."
    ls $dir_reads/*cram | grep -Ff $path_samples_txt > $scratch_dir/samples_cram.txt
    sample_count=$(cat $scratch_dir/samples_cram.txt | wc -l)
    
    # Calculate per-sample thread allocation
    local threads_per_sample=8  # Minimum threads per sample
    local parallel_samples=$(( threads / threads_per_sample ))
    [[ $parallel_samples -lt 1 ]] && parallel_samples=1
    [[ $parallel_samples -gt $sample_count ]] && parallel_samples=$sample_count
    threads_per_sample=$(( threads / parallel_samples ))
    
    echo "  Processing $sample_count samples ($parallel_samples samples in parallel, $threads_per_sample threads each)"
    
    # Process CRAM files in parallel
    cat $scratch_dir/samples_cram.txt | parallel --tmpdir $scratch_dir -j $parallel_samples \
        process_cram_file {} "$path_reference_cram" "$bed_file_path" "$threads_per_sample" \
        "$impg_dir" "$region" "$chrom" "$odgi_view_dir" "$odgi_paths_matrix_dir" "$clusters_dir" "$scratch_dir"
    
    # Copy results to output directory
    echo "  Copying results to output directory..."
    mkdir -p "$output_dir"
    cp -r "$scratch_dir/impg" "$output_dir/"
    cp -r "$scratch_dir/pggb" "$output_dir/"
    cp -r "$scratch_dir/odgi" "$output_dir/"
    cp -r "$scratch_dir/clusters" "$output_dir/"
    cp -r "$scratch_dir/cosigt" "$output_dir/"
    cp -r "$scratch_dir/bwa-mem2" "$output_dir/"
    cp -r "$scratch_dir/gfainject" "$output_dir/"
    cp -r "$scratch_dir/gafpack" "$output_dir/"
    cp "$log_file" "$output_dir/"
    rm -rf "$scratch_dir"

    echo "=== Completed processing of region $region at $(date) ==="
    echo "Results are in $output_dir"
}

# Execute main function with all arguments
main "$@"