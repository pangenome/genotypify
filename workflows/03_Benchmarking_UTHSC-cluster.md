# Benchmarking

## Paths

```shell
dir_base=/lizardfs/guarracino/genotypify

dir_pangenome=/lizardfs/guarracino/pangenomes/HPRCv2
dir_reads=$dir_base/data/HPRCv2/illumina

export PATH="/lizardfs/guarracino/tools/bedtools2/bin:$PATH"
export PATH="/lizardfs/guarracino/tools/samtools-1.21:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/impg/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/wfmash/build/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/seqwish/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/smoothxg/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/GFAffix/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/odgi/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/pggb:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/gfainject/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/gafpack/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/bwa-mem2-2.2.1_x64-linux:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/cosigt:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/agc-1.1_x64-linux:$PATH"
#export PATH="/lizardfs/guarracino/condatools/emboss/6.6.0/bin/stretcher:$PATH" 
```

## 200 HPRCv2's short read samples on regions from the same 200 HPRCv2's assemblies (the full pangenome is 466 haplotypes)

Samples for which we have short reads:

```shell
ls $dir_base/data/HPRCv2/illumina/*.cram | while read f; do echo $(basename $f .final.cram); done > $dir_base/data/HPRCv2/illumina/samples.txt
```

Pangenome vs reference alignment:

```shell
mkdir -p $dir_base/wfmash
cd $dir_base/wfmash

ls $dir_pangenome/*.fa.gz | while read fasta; do
    sample=$(basename $fasta .fa.gz)

    sbatch -p allnodes -c 24 --job-name $sample-vs-grch38 --wrap "hostname; cd /scratch; wfmash $dir_base/reference/GRCh38.fa.gz $fasta -s 10k -p 95 -t 24 > $sample-vs-grch38.aln.paf; mv $sample-vs-grch38.aln.paf $dir_base/wfmash/"
done
```

Genotyping:

```shell
mkdir -p $dir_base/genotyping/logs
cd $dir_base/genotyping/logs

ls $dir_base/data/HPRCv2/illumina/*.cram | while read f; do echo $(basename $f .final.cram); done > $dir_base/genotyping/samples-to-consider.txt

# Create regions file with all region info (sorted by length to start with the shortest ones)
for f in $dir_base/data/loci/for_testing/*.bed; do 
    sed '1d' $f | cut -f 1-4
done | awk -v OFS='\t' '{name=$4; gsub(/[^a-zA-Z0-9._-]/, "-", name); print $1,$2,$3,substr(name, 1, 32),$3-$2}' | \
    sort -k 5n > $dir_base/genotyping/regions.txt

# Count total regions
num_regions=$(wc -l < $dir_base/genotyping/regions.txt)
echo "Total regions to process: $num_regions"

# Submit job array, max 8 jobs at a time
sbatch --array=1-$num_regions%8 -c 48 -p allnodes --job-name genotype_array --output=$dir_base/genotyping/logs/genotype_%a.out --error=$dir_base/genotyping/logs/genotype_%a.err $dir_base/scripts/genotype_region.array-wrapper.sh $dir_base
```

Check results:

```shell
cat $dir_base/genotyping/regions.txt | head -n 3465 | awk '{print NR, $1"_"$2"_"$3}' | while read -r row_num d; do
    if [ ! -d "$dir_base/genotyping/$d" ]; then
        echo "Row $row_num: Directory $dir_base/genotyping/$d does not exist."
        #sbatch --array=$row_num-$row_num -c 48 -p allnodes --job-name genotype_array --output=$dir_base/genotyping/logs/genotype_%a.out --error=$dir_base/genotyping/logs/genotype_%a.err $dir_base/scripts/genotype_region.array-wrapper.sh $dir_base
    fi
done
```

Benchmarking:

```shell
# List completed genotyping results
mkdir -p /scratch/hprcv2-benchmark/
# cat "$dir_base/genotyping-no-sparse-map/regions.txt" | head -n 3465 | awk '{print $1, $2, $3, NR, $1"_"$2"_"$3}' | grep -Ff <(cut -f 5 -d ' ' /scratch/hprcv2-benchmark-no-sparse-map/completed_regions.txt  | rev | cut -f 1 -d '/' | rev) -vw | while read -r chrom start end row_num d; do
cat "$dir_base/genotyping/regions.txt" | head -n 3465 | awk '{print $1, $2, $3, NR, $1"_"$2"_"$3}' | while read -r chrom start end row_num d; do
    if [ -d "$dir_base/genotyping/$d" ]; then
        # If directory exists, attempt the grep
        if ! grep -q 'Results are in' "$dir_base/genotyping/logs/genotype_${row_num}.out"; then
            #echo "Row $row_num: Directory $dir_base/genotyping/$d exists but grep failed."
        else
            echo $chrom $start $end $row_num $dir_base/genotyping/$d
        fi
    else
        #echo "Row $row_num: Directory $dir_base/genotyping/$d does not exist."
    fi
done > /scratch/hprcv2-benchmark/completed_regions.txt

# Prepare for benchmarking
mkdir -p /scratch/hprcv2-benchmark/cosigt
mkdir -p /scratch/hprcv2-benchmark/clusters
mkdir -p /scratch/hprcv2-benchmark/odgi/dissimilarity
mkdir -p /scratch/hprcv2-benchmark/odgi
cat /scratch/hprcv2-benchmark/completed_regions.txt | while read -r chrom start end row_num out_dir; do
    echo "Preparing region: $chrom:$start-$end (row $row_num, dir: $out_dir)"

    cp -r $out_dir/cosigt/* /scratch/hprcv2-benchmark/cosigt
    cp -r $out_dir/clusters/* /scratch/hprcv2-benchmark/clusters

    cp -r $out_dir/odgi/dissimilarity/* /scratch/hprcv2-benchmark/odgi/dissimilarity
    find /scratch/hprcv2-benchmark/odgi/dissimilarity -type f -name "*.tsv" -exec pigz -f -9 -p 4 {} +

    chr_name=$(basename $out_dir | cut -d'_' -f1)
    region=$(basename $out_dir)
    pggb_dest_dir=/scratch/hprcv2-benchmark/pggb/$chr_name/$region
    mkdir -p $pggb_dest_dir
    cp -r $out_dir/pggb/*/*/*.final.og $pggb_dest_dir

    odgi_dest_dir=/scratch/hprcv2-benchmark/odgi/view/$chr_name
    mkdir -p $odgi_dest_dir
    cp -r $out_dir/odgi/view/*/*.gfa $odgi_dest_dir
done

# Benchmark
echo "Starting benchmarking..."
cd /scratch/hprcv2-benchmark
mkdir -p /scratch/hprcv2-benchmark/logs
cat /scratch/hprcv2-benchmark/completed_regions.txt | parallel --colsep ' ' -j 46 "$dir_base/scripts/benchmark_region.sh {1} {2} {3} > logs/{1}_{2}_{3}.log 2>&1"

# Final step: Plot TPR results
echo "Plotting TPR results..."
join <(cat $dir_base/genotyping/regions.txt | awk '{print($1"_"$2"_"$3,$0)}' | sort) <(cat completed_regions.txt | awk '{print($1"_"$2"_"$3,$0)}' | sort) | cut -f 2-5,10 -d ' ' | sort -k 5,5n > completed_regions_with_names.txt
Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/plot_tpr.r \
    benchmark/tpr \
    completed_regions_with_names.txt \
    benchmark/*/*/tpr_qv.tsv

# Make a summary of the results
echo "Creating summary of results..."
(echo -e "name\tchrom\tstart\tend\tnum.haps\tref.length\tmin.hap.len\tmax.hap.len\tavg.hap.len\tsd.hap.length\tmin.jaccard.dist\tmax.jaccard.dist\tavg.jaccard.dist\tsd.jaccard.dist\tnum.clusters\tperc.haps.QV.high\tperc.haps.QV.mid\tperc.haps.QV.low\tperc.haps.QV.verylow\tnum.haps.QV.high\tnum.haps.QV.mid\tnum.haps.QV.low\tnum.haps.QV.verylow\t"; \
cat completed_regions_with_names.txt | while read -r chrom start end name num_row; do
    length=$(awk -v chrom="$chrom" -v start="$start" -v end="$end" 'BEGIN {print end - start}')
    region="${chrom}_${start}_${end}"

    # Calculate all statistics
    hap_lengths=$(odgi paths -i /scratch/hprcv2-benchmark-no-sparse-map/pggb/$chrom/$region/$region.final.og -Ll | awk '{print $3}')
    read -r min_hap_len max_hap_len avg_hap_len sd_hap_len <<< $(echo "$hap_lengths" | \
    awk '
    {
        len[NR] = $1
        sum += $1
        
        if (NR == 1 || $1 < min) min = $1
        if (NR == 1 || $1 > max) max = $1
    }
    END {
        avg = sum / NR
        
        # Calculate standard deviation
        for (i = 1; i <= NR; i++) {
            diff = len[i] - avg
            sum_sq_diff += diff * diff
        }
        sd = sqrt(sum_sq_diff / NR)
        
        print min, max, avg, sd
    }')

    num_haps=$(echo "$hap_lengths" | wc -l)

    jaccards=$(cut -f 6 $dir_base/genotyping/$region/odgi/dissimilarity/$chrom/$region.tsv | sed 1d)
    read -r min_j max_j avg_j sd_j <<< $(echo "$jaccards" | \
    awk '
    {
        len[NR] = $1
        sum += $1
        
        if (NR == 1 || $1 < min) min = $1
        if (NR == 1 || $1 > max) max = $1
    }
    END {
        avg = sum / NR
        
        # Calculate standard deviation
        for (i = 1; i <= NR; i++) {
            diff = len[i] - avg
            sum_sq_diff += diff * diff
        }
        sd = sqrt(sum_sq_diff / NR)
        
        print min, max, avg, sd
    }')

    # Calculate and assign all 8 variables at once
    eval $(cut -f 16,17 /scratch/hprcv2-benchmark-no-sparse-map/benchmark/$chrom/$region/tpr_qv.tsv | sed 1d | tr '\t' '\n' | \
    awk '
    {
        total++
        
        if ($1 > 33) {
            count_gt33++
        } else if ($1 > 23) {
            count_23to33++
        } else if ($1 > 17) {
            count_17to23++
        } else {
            count_lte17++
        }
    }
    END {
        # Output variable assignments for counts
        print "num_hap_gt33=" count_gt33+0
        print "num_hap_23to33=" count_23to33+0
        print "num_hap_17to23=" count_17to23+0
        print "num_hap_lte17=" count_lte17+0
        
        # Output variable assignments for percentages
        printf "pct_hap_gt33=%.4f\n", (count_gt33/total)*100
        printf "pct_hap_23to33=%.4f\n", (count_23to33/total)*100
        printf "pct_hap_17to23=%.4f\n", (count_17to23/total)*100
        printf "pct_hap_lte17=%.4f\n", (count_lte17/total)*100
    }')

    num_clusters=$(sed 1d /scratch/hprcv2-benchmark-no-sparse-map/clusters/$chrom/$region.clusters.tsv | cut -f 2 | sort | uniq | wc -l)

    echo -e "$name\t$chrom\t$start\t$end\t$num_haps\t$length\t$min_hap_len\t$max_hap_len\t$avg_hap_len\t$sd_hap_len\t$min_j\t$max_j\t$avg_j\t$sd_j\t$num_clusters\t$pct_hap_gt33\t$pct_hap_23to33\t$pct_hap_17to23\t$pct_hap_lte17\t$num_hap_gt33\t$num_hap_23to33\t$num_hap_17to23\t$num_hap_lte17"
done) > summary.tsv


###########################################################################################################################################
# OK
(echo -e "region\tnum_samples\tstep"; cat $dir_base/genotyping/regions.txt | while read -r chrom start end name len; do
    region="${chrom}_${start}_${end}"
    
    # Count samples for each step in this specific region
    #raw_count=$(find "$dir_base/genotyping/$region/impg/" -name "*.projected.bedpe" 2>/dev/null | wc -l)
    #span_count=$(find "$dir_base/genotyping/$region/impg/" -name "*.projected.filtered.bedpe" 2>/dev/null | wc -l)
    #flagger_count=$(find "$dir_base/genotyping/$region/impg/" -name "*.merged.filtered.bed" 2>/dev/null | wc -l)
    raw_count=$(find "/scratch/HPRCv2-gt/impg/$chrom/$region" -name "*.projected.bedpe" 2>/dev/null | wc -l)
    span_count=$(find "/scratch/HPRCv2-gt/impg/$chrom/$region" -name "*.projected.filtered.bedpe" 2>/dev/null | wc -l)
    flagger_count=$(find "/scratch/HPRCv2-gt/impg/$chrom/$region" -name "*.merged.filtered.bed" 2>/dev/null | wc -l)

    # Output results for this region
    echo -e "${region}_${name}\t$raw_count\tRAW"
    echo -e "${region}_${name}\t$span_count\tSPAN"
    echo -e "${region}_${name}\t$flagger_count\tFLAGGER"
done) > filtering.tsv

# (echo region haplotype step; find impg/*/*/ -name "*.projected.bedpe" | while read bedpe; do
#     sample=$(basename $bedpe .projected.bedpe)
#     chrom=$(echo $bedpe | cut -d'/' -f2)
#     region=$(echo $bedpe | cut -d'/' -f3)
#     echo $region $sample RAW
# done;
# find impg/*/*/ -name "*.projected.bedpe" | while read bedpe; do
#     sample=$(basename $bedpe .projected.filtered.bedpe)
#     chrom=$(echo $bedpe | cut -d'/' -f2)
#     region=$(echo $bedpe | cut -d'/' -f3)
#     echo $region $sample SPAN
# done;
# find impg/*/*/ -name "*.merged.filtered.bed" | while read bed; do
#     sample=$(basename $bed .merged.filtered.bed)
#     chrom=$(echo $bedpe | cut -d'/' -f2)
#     region=$(echo $bedpe | cut -d'/' -f3)

#     echo $region $sample FLAGGER
# done) | tr ' ' '\t' > filtering.tsv


chrom=chr11; start=69809692; end=69819692; name=FGF3; len=10000
chrom=chr6; start=31891045; end=32123783; name=C4A,C4B; len=232738

chrom=chr10; start=104572217; end=104577217; name=HARsv2_0295-LOC101927523-SORCS3; len=5000

mkdir -p /scratch/HPRCv2-gt
cd /scratch/HPRCv2-gt
mkdir -p regions_of_interest
#ls $dir_base/data/HPRCv2/illumina/*.cram | while read f; do echo $(basename $f .final.cram); done > samples-to-consider.txt

dir_flagger="$dir_base/data/HPRCv2/flagger"
cat $dir_base/genotyping/regions2.txt | while read chrom start end name len; do
    echo "Processing region: $chrom:$start-$end ($name, len: $len)"

    region=${chrom}_${start}_${end}
    
    echo -e "$chrom\t$start\t$end\t$name" > regions_of_interest/$region.bed

    echo "  Projecting alignments..."
    mkdir -p impg/$chrom/$region
    #ls $dir_base/wfmash/*-vs-grch38.aln.paf | grep -Ff $dir_base/genotyping/samples-to-consider.txt | while read paf; do
    ls wfmash/*-vs-grch38.aln.paf | grep -Ff $dir_base/genotyping/samples-to-consider.txt | while read paf; do
        sample=$(basename $paf -vs-grch38.aln.paf)
        
        impg query \
            -p $paf \
            -b regions_of_interest/$region.bed \
            > impg/$chrom/$region/$sample.projected.bedpe
    done

    echo "  Filtering and merging sequences..."
    # What we want at the end for each haplotype:
    # - 1 contig
    # - 1 contig spanning 1000bp at the beginning/end of the locus
    # - 1 contig after merging (-d 200k)
    # - 1 contig with <= 5% FLAGGER-flagged regions
    ls impg/$chrom/$region/*.projected.bedpe | while read bedpe; do
        sample=$(basename $bedpe .projected.bedpe)
        pass=$(python3 $dir_base/scripts/check_locus_span.py -b "$bedpe" -l "${chrom}:${start}-${end}" -bp 1000 -d 200000)
        if $pass; then
            cp $bedpe impg/$chrom/$region/$sample.projected.filtered.bedpe
            # Now perform the merge on the filtered bedpe
            bedtools sort -i impg/$chrom/$region/$sample.projected.filtered.bedpe | bedtools merge -d 200000 > impg/$chrom/$region/$sample.merged.bed
        else
            echo "  Filtering $sample for not fully spanning the locus with a single contig"
        fi
    done

    echo "  Applying FLAGGER filtering..."
    ls impg/$chrom/$region/*.merged.bed | while read bed; do
        sample=$(basename $bed .merged.bed)

        # Get the total length of the merged interval
        merged_length=$(awk '{sum += $3-$2} END {print sum}' "$bed")
        
        # Check if flagger file exists for this sample
        if [ -f "$dir_flagger/$sample.bed" ]; then
            # Get the total length of overlapping regions
            overlap_length=$(bedtools intersect -a "$bed" -b "$dir_flagger/$sample.bed" -wao | awk '{sum += $NF} END {print sum}')
            # Calculate percentage
            percentage=$(echo "scale=4; ($overlap_length / $merged_length) * 100" | bc)
            # Check if percentage is greater than 5%
            if (( $(echo "$percentage <= 5" | bc -l) )); then
                cp impg/$chrom/$region/$sample.merged.bed impg/$chrom/$region/$sample.merged.filtered.bed
            else
                echo "    Filtering $sample for having $percentage% >= 5% of regions flagged by FLAGGER"
            fi
        else
            # If no flagger file exists for this sample, keep the sample
            echo "    No FLAGGER file for $sample, keeping the merged bed"
            cp impg/$chrom/$region/$sample.merged.bed impg/$chrom/$region/$sample.merged.filtered.bed
        fi
    done
done &> xxx.log

    echo "  Collecting sequences..."
    bedtools getfasta -fi $dir_base/reference/GRCh38.fa.gz -bed regions_of_interest/$region.bed | sed 's/^>chr/>GRCh38#0#chr/g' > impg/$chrom/$region/$region.pangenome.fa

    ls impg/$chrom/$region/*.merged.filtered.bed | while read bed; do
        sample=$(basename $bed .merged.filtered.bed)
        fasta=$dir_pangenome/$sample.fa.gz

        bedtools getfasta -fi $fasta -bed $bed
    done >> impg/$chrom/$region/$region.pangenome.fa
    bgzip -@ 48 impg/$chrom/$region/$region.pangenome.fa
    samtools faidx impg/$chrom/$region/$region.pangenome.fa.gz


    echo "  Pangenome graph building..."
    mkdir -p pggb/$chrom
    pggb_param="-c 2"
    pggb -i impg/$chrom/$region/$region.pangenome.fa.gz -o pggb/$chrom/$region -t 48 -D /scratch $pggb_param
    mv pggb/$chrom/$region/*smooth.final.og pggb/$chrom/$region/$region.final.og # Rename the final ODGI graph


    echo "  Getting the path coverage matrix..."
    #mkdir -p odgi/chop/$chrom
    #odgi chop -i pggb/$chrom/$region/$region.final.og -c 32 -o odgi/chop/$chrom/$region.chopped.og
    mkdir -p odgi/paths/matrix/$chrom
    #odgi paths -i odgi/$chrom/$region.chopped.og -H | cut -f 1,4- | gzip > odgi/paths/matrix/$chrom/$region.paths_matrix.tsv.gz
    odgi paths -i pggb/$chrom/$region/$region.final.og -H | cut -f 1,4- | gzip > odgi/paths/matrix/$chrom/$region.paths_matrix.tsv.gz
    mkdir -p odgi/view/$chrom
    #odgi view -i odgi/chop/$chrom/$region.chopped.og -g > odgi/view/$chrom/$region.gfa
    odgi view -i pggb/$chrom/$region/$region.final.og -g > odgi/view/$chrom/$region.gfa


    echo "  Pangenome clustering..."
    mkdir -p odgi/dissimilarity/$chrom
    #odgi similarity -i odgi/$chrom/$region.chopped.og --distances -t 48 > odgi/dissimilarity/$chrom/$region.tsv
    odgi similarity -i pggb/$chrom/$region/$region.final.og --distances --all -t 48 > odgi/dissimilarity/$chrom/$region.tsv
    mkdir -p clusters/$chrom
    grep '^S' odgi/view/$chrom/$region.gfa | awk '{{print("node."$2,length($3))}}' OFS="\t" > odgi/view/$chrom/$region.node.length.tsv
    Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/filter.r \
        odgi/paths/matrix/$chrom/$region.paths_matrix.tsv.gz \
        odgi/view/$chrom/$region.node.length.tsv \
        no_filter \
        odgi/paths/matrix/$chrom/$region.shared.tsv
    region_similarity=$(cut -f 3 odgi/paths/matrix/$chrom/$region.shared.tsv | tail -1)
    Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/cluster.r \
        odgi/dissimilarity/$chrom/$region.tsv \
        clusters/$chrom/$region.clusters.json \
        automatic \
        $region_similarity


    echo "  Pangenome indexing..."
    bwa-mem2.avx index impg/$chrom/$region/$region.pangenome.fa.gz # "bwa-mem2 index" returns "ERROR: prefix is too long!"

    echo "  Alignment, injection, genotyping..."
    ls $dir_reads/*cram | grep -Ff samples-to-consider.txt | while read cram; do
        echo $cram

        sample=$(basename $cram .cram)

        mkdir -p bwa-mem2/$sample/$chrom
        samtools view \
            -T $dir_base/reference/GRCh38.fa.gz \
            -L $dir_base/regions_of_interest/$region.bed \
            -M \
            -b \
            -@ 3 \
            $cram | \
            samtools sort -n -T /scratch | \
            samtools fasta | \
            bwa-mem2.avx mem -t 42 -p -h 10000 impg/$chrom/$region/$region.pangenome.fa.gz - | \
            samtools view -b -F 4 -@ 3 - \
            > bwa-mem2/$sample/$chrom/$region.realigned.bam

        mkdir -p gfainject/$sample/$chrom
        gfainject \
            --gfa odgi/view/$chrom/$region.gfa \
            --bam bwa-mem2/$sample/$chrom/$region.realigned.bam \
            --alt-hits 10000 | \
            gzip > gfainject/$sample/$chrom/$region.gaf.gz

        mkdir -p gafpack/$sample/$chrom
        gafpack \
            --gfa odgi/view/$chrom/$region.gfa \
            --gaf gfainject/$sample/$chrom/$region.gaf.gz \
            --len-scale \
            --weight-queries | \
            gzip > gafpack/$sample/$chrom/$region.gafpack.gz

        mkdir -p cosigt/$sample/$chrom
        cosigt \
            -i $sample \
            -p odgi/paths/matrix/$chrom/$region.paths_matrix.tsv.gz \
            -g gafpack/$sample/$chrom/$region.gafpack.gz \
            -c clusters/$chrom/$region.clusters.json \
            -o cosigt/$sample/$chrom/$region
    done
done


for f in $dir_base/data/loci/*.bed; do sed '1d' $f | cut -f 1-4; done | awk -v OFS='\t' '{name=$4; gsub(/[^a-zA-Z0-9._-]/, "-", name); print $1,$2,$3,substr(name, 1, 32),$3-$2}' | sort -k 5n | head -n 150 | while read chrom start end name len; do
    echo "Processing region: $chrom:$start-$end ($name, len: $len)"

    region=${chrom}_${start}_${end}
    mkdir -p benchmark/$chrom/$region

    # Step 1: Make TPR table
    Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/calc_tpr.r \
        odgi/dissimilarity/$chrom/$region.tsv \
        clusters/$chrom/$region.clusters.json \
        clusters/$chrom/$region.clusters.hapdist.tsv \
        benchmark/$chrom/$region/tpr.tsv \
        cosigt/*/$chrom/$region/sorted_combos.tsv

    # Step 2: Flip PGGB graph
    odgi flip \
        -i pggb/$chrom/$region/$region.final.og \
        -o benchmark/$chrom/$region/$region.flip.og -P \
        --ref-flips <(grep '^P' odgi/view/$chrom/$region.gfa | cut -f 2 | grep "GRCh38#0")

    # Step 3: Convert OG to FASTA
    odgi paths \
        -i benchmark/$chrom/$region/$region.flip.og \
        -f | sed 's/_inv$//g' > benchmark/$chrom/$region/$region.flip.fasta

    # Step 4: Prepare combinations for QV
    mkdir -p benchmark/$chrom/$region/qv_prep
    bash /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/prepare_qv.sh \
        benchmark/$chrom/$region/tpr.tsv \
        benchmark/$chrom/$region/$region.flip.fasta \
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
done

# Final step: Plot TPR results
for f in $dir_base/data/loci/*.bed; do sed '1d' $f | cut -f 1-4; done > tmp.bed
Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/plot_tpr.r \
    benchmark/tpr \
    tmp.bed \
    benchmark/*/*/tpr_qv.tsv
rm tmp.bed
###########################################################################################################################################


# OPTIONAL (TO RE-CHECK)
mkdir -p annotation
bedtools intersect -wa -a $dir_base/annotation/gencode.v47.annotation.gtf.gz -b $dir_base/regions_of_interest/$region.bed > annotation/${region}_${chrom}_${start}_${end}.gtf
Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/annotate.r annotation/${region}_${chrom}_${start}_${end}.gtf chr6 annotation/${region}_${chrom}_${start}_${end}.bed
odgi procbed -i pggb/$region/$region.final.og -b <(sed 's/^chr/GRCh38#0#chr/g' annotation/${region}_${chrom}_${start}_${end}.bed) > annotation/${region}_${chrom}_${start}_${end}.procbed.bed
odgi inject -i pggb/$region/$region.final.og -b annotation/${region}_${chrom}_${start}_${end}.procbed.bed -o annotation/$region.final.inject.og -P
odgi flip -i annotation/$region.final.inject.og --ref-flips <(cut -f 1 annotation/${region}_${chrom}_${start}_${end}.procbed.bed | uniq) -o annotation/$region.final.inject.flip.og
odgi untangle -i annotation/$region.final.inject.flip.og -R <(cut -f 4 annotation/${region}_${chrom}_${start}_${end}.procbed.bed) -j 0.5 -g -t 8 > annotation/$region.gggenes.tsv
Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/plotgggenes.r annotation/$region.gggenes.tsv clusters/$region.clusters.json annotation/$region.gggenes.png
```
