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

    sbatch -p allnodes -c 24 --job-name $sample-vs-grch38 --wrap "hostname; cd /scratch; wfmash $dir_base/reference/GRCh38.fa.gz $fasta -s 10k -p 90 -t 24 > $sample-vs-grch38.aln.paf; mv $sample-vs-grch38.aln.paf $dir_base/wfmash/"
done
```

Genotyping:

```shell
LENGTH_THRESHOLD=3000000

mkdir -p $dir_base/genotyping
ls $dir_base/data/HPRCv2/illumina/*.cram | while read f; do echo $(basename $f .final.cram); done | head -n 20 > $dir_base/genotyping/samples-to-consider.txt
# Project the aligments to the region-of-interest (ROI)
# Sort by length to start with the shortest ones
cd $dir_base/genotyping
cat $dir_base/data/loci.bed $dir_base/data/HPRC_SV_gt_10000bp.protein_coding_genes.collapsed.100kb_slop.no-dup.bed | awk -v OFS='\t' '{name=$4; gsub(/[^a-zA-Z0-9._-]/, "-", name); print $1,$2,$3,substr(name, 1, 32),$3-$2}' | sort -k 5n | while read chrom start end name len; do
    echo "Processing region: $chrom:$start-$end ($name, len: $len)"
    
    region=${chrom}_${start}_${end}
    echo -e "$chrom\t$start\t$end\t$name" > $dir_base/genotyping/$region.bed

    threads=48
    if [ $len -gt $LENGTH_THRESHOLD ]; then
        threads=96
    fi

    sbatch -c $threads -p allnodes --job-name $name --wrap "hostname; $dir_base/scripts/genotype_region.sh \
        $dir_base/genotyping/$region.bed \
        $dir_base/wfmash \
        $dir_base/reference/GRCh38.fa.gz \
        $dir_base/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
        $dir_base/genotyping/samples-to-consider.txt \
        /lizardfs/guarracino/pangenomes/HPRCv2 \
        $dir_base/data/HPRCv2/illumina \
        $threads \
        /scratch \
        $dir_base/genotyping/$region"
done

# Prepare for benchmark
cd $dir_base/genotyping
mkdir -p /scratch/hprcv2-benchmark/cosigt
for dir in chr*/cosigt; do
    parent=$(dirname "$dir")
    cp -r "$dir"/* /scratch/hprcv2-benchmark/cosigt
done
mkdir -p /scratch/hprcv2-benchmark/clusters
for dir in chr*/clusters; do
  parent=$(dirname "$dir")
  cp -r "$dir"/* /scratch/hprcv2-benchmark/clusters
done
mkdir -p /scratch/hprcv2-benchmark/odgi/dissimilarity
cp -r chr*/odgi/dissimilarity/* /scratch/hprcv2-benchmark/odgi/dissimilarity
mkdir -p /scratch/hprcv2-benchmark/pggb
for og_file in chr*/pggb/*/*/*.final.og; do
  chr_name=$(echo $og_file | cut -d'/' -f3)
  region=$(echo $og_file | cut -d'/' -f4)
  dest_dir=/scratch/hprcv2-benchmark/pggb/$chr_name/$region
  mkdir -p $dest_dir
  filename=$(basename "$og_file")
  cp $og_file $dest_dir
done
mkdir -p /scratch/hprcv2-benchmark/odgi
for gfa in chr*/odgi/view/*/*.gfa; do
  chr_dir=$(echo $gfa | cut -d'/' -f1)
  chr_name=$(echo $gfa | cut -d'/' -f4)
  dest_dir=/scratch/hprcv2-benchmark/odgi/view/$chr_name
  mkdir -p $dest_dir
  filename=$(basename "$gfa")
  cp $gfa $dest_dir
done

# Benchmark
cd /scratch/hprcv2-benchmark
cat $dir_base/data/loci.bed $dir_base/data/HPRC_SV_gt_10000bp.protein_coding_genes.collapsed.100kb_slop.no-dup.bed | 
  awk -v OFS='\t' '{name=$4; gsub(/[^a-zA-Z0-9._-]/, "-", name); print $1,$2,$3,substr(name, 1, 32),$3-$2}' | 
  sort -k 5n | 
  head -n 665 > tmp.bed
mkdir -p logs
cat tmp.bed | parallel --colsep '\t' -j 44 "$dir_base/scripts/benchmark_region.sh {1} {2} {3} > logs/{1}_{2}_{3}.log 2>&1"
# Final step: Plot TPR results
Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/plot_tpr.r \
    benchmark/tpr \
    tmp.bed \
    benchmark/*/*/tpr_qv.tsv
rm tmp.bed


###########################################################################################################################################
# OK


(echo region haplotype step; find impg/*/*/ -name "*.projected.bedpe" | while read bedpe; do
    sample=$(basename $bedpe .projected.bedpe)
    chrom=$(echo $bedpe | cut -d'/' -f2)
    region=$(echo $bedpe | cut -d'/' -f3)

    echo $region $sample RAW
done;
find impg/*/*/ -name "*.projected.bedpe" | while read bedpe; do
    sample=$(basename $bedpe .projected.filtered.bedpe)
    chrom=$(echo $bedpe | cut -d'/' -f2)
    region=$(echo $bedpe | cut -d'/' -f3)

    echo $region $sample SPAN
done;
find impg/*/*/ -name "*.merged.filtered.bed" | while read bed; do
    sample=$(basename $bed .merged.filtered.bed)
    chrom=$(echo $bedpe | cut -d'/' -f2)
    region=$(echo $bedpe | cut -d'/' -f3)

    echo $region $sample FLAGGER
done) | tr ' ' '\t' > filtering.tsv


#!/bin/bash

(echo -e "region\tnum_samples\tstep"; cat "$dir_base/data/loci.bed" "$dir_base/data/HPRC_SV_gt_10000bp.protein_coding_genes.collapsed.100kb_slop.no-dup.bed" | 
awk -v OFS='\t' '{name=$4; gsub(/[^a-zA-Z0-9._-]/, "-", name); print $1,$2,$3,substr(name, 1, 32),$3-$2}' | 
sort -k 5n | 
while read -r chrom start end name len; do
    region="${chrom}_${start}_${end}"
    
    # Count samples for each step in this specific region
    raw_count=$(find "impg/$chrom/$region/" -name "*.projected.bedpe" 2>/dev/null | wc -l)
    span_count=$(find "impg/$chrom/$region/" -name "*.projected.filtered.bedpe" 2>/dev/null | wc -l)
    flagger_count=$(find "impg/$chrom/$region/" -name "*.merged.filtered.bed" 2>/dev/null | wc -l)
    
    # Output results for this region
    echo -e "$region\t$raw_count\tRAW"
    echo -e "$region\t$span_count\tSPAN"
    echo -e "$region\t$flagger_count\tFLAGGER"
done) > filtering.tsv


chrom=chr11; start=69809692; end=69819692; name=FGF3; len=10000
chrom=chr6; start=31891045; end=32123783; name=C4A,C4B; len=232738

chr6    31891045        32123783        C4A,C4B

mkdir -p /scratch/HPRCv2-gt
cd /scratch/HPRCv2-gt
ls $dir_base/data/HPRCv2/illumina/*.cram | while read f; do echo $(basename $f .final.cram); done > samples-to-consider.txt

cat $dir_base/data/loci.bed $dir_base/data/HPRC_SV_gt_10000bp.protein_coding_genes.collapsed.100kb_slop.no-dup.bed | awk -v OFS='\t' '{name=$4; gsub(/[^a-zA-Z0-9._-]/, "-", name); print $1,$2,$3,substr(name, 1, 32),$3-$2}' | sort -k 5n | while read chrom start end name len; do
    echo "Processing region: $chrom:$start-$end ($name, len: $len)"

    region=${chrom}_${start}_${end}
    mkdir -p regions_of_interest
    echo -e "$chrom\t$start\t$end\t$name" > regions_of_interest/$region.bed

    echo "  Projecting alignments..."
    mkdir -p impg/$chrom/$region
    ls $dir_base/wfmash.s10kp95/*-vs-grch38.aln.paf| while read paf; do
        sample=$(basename $paf -vs-grch38.aln.paf)
        
        impg query \
            -p $paf \
            -b regions_of_interest/$region.bed \
            > impg/$chrom/$region/$sample.projected.bedpe
    done


    echo "  Filtering and merging sequences..."
    # # By length: not sure, because it is okay but not super okay
    # cat impg/$chrom/$region/*.merged.bed > impg/$chrom/$region/ALL.tmp.bed
    # Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/outliers.r \
    #     impg/$chrom/$region/ALL.tmp.bed \
    #     impg/$chrom/$region/ALL.merged.filtered.bed
    # rm impg/$chrom/$region/ALL.tmp.bed

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
        else
            echo "  Filtering $sample for not fully spanning the locus with a single contig"
        fi
    done
    ls impg/$chrom/$region/*.projected.filtered.bedpe | while read bedpe; do
        sample=$(basename $bedpe .projected.filtered.bedpe)

        bedtools sort -i $bedpe | \
            bedtools merge -d 200000 > impg/$chrom/$region/$sample.merged.bed

        # Get the total length of the merged interval
        merged_length=$(awk '{sum += $3-$2} END {print sum}' impg/$chrom/$region/$sample.merged.bed)
        # Get the total length of overlapping regions
        overlap_length=$(bedtools intersect -a impg/$chrom/$region/$sample.merged.bed -b $dir_base/data/HPRCv2/flagger/$sample.bed -wao | awk '{sum += $NF} END {print sum}')
        # Calculate percentage
        percentage=$(echo "scale=2; ($overlap_length / $merged_length) * 100" | bc)
        # Check if percentage is greater than 5%
        if (( $(echo "$percentage <= 5" | bc -l) )); then
            cp impg/$chrom/$region/$sample.merged.bed impg/$chrom/$region/$sample.merged.filtered.bed
        else
            echo "  Filtering $sample for having $percentage% >= 5% of regions flagged by FLAGGER"
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


cat $dir_base/data/loci.bed $dir_base/data/HPRC_SV_gt_10000bp.protein_coding_genes.collapsed.100kb_slop.no-dup.bed | awk -v OFS='\t' '{name=$4; gsub(/[^a-zA-Z0-9._-]/, "-", name); print $1,$2,$3,substr(name, 1, 32),$3-$2}' | sort -k 5n | head -n 150 | while read chrom start end name len; do
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
cat $dir_base/data/loci.bed $dir_base/data/HPRC_SV_gt_10000bp.protein_coding_genes.collapsed.100kb_slop.no-dup.bed > tmp.bed
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
