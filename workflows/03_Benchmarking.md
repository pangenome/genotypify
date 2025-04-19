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
mkdir -p $dir_base/genotyping
ls $dir_base/data/HPRCv2/illumina/*.cram | while read f; do echo $(basename $f .final.cram); done | head -n 100 > $dir_base/genotyping/samples-to-consider.txt
mkdir -p $dir_base/regions_of_interest
# Project the aligments to the region-of-interest (ROI)
# Sort by length to start with the shortest ones
cat $dir_base/data/loci.bed $dir_base/data/HPRC_SV_gt_10000bp.protein_coding_genes.collapsed.100kb_slop.no-dup.bed | awk -v OFS='\t' '{name=$4; gsub(/[^a-zA-Z0-9._-]/, "-", name); print $1,$2,$3,substr(name, 1, 32),$3-$2}' | sort -k 5n | while read chrom start end name len; do
    echo "Processing region: $chrom:$start-$end ($name, len: $len)"

    echo -e "$chrom\t$start\t$end\t$name" > $dir_base/regions_of_interest/${chrom}_${start}_${end}_${name}.bed

    sbatch -c 48 -p workers --job-name $name --wrap "hostname; $dir_base/scripts/genotype_region.sh \
        $dir_base/regions_of_interest/${chrom}_${start}_${end}_${name}.bed \
        /$dir_base/wfmash \
        $dir_base/reference/GRCh38.fa.gz \
        $dir_base/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
        $dir_base/genotyping/samples-to-consider.txt \
        /lizardfs/guarracino/pangenomes/HPRCv2 \
        $dir_base/data/HPRCv2/illumina \
        48 \
        /scratch \
        $dir_base/genotyping/${chrom}_${start}_${end}_${name}"
done



# OK
mkdir -p /scratch/HPRCv2-gt
cd /scratch/HPRCv2-gt

ls $dir_base/regions_of_interest | sed 's/.bed//g' | while read region; do
    echo "Processing $region..."


    echo "  Projecting alignments..."
    mkdir -p impg/$region
    ls $dir_base/wfmash/*-vs-grch38.aln.paf | grep -Ff $dir_base/data/HPRCv2/illumina/samples.txt | while read paf; do
        sample=$(basename $paf -vs-grch38.aln.paf)
        
        impg query \
            -p $paf \
            -b $dir_base/regions_of_interest/$region.bed \
            > impg/$region/$sample.projected.bedpe

        bedtools sort -i impg/$region/$sample.projected.bedpe | \
            bedtools merge -d 100000 > impg/$region/$sample.merged.bed
    done


    echo "  Filtering sequences..."
    # QC TO DOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    # - BY LENGTH: NOT SURE, BECAUSE IT IS OKAY, BUT NOT SUPER OKAY
    # - FULLY SPAN THE LOCUS (5% on both sides?)
    # - FLAGGER-FLAGGED REGIONS (max 5% of bad regions?)
    # from impg/$region/$sample.merged.bed to impg/$region/$sample.merged.filtered.bed
    # QC TO DOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO


    echo "  Collecting sequences..."
    bedtools getfasta -fi $dir_base/reference/GRCh38.fa.gz -bed $dir_base/regions_of_interest/$region.bed | sed 's/^>chr/>GRCh38#0#chr/g' > impg/$region/$region.pangenome.fa

    ls impg/$region/*.merged.bed | while read bed; do
        sample=$(basename $bed .merged.bed)
        fasta=$dir_pangenome/$sample.fa.gz

        bedtools getfasta -fi $fasta -bed impg/$region/$sample.merged.bed
    done >> impg/$region/$region.pangenome.fa
    bgzip -@ 48 impg/$region/$region.pangenome.fa
    samtools faidx impg/$region/$region.pangenome.fa.gz


    echo "  Pangenome graph building..."
    mkdir -p pggb
    pggb -i impg/$region/$region.pangenome.fa.gz -o pggb/$region -t 48 -D /scratch -c 2
    mv pggb/$region/*smooth.final.og pggb/$region/$region.final.og # Rename the final ODGI graph


    echo "  Getting the path coverage matrix..."
    mkdir -p odgi
    odgi chop -i pggb/$region/$region.final.og -c 32 -o odgi/$region.chopped.og
    odgi paths -i odgi/$region.chopped.og -H | \
        cut -f 1,4- | \
        gzip > odgi/$region.paths_matrix.tsv.gz


    echo "  Pangenome clustering..."
    mkdir -p odgi/dissimilarity
    odgi similarity -i odgi/$region.chopped.og -d > odgi/dissimilarity/$region.tsv
    mkdir -p clusters
    Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/cluster.r odgi/dissimilarity/$region.tsv clusters/$region.clusters.json automatic 0


    echo "  Pangenome indexing..."
    bwa-mem2.avx index impg/$region/$region.pangenome.fa.gz # "bwa-mem2 index" returns "ERROR: prefix is too long!"

    echo "  Alignment, injection, genotyping..." 
    odgi view -i odgi/$region.chopped.og -g > odgi/$region.chopped.gfa
    mkdir -p alignments/$region
    ls $dir_reads/*cram | grep -Ff $dir_base/data/HPRCv2/illumina/samples.txt | while read cram; do
        echo $cram

        name=$(basename $cram .cram)

        samtools view \
            -T $dir_base/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
            -L $dir_base/regions_of_interest/$region.bed \
            -M \
            -b \
            $cram | \
            samtools sort -n | \
            samtools fasta | \
            bwa-mem2.avx mem -t 42 impg/$region/$region.pangenome.fa.gz - | \
            samtools view -b -F 4 -@ 6 - \
            > alignments/$region/$name.reads-vs-pangenome.bam

        gfainject \
            --gfa odgi/$region.chopped.gfa \
            --bam alignments/$region/$name.reads-vs-pangenome.bam \
            > alignments/$region/$name.reads-vs-graph.gaf

        gafpack \
            --gfa odgi/$region.chopped.gfa \
            --gaf alignments/$region/$name.reads-vs-graph.gaf \
            --len-scale | \
            gzip > alignments/$region/$name.coverage.gafpack.gz

        mkdir -p cosigt/$name
        cosigt \
            -i $name \
            -p odgi/$region.paths_matrix.tsv.gz \
            -g alignments/$region/$name.coverage.gafpack.gz \
            -c clusters/$region.clusters.json \
            -o cosigt/$name/$region
    done
done

mkdir -p benchmark
Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/plottpr.r cosigt clusters odgi/dissimilarity/ benchmark/benchmark.pdf

# OPTIONAL
mkdir -p annotation
bedtools intersect -wa -a $dir_base/annotation/gencode.v47.annotation.gtf.gz -b $dir_base/regions_of_interest/$region.bed > annotation/${region}_${chrom}_${start}_${end}.gtf
Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/annotate.r annotation/${region}_${chrom}_${start}_${end}.gtf chr6 annotation/${region}_${chrom}_${start}_${end}.bed
odgi procbed -i pggb/$region/$region.final.og -b <(sed 's/^chr/GRCh38#0#chr/g' annotation/${region}_${chrom}_${start}_${end}.bed) > annotation/${region}_${chrom}_${start}_${end}.procbed.bed
odgi inject -i pggb/$region/$region.final.og -b annotation/${region}_${chrom}_${start}_${end}.procbed.bed -o annotation/$region.final.inject.og -P
odgi flip -i annotation/$region.final.inject.og --ref-flips <(cut -f 1 annotation/${region}_${chrom}_${start}_${end}.procbed.bed | uniq) -o annotation/$region.final.inject.flip.og
odgi untangle -i annotation/$region.final.inject.flip.og -R <(cut -f 4 annotation/${region}_${chrom}_${start}_${end}.procbed.bed) -j 0.5 -g -t 8 > annotation/$region.gggenes.tsv
Rscript /lizardfs/guarracino/tools_for_genotyping/cosigt/cosigt_smk/workflow/scripts/plotgggenes.r annotation/$region.gggenes.tsv clusters/$region.clusters.json annotation/$region.gggenes.png
```
