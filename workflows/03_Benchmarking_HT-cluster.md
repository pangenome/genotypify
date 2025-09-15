# Benchmarking

## Paths

- Genome Analysis Unit guide: https://rstudio-cn.fht.org/GAU_user_guide/
- Project folder: `/project/ham`
- HPRCv2 assemblies: `/processing_data/reference_datasets/HPRC/2_0.6.1`

## Little test

```shell
conda activate smk7324app
cd /group/soranzo/andrea.guarracino/cosigt/cosigt_smk
python workflow/scripts/organize.py \
    -a /group/soranzo/andrea.guarracino/chr_map.tsv \
    -g /project/ham/cosigt_paper/real_data/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -r /project/ham/cosigt_paper/real_data/data/1000G/ \
    -b /group/soranzo/andrea.guarracino/loci.bed \
    -o /scratch/andrea.guarracino/uffa/ \
    --map /group/soranzo/andrea.guarracino/aln_map.tsv \
    --tmp /scratch/andrea.guarracino/ \
    --profile config/slurm/
sh cosigt_smk.sh
#    #--gtf /group/soranzo/chiara.paleni/references/gencode.v47.annotation.gtf.gz --proteins /group/soranzo/chiara.paleni/references/gencode.v47.pc_translations.fa

sh cosigt_smk.sh #cosigt|benchmark|refine
```

## 200 HPRCv2's short read samples on regions from the same 200 HPRCv2's assemblies (the full pangenome is 466 haplotypes)

### Partitioning by chromosome

Download chromosome annotations:

```shell
mkdir -p /project/ham/HPRCv2/by-chromosome/chromAlias
cd /project/ham/HPRCv2/by-chromosome/chromAlias

wget -c https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/chrom_assignment/chrom_alias_pre_release_v0.1.index.csv # 2025/08/13

#######################################
# download.chromAlias.sh
########################
#!/bin/sh
cd /project/ham/HPRCv2/by-chromosome/chromAlias
cut -f 4 chrom_alias_pre_release_v0.1.index.csv -d ',' | sed 1d | while read url; do
    aws s3 --no-sign-request cp "$url" .
done
#######################################

conda activate /ssu/gassu/conda_envs/awsclienv_12586 # for aws-cli
sbatch -p cpuq -c 1 -t 0-02:00:00 download.chromAlias.sh # 2 hours

# Rename TXT files to match the FASTA file names
mv HG06807_mat_v1.chromAlias.txt HG06807_mat_v1.0.genbank.chromAlias.txt
mv HG06807_pat_v1.chromAlias.txt HG06807_pat_v1.0.genbank.chromAlias.txt
```

Split contigs by chromosome:

```shell
cd /project/ham/HPRCv2/by-chromosome

# TO DELETE
# mkdir -p /project/ham/HPRCv2/by-chromosome/tmp
# cd /project/ham/HPRCv2/by-chromosome/tmp
# module load htslib-tools samtools/1.18 # bgzip and samtools
# sbatch -p cpuq -c 1 -t 0-02:00:00 --wrap "cp /processing_data/reference_datasets/HPRC/2_0.6.1/*.fa /project/ham/HPRCv2/by-chromosome/tmp; cd /project/ham/HPRCv2/by-chromosome/tmp; ls *.fa | while read f; do samtools faidx $f; done"
sbatch -p cpuq -c 1 -t 0-02:00:00 --wrap "cd /project/ham/HPRCv2/by-chromosome/tmp; ls *.fa | while read f; do echo $f; samtools faidx $f; done"

#######################################
# partition-by-chromosome.sh
############################
#!/bin/sh
cd /project/ham/HPRCv2/by-chromosome
(seq 1 22; echo X; echo Y; echo M) | while read c; do
    chr=chr$c
    echo $chr

    mkdir $chr
    #ls /processing_data/reference_datasets/HPRC/2_0.6.1/*.fa.gz | while read f; do
    #    sample=$(basename $f .fa)
    ls /project/ham/HPRCv2/by-chromosome/tmp/*.fa | while read f; do
        sample=$(basename $f .fa)

        if [[ $f == *"chm13v2.0"* ]]; then
            grep $chr $f.fai -w
        elif [[ $f == *"hg002v1.1"* ]]; then
            grep $chr $f.fai -w
        elif [[ $f == *"GCA_000001405.15"* ]]; then
            grep $chr $f.fai -w
        else
            grep -E "${chr}[[:blank:]]|${chr}_" /project/ham/HPRCv2/by-chromosome/chromAlias/$sample.chromAlias.txt
        fi | cut -f 1 > contigs.$chr.txt

        samtools faidx $f $(cat contigs.$chr.txt)
    done | bgzip -l 9 -@ 24 > $chr/$chr.fa.gz
    samtools faidx $chr/$chr.fa.gz
    rm contigs.$chr.txt
done
#######################################

module load htslib-tools samtools/1.18 # bgzip and samtools
sbatch -p cpuq -c 24 -t 2-00:00:00 partition-by-chromosome.sh # 2 days
```

Prepare chromosome map:

```shell
ls /project/ham/HPRCv2/by-chromosome/chr*/chr*.fa.gz | while read f; do
    chr=$(basename $f .fa.gz)
    echo -e "$chr\t$f"
done | sort -V > /group/soranzo/andrea.guarracino/chr_map.hprcv2.tsv
```

Prepare sample map:

```shell
ls /project/ham/HPRCv2/illumina/*.final.cram | while read f; do
    sample=$(basename $f .final.cram)
    echo -e "$f\t$sample"
done > /group/soranzo/andrea.guarracino/sample_map.hprcv2.tsv
```

Prepare loci:

```shell
cut -f 1-4 /group/soranzo/andrea.guarracino/HGSVC3.GRCh38.sv.slop50kb.merge50kb.anno1Mbp.bed > /group/soranzo/andrea.guarracino/loci.bed
```

Refine:

```shell
mkdir -p /group/soranzo/andrea.guarracino/hprcv2-refinement/
cd /group/soranzo/andrea.guarracino/hprcv2-refinement/

split -l 50 /group/soranzo/andrea.guarracino/loci.bed # to avoid too many jobs

conda activate smk7324app
cd /group/soranzo/andrea.guarracino/cosigt/cosigt_smk

rm -rf resources

python workflow/scripts/organize.py \
    -a /group/soranzo/andrea.guarracino/chr_map.hprcv2.tsv \
    -g /project/ham/cosigt_paper/real_data/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -r /project/ham/HPRCv2/illumina/ \
    -b /group/soranzo/andrea.guarracino/hprcv2-refinement/xaa \
    -o /scratch/andrea.guarracino/hprcv2-refinement/ \
    --map /group/soranzo/andrea.guarracino/sample_map.hprcv2.tsv \
    --tmp /scratch/andrea.guarracino/ \
    --profile config/slurm/
sed 's/cosigt/refine/g' cosigt_smk.sh -i
sh cosigt_smk.sh
```