# Genotyping (TO UPDATE)

## Small test for cosigt on UTHSC cluster

### Preparation

```shell
mkdir -p /scratch/small_test && cd /scratch/small_test

# Region-of-interest
mkdir roi && echo -e "grch38#chr6\t31972057\t32055418" > roi/roi.bed

# Cram
mkdir cram && cd cram
wget https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3988768/HG00438.final.cram
wget https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3988768/HG00438.final.cram.crai
cd ..

# Reference
mkdir reference && cd reference
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
cd ..

# Assemblies
mkdir assemblies && cd assemblies
wget -O HPRC-yr1.agc "https://zenodo.org/record/5826274/files/HPRC-yr1.agc?download=1"
rm chr6.y1.fa
while read -r line; do
    agc getctg HPRC-yr1.agc $line >> chr6.y1.fa
done < $dir_base/data/chr6.y1.txt
samtools faidx ../reference/GRCh38_full_analysis_set_plus_decoy_hla.fa chr6 | sed 's/^>chr6/>grch38#chr6/' > chr6.grhch38.fa
cat chr6.grhch38.fa >> chr6.y1.fa
samtools faidx chr6.y1.fa
cd ..

# Prepare pipeline and two OPTIONAL input files
git clone https://github.com/davidebolo1993/cosigt
cd cosigt
git checkout eb36f56f210be9de9859fbe7902a21879267a94a
cd cosigt_smk

# Disable stuff we don't care about for now
sed "/\/ava\.pdf/s/^/#/" workflow/Snakefile -i
sed "/\/pgrtk/s/^/#/" workflow/Snakefile -i
sed "/annotations/s/^/#/" workflow/Snakefile -i
sed "/\/untangle/s/^/#/" workflow/Snakefile -i

# Fix current issues
sed 's/reg[3]/reg[2]/g' workflow/scripts/annotate.r -i
sed 's/reg[4]/reg[3]/g' workflow/scripts/annotate.r -i

# OPTIONAL: create a map of alignment name to id
for s in $(ls ../../cram/*.cram); do cram=$(basename $s) && id=$(echo $cram | cut -d "." -f 1) && echo -e "$cram\t$id"; done > sample.map.tsv

# OPTIONAL: add some annotations for the region of interest
cd resources/annotations
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
cd ../..

# Create a temporary directory
mkdir -p /scratch/small_test/tmp

# Prepare the configuration files
conda activate /lizardfs/guarracino/condatools/cosigt

cd /scratch/small_test/cosigt/cosigt_smk
python workflow/scripts/organize.py \
    -a /scratch/small_test/cram \
    -r /scratch/small_test/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --assemblies /scratch/small_test/assemblies/chr6.y1.fa \
    --roi /scratch/small_test/roi/roi.bed \
    --wfmash_tmpdir /scratch/small_test/tmp \
    --pggb_tmpdir /scratch/small_test/tmp \
    --output /scratch/small_test/output \
    --samplemap /scratch/small_test/cosigt/cosigt_smk/sample.map.tsv \
    --annotation /scratch/small_test/cosigt/cosigt_smk/resources/annotations/gencode.v47.annotation.gtf.gz

conda deactivate
```

### Genotyping

On the UTHSC cluster, put the following in `~/.bashrc` (or `~/.zshrc`):

```shell
export PATH=$(echo $PATH | tr ':' '\n' | awk '!(/\/gnu\/store\// || /guix/)' | paste -sd ':') # Remove guix's path

export GUIX_PROFILE="/home/guarracino/.guix-profile"
 . "$GUIX_PROFILE/etc/profile"
#export PYTHONNOUSERSITE=1 # Tells Python not to use the user site-packages directory and clear the Python path
#export PATH="/home/guarracino/.guix-profile/bin:$PATH"
#export GUIX_LOCPATH="/home/guarracino/.guix-profile/lib/locale"

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
export PATH="/lizardfs/guarracino/tools_for_genotyping//lizardfs/guarracino/tools_for_genotyping/bwa-mem2-2.2.1_x64-linux:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/cosigt:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/agc-1.1_x64-linux:$PATH"
```

Run the `cosigt` pipeline:

```shell
# Activate conda environment
conda activate /lizardfs/guarracino/condatools/cosigt

cd /scratch/small_test/cosigt/cosigt_smk
snakemake cosigt --cores 8

conda deactivate
```

Results will be in `/scratch/small_test/output`.
