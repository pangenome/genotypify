# Tools

## Installation on UTHSC cluster

### cosigt (tool)

Install `go`:

```shell
cd /lizardfs/guarracino/tools
wget -c https://go.dev/dl/go1.23.1.linux-amd64.tar.gz
tar -C /lizardfs/guarracino/tools -xzf go1.23.1.linux-amd64.tar.gz && rm go1.23.1.linux-amd64.tar.gz
#Add 'export PATH="/lizardfs/guarracino/tools/go/bin:$PATH"' to ~/.zshrc
```

Build `cosigt`:

```shell
cd /lizardfs/guarracino/git
git clone https://github.com/davidebolo1993/cosigt
cd cosigt
git checkout eb36f56f210be9de9859fbe7902a21879267a94a
go mod init cosigt && go mod tidy && go build cosigt
#Add 'export PATH="/lizardfs/guarracino/git/cosigt:$PATH"' to ~/.zshrc
```

### cosigt (pipeline)

Create a `conda` environment for `cosigt` with all its dependencies:

```shell
conda create --prefix /lizardfs/guarracino/condatools/cosigt -c conda-forge -c bioconda -c anaconda -c vikky34v snakemake=7.32.4 cookiecutter=2.6.0 bwa-mem2=2.2.1 megadepth=1.2.0 python=3.11 pyyaml=6.0.2 pandas=2.2.3 numpy=2.2.1 time -y # r-data.table r-rjson=0.2.23 r-ggplot2=3.5.1 r-gggenes=0.5.1 bioconductor-rtracklayer=1.62.0 r-reshape2=1.4.4 r-nbclust=3.0.1 r-dendextend=1.18.1 r-randomcolor=1.1.0.1 r-base=4.3.3 r-tidyverse=2.0.0 r-dbscan=1.2.2 r-devtools r-pafr r-ggnewscale r-ggforce bioconductor-s4vectors bioconductor-bsgenome bioconductor-genomicranges bioconductor-biostrings bioconductor-iranges bioconductor-genomeinfodb bioconductor-rsamtools r-wesanderson r-abind bioconductor-rtracklayer r-xml r-v8 r-curl -y
```

It assumes that all `pggb` and its tools (`wfmash`, `seqwish`, `smoothxg`, `odgi`, `gfaffix`), `samtools`, `bedtools` are in `$PATH` are installed and included in system's `$PATH` environment variable so they can be executed from any directory.

Prepare the R environment:

```shell
guix install r-dbscan
```

### agc (to obtain assemblies)

```shell
cd /lizardfs/guarracino/tools
curl -L https://github.com/refresh-bio/agc/releases/download/v1.1/agc-1.1_x64-linux.tar.gz | tar -zxvf - agc-1.1_x64-linux/agc
#Add 'export PATH="/lizardfs/guarracino/tools/agc-1.1_x64-linux:$PATH"' to ~/.zshrc
```

### TO IGNORE: MONI

I've added `moni` on `bioconda` (https://github.com/bioconda/bioconda-recipes/pull/50925), so we can just:

```shell
conda create --prefix /lizardfs/guarracino/condatools/moni/0.2.2 -c conda-forge -c bioconda moni=0.2.2 -y
conda create --prefix /lizardfs/guarracino/condatools/bwa/0.7.18 -c conda-forge -c bioconda bwa=0.7.18 -y
```

### bmws (beta-mixture-with-spikes model)

```shell
conda create --prefix /lizardfs/guarracino/condatools/bmws/0.2.1/ -c conda-forge -c bioconda python==3.10 numpy scipy matplotlib #zlib==1.2.13 pandas pip gcc bcftools -y
conda activate /lizardfs/guarracino/condatools/bmws/0.2.1
pip install -e git+https://github.com/jthlab/bmws#egg=bmws
```

## Small test for cosigt

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
sed "/\/ava\.pdf/s/^/#/" workflow/Snakefile -i
sed "/\/pgrtk/s/^/#/" workflow/Snakefile -i
sed "/\/annotations/s/^/#/" workflow/Snakefile -i
sed "/\/untangle/s/^/#/" workflow/Snakefile -i

# OPTIONAL: create a map of alignment name to id
for s in $(ls ../../cram/*.cram); do cram=$(basename $s) && id=$(echo $cram | cut -d "." -f 1) && echo -e "$cram\t$id"; done > sample.map.tsv
# OPTIONAL: add some annotations for the region of interest

cd resources/annotations
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz

sed 's/reg[3]/reg[2]/g' workflow/scripts/annotate.r -i
sed 's/reg[4]/reg[3]/g' workflow/scripts/annotate.r -i


cd ../..

# Create a temporary directory
mkdir -p /scratch/small_test/tmp
```

### Genotyping

On the UTHSC cluster, put the following in `~/.bashrc` (or `~/.zshrc`):

```shell
export PATH=$(echo $PATH | tr ':' '\n' | awk '!(/\/gnu\/store\// || /guix/)' | paste -sd ':') # Remove guix's path
#export PYTHONNOUSERSITE=1 # Tells Python not to use the user site-packages directory and clear the Python path
#export PATH="/home/guarracino/.guix-profile/bin:$PATH"
#export PATH="/scratch/cosigt:$PATH"

export PATH="/lizardfs/guarracino/tools/bedtools2/bin:$PATH"
export PATH="/lizardfs/guarracino/tools/samtools-1.21:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/gafpack/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/GFAffix/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/impg/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/gfainject/target/release:$PATH"

export PATH="/lizardfs/guarracino/tools_for_cosigt/wfmash/build/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/seqwish/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/smoothxg/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/odgi/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_cosigt/pggb:$PATH"

# Added for go and cosigt installation
export PATH="/lizardfs/guarracino/tools/go/bin:$PATH"
export PATH="/lizardfs/guarracino/tools/agc-1.1_x64-linux:$PATH"
export PATH="/lizardfs/guarracino/git/cosigt:$PATH"
```

Prepare the configuration files and run the `cosigt` pipeline:

```shell
# Activate conda environment
conda activate /lizardfs/guarracino/condatools/cosigt

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

# Genotyping
snakemake cosigt --cores 16

conda deactivate
```
