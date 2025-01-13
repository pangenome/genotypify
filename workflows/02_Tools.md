# Tools

## Installation on UTHSC cluster

### cosigt

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
go mod init cosigt && go mod tidy && go build cosigt
#Add 'export PATH="/lizardfs/guarracino/git/cosigt:$PATH"' to ~/.zshrc
```

Create a `conda` environment for `cosigt` with all its dependencies:

```shell
conda create --prefix /lizardfs/guarracino/condatools/cosigt -c conda-forge -c bioconda -c anaconda -c vikky34v snakemake=7.32.4 cookiecutter=2.6.0 bwa-mem2=2.2.1 megadepth samtools=1.21 bedtools=2.31.1 python=3.9 pyyaml=6.0.2 pandas -y #r-base r-rjson=0.2.23 r-reshape2=1.4.4 r-nbclust=3.0.1 r-data.table r-ggplot2=3.5.1 r-dendextend=1.18.1 r-gggenes=0.5.1 bioconductor-rtracklayer time -y
```

It assumes that all `pggb` and its tools (`wfmash`, `seqwish`, `smoothxg`, `odgi`, `gfaffix`), `samtools`, `bedtools` are in `$PATH` are installed and included in system's `$PATH` environment variable so they can be executed from any directory.

### MONI

I've added `moni` on `bioconda` (https://github.com/bioconda/bioconda-recipes/pull/50925), so we can just:

```shell
conda create --prefix /lizardfs/guarracino/condatools/moni/0.2.2 -c conda-forge -c bioconda moni=0.2.2 -y
conda create --prefix /lizardfs/guarracino/condatools/bwa/0.7.18 -c conda-forge -c bioconda bwa=0.7.18 -y
```

### bmws (beta-mixture-with-spikes model)

```shell
conda create --prefix /lizardfs/guarracino/condatools/bmws/0.2.1/ -c conda-forge -c bioconda python==3.10 #zlib==1.2.13 pandas numpy pip gcc bcftools -y
conda activate /lizardfs/guarracino/condatools/bmws/0.2.1
pip install -e git+https://github.com/jthlab/bmws#egg=bmws
```

## Small test for cosigt

Preparation:

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
curl -L https://github.com/refresh-bio/agc/releases/download/v1.1/agc-1.1_x64-linux.tar.gz | tar -zxvf - agc-1.1_x64-linux/agc
mkdir assemblies && cd assemblies
wget -O HPRC-yr1.agc "https://zenodo.org/record/5826274/files/HPRC-yr1.agc?download=1"
rm chr6.y1.fa
while read -r line; do
    ../agc-1.1_x64-linux/agc getctg HPRC-yr1.agc $line >> chr6.y1.fa
done < $DIR_BASE/data/chr6.y1.txt
samtools faidx ../reference/GRCh38_full_analysis_set_plus_decoy_hla.fa chr6 | sed 's/^>chr6/>grch38#chr6/' > chr6.grhch38.fa
cat chr6.grhch38.fa >> chr6.y1.fa
samtools faidx chr6.y1.fa
cd ..

# PAF
mkdir paf && cd paf
wfmash ../assemblies/chr6.grhch38.fa ../assemblies/chr6.y1.fa -X -t 16 -s 10k -p 95 > chr6.y1.s10p95.paf
cd ..
```

Genotyping:

```shell
export PATH=$(echo $PATH | tr ':' '\n' | awk '!(/\/gnu\/store\// || /guix/)' | paste -sd ':') # Remove guix's path
export PYTHONNOUSERSITE=1 # Tells Python not to use the user site-packages directory and clear the Python path
export PATH="/home/guarracino/.guix-profile/bin:$PATH"
export PATH="/scratch/cosigt:$PATH"

# Activate conda environment
conda activate /lizardfs/guarracino/condatools/cosigt

# Prepare input files
cd /scratch/cosigt/cosigt_smk
python workflow/scripts/organize.py -a /scratch/small_test/cram -r /scratch/small_test/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa --assemblies /scratch/small_test/assemblies/chr6.y1.fa --roi /scratch/small_test/roi/roi.bed --pggb_tmpdir /scratch

# Genotyping
snakemake cosigt --cores 16

conda deactivate
```
