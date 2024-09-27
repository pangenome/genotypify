# Genotypify
Genotyping lots of samples with big pangenomes

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
conda create --prefix /lizardfs/guarracino/condatools/cosigt -c conda-forge -c bioconda -c anaconda snakemake cookiecutter python=3.9 pip pyyaml pandas pyfaidx bwa-mem2 matplotlib Pillow scipy scikit-learn -y
```

It assumes that all `pggb` and its tools (`wfmash`, `seqwish`, `smoothxg`, `odgi`, `gfaffix`), `samtools`, `bedtools` are in `$PATH` are installed and included in system's `$PATH` environment variable so they can be executed from any directory.

### MONI

I've added `moni` on `bioconda`, so we can just:

```shell
conda create --prefix /lizardfs/guarracino/condatools/moni/0.2.2 -c conda-forge -c bioconda moni=0.2.2 -y
```

## Small test

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/genotypify
```

Preparation:

```shell
mkdir -p $DIR_BASE/small_test && cd $DIR_BASE/small_test

# Prepare folders
cd $DIR_BASE/small_test
mkdir cram reference assemblies paf roi

# Region-of-interest
echo -e "grch38#chr6\t31972057\t32055418" > roi/roi.bed

# Cram
cd cram
wget https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3988768/HG00438.final.cram
wget https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3988768/HG00438.final.cram.crai
cd ..

# Reference
cd reference
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
cd ..

# Assemblies
curl -o HPRC-yr1.agc https://zenodo.org/record/5826274/files/HPRC-yr1.agc?download=1
curl -L https://github.com/refresh-bio/agc/releases/download/v1.1/agc-1.1_x64-linux.tar.gz|tar -zxvf - agc-1.1_x64-linux/agc
cd assemblies
while read -r line; do
    ../agc-1.1_x64-linux/agc getctg $DIR_BASE/data/HPRC-yr1.agc $line >> /scratch/chr6.y1.fa
done < ../chr6.y1.txt
samtools faidx ../reference/GRCh38_full_analysis_set_plus_decoy_hla.fa chr6 | sed 's/^>chr6/>grch38#chr6/' > chr6.grhch38.fa
cat chr6.grhch38.fa >> /scratch/chr6.y1.fa
samtools faidx /scratch/chr6.y1.fa
mv /scratch/chr6.y1.fa*$DIR_BASE/small_test/assemblies
cd ..

# Paf
cd paf
wfmash ../assemblies/chr6.grhch38.fa ../assemblies/chr6.y1.fa -X -t 16 -s 10k -p 95 > chr6.y1.s10p95.paf
cd ..
```

Genotyping:

```shell
# Activate conda environment
conda activate /lizardfs/guarracino/condatools/cosigt

# Prepare input files
cd /lizardfs/guarracino/git/cosigt/cosigt_smk
python workflow/scripts/organize.py -a $DIR_BASE/small_test/cram -r $DIR_BASE/small_test/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa --fasta $DIR_BASE/small_test/assemblies/chr6.y1.fa --paf $DIR_BASE/small_test/paf/chr6.y1.s10p95.paf --roi $DIR_BASE/small_test/roi/roi.bed --temp_dir /scratch

# Genotyping
snakemake cosigt --cores 16

conda deactivate
```


