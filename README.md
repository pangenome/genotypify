# Genotypify
Genotyping lots of samples with big pangenomes

## Variables

```shell
DIR_BASE=/lizardfs/guarracino/genotypify
```

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
conda create --prefix /lizardfs/guarracino/condatools/cosigt -c conda-forge -c bioconda -c anaconda snakemake=7.32.4 cookiecutter=2.6.0 bwa-mem2=2.2.1 samtools=1.21 bedtools=2.31.1 python=3.9 pyyaml=6.0.2 pandas r-base r-rjson=0.2.23 r-reshape2=1.4.4 r-nbclust=3.0.1 r-data.table r-ggplot2=3.5.1 r-dendextend=1.18.1 r-gggenes=0.5.1 bioconductor-rtracklayer time -y
```

It assumes that all `pggb` and its tools (`wfmash`, `seqwish`, `smoothxg`, `odgi`, `gfaffix`), `samtools`, `bedtools` are in `$PATH` are installed and included in system's `$PATH` environment variable so they can be executed from any directory.

### MONI

I've added `moni` on `bioconda`, so we can just:

```shell
conda create --prefix /lizardfs/guarracino/condatools/moni/0.2.2 -c conda-forge -c bioconda moni=0.2.2 -y
```

## Small test

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

# Paf
mkdir paf && cd paf
wfmash ../assemblies/chr6.grhch38.fa ../assemblies/chr6.y1.fa -X -t 16 -s 10k -p 95 > chr6.y1.s10p95.paf
cd ..
```

Genotyping:

```shell
export PATH=$(echo $PATH | tr ':' '\n' | awk '!(/\/gnu\/store\// || /guix/)' | paste -sd ':') # Remove guix's path
export PYTHONNOUSERSITE=1 # Tells Python not to use the user site-packages directory and clear the Python path

# Activate conda environment
conda activate /lizardfs/guarracino/condatools/cosigt
export PATH="/scratch/cosigt:$PATH"

# Prepare input files
cd /scratch/cosigt/cosigt_smk
python workflow/scripts/organize.py -a /scratch/small_test/cram -r /scratch/small_test/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa --assemblies /scratch/small_test/assemblies/chr6.y1.fa --roi /scratch/small_test/roi/roi.bed --pggb_tmpdir /scratch

# Genotyping
snakemake cosigt --cores 16

conda deactivate
```

## Data

### Ancient

15 samples from Marchi et al., 2022 (https://doi.org/10.1016/j.cell.2022.04.008).

```shell
mkdir -p $DIR_BASE/sequencing_data/ancient

# Download the submitted FASTQ files, discarding the XXX_realg.noclips_splitRG.bam files
mkdir -p /scratch/Marchi2022
cd /scratch/Marchi2022
grep -Ff <(sed '1d' $DIR_BASE/data/Marchi2022.TableS1.csv | cut -f 9 | sed 's/SL_MU_/SLMU/g') $DIR_BASE/data/filereport_read_run_PRJEB50857_tsv.txt | awk -F'\t' 'NR>1 {split($8, a, ";"); for (i in a) print "ftp://"a[i]}' | xargs -n 1 wget

cd /scratch
mv /scratch/Marchi2022 $DIR_BASE/sequencing_data/ancient

# Marchi2022.TableS1.csv to map `SAMPLE <-> FILE`
```

315 (317 - 2 samples with missing FASTQ files) samples from Allentoft et al., 2024 (https://doi.org/10.1038/s41586-023-06865-0):

```shell
# Counts of ancient samples they collected from other studies
sed '1d' Allentoft2024.SupplementaryData7.csv | cut -f 8 | sort | uniq -c | sort -k 1,1n

mkdir -p $DIR_BASE/sequencing_data/ancient/Allentoft2024
cd /scratch
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB64656_tsv.txt | cut -f 7 | grep ftp | awk -F'\t' 'NR>1 {split($1, a, ";"); for (i in a) print "ftp://"a[i]}' | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Allentoft2024'

# To map `SAMPLE <-> FILE`
(sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB64656_tsv.txt | cut -f 7,8 | grep 'NEO898\|NEO813\|ERR12075080' -v | awk -v FS='/' '{print($7,$12)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' -e 's/.sort.rmdup.realign.md.bam;//g' -e 's/.neo.merge.bam;//g' -e 's/.neo.clean.bam;//g' | sed 's/ //g'; echo "ERR12075080\tNEO962")
```

442 samples from Margaryan et al., 2020 (https://doi.org/10.1038/s41586-020-2688-8; https://www.ebi.ac.uk/ena/browser/view/PRJEB37976):

```shell
mkdir -p $DIR_BASE/sequencing_data/ancient/Margaryan2020
cd /scratch/
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB37976_tsv.txt | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Margaryan2020'

# To map `SAMPLE <-> FILE`
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB37976_tsv.txt | cut -f 7,8 | awk -v FS='/' '{print($7,$12)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' -e 's/.final.bam//g' | sed 's/ //g'
```

134 samples from Antonio et al., 2019 (https://www.ebi.ac.uk/ena/browser/view/PRJEB32566; https://www.ebi.ac.uk/ena/browser/view/PRJEB32566):

```shell
mkdir -p $DIR_BASE/sequencing_data/ancient/Antonio2019
cd /scratch/
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB32566_tsv.txt | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Antonio2019'

# To map `SAMPLE <-> FILE`
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB32566_tsv.txt | cut -f 7,8 | awk -v FS='/' '{print($7,$12)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' -e 's/.bam;//g' -e 's/.bam//g' | sed 's/ //g'
```

137 samples from Damgaard et al., 2018, Nature (https://doi.org/10.1038/s41586-018-0094-2; https://www.ebi.ac.uk/ena/browser/view/PRJEB20658):

```shell
mkdir -p $DIR_BASE/sequencing_data/ancient/Damgaard2018_Nature
cd /scratch
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB20658_tsv.txt | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Damgaard2018_Nature'

# To map `SAMPLE <-> FILE`
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB20658_tsv.txt | cut -f 7,8 | awk -v FS='/' '{print($7,$12)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' -e 's/.sort.rmdup.realign.md.bam;//g' | sed 's/ //g'
```

70 samples from Damgaard et al., 2018, Science (https://doi.org/10.1126/science.aar7711; https://www.ebi.ac.uk/ena/browser/view/PRJEB25389 and https://www.ebi.ac.uk/ena/browser/view/PRJEB26349):

```shell
# 70 samples, but the paper mentions 74 samples
mkdir -p $DIR_BASE/sequencing_data/ancient/Damgaard2018_Science
cd /scratch
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB26349_tsv.txt | grep -Ff <(sed '1d' $DIR_BASE/data/Damgaard2018.Science.SupplementaryTable1.csv | cut -f 1) | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Damgaard2018_Science'

# The 4 missed samples are in the Damgaard et al., 2018, Nature paper
sed '1d' $DIR_BASE/data/Damgaard2018.Science.SupplementaryTable1.csv | cut -f 1 | sort | grep -Ff <(sed '1d' data/filereport_read_run_PRJEB26349_tsv.txt | grep -Ff <(sed '1d' data/Damgaard2018.Science.SupplementaryTable1.csv | cut -f 1) | cut -f 8 | cut -f 6 -d '/' | cut -f 1 -d '.' | cut -f 1 -d '_' ) -v
# DA343
# DA353
# DA356
# DA361
```

102 samples from Allentoft et al., 2015 (https://doi.org/10.1038/nature14507; https://www.ebi.ac.uk/ena/browser/view/PRJEB9021).

```shell
mkdir -p $DIR_BASE/sequencing_data/ancient/Allentoft2015
cd /scratch
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB9021_tsv.txt | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Allentoft2015'

# To map `SAMPLE <-> FILE`
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB9021_tsv.txt | cut -f 7,8 | awk -v FS='/' '{print($6,$11)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' | sed -e 's/.hg19.flt.sort.rmdup.realign.md.bam//g' | sed 's/ //g'
```

58 samples from Brunel et al., 2020 (https://doi.org/10.1073/pnas.1918034117; https://www.ebi.ac.uk/ena/browser/view/PRJEB38152)

```shell
#sed '1d' filereport_read_run_PRJEB38152_tsv.txt | grep -i -Ff <(sed '1,2d' Brunel2018.SupplementaryTable1.csv | cut -f 1 | sed -e 's/Schw72-15/Sch72-15/g' -e 's/CRE8C/CRE8-C/g' -e 's/BIS388-/BIS388/g'; echo CRE20D)
mkdir -p $DIR_BASE/sequencing_data/ancient/Brunel2020
cd /scratch
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB38152_tsv.txt | grep bam | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Brunel2020'

# To map `SAMPLE <-> FILE`
sed '1d' filereport_read_run_PRJEB38152_tsv.txt | cut -f 7,8 | grep bam | awk -v FS='/' '{print($7,$12)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' | sed -e 's/.mappedHs37d5.trim2.bam//g' -e 's/.mappedHs37d5.trim2.bam//g' -e 's/.mappedHs37.trim2.bam//g' | sed 's/ //g'
```

91??? samples frm Scheib et al., 2018 (https://doi.org/10.1126/science.aar6851; https://www.ebi.ac.uk/ena/browser/view/PRJEB25445)

```shell
mkdir -p $DIR_BASE/sequencing_data/ancient/Scheib2018
cd /scratch
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB25445_tsv.txt | cut -f 7 | grep ftp | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Scheib2018'

# Check Scheib2018.SupplementaryTable1.tsv and Scheib2018.SupplementaryTable3.tsv (number of reads) to process downloaded samples
```

34 samples from Sikora et al.m 2019 (https://doi.org/10.1038/s41586-019-1279-z; https://www.ebi.ac.uk/ena/browser/view/PRJEB29700 and https://www.ebi.ac.uk/ena/browser/view/PRJEB26336):

```shell
mkdir -p $DIR_BASE/sequencing_data/ancient/Sikora2019
cd /scratch
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB29700_tsv.txt | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Sikora2019'

# To map `SAMPLE <-> FILE`
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB29700_tsv.txt | cut -f 7,8 | awk -v FS='/' '{print($7,$12)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' -e 's/.sort.rmdup.realign.md.bam;//g' -e 's/.realigned.calmd.readsadded.bam;//g' -e 's/.sort.rmdup.uniq.rg.realn.md.bam;//g' | sed 's/ //g'
```

35 samples from Krzewinska et al., 2018 (https://doi.org/10.1126/sciadv.aat4457; https://www.ebi.ac.uk/ena/browser/view/PRJEB27628):

```shell
mkdir -p $DIR_BASE/sequencing_data/ancient/Krzewinska2018
cd /scratch
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB27628_tsv.txt | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Krzewinska2018'

# To map `SAMPLE <-> FILE`
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB27628_tsv.txt | cut -f 7,8 | awk -v FS='/' '{print($7,$12)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' | cut -f 1 -d '_' | sed 's/ //g'
```

25 samples from McColl et al., 2018 (https://doi.org/10.1126/science.aat3628; https://www.ebi.ac.uk/ena/browser/view/PRJEB26721)

```shell
mkdir -p $DIR_BASE/sequencing_data/ancient/McColl2018
cd /scratch
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB26721_tsv.txt | cut -f 7 | awk -F'\t' 'NR>1 {split($1, a, ";"); for (i in a) print "ftp://"a[i]}' | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/McColl2018'
```

Ebenesersdottir et al.m 2018, Science (https://doi.org/10.1126/science.aar2625): I NEED ACCESS TO THE PUBLICATION

15/24?? samples from Schroeder et al., 2019 (https://doi.org/10.1073/pnas.1820210116; https://www.ebi.ac.uk/ena/browser/view/PRJEB28451):

```shell
mkdir -p $DIR_BASE/sequencing_data/ancient/Schroeder2019
cd /scratch
sed '1d' $DIR_BASE/data/filereport_read_run_PRJEB28451_tsv.txt | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Schroeder2019'
```

### Modern

1000 Genomes Project sample collection to 30x coverage (from https://www.internationalgenome.org/data-portal/data-collection/30x-grch38). Initially, the 2504 unrelated samples from the phase three panel from the 1000 Genomes Project were sequenced. Thereafter, an additional 698 samples, related to samples in the 2504 panel, were also sequenced.

```shell
mkdir -p $DIR_BASE/sequencing_data/modern/1000G/2504_high_coverage
cd $DIR_BASE/sequencing_data/modern/1000G/2504_high_coverage
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index
awk '!/^#/ {print $1}' 1000G_2504_high_coverage.sequence.index | xargs -n 1 wget

mkdir -p $DIR_BASE/sequencing_data/modern/1000G/additional_698_related
cd $DIR_BASE/sequencing_data/modern/1000G/additional_698_related
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index
awk '!/^#/ {print $1}' 1000G_698_related_high_coverage.sequence.index | xargs -n 1 wget
```

Simons Genome Diversity Project (from https://www.ebi.ac.uk/ena/browser/view/PRJEB9586):

```shell
mkdir -p $DIR_BASE/sequencing_data/modern/SGDP/
cd $DIR_BASE/sequencing_data/modern/SGDP/
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/hgdp_wgs.sequence.index
awk -F'\t' 'NR>1 {split($7, a, ";"); for (i in a) print "ftp://"a[i]}' $DIR_BASE/data/filereport_read_run_PRJEB9586_tsv.txt | xargs -n 1 wget
```
