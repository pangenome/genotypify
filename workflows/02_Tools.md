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
guix install r-data-table r-rjson r-ggplot2 r-reshape2 r-nbclust r-dendextend r-tidyverse r-dbscan r-ggnewscale r-ggforce r-wesanderson r-abind r-xml r-v8 r-curl
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
