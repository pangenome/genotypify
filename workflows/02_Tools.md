# Tools (UTHSC cluster)

Create folder:

```shell
mkdir -p /lizardfs/guarracino/tools_for_genotyping
```

## impg

```shell
cd /lizardfs/guarracino/tools_for_genotyping

git clone https://github.com/pangenome/impg \
    && cd GFAffix \
    && git pull \
    && git checkout 393d6bb87496381fe6e45af77088958d7731f604 \
    && cargo build --release
cd ..
```

## PGGB and its tools (wfmash, seqwish, smoothxg, GFAffix, odgi)

```shell
cd /lizardfs/guarracino/tools_for_genotyping

git clone --recursive https://github.com/waveygang/wfmash \
    && cd wfmash \
    && git pull \
    && git checkout v0.13.1 \
    && git submodule update --init --recursive \
    && cmake -H. -DCMAKE_BUILD_TYPE=Release -Bbuild && cmake --build build -- -j $(nproc)
cd ..

git clone --recursive https://github.com/ekg/seqwish \
    && cd seqwish \
    && git pull \
    && git checkout 0eb6468be0814ab5a0cda10d12aa38cb87d086f1 \
    && git submodule update --init --recursive \
    && cmake -H. -DCMAKE_BUILD_TYPE=Release -Bbuild && cmake --build build -- -j $(nproc)
cd ..

git clone --recursive https://github.com/pangenome/smoothxg \
    && cd smoothxg \
    && git pull \
    && git checkout v0.8.1 \
    && git submodule update --init --recursive \
    && cmake -H. -DCMAKE_BUILD_TYPE=Release -Bbuild && cmake --build build -- -j $(nproc)
cd ..

git clone https://github.com/marschall-lab/GFAffix \
    && cd GFAffix \
    && git pull \
    && git checkout v0.2.1 \
    && cargo build --release
cd ..

git clone --recursive https://github.com/pangenome/odgi \
    && cd odgi \
    && git pull \
    && git checkout v0.9.1 \
    && git submodule update --init --recursive \
    && cmake -H. -DCMAKE_BUILD_TYPE=Release -Bbuild && cmake --build build -- -j $(nproc)
cd ..

git clone https://github.com/pangenome/pggb.git
cd pggb
git checkout cfcfba4f27042d414a30cc1e8fabf8f958aa838a
cd ..
```

## gfainject, gafpack

```shell
cd /lizardfs/guarracino/tools_for_genotyping

git clone https://github.com/AndreaGuarracino/gfainject \
    && cd gfainject \
    && git pull \
    && git checkout 48c98e01b39d0cec71c9130898cbd8fa95bdac80 \
    && cargo build --release

git clone https://github.com/pangenome/gafpack \
    && cd gfainject \
    && git pull \
    && git checkout 72d6d7b72ccbb7dd519fcf29ae00d381d1571464 \
    && cargo build --release
```

## BWA-MEM2

```shell
cd /lizardfs/guarracino/tools_for_genotyping

curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 \
  | tar jxf -
```

## emboss (for stretcher)

```shell
conda create --prefix /lizardfs/guarracino/condatools/emboss/6.6.0 -c conda-forge -c bioconda emboss=6.6.0 -y
```

## agc (to obtain assemblies)

```shell
cd /lizardfs/guarracino/tools_for_genotyping

curl -L https://github.com/refresh-bio/agc/releases/download/v1.1/agc-1.1_x64-linux.tar.gz | tar -zxvf - agc-1.1_x64-linux/agc
```

## cosigt (tool)

Download `go`:

```shell
cd /lizardfs/guarracino/tools_for_genotyping

wget -c https://go.dev/dl/go1.24.2.linux-amd64.tar.gz
tar -xzf go1.24.2.linux-amd64.tar.gz && rm go1.24.2.linux-amd64.tar.gz
```

Build `cosigt`:

```shell
export PATH="/lizardfs/guarracino/tools_for_genotyping/go/bin:$PATH"

cd /lizardfs/guarracino/tools_for_genotyping

git clone https://github.com/davidebolo1993/cosigt
cd cosigt
git checkout 582aa15b152a5963c88bac4bd602b0710b89307c
go mod init cosigt && go mod tidy && go build cosigt
```

## cosigt (pipeline)

Create a `conda` environment for `cosigt` with all its dependencies:

```shell
conda create --prefix /lizardfs/guarracino/condatools/cosigt -c conda-forge -c bioconda -c anaconda -c vikky34v snakemake=7.32.4 cookiecutter=2.6.0 megadepth=1.2.0 python=3.11 pyyaml=6.0.2 pandas=2.2.3 numpy=2.2.1 time -y # r-data.table r-rjson=0.2.23 r-ggplot2=3.5.1 r-gggenes=0.5.1 bioconductor-rtracklayer=1.62.0 r-reshape2=1.4.4 r-nbclust=3.0.1 r-dendextend=1.18.1 r-randomcolor=1.1.0.1 r-base=4.3.3 r-tidyverse=2.0.0 r-dbscan=1.2.2 r-devtools r-pafr r-ggnewscale r-ggforce bioconductor-s4vectors bioconductor-bsgenome bioconductor-genomicranges bioconductor-biostrings bioconductor-iranges bioconductor-genomeinfodb bioconductor-rsamtools r-wesanderson r-abind bioconductor-rtracklayer r-xml r-v8 r-curl -y
```

Prepare the R environment:

```shell
guix install r-data-table r-rjson r-ggplot2 r-reshape2 r-nbclust r-dendextend r-tidyverse r-dbscan r-ggnewscale r-ggforce r-wesanderson r-abind r-xml r-v8 r-curl
```

<!-- ### TO IGNORE: MONI

I've added `moni` on `bioconda` (https://github.com/bioconda/bioconda-recipes/pull/50925), so we can just:

```shell
conda create --prefix /lizardfs/guarracino/condatools/moni/0.2.2 -c conda-forge -c bioconda moni=0.2.2 -y
conda create --prefix /lizardfs/guarracino/condatools/bwa/0.7.18 -c conda-forge -c bioconda bwa=0.7.18 -y
``` -->

## bmws (beta-mixture-with-spikes model)

```shell
conda create --prefix /lizardfs/guarracino/condatools/bmws/0.2.1/ -c conda-forge -c bioconda python==3.10 numpy scipy matplotlib #zlib==1.2.13 pandas pip gcc bcftools -y
conda activate /lizardfs/guarracino/condatools/bmws/0.2.1
pip install -e git+https://github.com/jthlab/bmws#egg=bmws
```
