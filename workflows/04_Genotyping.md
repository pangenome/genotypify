# Genotyping

## Paths

```shell
dir_base=/lizardfs/guarracino/genotypify

dir_pangenome=/lizardfs/guarracino/pangenomes/HPRCv2
dir_reads=$dir_base/data/HPRCv2/illumina

export PATH="/lizardfs/guarracino/tools/bedtools2/bin:$PATH"
export PATH="/lizardfs/guarracino/tools/samtools-1.22.1:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/impg/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/wfmash/build/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/seqwish/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/smoothxg/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/GFAffix/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/odgi/bin:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/pggb:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/gfainject/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/gafpack/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/panplexity/target/release:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/bwa-mem2-2.2.1_x64-linux:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/cosigt:$PATH"
export PATH="/lizardfs/guarracino/tools_for_genotyping/agc-1.1_x64-linux:$PATH"
#export PATH="/lizardfs/guarracino/condatools/emboss/6.6.0/bin/stretcher:$PATH" 
```

## Refinement of regions of interest

Samples for which we have short reads:

```shell
ls $dir_base/data/HPRCv2/illumina/*.cram | while read f; do echo $(basename $f .final.cram); done > $dir_base/data/HPRCv2/illumina/samples.txt
```

Pangenome vs reference alignment:

```shell
conda activate /lizardfs/guarracino/condatools/minimap2/2.30

mkdir -p $dir_base/minimap2
cd $dir_base/minimap2

ls $dir_pangenome/*.fa.gz | while read fasta; do
    sample=$(basename $fasta .fa.gz)

    sbatch -p allnodes -c 24 --job-name $sample-vs-grch38 --wrap "hostname; cd /scratch; minimap2 -x asm20 --eqx -c -t 24 $dir_base/reference/GRCh38.fa.gz $fasta > $sample-vs-grch38.paf; mv $sample-vs-grch38.paf $dir_base/minimap2/"
done
    
conda deactivate
```

