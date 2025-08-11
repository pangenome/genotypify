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
    -g /group/soranzo/davide.bolognini/working/dev/cosigt_paper/real_data/data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -r /group/soranzo/davide.bolognini/working/dev/cosigt_paper/real_data/data/1000G/ \
    -b /group/soranzo/andrea.guarracino/C4.bed \
    -o /scratch/andrea.guarracino/uffa/ \
    --map /group/soranzo/andrea.guarracino/aln_map.tsv \
    --tmp /scratch/andrea.guarracino/ \
    --profile config/slurm/
sh cosigt_smk.sh
#    #--gtf /group/soranzo/chiara.paleni/references/gencode.v47.annotation.gtf.gz --proteins /group/soranzo/chiara.paleni/references/gencode.v47.pc_translations.fa

sh cosigft_smk.sh
cosigt|benchmark|refine
```

## 200 HPRCv2's short read samples on regions from the same 200 HPRCv2's assemblies (the full pangenome is 466 haplotypes)

