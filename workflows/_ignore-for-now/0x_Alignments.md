
# Alignments  (TO IGNORE FOR NOW)

## Variables

```shell
dir_base=/lizardfs/guarracino/genotypify
```

## Ancient

15 samples from Marchi et al., 2022 (https://doi.org/10.1016/j.cell.2022.04.008).

```shell
mkdir -p $dir_base/alignment/ancient/Marchi2022
cd $dir_base/alignment/ancient/Marchi2022

conda activate /lizardfs/guarracino/condatools/cosigt # for bwa-mem2

# bwa aln -l 1024 -n 0.01 -o 2 from https://pmc.ncbi.nlm.nih.gov/articles/PMC8717315/

cut -f 1 $dir_base/data/Marchi2022.Table1.csv | sed 1d | sort | while read sample; do
    sex=$(cat $dir_base/data/Marchi2022.Table1.csv| cut -f 1,7 | grep $sample | cut -f 2)

    if [ $sex = 'M' ]; then
        ref=/lizardfs/guarracino/rob_evolution/assemblies/chm13v2.0_masked_DJ_5S_rDNA_PHR_PAR_wi_rCRS.fa
    else
        ref=/lizardfs/guarracino/rob_evolution/assemblies/chm13v2.0_masked_DJ_5S_rDNA_PHR_noY_wi_rCRS.fa
    fi

    echo $sample $sex $ref

    # Single-End runs
    grep -Ff <(cut -f 1 $dir_base/data/Marchi2022.Table1.csv) $dir_base/data/Marchi2022.TableS1.csv -w | grep "^$sample" -w | awk -F'\t' '{if ($13 !~ /PE/) print($9)}' | cut -f 9 | while read prefix; do
        ls $dir_base/sequencing_data/ancient/Marchi2022/${prefix}_*.fastq.gz | while read fastq; do
            name=$sample.$prefix
            sbatch -c 48 -p tux --job-name bwa-SE-$name --wrap "hostname; \time -v bwa aln -t 48 $ref $fastq -l 1024 -n 0.01 -o 2 > /scratch/${name}.SE.sai; \time -v bwa samse $ref /scratch/${name}.SE.sai $fastq > /scratch/${name}.SE.sam; mv /scratch/${name}.SE.sam $dir_base/alignment/ancient/Marchi2022; rm /scratch/${name}.SE.sai"
        done
    done

    # Paired-End runs
    grep -Ff <(cut -f 1 $dir_base/data/Marchi2022.Table1.csv) $dir_base/data/Marchi2022.TableS1.csv -w | grep "^$sample" -w | awk -F'\t' '{if ($13 ~ /PE/) print}' | cut -f 9 | while read prefix; do
        name=$sample.$prefix
        R1_fastq=$(eval echo "$dir_base/sequencing_data/ancient/Marchi2022/$prefix*_R1_*.fastq.gz")
        R2_fastq=$(eval echo "$dir_base/sequencing_data/ancient/Marchi2022/$prefix*_R2_*.fastq.gz")

        sbatch -c 96 -p tux --job-name bwa-PE-$name --wrap "hostname; \time -v bwa aln -t 96 $ref $R1_fastq -l 1024 -n 0.01 -o 2 > /scratch/${name}.R1.sai; \time -v bwa aln -t 96 $ref $R2_fastq -l 1024 -n 0.01 -o 2 > /scratch/${name}.R2.sai; \time -v bwa sampe $ref /scratch/${name}.R1.sai /scratch/${name}.R2.sai $R1_fastq $R2_fastq > /scratch/${name}.PE.sam; mv /scratch/${name}.PE.sam $dir_base/alignment/ancient/Marchi2022; rm /scratch/${name}.R1.sai /scratch/${name}.R2.sai"
    done
done
```

Fixmate in PE alignments:

```shell
cut -f 1 $dir_base/data/Marchi2022.Table1.csv | sed 1d | sort | while read sample; do
    sex=$(cat $dir_base/data/Marchi2022.Table1.csv| cut -f 1,7 | grep $sample | cut -f 2)

    if [ $sex = 'M' ]; then
        ref=/lizardfs/guarracino/rob_evolution/assemblies/chm13v2.0_masked_DJ_5S_rDNA_PHR_PAR_wi_rCRS.fa
    else
        ref=/lizardfs/guarracino/rob_evolution/assemblies/chm13v2.0_masked_DJ_5S_rDNA_PHR_noY_wi_rCRS.fa
    fi

    echo $sample $sex $ref

    # Single-End runs
    grep -Ff <(cut -f 1 $dir_base/data/Marchi2022.Table1.csv) $dir_base/data/Marchi2022.TableS1.csv -w | grep "^$sample" -w | awk -F'\t' '{if ($13 !~ /PE/) print($9)}' | cut -f 9 | while read prefix; do
        ls $dir_base/sequencing_data/ancient/Marchi2022/${prefix}_*.fastq.gz | while read fastq; do
            name=$sample.$prefix
            mv $dir_base/alignment/ancient/Marchi2022/${name}.SE.sam $dir_base/alignment/ancient/Marchi2022/${name}.SE.fix.sam
        done
    done

    # Paired-End runs
    grep -Ff <(cut -f 1 $dir_base/data/Marchi2022.Table1.csv) $dir_base/data/Marchi2022.TableS1.csv -w | grep "^$sample" -w | awk -F'\t' '{if ($13 ~ /PE/) print}' | cut -f 9 | while read prefix; do
        name=$sample.$prefix
        R1_fastq=$(eval echo "$dir_base/sequencing_data/ancient/Marchi2022/$prefix*_R1_*.fastq.gz")
        R2_fastq=$(eval echo "$dir_base/sequencing_data/ancient/Marchi2022/$prefix*_R2_*.fastq.gz")
        
        samtools sort -@ 48 -n -m8G -T /scratch/${NAME}.tmp $dir_base/alignment/ancient/Marchi2022/${name}.PE.sam > /scratch/${name}.PE.bam
        samtools fixmate -m -@ 48 /scratch/${name}.PE.bam /scratch/${name}.PE.fix.bam
        rm /scratch/${name}.PE.bam
        mv /scratch/${name}.PE.fix.bam $dir_base/alignment/ancient/Marchi2022
    done
done
```

Sort and dedup:

```shell

samtools sort -@ 48 -O BAM -o /scratch/${NAME}.bam -T /scratch/${NAME}.tmp /scratch/${NAME}.fix.bam
rm /scratch/${NAME}.fix.bam
samtools index /scratch/${NAME}.bam
samtools markdup -r -@ 48 /scratch/${NAME}.bam /scratch/${NAME}.dedup.bam
rm /scratch/${NAME}.ba*

# Filter out secondary alignments
mkdir -p /scratch/$NAME/
samtools view -@ 48 -F0x100 -hb --write-index -o /scratch/${NAME}.dedup.pri.bam /scratch/$NAME/${NAME}.dedup.bam
rm /scratch/${NAME}.dedup.bam
```