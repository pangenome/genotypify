# Data
Genotyping lots of samples with big pangenomes

## Variables

```shell
dir_base=/lizardfs/guarracino/genotypify
```

## Loci

Human accelerated regions (HARs) from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180714:

```shell
mkdir -p $dir_base/data/loci
cd $dir_base/data/loci
wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180714/suppl/GSE180714%5FHARs.bed.gz
gunzip GSE180714_HARs.bed.gz

# Make HARs at least 6kbp long to be able to map them with WFMASH
awk -v OFS='\t' '{print($1,$2,$3,$4"-"$15)}' GSE180714_HARs.bed > tmp
(head -n 1 tmp; sed 1d tmp | python3 $dir_base/scripts/extend_bed_if_short.py $dir_base/data/grch38.fa.fai 6000 | bedtools sort | bedtools merge -d 0 -c 4 -o distinct) > GSE180714_HARs.extended.bed
rm tmp
mv GSE180714_HARs.bed GSE180714_HARs.bed.bak
```

## Reference (no PanSN-spec)

```shell
mkdir -p $dir_base/reference
cd $dir_base/reference

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

# Remove alt contigs to avoid fragmented ROI-specific pangenome sequences
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa $(grep chr GRCh38_full_analysis_set_plus_decoy_hla.fa.fai  | grep '_' -v | cut -f 1) | bgzip -l 9 -@ 24 > GRCh38.fa.gz
samtools faidx GRCh38.fa.gz
```

## Annotation

```shell
mkdir -p $dir_base/annotation
cd $dir_base/annotation

wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
```

## HPRCv2 assemblies

Notes at https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/README.md

```shell
# Create directories
mkdir -p /scratch/HPRCv2/annotation/repeat_masker
mkdir -p /scratch/HPRCv2/annotation/censat

cd /scratch/HPRCv2

# Download the index file
wget https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/assemblies_pre_release_v0.6.1.index.csv

## S3 locations for assembly are stored in column 13 (2025 Feb 24)
ASSEMBLY_COLUMN_NUM=13
tail -n +2 assemblies_pre_release_v0.6.1.index.csv | awk -F',' -v col="$ASSEMBLY_COLUMN_NUM" '{print $col}' | while read -r assembly_file; do
    echo "Downloading $assembly_file..."
    aws s3 --no-sign-request cp "$assembly_file" .
    samtools fadix $assembly_file

    dir_parent="$(dirname $assembly_file)";
    name_v1=$(basename $assembly_file)
    name_v1=$(echo $name_v1 | cut -f 1 -d '.')

    aws s3 --no-sign-request cp  "$dir_parent/annotation/repeat_masker/$name_v1.RepeatMasker.bed" /scratch/HPRCv2/annotation/repeat_masker
    bgzip -l 9 -@ 24 /scratch/HPRCv2/annotation/repeat_masker/$name_v1.RepeatMasker.bed

    aws s3 --no-sign-request cp "$dir_parent/annotation/censat/$name_v1.cenSat.bed" /scratch/HPRCv2/annotation/censat
    bgzip -l 9 -@ 24 /scratch/HPRCv2/annotation/censat/$name_v1.cenSat.bed
done

mkdir -p /lizardfs/guarracino/pangenomes/HPRCv2
mv /scratch/HPRCv2 /lizardfs/guarracino/pangenomes
```

## HPRCv2 short reads

Notes at https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/sequencing_data/README.md

```shell
mkdir -p /scratch/HPRCv2
cd /scratch/HPRCv2

# Download the index file
wget https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/sequencing_data/data_illumina_pre_release.index.csv

READS_COLUMN_NUM=22
# Manage that the ',' can be present as a column value, like in "NA20302,NA20313"
tail -n +2 data_illumina_pre_release.index.csv | awk -F, -v FPAT='([^,]+)|("[^"]+")' -v col="$READS_COLUMN_NUM" '{print $col}' | 
tr -d '"' | while read -r reads_file; do
    echo "Downloading $reads_file..."
    aws s3 --no-sign-request cp "$reads_file" .
done

ls *.cram | while read cram; do samtools index $cram; done

mkdir -p $dir_base/data/HPRCv2/illumina
mv /scratch/HPRCv2/*cram $dir_base/data/HPRCv2/illumina

rm -rf /scratch/HPRCv2
```

## FLAGGER

Notes at https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/data_tables/assembly_qc

```shell
mkdir -p /scratch/HPRCv2
cd /scratch/HPRCv2

# Download the index file (2025/04/24)
wget https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/assembly_qc/flagger/flagger_hifi_v0.1.csv

FLAGGER_COLUMN_NUM=4
tail -n +2 flagger_hifi_v0.1.csv | awk -F',' -v col="$FLAGGER_COLUMN_NUM" '{print $3, $col}' | while read -r haplotype flagger_file; do
    echo "Downloading $flagger_file for $haplotype..."
    name=$(basename $flagger_file)
    aws s3 --no-sign-request cp "$flagger_file" .
    mv $name $haplotype.bed
done

mkdir -p $dir_base/data/HPRCv2/flagger
mv /scratch/HPRCv2/*bed $dir_base/data/HPRCv2/flagger

rm -rf /scratch/HPRCv2

# Rename BED files to match the FASTA file names
mv HG06807_mat_v1.bed HG06807_mat_v1.0.genbank.bed
mv HG06807_pat_v1.bed HG06807_pat_v1.0.genbank.bed
```

## Ancient samples (TO DO)

<!-- Create folder:

```shell
mkdir -p $dir_base/sequencing_data/ancient
```

## TO IGNORE FOR NOW

15 samples from Marchi et al., 2022 (https://doi.org/10.1016/j.cell.2022.04.008).

```shell
sed 's/SL_MU_/SLMU/g' $dir_base/data/Marchi2022.TableS1.csv -i

cat $dir_base/data/Marchi2022.TableS1.csv | awk -F'\t' '{if ($6 == $9) $9 = $9"_1"; print}' OFS='\t' > x
mv x $dir_base/data/Marchi2022.TableS1.csv

# Download the submitted FASTQ files, discarding the XXX_realg.noclips_splitRG.bam files
mkdir -p /scratch/Marchi2022
cd /scratch/Marchi2022
grep -Ff <(sed '1d' $dir_base/data/Marchi2022.TableS1.csv | cut -f 9) $dir_base/data/filereport_tables/filereport_read_run_PRJEB50857_tsv.txt | awk -F'\t' '{split($8, a, ";"); for (i in a) print "ftp://"a[i]}' | xargs -n 1 wget

cd /scratch
mv /scratch/Marchi2022 $dir_base/sequencing_data/ancient

sed -e 's/VC3-2/VC3/g' -e 's/Bar25/BAR25/g' $dir_base/data/Marchi2022.Table1.csv -i
# Marchi2022.Table1.csv to get the 15 samples
# Marchi2022.TableS1.csv to map `SAMPLE <-> FILE`
grep -Ff <(cut -f 1 $dir_base/data/Marchi2022.Table1.csv) $dir_base/data/Marchi2022.TableS1.csv -w
```

### Likely not raw reads

315 (317 - 2 samples with missing FASTQ files) samples from Allentoft et al., 2024 (https://doi.org/10.1038/s41586-023-06865-0):

```shell
# Counts of ancient samples they collected from other studies
sed '1d' Allentoft2024.SupplementaryData7.csv | cut -f 8 | sort | uniq -c | sort -k 1,1n

mkdir -p $dir_base/sequencing_data/ancient/Allentoft2024
cd /scratch
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB64656_tsv.txt | cut -f 7 | grep ftp | awk -F'\t' 'NR>1 {split($1, a, ";"); for (i in a) print "ftp://"a[i]}' | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Allentoft2024'

# To map `SAMPLE <-> FILE`
(sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB64656_tsv.txt | cut -f 7,8 | grep 'NEO898\|NEO813\|ERR12075080' -v | awk -v FS='/' '{print($7,$12)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' -e 's/.sort.rmdup.realign.md.bam;//g' -e 's/.neo.merge.bam;//g' -e 's/.neo.clean.bam;//g' | sed 's/ //g'; echo "ERR12075080\tNEO962")
```

442 samples from Margaryan et al., 2020 (https://doi.org/10.1038/s41586-020-2688-8; https://www.ebi.ac.uk/ena/browser/view/PRJEB37976):

```shell
mkdir -p $dir_base/sequencing_data/ancient/Margaryan2020
cd /scratch/
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB37976_tsv.txt | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Margaryan2020'

# To map `SAMPLE <-> FILE`
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB37976_tsv.txt | cut -f 7,8 | awk -v FS='/' '{print($7,$12)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' -e 's/.final.bam//g' | sed 's/ //g'
```

134 samples from Antonio et al., 2019 (https://www.ebi.ac.uk/ena/browser/view/PRJEB32566; https://www.ebi.ac.uk/ena/browser/view/PRJEB32566):

```shell
mkdir -p $dir_base/sequencing_data/ancient/Antonio2019
cd /scratch/
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB32566_tsv.txt | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Antonio2019'

# To map `SAMPLE <-> FILE`
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB32566_tsv.txt | cut -f 7,8 | awk -v FS='/' '{print($7,$12)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' -e 's/.bam;//g' -e 's/.bam//g' | sed 's/ //g'
```

137 samples from Damgaard et al., 2018, Nature (https://doi.org/10.1038/s41586-018-0094-2; https://www.ebi.ac.uk/ena/browser/view/PRJEB20658):

```shell
mkdir -p $dir_base/sequencing_data/ancient/Damgaard2018_Nature
cd /scratch
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB20658_tsv.txt | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Damgaard2018_Nature'

# To map `SAMPLE <-> FILE`
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB20658_tsv.txt | cut -f 7,8 | awk -v FS='/' '{print($7,$12)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' -e 's/.sort.rmdup.realign.md.bam;//g' | sed 's/ //g'
```

70 samples from Damgaard et al., 2018, Science (https://doi.org/10.1126/science.aar7711; https://www.ebi.ac.uk/ena/browser/view/PRJEB25389 and https://www.ebi.ac.uk/ena/browser/view/PRJEB26349):

```shell
# 70 samples, but the paper mentions 74 samples
mkdir -p $dir_base/sequencing_data/ancient/Damgaard2018_Science
cd /scratch
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB26349_tsv.txt | grep -Ff <(sed '1d' $dir_base/data/Damgaard2018.Science.SupplementaryTable1.csv | cut -f 1) | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Damgaard2018_Science'

# The 4 missed samples are in the Damgaard et al., 2018, Nature paper
sed '1d' $dir_base/data/Damgaard2018.Science.SupplementaryTable1.csv | cut -f 1 | sort | grep -Ff <(sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB26349_tsv.txt | grep -Ff <(sed '1d' data/Damgaard2018.Science.SupplementaryTable1.csv | cut -f 1) | cut -f 8 | cut -f 6 -d '/' | cut -f 1 -d '.' | cut -f 1 -d '_' ) -v
# DA343
# DA353
# DA356
# DA361
```

102 samples from Allentoft et al., 2015 (https://doi.org/10.1038/nature14507; https://www.ebi.ac.uk/ena/browser/view/PRJEB9021).

```shell
mkdir -p $dir_base/sequencing_data/ancient/Allentoft2015
cd /scratch
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB9021_tsv.txt | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Allentoft2015'

# To map `SAMPLE <-> FILE`
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB9021_tsv.txt | cut -f 7,8 | awk -v FS='/' '{print($6,$11)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' | sed -e 's/.hg19.flt.sort.rmdup.realign.md.bam//g' | sed 's/ //g'
```

58 samples from Brunel et al., 2020 (https://doi.org/10.1073/pnas.1918034117; https://www.ebi.ac.uk/ena/browser/view/PRJEB38152)

```shell
#sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB38152_tsv.txt | grep -i -Ff <(sed '1,2d' Brunel2018.SupplementaryTable1.csv | cut -f 1 | sed -e 's/Schw72-15/Sch72-15/g' -e 's/CRE8C/CRE8-C/g' -e 's/BIS388-/BIS388/g'; echo CRE20D)
mkdir -p $dir_base/sequencing_data/ancient/Brunel2020
cd /scratch
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB38152_tsv.txt | grep bam | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Brunel2020'

# To map `SAMPLE <-> FILE`
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB38152_tsv.txt | cut -f 7,8 | grep bam | awk -v FS='/' '{print($7,$12)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' | sed -e 's/.mappedHs37d5.trim2.bam//g' -e 's/.mappedHs37d5.trim2.bam//g' -e 's/.mappedHs37.trim2.bam//g' | sed 's/ //g'
```

91??? samples frm Scheib et al., 2018 (https://doi.org/10.1126/science.aar6851; https://www.ebi.ac.uk/ena/browser/view/PRJEB25445)

```shell
mkdir -p $dir_base/sequencing_data/ancient/Scheib2018
cd /scratch
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB25445_tsv.txt | cut -f 7 | grep ftp | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Scheib2018'

# Check Scheib2018.SupplementaryTable1.tsv and Scheib2018.SupplementaryTable3.tsv (number of reads) to process downloaded samples
```

34 samples from Sikora et al.m 2019 (https://doi.org/10.1038/s41586-019-1279-z; https://www.ebi.ac.uk/ena/browser/view/PRJEB29700 and https://www.ebi.ac.uk/ena/browser/view/PRJEB26336):

```shell
mkdir -p $dir_base/sequencing_data/ancient/Sikora2019
cd /scratch
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB29700_tsv.txt | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Sikora2019'

# To map `SAMPLE <-> FILE`
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB29700_tsv.txt | cut -f 7,8 | awk -v FS='/' '{print($7,$12)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' -e 's/.sort.rmdup.realign.md.bam;//g' -e 's/.realigned.calmd.readsadded.bam;//g' -e 's/.sort.rmdup.uniq.rg.realn.md.bam;//g' | sed 's/ //g'
```

35 samples from Krzewinska et al., 2018 (https://doi.org/10.1126/sciadv.aat4457; https://www.ebi.ac.uk/ena/browser/view/PRJEB27628):

```shell
mkdir -p $dir_base/sequencing_data/ancient/Krzewinska2018
cd /scratch
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB27628_tsv.txt | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Krzewinska2018'

# To map `SAMPLE <-> FILE`
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB27628_tsv.txt | cut -f 7,8 | awk -v FS='/' '{print($7,$12)}' | sed -e 's/ftp.sra.ebi.ac.uk//g' -e 's/.fastq.gz//g' | cut -f 1 -d '_' | sed 's/ //g'
```

25 samples from McColl et al., 2018 (https://doi.org/10.1126/science.aat3628; https://www.ebi.ac.uk/ena/browser/view/PRJEB26721)

```shell
mkdir -p $dir_base/sequencing_data/ancient/McColl2018
cd /scratch
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB26721_tsv.txt | cut -f 7 | awk -F'\t' 'NR>1 {split($1, a, ";"); for (i in a) print "ftp://"a[i]}' | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/McColl2018'
```

Ebenesersdottir et al.m 2018, Science (https://doi.org/10.1126/science.aar2625): I NEED ACCESS TO THE PUBLICATION

15/24?? samples from Schroeder et al., 2019 (https://doi.org/10.1073/pnas.1820210116; https://www.ebi.ac.uk/ena/browser/view/PRJEB28451):

```shell
mkdir -p $dir_base/sequencing_data/ancient/Schroeder2019
cd /scratch
sed '1d' $dir_base/data/filereport_tables/filereport_read_run_PRJEB28451_tsv.txt | cut -f 7 | xargs -n 1 -I {} sh -c 'wget {}; mv $(basename {}) /lizardfs/guarracino/genotypify/sequencing_data/ancient/Schroeder2019'
```

### Check if reads are sanitized with hg38

With minimap2:

```shell
cd /lizardfs/guarracino/genotypify/sequencing_data/ancient

# Map
conda activate /lizardfs/guarracino/condatools/minimap2/2.28/

mkdir -p /lizardfs/guarracino/genotypify/sequencing_data/ancient/reads-vs-chm13
cd /lizardfs/guarracino/genotypify/sequencing_data/ancient/reads-vs-chm13
REF_FASTA=/lizardfs/guarracino/pangenomes/references/chm13v2.fa.gz
ls /lizardfs/guarracino/genotypify/sequencing_data/ancient | while read DATASET; do
    ls /lizardfs/guarracino/genotypify/sequencing_data/ancient/$DATASET/*fastq.gz | while read FASTQ; do
        NAME=$(basename $FASTQ .fastq.gz)
        echo $DATASET $FASTQ $NAME
        sbatch -c 48 -p tux --job-name "$DATASET-$NAME" --wrap "hostname; cd /scratch; \time -v minimap2 $REF_FASTA $FASTQ -x sr -t 48 | pigz -9 > $DATASET.$NAME.mm2.map.paf.gz; mv $DATASET.$NAME.mm2.map.paf.gz /lizardfs/guarracino/genotypify/sequencing_data/ancient/reads-vs-chm13/"
    done
done

# Coverage
cut -f 1,2 /lizardfs/guarracino/pangenomes/references/chm13v2.fa.gz.fai > /lizardfs/guarracino/genotypify/sequencing_data/ancient/reads-vs-chm13/chm13.genome
bedtools makewindows -g /lizardfs/guarracino/genotypify/sequencing_data/ancient/reads-vs-chm13/chm13.genome -w 200000 | bedtools sort > /lizardfs/guarracino/genotypify/sequencing_data/ancient/reads-vs-chm13/chm13v2.windows.200kbp.bed

ls /lizardfs/guarracino/genotypify/sequencing_data/ancient | while read DATASET; do
    ls /lizardfs/guarracino/genotypify/sequencing_data/ancient/$DATASET/*fastq.gz | while read FASTQ; do
        NAME=$(basename $FASTQ .fastq.gz)
        echo $DATASET $FASTQ $NAME

        PAF=/lizardfs/guarracino/genotypify/sequencing_data/ancient/reads-vs-chm13/$DATASET.$NAME.mm2.map.paf.gz
        NAME=$(basename $PAF .paf.gz)

        # For base-by-base coverage depth
        sbatch -c 8 -p tux --job-name "$DATASET-$NAME" --wrap "hostname; cd /scratch; cp $PAF /scratch/; zcat $NAME.paf.gz | cut -f 6,8,9 | sort -k1,1 -k2,2n -T /scratch > $NAME.bed; bedtools genomecov -i $NAME.bed -g /lizardfs/guarracino/genotypify/sequencing_data/ancient/reads-vs-chm13/chm13.genome -bga > $NAME.depth.bedgraph; bedtools map -a /lizardfs/guarracino/genotypify/sequencing_data/ancient/reads-vs-chm13/chm13v2.windows.200kbp.bed -b $NAME.depth.bedgraph -c 4 -o mean > $NAME.depth.windows.200kbp.bed; pigz -9 $NAME.bed; pigz -9 $NAME.depth.bedgraph; rm $NAME.paf.gz; mv $NAME.* /lizardfs/guarracino/genotypify/sequencing_data/ancient/reads-vs-chm13/"
    done
done

# # Uncalled regions in hg38
# wget -c https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
# zcat gap.txt.gz | cut -f 2,3,4,8 > hg38.gaps.bed
# rm gap.txt.gz

# Acro short arms vs the rest
wget -c https://raw.githubusercontent.com/pangenome/chromosome_communities/refs/heads/main/data/annotation/chm13.p_arms.approximate.acros.bed
sed 's/chm13#/chm13#1#/g' chm13.p_arms.approximate.acros.bed -i

cd /lizardfs/guarracino/genotypify/sequencing_data/ancient/reads-vs-chm13
ls /lizardfs/guarracino/genotypify/sequencing_data/ancient | while read DATASET; do
    ls /lizardfs/guarracino/genotypify/sequencing_data/ancient/$DATASET/*fastq.gz | while read FASTQ; do
        NAME=$(basename $FASTQ .fastq.gz)
        echo $DATASET $FASTQ $NAME

        PAF=/lizardfs/guarracino/genotypify/sequencing_data/ancient/reads-vs-chm13/$DATASET.$NAME.mm2.map.paf.gz
        NAME=$(basename $PAF .paf.gz)

        COVERAGE_BED=$NAME.depth.windows.200kbp.bed

        # We'll mark everything as "non-acro-parm" first
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,"non-acro-parm"}' $COVERAGE_BED > /scratch/$COVERAGE_BED.tmp
        # Then use bedtools intersect to identify and mark the acro-p-arm regions
        bedtools intersect -a /scratch/$COVERAGE_BED.tmp -b chm13.p_arms.approximate.acros.bed -wa -wb | \
            awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,"acro-p-arm"}' > /scratch/$COVERAGE_BED.acro.bed
        # Finally, use bedtools subtract to keep the non-acro regions with their original marking
        bedtools subtract -a /scratch/$COVERAGE_BED.tmp -b chm13.p_arms.approximate.acros.bed > /scratch/$COVERAGE_BED.non-acro.bed
        # Combine the results
        cat /scratch/$COVERAGE_BED.acro.bed /scratch/$COVERAGE_BED.non-acro.bed | sort -k1,1 -k2,2n -T /scratch > $(basename $COVERAGE_BED .bed).classified.bed
        rm /scratch/$COVERAGE_BED.tmp /scratch/$COVERAGE_BED.acro.bed /scratch/$COVERAGE_BED.non-acro.bed
    done
done

# Take info from the top 10 biggest files in each dataset
cd /lizardfs/guarracino/genotypify/sequencing_data/ancient/reads-vs-chm13
ls /lizardfs/guarracino/genotypify/sequencing_data/ancient | while read DATASET; do
    ls /lizardfs/guarracino/genotypify/sequencing_data/ancient/$DATASET/*fastq.gz -S | head -n 200 | while read FASTQ; do
        NAME=$(basename $FASTQ .fastq.gz)
        #echo $DATASET $FASTQ $NAME

        PAF=/lizardfs/guarracino/genotypify/sequencing_data/ancient/reads-vs-chm13/$DATASET.$NAME.mm2.map.paf.gz
        NAME=$(basename $PAF .paf.gz)

        COVERAGE_BED=$NAME.depth.windows.200kbp.bed

        NAME2=$(basename $PAF .paf.gz | cut -f 2 -d '.')
        cat $(basename $COVERAGE_BED .bed).classified.bed | awk -v OFS='\t' -v dataset=$DATASET -v name=$NAME2 '{print($0,dataset,name)}'
    done
done | pigz -9 > /scratch/ancientHead200.depth.windows.200kbp.bed.gz && mv /scratch/ancientHead200.depth.windows.200kbp.bed.gz .

#zcat $PAF | awk -v OFS="\t" '{print $6, $8, $9, $1, ".", $5}' | sort -k1,1 -k2,2n -T /scratch > $(basename $PAF .paf.gz).bed
# Compute the coverage of sequence alignments (file B) across XX kilobase windows (file A) tiling a genome of interest
#bedtools coverage -a /scratch/chm13v2.windows.200kbp.bed -b $(basename $PAF .paf.gz).bed > coverage.bed
``` -->

## Modern samples (INCOMPLETE)

1000 Genomes Project sample collection to 30x coverage (from <https://www.internationalgenome.org/data-portal/data-collection/30x-grch38>). Initially, the 2504 unrelated samples from the phase three panel from the 1000 Genomes Project were sequenced. Thereafter, an additional 698 samples, related to samples in the 2504 panel, were also sequenced.

```shell
mkdir -p $dir_base/sequencing_data/modern/1000G/2504_high_coverage
cd $dir_base/sequencing_data/modern/1000G/2504_high_coverage
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index
awk '!/^#/ {print $1}' 1000G_2504_high_coverage.sequence.index | xargs -n 1 wget

mkdir -p $dir_base/sequencing_data/modern/1000G/additional_698_related
cd $dir_base/sequencing_data/modern/1000G/additional_698_related
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index
awk '!/^#/ {print $1}' 1000G_698_related_high_coverage.sequence.index | xargs -n 1 wget
```

Simons Genome Diversity Project (from <https://www.ebi.ac.uk/ena/browser/view/PRJEB9586>):

```shell
mkdir -p $dir_base/sequencing_data/modern/SGDP/
cd $dir_base/sequencing_data/modern/SGDP/
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/hgdp_wgs.sequence.index
awk -F'\t' 'NR>1 {split($7, a, ";"); for (i in a) print "ftp://"a[i]}' $dir_base/data/filereport_tables/filereport_read_run_PRJEB9586_tsv.txt | xargs -n 1 wget
```
