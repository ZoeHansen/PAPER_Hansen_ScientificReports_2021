#####################################################

# Campylobacter Metagenomics Pipeline - AmrPlusPlus

#####################################################

# Note: Our samples are designated by an "ER" prefix, followed by a unique identifier (e.g. ER0043, ER0073, etc.)

### PATH to Metagenomes of Cases and Controls (Campylobacter)
/mnt/research/manninglab/Sequencing_Data/misc_sequencing/20161123_B_DNASeq_PE/ER*_NNNNNNNN-NNNNNNNN_L00*_R*.fastq.gz
/mnt/research/manninglab/Sequencing_Data/misc_sequencing/20150123_DNASeq_PE_metagenomes/ER*_NNNNNNNN-NNNNNNNN_L00*_R*.fastq.gz
/mnt/research/manninglab/Sequencing_Data/misc_sequencing/20170206_B_DNASeq_PE/ER*_NNNNNNNN-NNNNNNNN_L00*_R*.fastq.gz

### Merge paired-end reads
sample_names=(ER0* .. ER1*)

for sample in ${sample_names[@]}
do
cat ./${sample}_*_L001_R1_001.fastq.gz ${sample}_*_L002_R1_001.fastq.gz > /mnt/scratch/hansenzo/${sample}_R1.fastq.gz;   # merging forward reads
cat ./${sample}_*_L001_R2_001.fastq.gz ${sample}_*_L002_R2_001.fastq.gz > /mnt/scratch/hansenzo/${sample}_R2.fastq.gz    # merging reverse reads
done


### Run AmrPlusPlus on our samples

# Note: the AmrPlusPlus pipeline includes adapter trimming, quality filtering, non-host genome removal, alignment to MEGARes 2.0, and annotation
# The below script also includes commands to integrate haplotype discovery using the Resistance Gene Identifier (RGI) against the CARD database

# Submit a batch script to MSU's High Perfomance Computing Center (HPCC):

module purge
module load GCC/8.3.0
module load AmrPlusPlus/2.0.2

nextflow run $EBROOTAMRPLUSPLUS/main_AmrPlusPlus_v2_withRGI.nf --reads "/mnt/scratch/hansenzo/amrplusplus/*_R{1,2}.fastq.gz" --host $HOME/Database/HumanRef/GRCh38_latest_genomic.fna.gz --adapters $HOME/Database/adapters.fa --annotation "$EBROOTAMRPLUSPLUS/data/amr/megares_annotations_v1.02.csv" --amr "$EBROOTAMRPLUSPLUS/data/amr/megares_database_v1.02.fasta" --card_db $HOME/amrplusplus/amrplusplus_v2/card.json --output /mnt/scratch/hansenzo/amrplusplus_results -w /mnt/scratch/hansenzo/amrplusplus_work -profile singularity 
