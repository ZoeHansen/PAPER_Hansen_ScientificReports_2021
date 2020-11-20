################################################

# Using Nonpareil to estimate Sequencing Coverage

################################################

# Information on the installation and use of Nonpareil can be found here: https://nonpareil.readthedocs.io/en/latest/index.html

### Run Nonpareil on the trimmed, human-genome-removed reads

gzip -d ~/ERIN_trim.f.pe.hre.fq.gz;
~/nonpareil/nonpareil -s ~/ERIN_trim.f.pe.hre.fq -T kmer -f fastq -b nonpareil/nonpareil_output;
gzip ~/ERIN_trim.f.pe.hre.fq

### Output from Nonpareil was exported and merged using Python
# Note: see script "Nonpareil_Output_MergeFiles.py"

### The nonpareil file paths were used to generate Nonpareil curves in R
# See 'SequencingMetrics_Rcode.R' 