################################################

# Using MicrobeCensus to estimate Average Genome Size

################################################

# Information on the installation and use of MicrobeCensus can be found here: https://github.com/snayfach/MicrobeCensus

### Run MicrobeCensus using Python
# Python Version 3.6.4 was used

python /mnt/home/hansenzo/MicrobeCensus/microbe_census_script.py

### The contents of 'microbe_census_script.py' are shown below:

  #! /usr/bin/env python

  from microbe_census import microbe_census
  import os

  args = {'seqfiles':['ERIN_trim.f.pe.hre.fq.gz','ERIN_trim.r.pe.hre.fq.gz'], 'threads':2, 'verbose':True}

  average_genome_size, args = microbe_census.run_pipeline(args)

  count_bases = microbe_census.count_bases(args)
  genome_equivalents = count_bases/average_genome_size

  print('Average genome size:', average_genome_size)
  print('Base Count:', count_bases)
  print('Genome Equivalents:', genome_equivalents)

  f = open('/mnt/home/hansenzo/MicrobeCensus/mc_output_SR1.txt','a+', encoding = 'utf-8')
  #with open('/mnt/home/hansenzo/MicrobeCensus/mc_output_test1.txt', encoding = 'utf-8') as f:
  sample = str(os.getcwd())
  sample1 = os.path.split(sample)[1]
  vals = str(sample1) + ',' + str(average_genome_size) + ',' + str(count_bases) + ',' + str(genome_equivalents) + '\n'
  f.write(vals)
  f.close()
  
### The contents of the .txt file were used to normalize read abundances. Abundances were normalized by the number of estimated genome equivalents (GEs) per sample. 
### Normalization was performed in R by dividing the number of reads in each sample by GEs


