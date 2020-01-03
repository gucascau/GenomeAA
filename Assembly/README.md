# Hybrid assembly of Platygyra daedalea genome

## 1. de novo Assembly using FALCON

source env.sh
fc_run.py fc_run_new.cfg

## 2. de novo Assembly using HGAP3

2.1 referenceUploader (prepare reference repository)
2.2 prepare input.fofn. one-per-line, each “bax.h5” in input data set
2.3 convert “input.fofn” to “input.xml” file :
  $ fofnToSmrtpipeInput.py input.fofn >input.xml
2.4. prepare “params.xml” Prepare your "params.xml" file. Here is a params.xml template you can use; you should just need to edit the reference path.
How to set the parameter of params.xml:
https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/HGAP-Whitelisting-Tutorial


2.5. active your SMRTanalysis environment, and invoke smrtpipe
source <SMRT Analysis>/etc/setup.sh
smrtpipe.py –distribute –params=params.xml xml:input.xml
2.6 After successful execution is complete, the results should be available as data/consensus.fast[aq].gz and data/variants.gff.gz, etc.

## 3. Using Matepair and Pairend reads for assembly

### 3.1 Raw reads filtering 
Remove adaptor, low quality sequences ((< Q30 and < 30 bp) and duplicates reads using Trimmomatic v0.32

java -jar trimmomatic.jar PE -threads 40  -trimlog trimmoatic.log  Raw_final_R1.fastq Filter_R1.fastq unpaired_R1.fastq Raw_final_R2.fastq Filter_R2.fastq unpaired_R2.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:30

### 3.2 Mapping to S. minutum and S. microadriaticum
# Removing potential contamination from two Symbiodium sources using Bowtie2 v2.1.0.

### 3.3 Combine two Symbiodium genomes and build the index
bowtie2-build –f genome.fa genome_index
### 3.4 Map different libraries to two genomes
bowtie2 –x genome_index -1 Filtered_R1.fastq -2 Filtered_R2.fastq -q -t --un mapresult_un --al mapresult_al --un-conc mapresult_un_conc --al-conc mapresult_al_conc -p 10 –S mapping.sam

### 3.5 Genome size estimate
#Estimate genome size using KmerFreq_AR under different K-mer
KmerFreq_AR -k 17 -t 4 -c -1 -p test test_read.lst >kmerfreq.cout 2>kmerfreq.cerr

Digital normalization (Khmer v1.4)
### 3.6 Normalize everything to the coverage of 20 for pair end reads
interleave-reads.py Filter.1.fastq Filter.2.fastq > Lall_B_filter.fastq
normalize-by-median.py -k 20 -C 20 -N 4 -x 4e9 -p -s Lall_B_filter_norm_hash.kh Lall_B_filter.fastq
### 3.7 Trim off any k-mer that are abundance in high-coverage reads
filter-abund.py Lall_B_filter_norm_hash.kh Lall_B_filter.fastq.keep
extract-paired-reads.py Lall_B_filter.fastq.keep.abundfilt
### 3.8 Normalize down to C=10
normalize-by-median.py -k 20 -C 10 -N 4 -x 4e9 -p Lall_B_filter.fastq.keep.abundfilt.pe
extract-paired-reads.py Lall_B_filter.fastq.keep.abundfilt.pe.keep
split-paired-reads.py Lall_B_filter.fastq.keep.abundfilt.pe.keep.pe

ALLPATH-LG assembly 
### 3.9 Prepare data for ALLPATHS
PrepareAllPathsInputs.pl DATA_DIR=~/Allpath_assemble/data PLOIDY=2 IN_GROUPS_CSV=in_groups.csv IN_LIBS_CSV=in_libs.csv GENOME_SIZE=400000000 OVERWRITE=True | tee prepare2.out
### 3.10 Run ALLPATH-LG using HAPLOIDIFY and OVERWRITE parameters
RunAllPathsLG \
 PRE=$PWD\
 REFERENCE_NAME=norm\
 DATA_SUBDIR=data\
 RUN=run\
 SUBDIR=FINAL\
 OVERWRITE=True\
 TARGETS=standard\
 HAPLOIDIFY=True  | tee -a assemble2.out

## 4. Gapfilling
Close Gap using Gapcloser v1.12
GapCloser –a Final.corrected.fasta –b library.txt –o Gapcloser.fasta –l 125 –p 31



## 5. Use PBjelly for scaffolding 

source ./setup.sh
summarizeAssembly.py Braincoral_gapcloser.fasta
Jelly.py setup Protocol.xml
Jelly.py mapping Protocol.xml
Jelly.py support Protocol.xml
Jelly.py extraction Protocol.xml

## 6. Use SSPace for scaffolding

nohup perl /home/xinw/software/SSPACE-BASIC-2.0_linux-x86_64/SSPACE_Basic_v2.0.pl -l library.txt  -s Brain_coral_pbjelly.fasta -b Brain_coral_pbjelly_sspace  -T 60 >log.sspace 

## 7. Final round of Gapfilling

nohup /home/xinw/software/Gapcloser/GapCloser  -a Brain_coral_pbjelly_sspace.final.scaffolds.fasta -b Braincoral_config.txt -o Brain_coral_pbjelly_sspace.final.scaffolds.gapcloser.fast -t 60 >log.gapcloser &




