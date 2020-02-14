# Annotation of Platygyra daedalea

## 1. Prepared the datasets
Download the current replease repeatbase:

▪	EMBL format (87.86 MB) 01-24-2017:
▪	RepBase22.01.embl.tar.gz 
▪	FASTA format (42.11 MB) 01-24-2017:
▪	RepBase22.01.fasta.tar.gz  
▪	Repbase-derived RepeatMasker libraries:
▪	RepBaseRepeatMaskerEdition-20170127.tar.gz (48.84 MB)
▪	REPET edition:
▪	RepBase20.05_REPET.embl.tar.gz (42.84 MB)
▪	Repbase source data for DFAM:
▪	dfamrepref.embl.tgz (0.75 MB)

## 2. Mask the assembled genome
/home/xinw/software/RepeatModeler/BuildDatabase -name Braincoral_genome_repeat Braincoral_assembly_novirus_v1.fasta -engine ncbi
/home/xinw/software/RepeatModeler/RepeatModeler -database Braincoral_genome_repeat -enginie ncbi -pa 25 -recoverDir Braincoral_failed >log.repeatmodel &
/home/xinw/software/RepeatMasker/RepeatMasker -e ncbi -pa 30 -s -lib /home/xinw/database/Libraries/RMRBSeqs.embl  -gff Braincoral_assembly_novirus_v1.fasta.masked >log.repbase &

## 3. Assembly the transcriptome

3.1 quality control
java -jar /home/xinw/software/trimmomatic/classes/trimmomatic.jar PE  -threads 40  -trimlog trimmoatic.log  R1.fastq R2.fastq  Braincoral_rna_filter_1.fastq Braincoral_rna_filter_2.fastq Braincoral_rna_unpaired_1.fastq Braincoral_rna_unpaired_2.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
TrimmomaticPE: Started with arguments: -threads 40 -trimlog trimmoatic.log R1.fastq R2.fastq Braincoral_rna_filter_1.fastq Braincoral_rna_filter_2.fastq Braincoral_rna_unpaired_1.fastq Braincoral_rna_unpaired_2.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

3.2 Use hisat2-strintie to create the trascript
hisat2-build Braincoral_assembly_novirus_v1.fasta Braincoral_genome_index -p 30
hisat2 -p 60 --mp 8,1 --rdg 6,2 --rfg 6,2 --dta -x ../Braincoral_genome_index -1 Braincoral_rna_filter_R1.fastq -2 Braincoral_rna_filter_R2.fastq -S Braincoral_rna_map_v2.sam --un Braincoral_rna_filter_unpaired --un-conc Braincoral_rna_filter_paired_unconc --al-conc Braincoral_rna_filter_al_conc 



## 4. Evaluation of transcriptome

Formatdb –p F –I Braincoral_assembly_v1.fasta
python3 /home/xinw/software/BUSCO_v1.1b1/BUSCO_v1.1b1.py  -o Braincoral_trancript_evalution -in Pdae_v1.fasta -l /home/xinw/software/BUSCO_v1.1b1/metazoa/ -m tran -c 40 

## 5. using Scipio to create a trainingset of gene structures
ln -s ~/project/Brain_coral_new/Gene_prediction/database/Braincoral_assembly_novirus_v1.fasta genome.fa
ln -s ~/project/Brain_coral_new/Gene_prediction/database/closely_related_protein.fa proteins.aa

perl /home/xinw/software/scipio-1.4/scipio.1.4.1.pl --blat_output=prot.vs.genome.psl genome.fa proteins.aa > scipio.yaml
scipio.1.4.pl --blat_output=prot.vs.genome.psl genome.fa proteins.aa > scipio.yaml
cat scipio.yaml | yaml2gff.1.4.pl > scipio.scipiogff
scipiogff2gff.pl --in=scipio.scipiogff --out=scipio.gff
cat scipio.yaml |/home/xinw/software/scipio-1.4/yaml2log.1.4.pl >scipio.log 
gff2gbSmallDNA.pl scipio.gff genome.fa 1000 genes.raw.gb

### 6. Using the scripio results to train AUGUSTUS
6.1 build the trainingset
etraining --species=platygyra --stopCodonExcludedFromCDS=true genes.raw.gb 2>train.error
6.2 filter out the problematic genes
cat train.err | perl -pe 's/.*in sequence (\S+): .*/$1/' > badgenes.lst
filterGenes.pl badgenes.lst genes.raw.gb > genes.gb
grep -c "LOCUS" genes.raw.gb genes.gb

### 7. improve the accurate of trainingsets
SPlitgene structure set into training and test set
randomSplit.pl genes.gb 100
etraining --species=platygyra --stopCodonExcludedFromCDS=true genes.gb.train
make a first try and predict the genes in gene.gb.train ab initio
augustus --species=bug genes.gb.test | tee firsttest.out # takes ~1m

### 8. Use cDNA to create a hits file
blat -minIdentity=92 Braincoral_assembly_novirus_v1.fasta Pdae_v1.fasta cdna.psl Braincoral_cDNA.psl
filter the cDNA alignments and report only the highest-scoring spliced alignments for each cDNA.
/home/xinw/software/ucsc_tools/ucsc_tools/executables/pslCDnaFilter -maxAligns=1 Braincoral_cDNA.psl  Braincoral_cDNA_uniq.psl
create hints file:
blat2hints.pl --in=Braincoral_cDNA_uniq.psl --out=hits.gff
Furthe traing with augustus
etraining --species=platygyra --stopCodonExcludedFromCDS=true genes.raw.gb 2>train.error

### 9. Annotated with exonerate using closest proteins
exonerate  --model protein2genome --percent 70 --softmasktarget yes  -n 1 --showtargetgff -t $file -q closely_related_protein.fa 

### 10. Map the cDNA to genome
gmap_build -d Braincoral -D /home/xinw/project/Brain_coral_new/Gene_prediction/Gmap/GMAP_index Braincoral_assembly_v1.fasta
/home/xinw/software/gmap-2017-01-14/bin/gmap -D /home/xinw/project/Brain_coral_new/Gene_prediction/Gmap/GMAP_index/ -d Braincoral Pdae_v1.fasta --min-trimmed-coverage=75 --min-identity=90 -f gff3_gene >Braincoral.gff3 2>log.txt &

#### 11. Find candidated coding regions from Transdcirpts
/home/xinw/software/TransDecoder/TransDecoder.LongOrfs -t Pdae_v1.fasta
/home/xinw/software/TransDecoder/TransDecoder.Predict -t Pdae_v1.fasta 

### 12. Use Portcullis to identify junction
portcullis prep Braincoral_assembly_novirus_v1.fasta Braincoral_sort.bam -o Brain_portcullis_prep
portcullis junc -t 30 -o Brain_portcullis_junc Brain_portcullis_prep >log.junc &
portcullis filter Brain_portcullis_prep/ portcullis_junc/portcullis.junctions.tab
portcullis filter --min_cov 2 Brain_portcullis_prep/ portcullis_junc/portcullis.junctions.tab

### 13. Use Mikado to combine different data

mikado configure --list configure.txt --reference Braincoral_assembly_novirus_v1.fasta --mode permissive --scoring worm.yaml  --copy-scoring worm.yaml --junctions junctions.bed -bt uniprot_sprot.fasta configuration.yaml
mikado prepare --json-conf configuration.yaml
makeblastdb -in uniprot_sprot.fasta  -dbtype prot -parse_seqids > blast_prepare.log
blastx -max_target_seqs 5 -num_threads 30 -query mikado_prepared.fasta -dd uniprot_sprot.fasta -evalue 0.000001 2> blast.log | sed '/^$/d' | gzip -c - > mikado.blast.xml.gz
mikado serialise --json-conf configuration.yaml --xml mikado.blast.xml.gz --orfs orfs.bed --blast_targets uniprot_sprot.fasta
mikado pick --json-conf configuration.yaml --subloci_out mikado.subloci.gff3
Statistic the gff files:
mikado util stats mikado.loci.gff3 mikado.loci.stats

### 14. Create hints from proteins

exonerate2hints.pl --in=Braincoral_exonerate_all.out --maxintronlen=20000 --source=P --minintron=31 --out=hints.protein.gff

### 15. Create hints from RepeatMasker 

perl -pe 's/^\s*$//' | perl -ne 'chomp; s/^\s+//; @t = split(/\s+/); print $t[4]."\t"."repmask\tnonexonpart\t".$t[5]."\t".$t[6]."\t0\t.\t.\tsrc=RM\n";' > hints.repeats.gff

### 16. Create hints from RNA-seq

bam2hints --in=Braincoral_sort.bam --out=hints.intron.gff --maxgenelen=30000 --intronsonly
Get exon hints:
bam2wig Braincoral_sort.bam >Braincoral_hisat.wig
cat Braincoral_hisat.wig | wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --UCSC=unstranded.track --radius=4.5 --pri=4 --strand="." > hints.exonpart.gff

### 17 Run PASA for trainingset

17.1 seqclean to clean transcript sequences
seqclean transcripts.fasta -v vectors.fasta
17.2 Run PASA
cp $PASA_HOME/pasa_conf/pasa.alignAssembly.Template.txt alignAssembly.config
perl -p -i -e 's/MYSQLDB=.*/MYSQLDB=sample_mydb_pasa/' alignAssembly.config
$PASA_HOME/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g genome_sample.fasta -t all_transcripts.fasta.clean -T -u all_transcripts.fasta -f FL_accs.txt --ALIGNERS blat,gmap --CPU 20
17.3 Filter PASA for best training sets
~/PASApipeline/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta Braincoral_mydb_pasa.pasa_assemblies.fasta --pasa_transcripts_gff3 Braincoral_mydb_pasa.pasa_assemblies.gff3
~/PASApipeline/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta ../Braincoral_mydb_pasa.pasa_assemblies.fasta --pasa_transcripts_gff3 ../Braincoral_mydb_pasa.pasa_assemblies.gff3 Braincoral_pasa_assemble.gff3

17.4 Check for complete ORF with at least 2 exons
perl ../selectComplete.pl  Braincoral_mydb_pasa.pasa_assemblies.fasta.transdecoder.genome.gff3 Braincoral_mydb_pasa.pasa_assemblies.fasta.transdecoder.gff3 >selected_models_2exon.gff3

17.5 blast against Uniprot to remove unfunction genome models
blastp -db uniprot_sprot.fasta -query select.protein_correctedID.fasta  -out select.protein_correctedID.blast.uniprot.out -evalue 1e-5 -outfmt 6 -num_threads 60
perl retain_uniprot4gff.pl  select.protein_correctedID.blast.uniprot.out selected_models_2exon_noredu.gff3 selected_models_2exon_noredu_uniprotID.gff3 >log.uniprot

pretraining with the trainingsets (476)
gff2gbSmallDNA.pl  selected_models_2exon_noredu_uniprotID.gff3 ../Braincoral_assembly_novirus_v1.fasta 100 trainingSetComple.gb
new_species.pl --species=platygyria_test
randomSplit.pl trainingSetComple.gb 50
etraining --species=platygyria_test trainingSetComple.gb.train
augustus --species=platygyria_test trainingSetComple.gb.test |tee firsttest.out 


### 18 use Glimmer to predict
Obtain exon file:
perl -ne '{chomp;my @array=split/\s+/,$_; if($array[2] eq "gene"){print "\n"}elsif ($array[2] eq "exon") {if ($array[6] eq "+"){print "$array[0] $array[3] $array[4]\n"}else{print "$array[0] $array[4] $array[3]\n"}}}' selected_models_2exon_noredu_uniprotID.gff3 >Braincoral_glimmer_train.mmtrain
training for glimmer:
~/GlimmerHMM/train/trainGlimmerHMM Braincoral_assembly_novirus_v1.fasta Braincoral_glimmer_train.mmtrain -d Glimmer_train
perl ~/script/auto_glimmerhmm.pl -g Braincoral_assembly_novirus_v1.fasta -d Glimmer_train -w /home/xinw/project/Briancoral/gene_prediction/Glimmer/ -o Braincoral_glimmer -t 35 -f -gff > log.glimmer &

### 19 Convert into EVM format:
Exonerate:
perl ~/exonerate2evm.pl -i Braincoral_exonerate_all.out -o exonerate.gff3
sed 's/exonerate:protein2genome:local/exonerate/' exonerate.gff3 -i
Mikado:
perl ~/script/mikado_2evm.pl -i ../mikado.loci.gff3  -o mikado.gff3 -n mikado.ncRNA.gff3
Cufflink:
perl ~/script/cufflinks_2evm.pl -i  ../Brain_cufflinks.gtf -o cufflinks.gff3
Stringtie:
perl ~/script/cufflinks_2evm.pl -i ../Braincoral_hisat_stringtie.gtf -o stringtie.gff3
Augustus:
perl ~/script/augustus_add_order.pl -i Braincoral_augustus.out -o Braincoral_augustus.gff3
grep "^#" ../Braincoral_augustus.gff3 -v |grep "gene|CDS|transcript" -E >augustus.gff3

Trinity:
grep "\bexon\b" ../Braincoral_mydb_pasa.pasa_assemblies.fasta.transdecoder.genome.gff3 >trinityGmap.gff3


sed 's/|size/_size/g' transcript_alignments.gff3 –i


### 20 Final intergrated with EVM
~/software/EVidenceModeler/EvmUtils/partition_EVM_inputs.pl --genome genome.fasta --gene_predictions gene_prediction.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcript_alignments.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out]
$EVM_HOME/EvmUtils/write_EVM_commands.pl --genome genome.fasta --weights `pwd`/weights.txt \
      --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 \
      --transcript_alignments transcript_alignments.gff3 \
      --output_file_name evm.out  --partitions partitions_list.out >  commands.list

$EVM_HOME/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log











