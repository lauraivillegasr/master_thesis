**Master thesis**
Commands implemented during my master thesis on the topic understanding genomic features of parthenogenetic triploid nematodes

**POPULATION ANALYSIS **

Part 1. **PRE-PROCESSING of data**
commands implemented on re-sequencing data from Illumina sequencing. 

make separate CHEOPS software accessible:
```module use /opt/rrzk/modules/experimental```

1. Trimming using fastp: 

```fastp -i forward.fq.gz -I reverse.fq.gz -o out.forward.fastq -O out.reverse.fastq```


2. Mapping to reference genome via bwa mem

	2.1 Create index for mapping: 
```bwamem2 index ref.fa```

	2.2 Mapping: 
```bwamem2 mem -M -t 30 -R “@RG\tID:sample-id\tSM:sample\tPL:ILLUMINA\tPU:1” /path/to/reference-genome.fasta out.forward.fastq out.reverse.fastq > sample-1_bwamem.sam``` 


3. Creating list for looping further analysis
  
 #create file list
```ls -1 | sed 's/_bwamem.sam//g' > list-XX```

4. Convert sam to bam and sort files


```ls -1 | sed 's/_bwamem.sam//g' > list-XX```
```while read f; do samtools view -b $f"_bwamem.sam" > $f".bam" ;done < list-XX```


```cat list-XX | parallel -j 12 'samtools sort -@ 4 -o {}bwamem.sort.bam {}.sam'```


5. Indexing sorted bam files 

 ```ls *.sort.bam | parallel samtools index '{}'```
	
6. quality statistics of your mappings via qualimap. Info on the tool  http://qualimap.conesalab.org/
	#create data-description-file for multi-bamqc
```ls -d -1 $PWD/** | grep sort.bam$ > qualimap.list  ### requires additional manual editing - I used nano editor```

	#example:

|NAME        | PATH	|GROUP |  
-------------|----------|-------|
|c12_JU_100|	/scratch/lvilleg1/MAL/c12_JU_100bwamem.sort_stats |      c12_JU_100|
|c12_JU_47|	/scratch/lvilleg1/MAL/c12_JU_47bwamem.sort_stats   |     c12_JU_47|
|c12_JU_60|	/scratch/lvilleg1/MAL/c12_JU_60bwamem.sort_stats    |    c12_JU_60|



	qualimap multi-bamqc -r -d qualimap.list --java-mem-size=2400M #mem size needs to be edited as default 1220M was not enough for the data used

7. Removing duplicates with picard tools. Info on the tool: https://broadinstitute.github.io/picard/

		while read f; do java -jar picard.jar MarkDuplicatesWithMateCigar I=$f".sort.bam" O=$f".sort.rmd.bam" M=$f"sort.bam.metrics" VALIDATION_STRINGENCY=SILENT MINIMUM_DISTANCE=300 REMOVE_DUPLICATES=true & done < list


8. Filtering the data to keep just quality higher than 30:

```samtools view -q 30 -b in.bam > aligned_reads.q30.bam```

9.A. Gatk haplotypecaller. Info on the tool: https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller


9.A.1 Indexing files after filtering
```ls *.bam | parallel samtools index '{}'```

9.A.2 Haplotype calling 
```gatk HaplotypeCaller -R reference_genomes/panagrolaimus_es5.PRJEB32708.WBPS15.genomic.fa -I filesq30/P_bromber.sort.rmd.q30.bam -O P-bromber.vcf```


**ADDITIONAL INFORMATION ON PRE-PROCESSING:**

1. For p_superbus I used NGM - reads shorter than 70bp. Info on the tool: https://github.com/Cibiv/NextGenMap/wiki

	1.1 Create index for mapping:
```ngm -r reference-sequence.fa```

 	1.2 Mapping: 
```ngm -r reference-sequence.fa -1 forward_reads.fq -2 reverse_reads.fq -o results.sam ```


2. For files with insert size smaller than double read lenght I used pear. Info on the tool: https://cme.h-its.org/exelixis/web/software/pear/doc.html#cl-main

```pear -f forwardread.sanfastq.gz -r reverseread.sanfastq.gz -o pearoutputname```
				
pear output gives: file.unassembled.forward.fastq, file.unassembled.reverse.fastq, file.assembled.forward.fastq file.assembled.reverse.fastq and discarded reads

Mapping for this files was done twice: once for assembly and once for unassembled reads and then they were merged. 

```bwamem2 mem -M -t 30 -R "@RG\tID:ASEX\tSM:PS1806\tPL:ILLUMINA\tPU:1" /path/to/reference-genome.fasta file.unassembled.forward.fastq file.unassembled.reverse.fastq > mapped_unassembledpear.sam```  #no extra indexing required, as the same index from previous mapping could be used

After re-mapping in either case, all steps from pre-processing were followed again. 

Part 2. **POPULATION ANALYSIS USING POPOOLATION:**

remember to specify -fastq-type!

The sorted final bam files are used as initial input for the tool. 

A. Creating pileup file for (A) individual files (further estimation of pi, theta) and (B) merged files (for Fst estimation)

	samtools mpileup P_bornheim.sort.rmd.q30.bam > P_bornheim.pileup
	samtools mpileup -B -b list-samtoolpileup_sex > sexpop.mpileup

B. Creating syncronized files for further estimations 

	java -ea -Xmx7g -jar popoolation2_1201/mpileup2sync.jar --input Sex_network/asexpop.mpileup --output Sex_network/asexpop_java.sync --fastq-type sanger —min-qual 30 --threads 4

B.1. Fst estimated on sliding window non overlapping —> to compare populations amongst each other

	perl popoolation2_1201/fst-sliding.pl --input Sex_network/asexpop_java.sync  --output Sex_network/asexpop_w1kb_corrected.fst  --suppress-noninformative --min-count 4 --min-coverage 10 --max-coverage 80 --min-covered-fraction 0,5 --window-size 1000 --step-size 1000 --pool-size 3000
	
About coverage ranges: for asex populations 10 - 80 and for sex pops 10-56 (meaning 23 is the covergae per gene copy, ase population are triploid whereas sexual ones are diploid)

B.1.2 Extracting values of Fst from all columns - needs to be done one at a time

	awk '{ gsub(/1:2=/,"", $6); print } ' filename > newfilename

1:2 means fst of pop1 vs pop2. The order of populations is the same as provided on 1(B)


B.2. Fst estimation on single positions to obtain common positions between all populations (asex as a group and sex as another)

	perl popoolation2_1201/fst-sliding.pl --input Sex_network/asexpop_java.sync  --output Sex_network/asexpop_w1kb_corrected.fst  --suppress-noninformative --min-count 4 --min-coverage 10 --max-coverage 80 --min-covered-fraction 0,5 --window-size 1 --step-size 1 --pool-size 3000

B.2.1 From fst per portion, grab only the first two columns that have the information on contain and position present in all pops

	awk '{print $1, $2 }' file > positions_asex
			
C.1. With he positions file obtained in B.2.1, common positions between all populations were extracted on the individual pileup files using: 

	awk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' positions_asex PS1159.pileup > PS1159.corrected.pileup

Total positions asex: 48704  
Total positions sex: 122138 

C.1.2 Watterson theta and pi estimation using corrected pileup files

	perl popoolation_1.2.2/Variance-sliding.pl --measure theta --input Sex_network/p_davidi.corrected.pileup --output Sex_network/p_davidi_WT.file --pool-size 3000 --min-count 2 --min-coverage 10 --max-coverage 80 --window-size 1 --step-size 1 --fastq-type sanger
 
	perl popoolation_1.2.2/Variance-sliding.pl --measure pi --input Sex_network/p_davidi.corrected.pileup --output Sex_network/P_davidi_pi.file --pool-size 3000 --min-count 2 --min-coverage 10 --max-coverage 80 --window-size 1 --step-size 1 --fastq-type sanger

C.1.3 Remove all “fake/empty” positions, created by population, to only have meaningful results

	awk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' positions_asex pi_PS1579 > PS1159.corrected.pi

C.1.4 Obtaining average values for watterson theta and pi 

	awk '{ total += $5; count++ } END { print total/count }' PS1159.pi


For visualisation of results as violin plots, maps or fst distribution, codes were generated using R (see files: violin_R_WT.R, w1kb_pi_violin_plot.R, pi_box_plot.R. pi_map_plot.R and coverage script.R)

Alternative theta estimation using Tetmer - This tool estimates theta on only homologous copies in triloid genomes, which allows us to see how the third copy porivdes more genetic diversity.

	kat hist -o DL137G2 Panagrolaimus_rawreads/DL137G2/SN7640087_3181_DL137G2_1_sequence.fq.gz Panagrolaimus_rawreads/DL137G2/SN7640087_3181_DL137G2_2_sequence.fq.gz 

Provide histogram to tetmer and obtain per-k-mer theta (if k=21, then you'd have to multiply your result times 21).

PArt 3. **GENE NETWORK USING BUSCO GENES**

 **ESTIMATION OF MUTATION RATES**

Installing accuMulate gave some problems. Remember to add the path of bamtools folder to my path using `export PATH=$PATH:/your/new/path/here`

For creating ref pool of PS1159 more reads were added, they had solexa encoding and were transformed to Sanger encoding using: 

	cat 110815_0133_D0738ACXX_5_SA-PE-005_1.solfastq | seqret -filter -auto -sformat fastq-solexa -osformat fastq-sanger -out PS1159_5_1.fastq

All reads coming from PS1159 where merged

We checked through fastqc to estimate/know wether the conversion worked: the result showed: SANGER :D

Pre-processing as done on part 1 from Population analysis. 
Trimming

	./fastp -i /scratch/a200302/L19G31/SN7640087_3184_L19G31_1_sequence.fq.gz  -I /scratch/a200302/L19G31/SN7640087_3184_L19G31_2_sequence.fq.gz -o /scratch/lvilleg1/L19G31_1 -O /scratch/lvilleg1/L19G31_1 -h /scratch/lvilleg1/report_L19G31

Mapping

	bwamem2 mem -M -t 30 -R "@RG\tID:ASEX\tSM:PS1159_c12_PS_8\tPL:ILLUMINA\tPU:1" /home/lvilleg1/reference_genomes/panagrolaimus_ps1159.PRJEB32708.WBPS15.genomic.fa.gz /scratch/lvilleg1/MAL_fastp/c12_PS_8_1 /scratch/lvilleg1/MAL_fastp/c12_PS_8_2 > /scratch/lvilleg1/c12_PS_8_bwamem.sam

Convert sam to bam

	ls -1 | sed 's/_bwamem.sam//g' > list-XX
	while read f; do samtools view -b $f"_bwamem.sam" > $f".bam" ;done < list-XX

Sort files 

	samtools sort -o /scratch/lvilleg1/MAL/L19G31.sort.bam scratch/lvilleg1/MAL/L19G31.bam #had to be run independently on each file, after a whole day command with parallel was not running

Learn about coverage from sorted files

	ls *.sort.bam | parallel -j 5 'samtools depth {} > {}.sort.bam.depth' 

Remove duplicates using picard


module purge #remove the current default java version
	module load openjdk/1.8.0_202 #need to be loaded, otherwise it tries to run with the incorrect java version and won’t work 


	java -jar /home/lvilleg1/picard/build/libs/picard.jar MarkDuplicatesWithMateCigar I=/scratch/lvilleg1/MAL/c12_PS_86.sort.bam O=c12_PS_86.sort.rmd.bam M=c12_PS_86.sort.bam.metrics VALIDATION_STRINGENCY=SILENT MINIMUM_DISTANCE=300 REMOVE_DUPLICATES=true

Remove low quality reads 

ls *.sort.rmd.bam | parallel 'samtools view -q 30 -b {} > {.}.q30’


Samtools flagstat can be used to check quality of mapping 
Indexing files for haplotype caller - maybe ignore

	ls *.sort.rmd.q30.bam | parallel samtools index '{}'

Creating dictionary of reference genome - maybe ignore

	/home/lvilleg1/gatk-4.2.0.0/gatk CreateSequenceDictionary -R /home/lvilleg1/reference_genomes/panagrolaimus_ps1159.PRJEB32708.WBPS15.genomic.fa


To check that the header is correctly stablished

	samtools view -H c12_JU_100.sort.rmd.bam | grep '^@RG'
	for f in *q30.bam ; do samtools view -H $f | grep '^@RG'; done



Merging files for accuMUlate 

	samtools merge -r partheno_merged.bamL19G31.sort.rmd.q30.bam PS83Q.sort.rmd.q30.bam c12_PS_22.sort.rmd.q30.bam c12_PS_8.sort.rmd.q30.bam c12_PS_84.sort.rmd.q30.bam c12_PS_86.sort.rmd.q30.bam

	c12_JU_100.sort.rmd.q30.bam c12_JU_47.sort.rmd.q30.bam c12_JU_60.sort.rmd.q30.bam c12_JU_71.sort.rmd.q30.bam c12_JU_73.sort.rmd.q30.bam c12_JU_88.sort.rmd.q30.bam


Prepare data for accumulate, obtain ini file and obtain GC content using accumulate tools (pre-requisite: pip3.6 install biopython) ALL THINGS THAT NEEDED PYTHON WHERE SUBMITTED TO CHEOPS0 - note on installing boost for accuMUlate: version 1.73 wouldn’t work, I used 1.62

	module load samtools
	module load python/3.6.8
	cd /scratch/lvilleg1/MAL/Final_preprocessing

	samtools view -H parthenoPS1159_merged.bam | python3 /home/lvilleg1/accuMulate-tools/extract_samples.py PS1159_refpool - >> params.PS1159.ini



	python3 /home/lvilleg1/accuMulate-tools/GC_content.py /home/lvilleg1/reference_genomes/PS1159_reference_genome >> params.PS1159.ini


NOTE: CHEOPS ONLY allows to install Biopython on python 3 using pip, however, the scripts of accu-tools are written for python 2. Had to do some editing on printing statement (adapt it to python3 - was written in python2):

Before: 
 	print "{}\t{}".format(*pair)
After:
 	print ("{}\t{}".format(*pair))

	python3 /home/lvilleg1/accuMulate-tools/dictionary_converter.py /home/lvilleg1/reference_genomes/panagrolaimus_ps1159.PRJEB32708.WBPS15.genomic.fa  > /home/lvilleg1/reference_genomes/panagrolaimus_ps1159.PRJEB32708.WBPS15.genomic.dict



FOLLOWING STEP ON CHEOPS1 WHERE BEDTOOLS IS AVAILABLE

	module load /opt/rrzk/modules/experimental/bedtools/2.29.2 (bedtools on cheops1)

	mkdir -p /scratch/lvilleg1/MAL/Final_preprocessing/tmp


	bedtools makewindows -g /home/lvilleg1/reference_genomes/panagrolaimus_ps1159.PRJEB32708.WBPS15.genomic.accu.dict -w 100000 | split -l 1 - /scratch/lvilleg1/MAL/Final_preprocessing/tmp


When using n=3 ```terminate called after throwing an instance of 'boost::program_options::invalid_option_value'
  what():  the argument ('accuMUlate can't only deal with haploid or diploid ancestral samples') for option is invalid```

Running accumulate: 

A change had to be done on the parsers.cc file from accuMUlate. 

1. As the error was coming from the condition of the if defined in line XXX not being fulfilled, we added a print statement that would show exactly what the error was: 

```std::cout << "HOLA SOY LAURA***********: " << start_index; -> prints the start_index that is problematic```
2. We add a statement that tells the script to ignore this specific index so the if condition can be fulfilled. 

```if (start_index != std::string::npos || start_index == 18446744073709551615) {```

What is 18446744073709551615? Is probably a value defined as a maximum by boost (when not specified differently), one of the tools used by accuMUlate. I think it is plausible that the error is this definition of maximum and not on the data itself as 4 different data sets where tested and the error persisted the same ´18446744073709551615´ 

AccuMUlate was then compiled again with the new “version” of the parsers.cc file.

The command implemented for obtaining candidate mutations was: 
On CHEOPS1

	parallel -j 6 home/lvilleg1/accuMUlate-0.2.1/build/accuMUlate -c /scratch/lvilleg1/accu_JU/params.JU765.ini -b /scratch/lvilleg1/accu_JU/hermaphroJU765_merged.sort.bam -r /scratch/lvilleg1/accu_JU/propanagrolaimus_ju765.PRJEB32708.WBPS15.genomic.fa -i {} ::: /scratch/lvilleg1/accu_JU/tmp/* > /scratch/lvilleg1/JU765_mutationcandidates #if -o is used to write the output, it goes to the .out file defined on the job batch script. 


After obtaining the candidate mutations, filtering stars: 

1. Define coverage range: 

With samtools depth we obtain a file with coverage at several positions, on R we can get the summary statistics for knowing the lower and upper range. We used 2 times the standard deviation of the ref pool for its upper limit: sd(file$V3)

Filtering was done following the code: 

The value for $11 changed according to the coverage range defined for the specific data set. 

	cat PS1159_mutationcandidates | awk '{if ($11 >=332 && $11 <=575 && $15 ==0 && $7 >=0.9 && $8 >=0.9 && $9 >=0.9 && $16 <=1.96 && $17 <=1.96 && $18 >=0.05 && $19 >=0.05) print $0}' > PS1159_mutationcandidates.filter-A.bed

	cat PS1159_mutationcandidates | awk '{if ($11 >=332 && $11 <=575) print $0}' | wc -l 

	cat PS1159_mutationcandidates | awk '{if ($11 >=332 && $11 <=575 && $15 ==0) print $0}' | wc -l

	cat PS1159_mutationcandidates | awk '{if ($11 >=332 && $11 <=575 && $15 ==0 && $7 >=0.9 && $8 >=0.9 && $9 >=0.9) print $0}' | wc -l

	cat PS1159_mutationcandidates | awk '{if ($11 >=332 && $11 <=575 && $15 ==0 && $7 >=0.9 && $8 >=0.9 && $9 >=0.9 && $16 <=1.96 && $17 <=1.96) print $0}' | wc -l


2. We obtained the number of callable sites for each of the lines. We used the already obtained depth files from samtools depth filename.bam > filename.depth

```cat L19G31.sort.bam.sort.bam.depth | awk '{if ($3 >=10 && $3 <=50) print $0}' | wc -l ```


Obtaining mutation rates and confidence intervals: 

1. For mutation rates the following equation was used for each of the mutation lines: 

	μ=(called mutations)/(generations∗callable sites)

2. Average of μ for each of the strains was calculated. (Can’t insert equation, basically all μ divided the total number of μ for the strain). 

3. Estimation of confidence intervals: 

Downloading Bayesian first aid on R

	install.packages("devtools")
	devtools::install_github("rasmusab/bayesian_first_aid")

NOTE: as I was working on a Mac computer, an error on installation or jags occurred (package required for bayesian_first_aid). I directly downloaded the package from https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/ and installed it before installing bayesian_first_aid. 


	JU_sex = c(1.62637E-09, 6.4472E-10, 1.60744E-09, 5.93941E-10, 1.17891E-09)
	PS_asex =c(9.64506E-10, 6.94089E-10, 4.08937E-10, 7.2004E-10, 7.9846E-10)
	comparing= c(9.41896E-10, 5.67965E-10)
	sites = c(278665301,374228169)
	bayes.poisson.test(comparing, sites)
	plot(bayes.poisson.test(comparing, sites))

Sites refers to the number of callable sites for each reproduction mode. The result of this analysis tells us how different our groups are and provides confidence intervals of the values: a lower limit and an upper limit.





