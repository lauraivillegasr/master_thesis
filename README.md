**Master thesis**
Commands implemented during my master thesis on the topic understanding genomic features of parthenogenetic triploid nematodes

POPULATION ANALYSIS 

A. **PRE-PROCESSING of data**
commands implemented on re-sequencing data from Illumina sequencing. 

make separate CHEOPS software accessible:
```module use /opt/rrzk/modules/experimental```

1. Trimming using fastp: 

```fastp -i forward.fq.gz -I reverse.fq.gz -o out.forward.fastq -O out.reverse.fastq```


2. Mapping to reference genome via bwa mem

	2.1 Create index for mapping: 
```bwa index ref.fa```

	2.2 Mapping: 
```bwa mem -M -t 30 -R “@RG\tID:sample-id\tSM:sample\tPL:ILLUMINA\tPU:1” /path/to/reference-genome.fasta out.forward.fastq out.reverse.fastq > sample-1_bwamem.sam``` 


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

```bwa mem -M -t 30 -R "@RG\tID:ASEX\tSM:PS1806\tPL:ILLUMINA\tPU:1" /path/to/reference-genome.fasta file.unassembled.forward.fastq file.unassembled.reverse.fastq > mapped_unassembledpear.sam```  #no extra indexing required, as the same index from previous mapping could be used

After re-mapping in either case, all steps from pre-processing were followed again. 

B. **POPULATION ANALYSIS USING POPOOLATION:**

remember to specify -fastq-type!

The sorted final bam files are used as initial input for the tool. 

1. Creating pileup file for (A) individual files (further estimation of pi, theta) and (B) merged files (for Fst estimation)

	(A)samtools mpileup P_bornheim.sort.rmd.q30.bam > P_bornheim.pileup
	(B)samtools mpileup -B -b list-samtoolpileup_sex > sexpop.mpileup

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
			
A.1. With he positions file obtained in B.2.1, common positions between all populations were extracted on the individual pileup files using: 

	awk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' positions_asex PS1159.pileup > PS1159.corrected.pileup

Total positions asex: 48704  
Total positions sex: 122138 

A.1.2 Watterson theta and pi estimation using corrected pileup files

	perl popoolation_1.2.2/Variance-sliding.pl --measure theta --input Sex_network/p_davidi.corrected.pileup --output Sex_network/p_davidi_WT.file --pool-size 3000 --min-count 2 --min-coverage 10 --max-coverage 80 --window-size 1 --step-size 1 --fastq-type sanger
 
	perl popoolation_1.2.2/Variance-sliding.pl --measure pi --input Sex_network/p_davidi.corrected.pileup --output Sex_network/P_davidi_pi.file --pool-size 3000 --min-count 2 --min-coverage 10 --max-coverage 80 --window-size 1 --step-size 1 --fastq-type sanger

A.1.3 Remove all “fake/empty” positions, created by population, to only have meaningful results

	awk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' positions_asex pi_PS1579 > PS1159.corrected.pi

A.1.4 Obtaining average values for watterson theta and pi 

	awk '{ total += $5; count++ } END { print total/count }' PS1159.pi


For visualisation of results as violin plots, maps or fst distribution, codes were generated using R (see files: violin_R_WT.R, w1kb_pi_violin_plot.R, pi_box_plot.R. pi_map_plot.R and coverage script.R)

Alternative theta estimation using Tetmer - This tool estimates theta on only homologous copies in triloid genomes, which allows us to see how the third copy porivdes more genetic diversity.

	kat hist -o DL137G2 Panagrolaimus_rawreads/DL137G2/SN7640087_3181_DL137G2_1_sequence.fq.gz Panagrolaimus_rawreads/DL137G2/SN7640087_3181_DL137G2_2_sequence.fq.gz 

Provide histogram to tetmer and obtain per-k-mer theta (if k=21, then you'd have to multiply your result times 21).



