**Master thesis**
Commands implemented during my master thesis on the topic understanding genomic features of parthenogenetic triploid nematodes

**PRE-PROCESSING of data**
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


