#2. The Genome
#All living organisms contain the instructions for life in their genome, which is present in the nuclei of their cells. The genome is comprised of double-stranded DNA divided into chromosomes; for humans there are 23, but different organisms will have differing numbers of chromosomes. The building blocks of our DNA are called nucleotides, and there are four different nucleotide bases in DNA: guanine, adenine, cytosine, and thymine. We will refer to these nucleotides as G, A, C, and T.

#3. Nucleotides
#The double-stranded DNA forms a helix with a sugar-phosphate backbone, and within this helix, A nucleotides pair with T and G nucleotides pair with C. The order of these nucleotides is called the DNA sequence.

#4. Genes
#Within this sequence are regions called genes. Genes provide instructions to make proteins, which perform some function within the cell. To make proteins, the DNA is transcribed into messenger RNA, or mRNA, which is translated by the ribosome into protein. Some genes encode RNA that does not get translated into protein; these RNAs are called non-coding RNAs, or ncRNAs. Often these RNAs have a function in and of themselves and include rRNAs, tRNAs, and siRNAs, among others. All RNAs transcribed from genes are called transcripts.

#5. RNA processing
#To be translated into proteins, mRNA must undergo processing. In this figure, the top strand in the image represents a gene in the DNA, comprised of the untranslated regions (UTRs), highlighted in blue, and the open read frame, highlighted in red. Genes are transcribed into pre-mRNA, which still contains the intronic sequences. Transcription represents the blue portion of the image. After post-transcriptional processing, shown in the grey section of the image, the introns are spliced out and a polyA tail and 5' cap are added to yield mature mRNA transcripts. The mature mRNA transcripts can be translated into protein, shown in the red portion of the image. While mRNA transcripts have a polyA tail, which is a sequence of As at the end of the transcript, many of the non-coding RNAs do not.

#6. Gene expression in cells
#Although all cells contain the same DNA sequence, muscle cells are different from nerve cells and other types of cells because of the different genes that are turned on in these cells and the different RNAs and proteins produced.

#7. Gene expression in disease
#Similarly, a disease-causing mutation can lead to differences in what genes are turned on, or expressed, and which genes are turned off. A mutation can affect the type and quantity of RNAs and proteins produced. To explore the gene expression changes that occur in disease or between different conditions, it can be useful to measure the quantity of RNA expressed by all genes using RNA-Seq. Then, differential expression analysis of RNA-Seq data can be used to determine whether there are significant differences in gene expression between conditions.

#8. RNA-Seq questions
#Using differential expression analyses, we can ask various questions, including: Which genes are differentially expressed between sample groups? Are there any trends in gene expression over time or across conditions? Which groups of genes change similarly over time or across conditions? What processes or pathways are enriched for my condition of interest?


# Load library for DESeq2
install.packages("Deseq2")
BiocManager::install("DESeq2") 
library(DEseq2)

# Load library for RColorBrewer
install.packages("RColorBrewer")
library(RColorBrewer)

# Load library for pheatmap
install.packages("pheatmap")
library(pheatmap)

# Load library for tidyverse
install.packages("tidyverse")
library(tidyverse)

#1. RNA-Seq Workflow
#Now that you know a bit about the types of questions that RNA-Seq experiments can address, and how we use this technique to understand more about the genes important to a particular disease or condition, let's explore the steps required for the analysis workflow.

#2. RNA-Seq Workflow: RNA-Seq Experimental Design
#Prior to starting the RNA-Seq workflow, planning is essential. This step in the analysis is crucial for good results, as there is often no saving a poorly designed experiment. There are a couple of important considerations during planning, including replicates, batch effects, and confounding: For RNA-Seq experiments there is generally low technical variation, so invest in biological replicates instead. The more biological replicates you have, the better the estimates are for mean expression and variation, leading to more robust analyses; be sure to have at least 3. Also, an experiment performed as different batches can confound your analysis. As much as possible try to perform experimental steps across all conditions at the same time, and if you cannot avoid batches, distribute the samples from each sample group into each batch. Finally, avoid confounding your experiment with major sources of variation. For example, if your animals are of different sexes, don't have all male mice as control and all female mice as treatment, as you won't be able to differentiate the treatment effect from the effect of sex.

#3. RNA-Seq Workflow: Sample prep
#The first step in a successful RNA-Seq DE analysis is a well-planned experiment. A well-planned experiment should avoid batch effects, should divide known major sources of variation, such as different sexes or ages, equally between sample groups, and should have a good number of biological replicates, preferably more than 3. The more biological replicates we have, the better our ability to detect DE genes with better estimates of mean expression and variation. Also, if we need to remove an outlier sample, we will still have biological replicates for the analysis.

#4. RNA-Seq Workflow: Sample prep
#When preparing RNA-Seq libraries, the samples are harvested, the RNA is isolated and DNA contamination is removed. The rRNA is removed or mature mRNAs are selected by their polyA tails. Then, the RNA is turned into cDNA, fragmented, size selected and adapters are added to generate the RNA-Seq libraries to be sequenced. The sequencing generates millions of nucleotide sequences called reads. The reads correspond to ends of the fragments sequenced. The sequence of each read is output into FASTQ files.

#5. RNA-Seq Workflow: Quality control
#With the sequenced reads in the FASTQ files, a series of analytical steps is performed on the command line, beginning with the assessment of the raw data quality.

#6. RNA-Seq Workflow: Quality control
#At this step, we ensure something didn't go wrong at the sequencing facility and explore the data for contamination, such as vector, adapter, or ribosomal.

#7. RNA-Seq Workflow: Alignment
#The next step is alignment or mapping of the reads to the genome to determine the location on the genome where the reads originated.

#8. RNA-Seq Workflow: Alignment
#Since mRNA contains only the exons needed to create the proteins, when the mRNA is aligned to the genome containing introns, some of the reads will be split across introns. Therefore, tools for aligning reads to the genome need to align across introns for RNA-seq. The output of alignment gives the genome coordinates for where the read most likely originated from in the genome and information about the quality of the mapping.

#9. RNA-Seq Workflow: Quantitation
#Following alignment, the reads aligning to the exons of each gene are quantified to yield a matrix of gene counts.

#10. RNA-Seq Workflow: Count matrix
#We can read into R the count matrix using the `read.csv()` function and specifying the file. The gene count matrix is arranged with the samples as columns and gene IDs as rows. The count values represent the number of reads or fragments aligning to the exons of each gene.

#11. RNA-Seq Workflow: Differential expression
#Once we have count data, differential expression analysis is performed. With differential expression analysis, our goal is to determine whether the gene counts between the sample groups are significantly different given the variation in the counts within the sample group.

#12. RNA-Seq Workflow: DE results
#The output of the statistical analysis includes the log2 fold changes of expression between conditions and the adjusted p-values for each gene. Genes that reach a threshold for significance can be subset to define a list of significant differentially expressed genes.

#13. RNA-Seq Workflow
#Now that we have a general understanding of the workflow and have the counts file loaded, we can get started.

# Explore the first six observations of smoc2_rawcounts
###DONT HAVE THIS DATASET
head(smoc2_rawcounts)

# Explore the structure of smoc2_rawcounts
str(smoc2_rawcounts)