# MiSeq SOP tutorial from https://mothur.org/wiki/miseq_sop/


###### NOTES ############
# Operational taxonomic unit or OTU is considered as the basic unit used in numerical taxonomy. These units may refer to an individual, species, genus, or class.
# OTUs are cluster of similar sequence variants of the 16S rDNA marker gene sequence. Each of these cluster is intended to represent a taxonomic unit of a bacteria species or genus depending on the sequence similarity threshold. 
# A OTU table contains the number of sequences that are observed for each taxonomic unit (OTUs) in each samples. Columns usually represent samples and rows represent genera or species specific taxonomic units (OTUs).

# Alpha diversity: Variation of microbes in a single sample
# Beta diversity: Variation of microbial communities between samples

# What I did:
#in bash
ln silva.bacteria/silva.bacteria.fasta ./ 
# in mothur THIS IS WRONG?
#make.contigs(file=stability.files, processors=16)
#summary.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table)
#screen.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table, maxambig=0, maxlength=275, maxhomop=8)
#summary.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table)
#get.current()
#make.contigs(file=stability.files, maxambig=0, maxlength=275, maxhomop=8)
#summary.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table)
#unique.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table)
#summary.seqs(count=stability.trim.contigs.count_table)

### REDUCING SEQUENCING AND PCR ERRORS: ###

# make contigs command combines our two sets of reads for each sample and then to combine the data from all of the samples. 
# This command will extract the sequence and quality score data from your fastq files, create the reverse complement of the reverse read and then join the reads into contigs.
make.contigs(file=stability.files, processors=16)
summary.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table)
# This tells us that we have 152360 sequences that for the most part vary between 248 and 253 bases. Interestingly, the longest read in the dataset is 502 bp. Be suspicious of this. 
# Recall that the reads are supposed to be 251 bp each. This read clearly didn’t assemble well (or at all). 
# Also, note that at least 2.5% of our sequences had some ambiguous base calls. Finally, when we’ve previously looked at V4 sequence data we rarely/never see good sequences with a stretch where the same nucleotide is repeated more than 8 times. 
# We can use the maxambig, maxlength, and maxhomop options in make.contigs to resolve these issue while assembling the reads:
make.contigs(file=stability.files, maxambig=0, maxlength=275, maxhomop=8)

# Alternatively, we can take care of these issues in a separate step by running the screen.seqs command.
# This implementation of the command will remove any sequences with ambiguous bases and anything longer than 275 bp. 
screen.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table, maxambig=0, maxlength=275, maxhomop=8)
summary.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table)


### PROCESSING IMPROVED SEQUENCES ###
# We anticipate that many of our sequences are duplicates of each other. Because it’s computationally wasteful to align the same thing a bazillion times, we’ll unique our sequences using the unique.seqs command:
# If two sequences have the same identical sequence, then they’re considered duplicates and will get merged. In the screen output there are two columns - the first is the number of sequences characterized and the second is the number of unique sequences remaining.
unique.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table)
summary.seqs(count=stability.trim.contigs.count_table)

# Now we need to align our sequences to the reference alignment. Again we can make our lives a bit easier by making a database customized to our region of interest using the pcr.seqs command. To run this command you need to have the reference database (silva.bacteria.fasta) and know where in that alignment your sequences start and end.
# To remove the leading and trailing dots we will set keepdots to false. 
pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F)
rename.file(input=silva.bacteria.pcr.fasta, new=silva.v4.fasta)

# Now we have a customized reference alignment to align our sequences to.
align.seqs(fasta=stability.trim.contigs.unique.fasta, reference=silva.v4.fasta)
summary.seqs(fasta=stability.trim.contigs.unique.align, count=stability.trim.contigs.count_table)

# You’ll see that the bulk of the sequences start at position 1969 and end at position 11551. Some sequences start at position 1251 or 1969 and end at 10694 or 13401. These deviants from the mode positions are likely due to an insertion or deletion at the terminal ends of the alignments.
# Sometimes you’ll see sequences that start and end at the same position indicating a very poor alignment, which is generally due to non-specific amplification. To make sure that everything overlaps the same region we’ll re-run screen.seqs to get sequences that start at or before position 1969 and end at or after position 11551.
# We’ll also set the maximum homopolymer length to 8 since there’s nothing in the database with a stretch of 9 or more of the same base in a row
screen.seqs(fasta=stability.trim.contigs.unique.align, count=stability.trim.contigs.count_table, start=1969, end=11551)
summary.seqs(fasta=current, count=current)

# Now we know our sequences overlap the same alignment coordinates, we want to make sure they only overlap that region. So we’ll filter the sequences to remove the overhangs at both ends
# 
filter.seqs(fasta=stability.trim.contigs.unique.good.align, vertical=T, trump=.)
unique.seqs(fasta=stability.trim.contigs.unique.good.filter.fasta, count=stability.trim.contigs.good.count_table)
pre.cluster(fasta=stability.trim.contigs.unique.good.filter.unique.fasta, count=stability.trim.contigs.unique.good.filter.count_table, diffs=2)
chimera.vsearch(fasta=stability.trim.contigs.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.unique.good.filter.unique.precluster.count_table, dereplicate=t)
summary.seqs(fasta=current, count=current)
classify.seqs(fasta=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.fasta, count=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.count_table, reference=trainset9_032012.pds.asta, taxonomy=trainset9_032012.pds.tax)
remove.lineage(fasta=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.fasta, count=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.count_table, taxonomy=stability.trim.contig.unique.good.filter.unique.precluster.denovo.vsearch.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
summary.tax(taxonomy=current, count=current)
get.groups(count=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, fasta=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.fasta, groups=Mock)
seq.error(fasta=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.fasta, count=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, reference=HMP_MOK.v35.fasta, aligned=F)
dist.seqs(fasta=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.fasta, cutoff=0.03)
cluster(column=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.dist, count=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
make.shared(list=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.list, count=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=.03)
rarefaction.single(shared=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.shared)
remove.groups(count=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, fasta=stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.fasta, taxonomy=stability.trm.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pds.wang.pick.taxonomy, groups=Mock)
rename.file(fasta=current, count=current, taxonomy=current, prefix=final)
dist.seqs(fasta=final.fasta, cutoff=0.03)
cluster(column=final.dist, count=final.count_table)
cluster.split(fasta=final.fasta, count=final.count_table, taxonomy=final.taxonomy, taxlevel=4, cutoff=0.03)
make.shared(list=final.opti_mcc.list, count=final.count_table, label=0.03)
classify.otu(list=final.opti_mcc.list, count=final.count_table, taxonomy=final.taxonomy, label=0.03)
make.shared(count=final.count_table)
classify.otu(list=final.asv.list, count=final.count_table, taxonomy=final.taxonomy, label=ASV)
phylotype(taxonomy=final.taxonomy)
make.shared(list=final.tx.list, count=final.count_table, label=1)
classify.otu(list=final.tx.list, count=final.count_table, taxonomy=final.taxonomy, label=1)
dist.seqs(fasta=final.fasta, output=lt)
clearcut(phylip=final.phylip.dist)
#Analysis
count.groups(shared=final.opti_mcc.shared)
sub.sample(shared=final.opti_mcc.shared, size=2403)

# generate files ending in *.rarefaction, which again can be plotted in your favorite graphing software package.
# rarefaction is not a measure of richness, but a measure of diversity. If you consider two communities with the same richness, 
#but different evenness then after sampling a large number of individuals their rarefaction curves will asymptote to the same value. Since they have different evennesses the shapes of the curves will differ. 
#Therefore, selecting a number of individuals to cutoff the rarefaction curve isn’t allowing a researcher to compare samples based on richness, but their diversity.
rarefaction.single(shared=final.opti_mcc.shared, calc=sobs, freq=100)


summary.single(shared=final.opti_mcc.shared, calc=nseqs-coverage-sobs-invsimpson, subsample=T)
dist.shared(shared=final.opti_mcc.shared, calc=thetayc-jclass, subsample=t)
pcoa(phylip=final.opti_mcc.thetayc.0.03.lt.ave.dist)
nmds(phylip=final.opti_mcc.thetayc.0.03.lt.ave.dist)

# We can plot the three dimensions of the NMDS data by plotting the contents of final.opti_mcc.subsample.pick.thetayc.0.03.lt.nmds.axes. 
nmds(phylip=final.opti_mcc.thetayc.0.03.lt.ave.dist, mindim=3, maxdim=3)

amova(phylip=final.opti_mcc.thetayc.0.03.lt.ave.dist, design=mouse.time.design)

# these two didn't work
# These data can be plotted in what’s known as a biplot where lines radiating from the origin (axis1=0, axis2=0, axis3=0) to the correlation values with each axis are mapped on top of the PCoA or NMDS plots.
homova(phylip=stability.opti_mcc.thetayc.0.03.lt.ave.dist, design=mouse.time.design)
corr.axes(axes=stability.opti_mcc.thetayc.0.03.lt.ave.pcoa.axes, shared=stability.opti_mcc.0.03.subsample.shared, method=spearman, numaxes=3)


get.communitytype(shared=final.opti_mcc.0.03.subsample.shared)
metastats(shared=final.opti_mcc.0.03.subsample.shared, design=mouse.time.design)
#something weird happend with this command. It exited out of mothur
lefse(shared=final.opti_mcc.0.03.subsample.shared, design=mouse.time.design)


