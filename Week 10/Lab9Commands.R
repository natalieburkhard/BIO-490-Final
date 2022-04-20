### ABOUT THE DATA 
# The Schloss lab is interested in understanding the effect of normal variation in the gut microbiome on host health. To that end we collected fresh feces from mice on a daily basis for 365 days post weaning
# During the first 150 days post weaning (dpw), nothing was done to our mice except allow them to eat, get fat, and be merry. 
# We were curious whether the rapid change in weight observed during the first 10 dpw affected the stability microbiome compared to the microbiome observed between days 140 and 150.
# you are given the flow files for one animal at 10 time points (5 early and 5 late).

# set up a new directory for this lab
# remember, your path might look different - modify the lines as needed
#mkdir ~/Desktop/BIO260/Lab9
# what does the next line do?
# cd ~/Desktop/BIO260/Lab9

# download and unzip the latest Mac version of Mothur from https://github.com/mothur/mothur/releases into the new folder
# you can manually drag and drop it from your Downloads into Lab9, or use the cp command in Terminal
# double-clicking will unzip the files, or you can use the unzip command in Terminal

# then proceed to obtain the Schloss lab data from Brightspace - MiSeqSOPData.zip
# unpack it (you can also double-click on the icon in Finder)

#unzip miseqsopdata.zip
# copies all files from ./MiSeq_SOP/ that start with F3D0
cp ./MiSeq_SOP/F3D0* ./
  
  # Download SILVA bacterial database from Brightspace - it's a starting alignment of sequences belonging to known bacteria
  #  unzip silva.bacteria.zip
  # makes a hard link to the same file to avoid wasting space
  ln silva.bacteria/silva.bacteria.fasta ./  
  
  
  # Download RDP training set formatted for mothur from Brightspace
  #  unzip trainset9_032012.pds.zip
  
  
  # Phylotyping - Identifying reads by comparing them to a database of known sequences
  
  # make a new file called stability.files to hold the names of our fastq files.
  touch stability.files

# put the following in the file using nano: F3D0	<name of first file>	<name of second file>
nano stability.files

# run the mothur command line interface 
mothur/mothur

# reduce error by creating contigs from your forward and reverse reads. the more they overlap the better.
# we want to combine our two sets of reads for each sample and then to combine the data from all of the samples.
# This command will extract the sequence and quality score data from your fastq files, create the reverse complement of the reverse read and then join the reads into contigs.
#  In the end it will tell you the number of sequences in each sample
make.contigs(file=stability.files)


# summarize the output contigs
summary.seqs(fasta=stability.trim.contigs.fasta)
# This command shows us that there are almost 8,000 sequences that were processed with the previous command, assembled into contigs of varying lengths, but with most being about 252 bases long
# This command also put all that in a file

# The next step is to trim sequences that don't match expected lengths
# The short 16s gene is expected to always be the same length (275), so anything a lot shorter or longer is porbably bad
screen.seqs(fasta=stability.trim.contigs.fasta, maxambig=0, maxlength=275)


# the output files from screen.seqs:
# We got two files, one with bad sequences (too short or too long) and good ones (the expected length)

stability.trim.contigs.good.fasta   # good contigs
stability.trim.contigs.bad.accnos   # accession numbers for bad contigs

# show current files that mothur has remebered 
get.current()

# it is computationally wasteful to work with identical sequences, so get rid of them (mothur will remember that they exist, b/c we do want to know the abundance of each type of bacteria)
unique.seqs(fasta=stability.trim.contigs.good.fasta)

# count up the sequences to make future calculations easier. "generate a table where the rows are the names of the unique sequences and the columns are the names of the groups. The table is then filled with the number of times each unique sequence shows up in each group."
count.seqs(name=stability.trim.contigs.good.names)

# summarize again
summary.seqs(count=stability.trim.contigs.good.count_table)

# we want a good database to align our sequences to, the smaller the better. this will output a trimmed database starting with a full-length 16S database. 
# if you have your primer sequences you can make your own database. here we use the start and end positions for the primers used by the Schloss lab.
pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F)

# rename the output database using the system command mv
system(mv silva.bacteria.pcr.fasta silva.v4.fasta)
# show the summary of the database
summary.seqs(fasta=silva.v4.fasta)

# let's align our sequences to the database
align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.v4.fasta)
# and show the summary
summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table)

# clean up 
# Some sequences start at a different position than they should 
# To make sure that everything overlaps the same region we’ll re-run screen.seqs to get sequences that start at or before position 1969 and end at or after position 11551.
# We’ll also set the maximum homopolymer length to 8 since there’s nothing in the database with a stretch of 9 or more of the same base in a row 
# Note that we need the count table so that we can update the table for the sequences we’re removing:
screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.contigs.good.unique.summary, start=1968, end=11550, maxhomop=8)
summary.seqs(fasta=current, count=current)

# filter sequences to remove any bases that don't match the reference database
filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)
# Because we’ve perhaps created some redundancy across our sequences by trimming the ends, we can re-run unique.seqs:
unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)

# precluster sequences to remove some sequencing error
# This command will split the sequences by group and then sort them by abundance and go from most abundant to least and identify sequences that are within 2 nt of each other. If they are then they get merged.
pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=2)

# remove chimeras
chimera.vsearch(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
# if no chimeras were found, this will not do a whole lot


# classify sequences using the RDP training set downloaded earlier. The cutoff of 80 selects only the best matches. If doing phylotyping prior to OTU clustering it is probably a good idea to set this much lower to catch all of the matches.
#### classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=trainset9_032012.pds/trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds/trainset9_032012.pds.tax, cutoff=80)
### if no chimeras were found:
classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, reference=trainset9_032012.pds/trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds/trainset9_032012.pds.tax, cutoff=80)


# make a tree using 
clearcut(fasta=<name of aligned fasta file>, DNA=T)
# instead of the looong name of the fasta file, you can also just say current
system(ls)
# the .tre file is the file to load into Mega or Figtree

##when ready to quit mothur, type 
quit


