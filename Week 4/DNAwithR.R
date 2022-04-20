setwd("/Users/natalieburkhard/bio-490")

# If you don't have the bioconductor package suite installed, run the following

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install()
# the above only if you want to install every package available to Bioconductor

BiocManager::install(c("GenomicFeatures", "AnnotationDbi", "Biostrings"))

# again, if you already have these, don't install them again:

install.packages("rentrez")
install.packages("seqinr")

library(BiocManager)

library("genbankr")

library("Biostrings")

library("seqinr")
library("rentrez")
# seqinr and Biostrings both have the function translate; this creates issues
# we'll load these packages later

myDNA = c("ATGTTGCATTCATTAATTAAGAACGACCCAATACA")
myDNASeq = DNAString("CTGATTT-GATGGTC-NAT")
myDNASeq2 = DNAString("ATGTTGCATTCATTAATTAAGAACGACCCAATACA")

length(myDNA)
length(myDNASeq2)
# what is the difference between the two?
# The length of myDNA is 1 and the length of myDNAseq2 is 35

# next, use the class() command to see what type of data myDNA and MyDNASeq are
# compose the line below and in comments report the outcome
class(myDNA)
class(myDNASeq2)
# myDNA is character and myDNAseq2 is a DNAString which is an attribute in the Biostrings package

# next, try some basic manipulations and describe what happens
# use complement(), reverse(), and reverseComplement() on the DNA objects you created
complement(myDNASeq)
complement(myDNASeq2)

reverse(myDNASeq)
reverse(myDNASeq2)

reverseComplement(myDNASeq)
reverseComplement(myDNASeq2)

# next, put the reverse complement of myDNASeq2 into a new object called myDNASeq3
# what type (class) of data is myDNASeq3?
# myDNASeq3 is a DNAString
myDNAseq3 = reverseComplement(myDNASeq2)
class(myDNAseq3)

# try the following for all of the DNAStrings objects you've made - what information do you get?
alphabetFrequency(myDNASeq)
alphabetFrequency(myDNASeq2)
alphabetFrequency(myDNASeq3)

# next we'll try to translate our DNA into protein
# of course, that's not how it works in cells, but here we're conveniently skipping transcription
# use translate() on your DNA objects from above
# what if we wanted to use a different genetic code?

GENETIC_CODE_TABLE
# what did this do?

# there are a couple of ways to translate a sequence x under different genetic codes
# you can use getGeneticCode(id_or_name2 = "") with the codes name or number in the quotes
# you can place your genetic code into a new object, let's call it y, then put that object into the following line:
# translate(x, genetic.code = y)
# that would take two lines. can you combine them to accomplish the same in a single line?


# try translating myDNASeq2 under a few genetic codes. Any differences?

# if this didn't work, you might have seqinr loaded
# unload it like you've learned in your homework, or unclick it in Rstudio in the Packages tab, lower right window

# another way to translate a simple string like myDNA is to load the seqinr package; read the next few lines before trying
# the translate function in seqinr works a little different, and doesn't handle DNAString objects
# it handles simple text like myDNA
# there is a way to convert string to characters - or s2c, which can be integrated into the translate command:
# translate(s2c(myDNA))
# or if the DNA object is already a list of separate nucleotides - "a" "c" "c" "t" etc, you can just do translate(dnaobject)
# seqinr also uses a different way to call up genetic codes, by their number: translate(dnaobject, numcode=15)
# feel free to try this later, but this is the reason why loading Biostrings and seqinr at the same time is not a good idea

# now we'll load seqinr and explore its functions
library(seqinr)
# note the warning about the object translate being masked from Biostrings. UGH!

# next, read in a fasta file (make sure you have it downloaded and modify the line below to make sure R can find it)
# or use the file.choose() function
dengue = read.fasta(file="DNA_with_R/NC_001477.fasta")
dengue
# what is in it?

dengueseq = dengue[[1]]
dengueseq

# what's the difference between dengue and dengueseq?
# you might need to google some

length(dengueseq)
denguetable = table(dengue)
# how long is the dengue genome? what are the units of length here?
# what's in denguetable?

# another useful thing is looking at just subsets of data, not the whole genome
# because dengue has the list of nucleotides as one item, it's most convenient to do this on dengueseq
dengueseq[1:50]
# what did this line do?

# translate the whole dengue genome
# then translate only the first 99 nucleotides
# there are a few ways to do this, depending on which packages you use and how you've prepared your data
# find at least two different ways
# finally, translate the first 99 nucleotides of the dengue genome under a different genetic code

### next up we'll do some brief summary statistics on the genome
# genomes are complex - there is not a simple way to describe them in a few numbers
# one way is to look at over- and underrepresented "words" in the genome
count(dengueseq, 2)
# what does this do?
# play with the number in the line for a bit

# we can for example calculate the over- or underrepresentedness of the GC "word" given the frequency of "G" and "C" in the genome
# 𝜌(GC) = fGC/(fG * fC)- the actual formula for the rho statistic, f is for frequency (count per length of sequence)
denguetable[["g"]]
# gets you the frequency of a particular nucleotide
# similarly you can pick out parts of count(dengueseq, x)
count(dengueseq, 2)[["gc"]]
# have a go at calculating rho for "GC", according to the formula above
# it might take a few steps the length function may come in handy

# if I calculated correctly, rhoGC should be about 0.86, so slightly underrepresented given the frequency of G and C


# GC content is an important and commonly used gene and genome descriptive statistic
# any idea why GC matters? (think biology, DNA structure, base pairing)

# let's look into it
denguegc = (denguetable[["g"]] + denguetable[["c"]])/length(dengueseq)
denguegc == GC(dengueseq)
# what does the last line mean? What result does it return?

# furthermore, GC can differ in different parts of the genome
GC(dengueseq[1:2000])
GC(dengueseq[2001:4000])
GC(dengueseq[4001:6000])
GC(dengueseq[6001:8000])
GC(dengueseq[8001:10000])
GC(dengueseq[10001:10735])
# the above approach is a rather labor-intensive way to scan the GC content along the genome
# this is a better way to create such a "sliding window"
# below also note that comments can be on the code lines, after the actual code - useful!

starts = seq(1, length(dengueseq)-2000, by = 2000)
starts
n = length(starts) # Find the length of the vector "starts"
for (i in 1:n) {
  chunk = dengueseq[starts[i]:(starts[i]+1999)]
  chunkGC = GC(chunk)
  print (chunkGC)
}
# this particular set of lines enclosed in curly brackets is called a for-loop - loops are SUPER useful
# loops allow us to perform the same set of tasks on a list of objects
# see how the loop goes through our list n
# and for every item in the list it makes a corresponding DNA chunk, which is 2000 nucleotides long,
# calculates GC content for that chunk
# and prints it - aka shows it in the console for us to see
# does GC content vary across the 2000nt windows along the dengue genome?


# it's been a while since we've graphed something
# let's plot the GC content along the dengue virus genome

n = length(starts) # Find the length of the vector "starts"
# what is the length of starts? (check your answer and make sure it makes sense to you)
chunkGCs <- numeric(n) # Make a vector of the same length as vector "starts", but just containing zeroes
for (i in 1:n) {
  chunk = dengueseq[starts[i]:(starts[i]+1999)]
  chunkGC = GC(chunk)
  print(chunkGC)
  chunkGCs[i] = chunkGC # a new element - saves the GC content numbers into a new object, which we'll use below
}

# similar for loop as before, but now we'll graph the output
plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")
# remember what all these plot parts mean? annotate your code with your explanation
# an easy way to find out what the type="b" argument means is to remove it and plot again - what changed?

# now let's create a FUNCTION - a longer, more complex script that can again do repetitive tasks for us
# our function will be called slidingwindowplot - but you could name it anything you like
# this way you can call up a set of operations by just calling slidingwindowplot, without having to re-type everything in it
# our function slidingwindowplot will take in two inputs - windowsize and inputseq
# note the elements we're taking from the for loop above

slidingwindowplot = function(windowsize, inputseq)
{
  starts = seq(1, length(inputseq)-windowsize, by = windowsize)
  n = length(starts) # Find the length of the vector "starts"
  chunkGCs <- numeric(n) # Make a vector of the same length as vector "starts", but just containing zeroes
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
    chunkGC <- GC(chunk)
    print(chunkGC)
    chunkGCs[i] <- chunkGC}
  plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")
}

# let's test it:
slidingwindowplot(500, dengueseq)
# what does the above line do? Describe it in a sentence or two

# describe how/why this function is useful
# try a few more window sizes

# final challenge - get the fasta sequence of the Zika genome from NCBI (accession AY632535.2), save it into a fasta file
# to read it in without issues, it's best to download using the Send to - file option in NCBI
# which genome is longer - zika or dengue?
# compute the overall GC content, translate the genome, and use your slidingwindowplot function to plot GC by 200nt