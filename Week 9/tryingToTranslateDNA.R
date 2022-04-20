# What is a genetic code?
#   From the National Human Genome Research Institute: 
#   Genetic code is the term we use for the way that the four bases of DNA--the A, C, G, and Ts--are strung together in a way that the cellular machinery, the ribosome, can read them and turn them into a protein. 
#   In the genetic code, each three nucleotides in a row count as a triplet and code for a single amino acid. So each sequence of three codes for an amino acid. And proteins are made up of sometimes hundreds of amino acids. 
#   So the code that would make one protein could have hundreds, sometimes even thousands, of triplets contained in it.

# search for this accession: KJ806268 in the NCBI (National Center for Biotechnology Information)
# Find a protein coding gene in the genome
#        exon            7733..7990
#                        /number=1
#        intron          7991..9409
#                        /number=1
#        gene            8082..9134
#                        /gene="cox1"
#                        /locus_tag="EL90_g03"
#        CDS             8082..9134
#                        /gene="cox1"
#                        /locus_tag="EL90_g03"
#                        /codon_start=1
#                        /transl_table=22
#                        /product="hypothetical LAGLIDADG homing endonuclease"
#                        /protein_id="AIK29129.1"
#                        /translation="MNKQLLLCGTRKGTPIKVLENGDHNKPKKVTNLRNTGINQYTKP
#                        GSVTLNASMDHLSVHNQELLGSYLADLVEADGCIICPKTFSKQAHIDICFCIADKPFA
#                        EFLQTLLGGNIVIGASKNFVKLKINQQHKLLRLCQLMNDFFRTPKIDVLHNLIKYFNE
#                        KYNANLMMKGLNTSPLCSNAWLAGFSDGDANFLLEISKPRSGRVYVRIQYRLRLAQIY
#                        AKCKSLQLECNSMYGICFAIAILFDCNLYSYTHQRALGHNKSLLKAYHSYSVVTTNIK
#                        SDHLVCSYFDNYPLFSSKRLNYLDWRCVYELRNFKQHYTPQGVAKISEIKARFNTNRT
#                        VFTWDHLSTFYAKTRE"
# then find its genetic code
#   22. Scenedesmus obliquus Mitochondrial Code (transl_table=22)
   AAs  = "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
#   Starts = ------*---*---*--------------------M----------------------------
   Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
   Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
   Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

DNA <- DNAStringSet(x=character(Base1), start=1, end=length(Base1), width=length(Base1), use.names=TRUE)
translate(DNA, genetic.code=GENETIC_CODE, if.fuzzy.codon="error")

#Translating the DNAString
dnastring = DNAString("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG")
translate(dnastring, genetic.code=GENETIC_CODE, if.fuzzy.codon="error")

#show the frequency of each character
alphabetFrequency(dnastring)

#reverse complement of the DNAString
reverseComplement(dnastring)

#Access the second element of the DNAString
dnastring[2]

#Subseting the DNAString
dnastring[7:12]
#A better way to subset
subseq(dnastring, 7, 12)

#Convert DNAString to a character vector
toString(dnastring)








#look at the Escherichia coli APEC O1 genome (NC 008563). This can be found in the BSgenome package (along with genomes from other model organisms). Download and install the package data:
install.packages("biocLite")
BiocManager::install("biocLite")
library("biocLite")
biocLite("BSgenome.Ecoli.NCBI.20080805")

