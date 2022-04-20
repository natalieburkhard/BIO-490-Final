# Genomic annotation in Bioconductor 
# tutorial from http://genomicsclass.github.io/book/pages/bioc1_annoOverview.html
# goal: explore the Bioconductor package in R

# Background info: Bioconductor includes many different types of genomic annotation. 
# We can think of these annotation resources in a hierarchical structure:
#   At the base is the reference genomic sequence for an organism. 
#   This is always arranged into chromosomes, specified by linear sequences of nucleotides.

#   Above this is the organization of chromosomal sequence into regions of interest. 
#   The most prominent regions of interest are genes, but other structures like SNPs or CpG sites are annotated as well.
#   Genes have internal structure, with parts that are transcribed and parts that are not, and “gene models” define the ways in which these structures are labeled and laid out in genomic coordinates.

#       Within this concept of regions of interest we also identify platform-oriented annotation. 
#       This type of annotation is typically provided first by the manufacturer of an assay, but then refined as research identifies ambiguities or updates to initially declared roles for assay probe elements. 

#   Above this is the organization of regions (most often genes or gene products) into groups with shared structural or functional properties. 
#   Examples include pathways, groups of genes found together in cells, or identified as cooperating in biological processes.

# installing bioconductor and BSgenome

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("bioconductor")

BiocManager::install("BSgenome")

# Doesn't work
#asource("http://bioconductor.org/biocLite.R")
#biocLite()

# Bioconductor’s collection of annotation packages brings all elements of this hierarchy into a programmable environment.
# Reference genomic sequences are managed using the infrastructure of the Biostrings and BSgenome packages, 
# and the available.genomes function lists the reference genome build for humans and various model organisms now available.
library(Biostrings)
library(bioconductor)
library(BSgenome)
ag = available.genomes()
length(ag)

# we can see there is 106 genomes available 
# lets display the first couple using head()
head(ag)

# A reference genomic sequence for H. sapiens
# The reference sequence for Homo sapiens is acquired by installing and attaching a single package. 
# This is in contrast to downloading and parsing FASTA files. 
# The package defines an object Hsapiens that is the source of chromosomal sequence, but when evaluated on its own provides a report of the origins of the sequence data that it contains.
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens
head(genome(Hsapiens))

# Get the chromosome's sequence using the $ operator 
Hsapiens$chr17


# The transcripts and genes for a reference sequence

# UCSC annotation
# The TxDb family of packages and data objects manages information on transcripts and gene models. 
# We consider those derived from annotation tables prepared for the UCSC genome browser.
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# abbreviate
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene 
txdb

# We can use genes() to get the addresses of genes using Entrez Gene IDs.
ghs = genes(txdb)
ghs

# Filtering is permitted, with suitable identifiers. 
# Here we select all exons identified for two different genes, identified by their Entrez Gene ids:
exons(txdb, columns=c("EXONID", "TXNAME", "GENEID"),
      filter=list(gene_id=c(100, 101)))


# ENSEMBL annotation
# Ensembl creates, integrates and distributes reference datasets and analysis tools that enable genomics
# The ensembldb package provides functions to create and use transcript centric annotation databases/packages.
# The functionality and data is similar to that of the TxDb packages from the GenomicFeatures package, but, 
# in addition to retrieve all gene/transcript models and annotations from the database, 
# the ensembldb package provides also a filter framework allowing to retrieve annotations for specific entries like genes encoded on a chromosome region or transcript models of lincRNA genes. 
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v75")
library(ensembldb)
library(EnsDb.Hsapiens.v75)
names(listTables(EnsDb.Hsapiens.v75))

edb = EnsDb.Hsapiens.v75  # abbreviate
txs <- transcripts(edb, filter = GeneNameFilter("ZBTB16"),
                   columns = c("protein_id", "uniprot_id", "tx_biotype"))
txs

# import/export
# Your data will be someone else’s annotation
# As an example, we consider estrogen receptor (ER) binding data, published by ENCODE as narrowPeak files. 
# This is ascii text at its base, so can be imported as a set of textual lines with no difficulty. 
# If there is sufficient regularity to the record fields, the file could be imported as a table.
# We want to go beyond this, so that the import is usable as a computable object as rapidly as possible. 
# Recognizing the connection between the narrowPeak and bedGraph formats, we can import immediately to a GRanges.
# To illustrate this, we find the path to the narrowPeak raw data file in the ERBS package.
# !!! Find thar file "exdata" and then this will work (hopefully) !!!
f1 = dir(system.file("extdata",package="ERBS"), full=TRUE)[1]
readLines(f1, 4) # look at a few lines

# import command
BiocManager::install("rtracklayer")
library(rtracklayer)
imp = import(f1, format="bedGraph")
imp

genome(imp)  # genome identifier tag not set, but you should set it

#Exporting
# We obtain a GRanges in one stroke. There are some additional fields in the metadata columns whose names should be specified, 
# but if we are interested only in the ranges, we are done, with the exception of adding the genome metadata to protect against illegitimate combination with data recorded in an incompatible coordinate system.
# For communicating with other scientists or systems we have two main options. We can save the GRanges as an “RData” object, easily transmitted to another R user for immediate use. 
# Or we can export in another standard format. For example, if we are interested only in interval addresses and the binding scores, it is sufficient to save in “bed” format.
export(imp, "demoex.bed")  # implicit format choice
cat(readLines("demoex.bed", n=5), sep="\n")

# We have carried out a “round trip” of importing, modeling, and exporting experimental data that can be integrated with other data to advance biological understanding.


# Annotation hub
# The AnnotationHub package can be used to obtain GRanges or other suitably designed containers for institutionally curated annotation.
BiocManager::install("AnnotationHub")
install.packages("httpuv")
library(AnnotationHub)
ah = AnnotationHub()
ah

# There are a number of experimental data objects related to the HepG2 cell line available through AnnotationHub.
query(ah, "HepG2")

# The query method can take a vector of filtering strings. To limit response to annotation resources addressing the histone H4K5, simply add that tag:
query(ah, c("HepG2", "H4K5"))



# The OrgDb Gene annotation maps
# Packages named org.*.eg.db collect information at the gene level with links to location, protein product identifiers, KEGG pathway and GO terms, PMIDs of papers mentioning genes, and to identifiers for other annotation resources.
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db) # columns() gives same answer
head(select(org.Hs.eg.db, keys="ORMDL3", keytype="SYMBOL", 
            columns="PMID"))



# Resources for gene sets and pathways
# Gene Ontology (GO) is a widely used structured vocabulary that organizes terms relevant to the roles of genes and gene products in
#   biological processes,
#   molecular functions, and
#   cellular components. 
# The vocabulary itself is intended to be relevant for all organisms. It takes the form of a directed acyclic graph, with terms as nodes and ‘is-a’ and ‘part-of’ relationships comprising most of the links.
# The annotation that links organism-specific genes to terms in gene ontology is separate from the vocabulary itself, and involves different types of evidence. These are recorded in Bioconductor annotation packages.
# We have immediate access to the GO vocabulary with the GO.db package.
BiocManager::install("GO.db")
library(GO.db)
GO.db # metadata

# The keys/columns/select functionality of AnnotationDbi is easy to use for mappings between ids, terms and definitions.
k5 = keys(GO.db)[1:5]
cgo = columns(GO.db)
select(GO.db, keys=k5, columns=cgo[1:3])

# The graphical structure of the vocabulary is encoded in tables in a SQLite database. We can query this using the RSQLite interface.
install.packages("DBI")
library(DBI)
con = GO_dbconn()
dbListTables(con)

# The following query reveals some internal identifiers:
dbGetQuery(con, "select _id, go_id, term from go_term limit 5")

# We can trace the mitochondrion inheritance term to parent and grandparent terms:
dbGetQuery(con, "select * from go_bp_parents where _id=30")

dbGetQuery(con, "select _id, go_id, term from go_term where _id=26616")
dbGetQuery(con, "select * from go_bp_parents where _id=26616")
dbGetQuery(con, "select _id, go_id, term from go_term where _id=5932")

# It makes sense to regard “mitochondrion inheritance” as a conceptual refinement of processes “mitochondrion distribution”, and “organelle inheritance”, the two terms that are regarded as parents in this database scheme.
# The entire database schema can be viewed with GO_dbschema().

# KEGG: Kyoto Encyclopedia of Genes and Genomes
# The KEGG annotation system has been available in Bioconductor since the latter’s inception, but licensing of the database has changed. When we attach KEGG.db we see
library(KEGG.db)

# Therefore we focus on KEGGREST, which requires active internet connection. A very useful query resolution facility is based on Entrez identifiers. The Entrez ID for BRCA2 is 675. We’ll perform a general query.
BiocManager::install("KEGGREST")
library(KEGGREST)
brca2K = keggGet("hsa:675")
names(brca2K[[1]])

# The list of genes making up a pathway model can be obtained with another keggGet:
brpat = keggGet("path:hsa05212")
names(brpat[[1]])

brpat[[1]]$GENE[seq(1,132,2)] # entrez gene ids

# we can acquire a static image of the (human) pancreatic cancer pathway in which BRCA2 is implicated.
BiocManager::install("png")
BiocManager::install("grid")
library(png)
library(grid)
brpng = keggGet("hsa05212", "image")
grid.raster(brpng)


# The rols package interfaces to the EMBL-EBI Ontology Lookup Service.
BiocManager::install("rols")
library(rols)
oo = Ontologies()
oo
oo[[1]]

# To control the amount of network traffic involved in query retrieval, there are stages of search.
glis = OlsSearch("glioblastoma")
glis

res = olsSearch(glis)
dim(res)

resdf = as(res, "data.frame") # get content
resdf[1:4,1:4]

resdf[1,5]  # full description for one instance

# The GSEABase package has excellent infrastructure for managing gene sets and collections thereof. 
#We illustrate by importing a glioblastoma-related gene set from MSigDb.
BiocManager::install("GSEABase")
library(GSEABase)
glioG = getGmt(system.file("gmt/glioSets.gmt", package="ph525x"))
glioG
head(geneIds(glioG[[1]]))
# Somethings  not right here...

# The OrganismDb packages simplify access to annotation. Queries that succeed against TxDb, and org.[Nn].eg.db can be directed at the OrganismDb object.
BiocManager::install("Homo.sapiens")
library(Homo.sapiens)
class(Homo.sapiens)
Homo.sapiens
tx = transcripts(Homo.sapiens)
keytypes(Homo.sapiens)
columns(Homo.sapiens)


# By sorting information on the GPL information page at NCBI GEO, we see that the most commonly use oligonucleotide array platform (with 4760 series in the archive) is the Affy Human Genome U133 plus 2.0 array (GPL 570). 
# We can work with the annotation for this array with the hgu133plus2.db package.
BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
hgu133plus2.db

# The basic purpose of this resource (and of all the instances of the ChipDb class) is to map between probeset identifiers and higher-level genomic annotation.
# The detailed information on probes, the constituents of probe sets, is given in packages with names ending in ‘probe’.
BiocManager::install("hgu133plus2probe")
library(hgu133plus2probe)
head(hgu133plus2probe)
dim(hgu133plus2probe)

# Mapping out a probe set identifier to gene-level information can reveal some interesting ambiguities:
select(hgu133plus2.db, keytype="PROBEID", 
       columns=c("SYMBOL", "GENENAME", "PATH", "MAP"), keys="1007_s_at")

# Apparently this probe set includes sequence that may be used to quantify abundance of both an mRNA and a micro RNA. 
#As a sanity check we see that the distinct symbols map to the same cytoband.







