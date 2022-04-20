## from https://riffomonas.org/minimalR/07_line_plots.html

setwd = "/Users/natalieburkhard/bio-490/independentStudy/MiSeqSOP"
rarefy <- read_tsv(file="~/bio-490/independentStudy/MiSeqSOP/MiSeq_SOP/stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.groups.rarefaction")
rarefy
# In rarefy it has a column called numsampled. numsampled has the number of sequences that have been sampled
#A The columns are displayed in sets of threes. For example, 0.03-2003650 lci-2003650 hci-2003650. 
#The first column in the triplet is the average number of OTUs observed for that sample (e.g. 2003650) at the specified number of sequences sampled by the value in the numsampled column. 
#The second and third columns in the triplet represent the lower (lci) and higher (hci) confidence interval. 

rarefyTidy <- rarefy %>% select(contains("lci-"), contains("hci-"))
rarefyTidy <- rarefy %>% pivot_longer(cols=c(numsampled), names_to='sample', values_to='sobs')


ggplot(rarefyTidy, aes(x="numsampled", y=sobs, group=sample)) +
  geom_line()
ggplot(rarefyTidy, aes(x=sample, y=sobs, group=sample)) +
  geom_line()

ggplot(rarefy, aes(x="0.03-", y="numsampled", group=sample)) +
  geom_line()
