library(vegan)
library(tidyverse)
setwd("~/bio-490/independentStudy/MiSeqSOP/MiSeq_SOP/rarefy")

#set the random generator seed
set.seed(19760620)

days_wanted <- c(0:9, 141:150)

#read in the shared data filed
# rows are different samples and columns are the different OTU's
# as we read it in we make it have a column for the smaple, a column for the OTu and a column for the number of times that OTU was seen 
shared <- read_tsv("mice.shared") %>%
  select(Group, starts_with("Otu")) %>%
  mutate(day = str_replace(Group, ".*D", "")) %>%
  filter(day %in% days_wanted) %>%
  select(-day) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total = sum(value)) %>%
  filter(total > 1800) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total)

# takes all the data from shared and pool all the counts while keeping each OTU separate. 
# It kind of creates a "meta community" that we can sample from to regenerate the distribution of the samples 
# basically a bunch of random samplings from the same larger community but while keeping the same orginal number of sequences that were in each sample
rand <- shared %>%
  uncount(value) %>%
  mutate(name = sample(name)) %>%
  count(Group, name, name="value")

#Goal is to calculate a variety of alpha diversity metrics for each sample including:
#   The richness (totaal number of  OTU's in each sample)
#   The shannon weaver index
# relative abundance times the log of the relative abundance so it gives a more balanced description of the OTU's in the sample
#   The Simpson index (squares the relative abundance)
# so as your sample more, you're emphasizing the rare OTU's
#   The inverse simpson index
#   The total number of individuals in each sample

# To do this write 3 functions to calculate richness, Shannon, and Simpson index

# these functions take in one parameter, x, which is the counts ( THE VALUE COLUMN FROM THE SHARED DF)

# This= function returns the number of time x was greater than 0 
#commented out part is more verbose, but still works
richness <- function(x){
  
  # r <- sum(x > 0)
  # return(r)
  
  sum(x>0)
}

shannon <- function(x){
  
  #rabund gives the relative abundance 
  rabund <- x[x>0]/sum(x) 
  -sum(rabund * log(rabund))
  
}

simpson <- function(x){
  
  # start by calculating the total number of individuals in x
  n <- sum(x)
  
  # find the sum of the calculated values
  # sum(x * (x-1) / (n * (n-1)))
  1 - sum((x/n)^2)
}

#PLOT
# Y axis = diversity measure value 
# x axis total number of sequences for each  sample
# each point represents a different sample
# remeber rand is the random samplig of the pooled data
rand %>%
  group_by(Group) %>%
  summarize(sobs = richness(value),
            shannon = shannon(value),
            simpson = simpson(value),
            invsimpson = 1/simpson,
            #number of individuals in each sample
            n = sum(value)) %>%
  pivot_longer(cols=c(sobs, shannon, invsimpson, simpson),
               names_to="metric") %>%
  ggplot(aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~metric, nrow=4, scales="free_y")

# Does the same thing, but using vegan's fucntions instead of the user defined functions for each diversity measure
# The only difference between the two is in the calculations for the Simpson calculation
rand %>%
  group_by(Group) %>%
  summarize(sobs = specnumber(value),
            shannon = diversity(value, index="shannon"),
            simpson = diversity(value, index="simpson"),
            invsimpson = 1/simpson,
            n = sum(value)) %>%
  pivot_longer(cols=c(sobs, shannon, invsimpson, simpson),
               names_to="metric") %>%
  ggplot(aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~metric, nrow=4, scales="free_y")

