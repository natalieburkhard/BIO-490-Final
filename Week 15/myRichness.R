library(vegan)
# BiocManager::install("phyloseq")
#package used to estimate richness
install.packages("breakaway")
library(breakaway)
library(tidyverse)
setwd("~/bio-490/independentStudy/MiSeqSOP/MiSeq_SOP")

set.seed(19760620)

days_wanted <- c(0:9, 141:150)

# Load in data
shared <- read_tsv("final.opti_mcc.shared") %>%
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

# create the rand data frame which has the columns group, name, and value
rand <- shared %>%
  uncount(value) %>%
  mutate(name = sample(name)) %>%
  count(Group, name, name="value")

#######################################################################
# FROM THE VIDEO (SHOWS THE PROCESS BEHIND THE FINAL SCRIPT)
# breakaway needs the data in a different format
# It wants a column for the value and the number of times each of those values occurs
ba <- rand %>%
  filter(Group == "F3D0") %>%
  count(Group, value) %>%
  select(-Group) %>%
  breakaway()

# At this point we have a data frame made from the frequency counts of each group 
# We want to make a data frame for each group
# We have a column called data which contains all of the data except for the group column
# The final result is a two column tibble wiht the first column being the group id 
# and the second column (data) which is a tibble in each row indicating a frewuncy and a frequency of those frequencies
# Now we can feed this into mutate using map functions
# From this code we get a data frame with a column for the group, the data, the number of sequences and the observed richness, 
#and a column called ba with a list of alpha diversity objects (not helpful)
rand %>% 
  count(Group, value) %>%
  nest(data = -Group) %>%
  mutate(n_seqs = map_dbl(data, ~sum(.x$value * .x$n)),
          sobs = map_dbl(data, ~sum(.x$n)),
          ba = map(data, ~breakaway(.x)))


##################################################################################################################################

# function that takes in a data frame (x) and outputs a tibble 
get_breakaway <- function(x){
  
  ba <- breakaway(x)
  tibble(est= ba$estimate,
         lci=ba$interval[1], uci=ba$interval[2],
         model=ba$model)
  
}

# run this again and we see that ba column is a list of tibbles instead of alpha diversity objects
b_analysis2 <- rand %>% 
  count(Group, value) %>%
  nest(data = -Group) %>%
  mutate(n_seqs = map_dbl(data, ~sum(.x$value * .x$n)),
        sobs = map_dbl(data, ~sum(.x$n)),
        ba = map(data, ~get_breakaway(.x))) %>%
  # and then remove the data column and unnest the ba column
  select(-data) %>%
  unnest(ba)

# plot this ^
# want to see if the richness estimate increases as the number of sequences increases
b_analysis2 %>% 
  ggplot(aes(x=n_seqs, y=est, color=model)) + 
  geom_point() + coord_cartesian(ylim=c(0, 2000))

# better graph
b_analysis2 %>% 
  ggplot(aes(x=n_seqs, y=est)) + 
  geom_smooth() + 
  geom_point() + coord_cartesian(ylim=c(0, 2000))

# lets look at rarefied data, the observed data, and another estimae (chao 1)
############################################################################
#FROM THE SCRIPT
# sobs = observed richness
get_chao <- function(x){
  
  sobs <- sum(x$n)
  sing <- x[x$value == 1, "n"] %>% pull(n)
  doub <- x[x$value == 2, "n"] %>% pull(n)
  
  sobs + sing^2 / (2*doub)
  
}

b_analysis <- rand %>%
  count(Group, value) %>%
  nest(data = -Group) %>%
  mutate(n_seqs = map_dbl(data, ~sum(.x$value * .x$n)),
         sobs = map_dbl(data, ~sum(.x$n)),
         #get rarefied data. The input is just the frequencies (where as with breakawy we were using frequencies of frequencies (so we need to kinda unpack that ))
         # I changed min n_seqs from 1828 to 2403
         rare = map_dbl(data, ~rarefy(rep(.x$value, .x$n), sample=2403)),
         # get the chao data
         chao = map_dbl(data, ~get_chao(.x)),
         ba = map(data, ~get_breakaway(.x))) %>%
  select(-data) %>%
  unnest(ba)

# make a scatter plot colored by the three different methods
# rarefy is flat
# sobs is the observed richness (increases with the sampling effort)
# est is the estimated richness (also increases with the number of sequences)
# There is definitely an influence of sampling effort on richness estimate
b_analysis %>%
  select(Group, n_seqs, sobs, chao, rare, est) %>%
  pivot_longer(-c(Group, n_seqs)) %>%
  ggplot(aes(x=n_seqs, y=value, color=name)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  coord_cartesian(ylim=c(0, 2000)) +
  theme_classic() +
  ggtitle("My title")

