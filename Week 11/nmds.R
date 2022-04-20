setwd("~/bio-490/independentStudy/MiSeqSOP/")
library(tidyverse)
library(vegan)
source("readMatrix.R")

### NOTES ###
# https://www.youtube.com/watch?v=h7OrVmT7Ja8
# ordination = visualization tool to show clusters in the data. 
# NMDS (non-metric dimensional scaling) is a form of ordination. alternative to PCOA
# In NMDS you get to tell it how many dimensions youy want and the algorithm will figure out the best way to plot the data in a way that represents the variation of the data
# repeats the optimization many times to get the best placement
# https://sites.ualberta.ca/~lkgray/uploads/7/3/6/2/7362679/slides-nmds.pdf
# The goal of NMDS is to represent the original position of data in multidimensional space as accurately as possible using a reduced number of
# dimensions that can be easily plotted and visualized (like PCA). NMDS relies on rank orders (distances) for ordination (i.e non-metric)
# NMDS is an iterative procedure which takes place over several steps:
# 1. Define the original data point positions in multidimensional space
# 2. Specify the number of reduced dimensions you want (typically 2)
# 3. Construct an initial configuration of the data in 2-dimensions
# 4. Compare distances in this initial 2D configuration against the calculated distances
# 5. Determine the stress on data points
# 6. Correct the position of the points in 2D to optimize the stress for all points
# When we compress our 3D image to 2D we cannot accurately plot the true distances
# The difference between the data point position in 2D (or # of dimensions we consider with NMDS) and the distance calculations (based on multivariate) is the STRESS we are trying to optimize
# Think of optimizing stress as: “Pulling on all points a little bit so no single point is completely wrong, all points are a little off compared to distances” Ideally we want as little stress as possible 
# Scores – these are the data point outputs that have be pulled to optimize the stress from multi dimensions in 2D (or the # of dimensions considered) These are the values we plot to look at which data points group together

#reads in a lower triangule distance matrix
dist_matrix <- read_matrix("MiSeq_SOP/final.opti_mcc.jclass.0.03.lt.ave.dist")

dist_tbl <- as_tibble(dist_matrix, rownames="samples")

sample_lookup <- dist_tbl %>% 
  select(samples) %>%
  mutate(delimited = str_replace(samples,
                                 "^(([FM])\\d+)D(\\d+)$",
                                 "\\2-\\1-\\3")) %>%
  separate(col=delimited,
           into=c("sex", "animal", "day"), sep="-",
           convert=TRUE)

# looks at samples from the first 10 days post weaning and the from 141-150 days post weaning
days_wanted <- c(0:9, 141:150)

#outputs a distance matrix (dist_matrix) output and a sample_lookup dataframe which contains information about each of the samples (sex, animal identifier, and number of days post weaning)
dist_matrix <- dist_tbl %>%
  pivot_longer(cols=-samples, names_to="b", values_to="distances") %>%
  inner_join(., sample_lookup, by="samples") %>%
  inner_join(., sample_lookup, by=c("b" = "samples")) %>%
  filter(day.x %in% days_wanted & day.y %in% days_wanted) %>%
  select(samples, b, distances) %>%
  pivot_wider(names_from="b", values_from="distances") %>%
  select(-samples) %>%
  as.dist()

#set seed of the random number generator used my metaNMDS to get convergence
set.seed(10)
#runs the ordination
# nmds turns out to be a list of a bunch of objects. We are interested in the object called points which is a matrix which contains the poistions of MDS1 and MDS2
nmds <- metaMDS(dist_matrix)
# no convergence using file MiSeq_SOP/final.opti_mcc.jclass.0.03.lt.std.dist

# stress value = the amount of distortion that happens when you take highly dimensional and squish it into two dimensions
stress <- nmds$stress

#Access values (points) of the list (nmds) gi es the matrix which contains all the points data
nmds$points
#function from vegan that gives you the position of all the points in the two different axes
scores(nmds)
#generates the  plot of the ordination 
plot(nmds)
# this plots the ordination but with ggplot so you can get fancy
# %>% means pipe into. So all of that beneath scores(nmds) is being piped into scores(nmds)
scores(nmds) %>%
#tibbles don't have rownames so you kind of have to put it in a new column (c alled samples)
  as_tibble(rownames = "samples") %>%
  #have to get join the current data (which does not have the days post weaning) and the data from sample_lookup dataframe to get the days post weaning
  inner_join(., sample_lookup, by="samples") %>%
  #defining early or late
  mutate(period = if_else(day < 10, 'early', 'late')) %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=period)) +
  geom_point()

