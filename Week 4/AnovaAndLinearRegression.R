# most of the following script can be pasted directly into Rstudio
# but the script is intentionally incomplete - follow directions and think about what each line means
# comments are super important to help you remember what you did and why you did it when coding
# comment and annotate often - it is a gift to your future self who will not remember your present thought process
# your script should not only work in Rstudio, it should also have clear comment annotations

# first, if you haven't already, set the working directory - CHANGE THE FOLLOWING LINE to navigate to your own working directory
setwd("/Users/natalieburkhard/bio-490")

# you will need to set your working directory each time you start Rstudio
# so all of your scripts should probably start with that line
# use getwd() to check your working directory if you like

# read in data from a file in that directory - make sure to download today's data from Brightspace first
# there are other ways but this one is very intuitive and visual
# next we will try a new test - one comparing more than two samples

sward2 = read.csv(file = file.choose())

names(sward2)
# this just shows you the names of your variables, so that you can correctly enter them into the analysis script below
# view sward2 and note the scientific recording format

# visualize your data by REPLACING the variable names with their actual names and the data object with its actual name
boxplot(Height ~ Site, data=sward2)
# write down your summary of the graph

# try to color it up
boxplot(Height ~ Site, data=sward2, col=c("orange", "navy", "turquoise"))

# perform a one-way ANOVA and see its results
sward2.aov = aov(Height ~ Site, data=sward2)
summary(sward2.aov)

# if ANOVA result is significant, a post-hoc test is in order
# we will use Tukey's test
# in your post-lab, summarize the results of the test
# focus on the significant differences - between which samples?
TukeyHSD(sward2.aov, conf.level = 0.95)


# now for a different test, a more complex one
# two-way ANOVA has more than one experimental variable
# in the following example it's amount of water and plant species
# the response variable is still plant growth
# sward3 = read.csv(file = file.choose())
# instead of the interactive file.choose above, let's try this - MODIFY the line below

sward3 = read.csv("/Users/natalieburkhard/bio-490/sward3.csv")
names(sward3)

# the advantage is that you can run it as part of a code without manual input
# view sward3 and note how the format differs from sward2 - it is the same basic format but now has two factors

# boxplot is a little different with two variables
# growth is still the response variable
boxplot(growth~species*Water, data=sward3)

# your next challenge: make the boxes three colors - one color for each watering treatment
# paste the boxplot in your post-lab
boxplot(growth ~ species*Water, data=sward3, col=c("orange", "navy", "turquoise"))

# now the actual ANOVA - set it up similarly to setting up the box plot
sward3.aov = aov(growth~species*Water, data=sward3)
summary(sward3.aov)
TukeyHSD(sward3.aov, conf.level = 0.95)
# we can reorganize those same results to just focus on interactions of watering and species
TukeyHSD(sward3.aov, ordered = TRUE, which = "species:Water", conf.level = 0.95)
# interpret these results briefly in your post-lab

#######

# your next challenge is a correlation analysis - a simple linear regression
#before you start though - you will need to get your data from the "invert" tab in the shared spreadsheet into the right format for R
#examine the csv files (sward 2 and 3, or any previous data)
# save the invert data in csv format and make sure the data are in the first rows and columns, not in the middle of the spreadsheet
# name your file mayfly.csv or something similar; load the data into R and store in object "mayfly"
mayfly = read.csv("/Users/natalieburkhard/bio-490/mayfly.csv")
# make a histogram of each variable - if normally distributed, Pearson correlation is appropriate. Spearman correlation is for non-normal data.
hist(mayfly$abund)
hist(mayfly$flow)
# testing for correlation is simple; the method is either "pearson" or "spearman"
cor.test(mayfly$abund, mayfly$flow, method = "pearson")
cor.test(mayfly$abund, mayfly$flow, method = "spearman")

# now to fit a linear model to the data you simply use (MODIFIED)

lm(mayfly$abund ~ mayfly$flow)

# in the parentheses put (response_variable ~ predictor_variable) - you'll need to figure out how to tell R where to find them

# and similarly you can plot them by MODIFYING

plot(mayfly$abund~mayfly$flow)

# add a trendline, why don't you (MODIFY)

abline()

# want the formula on the graph? no problem (MODIFY), you'll need to click somewhere on the graph afterwards to place the text

text(locator(1), ”your_formula")

# that's it! (of course, straight lines don't fit every data set, but we'll talk about that some other time...)