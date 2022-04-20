## Principal Component Analysis (PCA) 101, using R
# Improving predictability and classification one dimension at a time! “Visualize” 30 dimensions using a 2D-plot!
## from https://towardsdatascience.com/principal-component-analysis-pca-101-using-r-361f4c53a9ff

#load the data and name all 32 variables. The ID, diagnosis and ten distinct (30) features
wdbc <- read.csv("wdbc.data", header = F)
features <- c("radius", "texture", "perimeter", "area", "smoothness", "compactness", "concavity", "concave_points", "symmetry", "fractal_dimension")
names(wdbc) <- c("id", "diagnosis", paste0(features,"_mean"), paste0(features,"_se"), paste0(features,"_worst"))

# Why PCA? capture 63.3% (Dim1 44.3% + Dim2 19%) of variance in the entire dataset by just using those two principal components, 
# pretty good when taking into consideration that the original data consisted of 30 features which would be impossible to plot in any meaningful way.
# Steps of how PCS works:
# Standardize the data (Center and scale).
# Calculate the Eigenvectors and Eigenvalues from the covariance matrix or correlation matrix 
# Sort the Eigenvalues in descending order and choose the K largest Eigenvectors (Where K is the desired number of dimensions of the new feature subspace k ≤ d).
# Construct the projection matrix W from the selected K Eigenvectors.
# Transform the original dataset X via W to obtain a K-dimensional feature subspace Y.


# the ‘prcomp’ function runs PCA on the data we supply it, in our case that’s ‘wdbc[c(3:32)]’ which is our data excluding the ID and diagnosis variables
# then we tell R to center and scale our data (thus standardizing the data). 
wdbc.pr <- prcomp(wdbc[c(3:32)], center = TRUE, scale = TRUE)

# Finally we call for a summary:
# Standard deviation: This is simply the eigenvalues in our case since the data has been centered and scaled (standardized)
# Proportion of Variance: This is the amount of variance the component accounts for in the data, ie. PC1 accounts for >44% of total variance in the data alone!
# Cumulative Proportion: This is simply the accumulated amount of explained variance, ie. if we used the first 10 components we would be able to account for >95% of total variance in the data.
summary(wdbc.pr)

# so how many components do we want? We obviously want to be able to explain as much variance as possible but to do that we would need all 30 components, 
# at the same time we want to reduce the number of dimensions so we definitely want less than 30!
# Since we standardized our data and we now have the corresponding eigenvalues of each PC we can actually use these to draw a boundary for us. 
# Since an eigenvalues <1 would mean that the component actually explains less than a single explanatory variable we would like to discard those. 
# If our data is well suited for PCA we should be able to discard these components while retaining at least 70–80% of cumulative variance. 
#Lets plot and see:
screeplot(wdbc.pr, type = "l", npcs = 15, main = "Screeplot of the first 10 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
# We notice is that the first 6 components has an Eigenvalue >1 and explains almost 90% of variance, this is great! 
#We can effectively reduce dimensionality from 30 to 6 while only “loosing” about 10% of variance!

cumpro <- cumsum(wdbc.pr$sdev^2 / sum(wdbc.pr$sdev^2))
plot(cumpro[0:15], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 6, col="blue", lty=5)
abline(h = 0.88759, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC6"),
       col=c("blue"), lty=5, cex=0.6)
# We also notice that we can actually explain more than 60% of variance with just the first two components. Let’s try plotting these:
plot(wdbc.pr$x[,1],wdbc.pr$x[,2], xlab="PC1 (44.3%)", ylab = "PC2 (19%)", main = "PC1 / PC2 - plot")

# There’s some clustering going on in the upper/middle-right. Lets also consider for a moment what the goal of this analysis actually is. 
# We want to explain difference between malignant and benign tumors. Let’s actually add the response variable (diagnosis) to the plot and see if we can make better sense of it:
install.packages("factoextra")
library("factoextra")
fviz_pca_ind(wdbc.pr, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = wdbc$diagnosis, 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Diagnosis") +
  ggtitle("2D PCA-plot from 30 feature dataset") +
  theme(plot.title = element_text(hjust = 0.5))
