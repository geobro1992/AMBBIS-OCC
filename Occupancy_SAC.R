library(AICcmodavg)
library(cvAUC)
library(ISLR)
library(boot)
library(ncf)
library(spdep)
library(codep)
library(dismo)
library(raster)
library(gbm)

#------------------
# SAC analyses
#------------------

dd <- read.csv("5894DAT.csv", sep = ",", header = T, row.names = NULL)

#------------
# TSA
#------------
# get polynomials 
xy = as.matrix(cbind(dd[,9],dd[,8]))

poly3 = poly(xy, degree = 3, raw = T)
colnames(poly3) = c("X","X2","X3","Y1","XY","X2Y","Y2","XY2","Y3")

hab = cbind(dd, poly3)
# TSA
poly1fit = glm(Y ~ X + Y1, data = hab, family = binomial)
AICc(poly1fit)

# leave one out cross validation
ci.cvAUC(poly1fit$fitted.values, poly1fit$y, label.ordering = NULL, folds = NULL, confidence = 0.95)
# using residulas did not change predictions

#------------
# Correlogram
#------------
par(mfrow = c(2, 4))
count = 0
dis = list()
inte = vector()
for (i in seq(from = 0.3, to = 2, length.out = 8)) {
  count = count + 1
  AB.corr = correlog(dd$LON, dd$LAT, dd$Y, w = NULL, i, resamp = 1000, latlon = TRUE, na.rm = TRUE, quiet = TRUE)
  signif = which(AB.corr[[5]] < 0.05)
  dis[[count]] = AB.corr[2]$mean.of.class[signif]
  inte[count] = AB.corr$x.intercept
  
  plot(AB.corr)
  abline(h = 0, col = "red")
}

seq(from = 0.3, to = 2, by = 0.05)
AB.corr = correlog(dd$LON, dd$LAT, dd$Y, w = NULL, 0.75, resamp = 1000, latlon = TRUE, na.rm = TRUE, quiet = TRUE)
AB.corr[2]$mean.of.class
#----------------------------
# VARIANCE PARTITIONING
#----------------------------

# extract latlons
xy = as.matrix(cbind(dd[, 9], dd[, 8]))

# create neighbourhood matrix <4km based on GPS locations
nbnear4 = dnearneigh(xy, 0, 4, row.names = dd$ID, longlat = TRUE) 

par(mfrow = c(1,1))
plot(nbnear4, xy, col = "red", pch = 20, cex = 1)
title(main = "neighbours if 0<d<4km")

###We now use the multiscale co-dependence analysis "codep" package to perform the more general Moran Eigenvector Maps (MEM) analysis
mem = eigenmap(xy, weighting = Wf.Drayf1, boundaries = c(0, 4))#performing MEM with the 1-dij/max(dij) weighting function and threshold = 4km
names(mem)#Of immediate interest are lambda and U
par(mfrow = c(1, 1), mar = c(4, 4, 4, 4))
barplot(mem$lambda, ylab = "Lambda")#Note differences in the number of positive eigenvectors compared to the PCNM
plot(mem)#A scree plot of eigenvalues and spatial pattern (map) plots of spatial eigenvectors (i.e., the MEMs)
##Now take your large positive eigenvectors and work with them as you :please
myEigenvecs = mem$U[, 1:2]


library(vegan)
vapa = varpart(Y = dd$Y, myEigenvecs, ~ Hydro + Area + MaxDepth + Habitat + Habitat_prop, data = dd)
par(mfrow = c(1, 1))
plot(vapa, cutoff = -Inf, cex = 1.5, bg = 2:5)

#-----------------------
# SAC test & Correlogram
#-----------------------
##Compute pairwise Euclidean distances with the function nbdist and calculate spatial weight matrix 
#The Euclidean distances are reported in units of decimal degrees 
dist_ecase = nbdists(nbnear4, xy)
#Now, define weight as a function of distance (e.g. 1-dij/max(dij))
#You could use other definitions here but this should work for starters
fdist = lapply(dist_ecase, function(x) 1 - x / max(dist(xy)))
#Now define your weight matrix, with 0 for non-neighbors and inverse-weighted distance for neighbors 
listw_ecase = nb2listw(nbnear4, glist = fdist, style = "W")
#Creating a neighborhood matrix from the nblistw file
mynbmatrix = listw2mat(listw_ecase)


#---------------------------------------------------------
# neighbourhood matrices and Moran's I for 1st eigenvector
#---------------------------------------------------------
n = length(myEigenvecs[, 1])
W = sum(mynbmatrix)
I_obs = (1 / W * t(myEigenvecs[, 1] - mean(myEigenvecs[, 1])) %*% mynbmatrix %*% (myEigenvecs[, 1] - mean(myEigenvecs[, 1])))/(1 / n * sum((myEigenvecs[, 1] - mean(myEigenvecs[, 1])) ^ 2))#Calculating observed value of Moran's I



###Randomization test for first eigenvector

nreps = 1000         # Note: 10,000 replications
I_ran = vector()

for (i in 1:nreps) {
  randScore <- sample(myEigenvecs[, 1], n, replace = FALSE)
  I_ran[i] = (1 / W * t(randScore - mean(randScore)) %*% mynbmatrix %*% (randScore - mean(randScore)))/(1 / n * sum((randScore - mean(randScore)) ^ 2))#Calculating observed value of Moran's I
}

##We can plot the distribution of Morans and test it for normality
hist(I_ran)
qqnorm(I_ran)
qqline(I_ran)
shapiro.test(I_ran)
###We know that the expected value of Moran's I is -1/(n-1) = -0.2 in this case, but we can't proceed with a Z-test due to lack of normality###
##Our other option is a permutation test based on the empirical CDF and the observed value of Moran's I
#Let's plot the hostogram and shade the tail of values > I_obs
par(mfrow = c(1, 1))
histMorans = hist(I_ran, breaks = 50, plot = FALSE)
plot(histMorans, col = ifelse(histMorans$breaks > I_obs[, 1], "red", "grey50"), main = "Randomization Test for 1st Eigenvector", xlab = "Moran's I")
text(0.2, 40, "Moran's I = 0.98", cex = 1.5)
#The Test
perc.rank = ecdf(I_ran)#This establishes percentile values for the empirical CDF
perc.rank(I_obs)
p = 1 - (perc.rank(I_obs))
p
###So with a p-value of 0, there is spatial autocorrelation in our data

#---------------------------------
# repeat with 2nd eigenvector
#----------------------------------
##Moran's I and Randomization test for 2nd positive eigenvector
I_obs = (1 / W * t(myEigenvecs[, 2] - mean(myEigenvecs[, 2])) %*% mynbmatrix %*% (myEigenvecs[, 2] - mean(myEigenvecs[, 2]))) / (1 / n * sum((myEigenvecs[, 2] - mean(myEigenvecs[, 2])) ^ 2))#Calculating observed value of Moran's I

nreps = 1000         # Note: 10,000 replications
I_ran = vector()

for (i in 1:nreps) {
  randScore <- sample(myEigenvecs[, 2], n, replace = FALSE)
  I_ran[i] = (1 / W * t(randScore - mean(randScore)) %*% mynbmatrix %*% (randScore - mean(randScore))) / (1 / n * sum((randScore - mean(randScore)) ^ 2))#Calculating observed value of Moran's I
}

##We can plot the distribution of Morans and test it for normality
hist(I_ran)
qqnorm(I_ran)
qqline(I_ran)
shapiro.test(I_ran)
###We know that the expected value of Moran's I is -1/(n-1) = -0.2 in this case, but we can't proceed with a Z-test due to lack of normality###
##Our other option is a permutation test based on the empirical CDF and the observed value of Moran's I
#Let's plot the hostogram and shade the tail of values > I_obs
histMorans = hist(I_ran, breaks = 50, plot = FALSE)
plot(histMorans, col = ifelse(histMorans$breaks > I_obs[, 1], "red", "grey50"), main = "Randomization Test for 2th Eigenvector", xlab = "Moran's I")
text(0.2, 40, "Moran's I = 0.37", cex = 1.5)
#The Test
perc.rank = ecdf(I_ran)#This establishes percentile values for the empirical CDF
perc.rank(I_obs)
p = 1 - (perc.rank(I_obs))
p
###So with a p-value of 1, there is no spatial autocorrelation in our 2nd eigenvector



#-------------------------
# BRT without eigenvectors
#-------------------------
sub.eig = dd[, -c(8, 9)]

sub.eig.boost = gbm.step(data = sub.eig, gbm.x = 3:7, gbm.y = 2,  family = "bernoulli", tree.complexity = 3, max.trees = 100000, learning.rate = 0.0001, bag.fraction = 0.85)
#Model has some good outcomes for both training and CV. You can tweak some of the tuning parameters to see how they affect the model performance
#Variable importance & plot
par(mfrow = c(1, 1), mar = c(5, 8, 4, 2))
summary(sub.eig.boost, las = 1)

#Partial dependence plots
gbm.plot(sub.eig.boost, write.title = FALSE, plot.layout = c(2,2))

#-----------------------
# BRT with eigenvectors
#-----------------------
sub.eig = cbind(dd[,-c(8,9)], myEigenvecs)


sub.eig.boost = gbm.step(data = sub.eig,gbm.x = 3:9, gbm.y = 2,  family = "bernoulli", tree.complexity = 2, learning.rate = 0.001, bag.fraction = 0.85)
#Model has some good outcomes for both training and CV. You can tweak some of the tuning parameters to see how they affect the model performance
#Variable importance & plot
par(mfrow = c(1, 1), mar = c(5, 8, 4, 2))
summary(sub.eig.boost, las = 1)

#Partial dependence plots
par(mfrow = c(2,4), mar = c(5, 2, 2, 2))
gbm.plot(sub.eig.boost, write.title = FALSE, plot.layout = c(2, 4))




