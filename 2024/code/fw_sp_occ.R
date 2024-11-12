rm(list=ls())

#############################################################################################
# Spatially Explicit Dynamic Occupancy Model For the Flatwoods Salamander on Eglin Air Force Base

# libraries
library(jagsUI)
library(wiqid)
library(tidyverse)
library(geosphere)

# load in detection data and site data
dat = read.csv("fw_surveys.csv")
covs = read.csv("fw_site_dat.csv")


###################
# data cleaning
###################
# convert NA captures to zero captures
dat[is.na(dat[,3]),3] <- 0


# extract sites with LAT LONS
ids = unique(covs$ID)
dat = dat[dat$SiteName %in% ids,]

# extract year from dates and merge to data
dat$Year = as.numeric(format(as.Date(dat[,2], format = "%m/%d/%Y") + 180, "%Y"))

# filter
dat = dat[which(dat$Year > 2002), ]  # filter years
dat = dat[which(dat$SurveyTypeID_FK == 1 | dat$SurveyTypeID_FK == 6 | 
                dat$SurveyTypeID_FK == 15), ]  # filter survey type

# select relevant columns
dat = dat %>% select(SiteName, Year, AMBBIS)

# convert captures into binary presence-absence
dat[which(dat$AMBBIS > 0), "AMBBIS"] = 1

# only use sites with at least 20 surveys
ids = as.numeric(names(which(table(dat$SiteName) > 30)))
dat = dat[dat$SiteName %in% ids,]

#------------------------------------------------------------------
# run if you want to exclude translocation ponds
dat[which(dat$SiteName == 2 & dat$Year > 2019), "AMBBIS"] = 0  #except 2004 & 2006 (translocations began 2020)
dat[which(dat$SiteName == 30), "AMBBIS"] = 0
dat[which(dat$SiteName == 31), "AMBBIS"] = 0
#--------------------------------------------------------------------


###################
# data formatting
###################

# calculate distances between sites
dmat = distm(covs[,2:3],fun = distHaversine)

# store first and last years
year1 = min(dat$Year)
yearn = max(dat$Year)

# create column for number of surveys within years
dat = dat %>%
  group_by(Year, SiteName) %>%
  mutate(survey = 1:n())

# take only the first 6 surveys in each year
#dat = dat %>% filter(survey < 7)

# maximum number of surveys at a site within a year
n.surveys = max(dat$survey)

# number of years
n.years = length(unique(dat$Year))

# number of sites 
n.sites = length(unique(dat$SiteName))

# spread columns by year
dat = spread(dat, key = Year, value = AMBBIS)


# spread by surveys within years
dat = pivot_wider(
  dat, 
  id_cols = 'SiteName', 
  names_from = 'survey', 
  values_from = starts_with("20"), 
  names_glue = '{.value}.{survey}'
)


# filter LAT LONs to sites that have surveys
ids = unique(dat$SiteName)
covs = covs[covs$ID %in% ids,]

# Extract the detection histories
DH <- as.matrix(dat[, 2:length(dat)])
head(DH)


# convert to sites x occasions x year array
Y <- array(DH, dim=c(n.sites, n.surveys, n.years))
head(Y)

# Aggregate detection histories across occasions
y <- apply(Y, c(1, 3), sum, na.rm=TRUE)  # sites by years, number of detections
n <- apply(!is.na(Y), c(1, 3), sum)      # ... and no. of surveys

# specify known z's 
z <- (y > 0)*1
z[z == 0] <- NA

z.naive = colSums(z, na.rm = T)

# indicator variable for connectivity metric (don't include self)
indica = matrix(data = 1, nrow = n.sites, ncol = n.sites)

for(i in 1:n.sites){
  indica[i, i] = 0
}



#############
# JAGS model
#############

sink("fw_sp_occ_notranslocation.txt")
cat("

model {

     # Priors

  PSI1 ~ dbeta(1, 1)
  
  alpha ~ dgamma(1, 1)
  
  phi ~ dbeta(1, 1)


  for (k in 1:nYears) {
      p[k] ~ dbeta(1, 1)
    }


  # Initial State and Site Connectivities
  for (i in 1:nSites) {

    z[i, 1] ~ dbern(PSI1)

    for (j in 1:nSites) {
      con[i, j] <- exp(-alpha * dmat[i, j]) * indica[i, j]
    }
  }

  # Ecological submodel: Define state conditional on parameters
  for (k in 2:nYears) {
    for (i in 1:nSites) {
      for (j in 1:nSites) {
        disp[i, j, k-1] <- z[j, k-1] * con[i, j]
      } #j

      s[i, k-1] <- sum(disp[i, , k-1])
      gamma[i, k-1] <- 1 - exp(-s[i, k-1])

      z[i, k] ~ dbern( z[i, k-1]*phi + (1 - z[i, k-1])*gamma[i, k-1] )
    
    } #i
  } #k

  # Observation model

  for (i in 1:nSites) {
      for (k in 1:nYears) {
      y[i, k] ~ dbin(p[k] * z[i, k], n[i, k])
      } #k
  } #i
  
  # Derived variable
  for(k in 1:nYears) {
    N[k] <- sum(z[,k]) # no. sites occupied for each year
  }

} # model

",fill = TRUE)
sink()


# bundle data for jags
jdata <- list(nSites = n.sites, nYears = n.years, y = y,
              n = n, dmat = dmat/1000, indica = indica, z = z)

# parameters to monitor
wanted <- c("PSI1", "phi", "p", "alpha", "N", "gamma")

inits <- function(){ list(p = runif(n.years, 0.9,1), 
                          phi = 0.9, 
                          alpha = 0.3, 
                          PSI1 = 0.1)}

# run specs
nc = 3
ni = 50000
na = 10000
nb = 10000
nt = 400

# Run the model
jagsOut2 <- jags(data = jdata, inits = inits, parameters.to.save = wanted, model.file = "fw_sp_occ_notranslocation.txt",
                n.chains=nc, n.iter=ni, n.adapt=na, n.burnin = nb, n.thin = nt, DIC=FALSE, parallel=TRUE)


traceplot(jagsOut2)
hist(jagsOut2$sims.list$PSI1)
hist(jagsOut2$sims.list$phi)
hist(jagsOut2$sims.list$alpha)
hist(jagsOut2$sims.list$p)
hist(jagsOut2$sims.list$N)

# save output
#save(jagsOut2, file = "fw_sp_occ_10Oct2024.RData")
save(jagsOut2, file = "fw_sp_occ_NO_TRANSLOCATIONS_10Oct2024.RData")

jagsOut2

z.naive
                                                 #2007
z.naive_trans <- as.numeric(c("1", "6", "3", "4", "0",  "4",  "1",  "9",  "4",  "2",  
                              "2", "11",  "9", "11", "13",  "7",  "9", "18",  "7",
                              "2", "12", "16")) #added to account for 2007

##########
# OUTPUT

pdf(file = "occ_trends2.pdf", width = 10, height = 12)

op <- par(mfrow=c(2,1), cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), 
          cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)

#one with and other without translocations
b <- barplot(z.naive_trans, las = 1, xlab = " ", ylab = " ", col = "grey", cex.lab = 1.7, 
             cex.main = 1.5, axes = FALSE, ylim = c(0, 27))

n.years <- 22 #added to account for 2007

axis(1, seq(from = 0.7, by = 1.2, length.out = n.years), year1:yearn, cex.axis = 1.3)
axis(2, cex.axis = 1.3, las = 1)
par(las = 0)
mtext("Year", side = 1, line = 2.5, cex = 1.5)
mtext("Naive Occupancy", side = 2, line = 3, cex = 1.5)

a <- load("fw_sp_occ_NO_TRANSLOCATIONS_10Oct2024.RData")
a=jagsOut2
#occupied sites   (N is the predicted occupancy)
n = c("N[1]", "N[2]", "N[3]", "N[4]", "N[5]", "N[6]", "N[7]", "N[8]", "N[9]", "N[10]", "N[11]", "N[12]", "N[13]",  "N[14]", "N[15]", "N[16]",  "N[17]", "N[18]", "N[19]", "N[20]", "N[21]")
length(n)
# Number of Occupied Sites (summary of the mean # of occupied sites per year)
toplot <- jagsOut2$summary[n, c(1,3,7)]

plotsegraph <- function(loc, lower, upper, wiskwidth, color = "grey", linewidth = 2) {
  
  w <- wiskwidth/2
  segments(x0 = loc, x1 = loc, y0 = lower, y1 = upper, col = color, 
           lwd = linewidth)
  segments(x0 = loc - w, x1 = loc + w, y0 = upper, y1 = upper, 
           col = color, lwd = linewidth)  # upper whiskers
  segments(x0 = loc - w, x1 = loc + w, y0 = lower, y1 = lower, 
           col = color, lwd = linewidth)  # lower whiskers
}

newyears <- c(year1:yearn) #added to account for 2007
newyears <- subset(newyears, newyears != 2007)  #added to omit 2007
length(newyears) #21, now matches

     #Changed to "newyears" to account for 2007
plot(newyears, toplot[,1], col = "grey", pch = 21, bg = "grey", cex = 1.5,
     xlim = c(year1, yearn), ylim = c(0, 25), ylab = "", xlab = "", axes = FALSE)
lines(newyears, toplot[,1], lwd = 2, type = "c", col = "grey")
plotsegraph(newyears, toplot[,2], toplot[,3], 0.2, color = "grey")
axis(1, newyears)
axis(2) 

par(las = 0)
mtext("Year", side = 1, line = 2.5, cex = 1.5)
mtext("Predicted Occupancy", side = 2, line = 3.7, cex = 1.5)

lines(year1:yearn, rep(toplot[1,1], n.years), lty = "dashed")

#########################################
load("fw_sp_occ_10Oct2024.RData")
toplot <- jagsOut2$summary[n, c(1,3,7)]
points(newyears+0.1, toplot[,1], col = "black", pch = 21, bg = "black", cex = 1.5)
lines(newyears+0.1, toplot[,1], lwd = 2, type = "c")
plotsegraph(newyears+0.1, toplot[,2], toplot[,3], 0.2, color = "black")

points(2008, 60, pch = 19, lwd = 2, cex = 1.7, col = "grey")
text(2008.5, 60, "Without Translocations", cex = 1.5, font = 1, adj = 0)
points(2008, 55, pch = 19, lwd = 2, cex = 1.7)
text(2008.5, 55, "With Translocations", cex = 1.5, font = 1, adj = 0)


dev.off()

#################################################################################
# plot predicted relationship between distance and colonization probability
# makes more sense to use data without translocations to estimate natural dispersal

load("fw_sp_occ_NO_TRANSLOCATIONS.RData")

toplot <- jagsOut2$summary[c("alpha"), c(1,3,7)]

a.mu = toplot[1]
a.L = toplot[2] 
a.U = toplot[3]

dists = seq(0,1, length.out = 100)

con.mu <- exp(-a.mu[1] * dists)
con.L <- exp(-a.L[1] * dists)
con.U <- exp(-a.U[1] * dists)

# colonization with distance
gam.mu = (1 - exp(-con.mu))
gam.L = (1 - exp(-con.L))
gam.U = (1 - exp(-con.U))

pred.dat = data.frame(gam.mu, gam.L, gam.U, dists = dists)

g1 = ggplot(pred.dat, aes(x = dists, y = gam.mu))+
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = gam.L, ymax = gam.U), alpha = 0.3) +
  xlab("Distance (km)") + ylab("Colonization Probability")

g1 = g1 + theme_Publication()

ggsave(filename = "colonization.pdf", g1, width = 4, height = 4, units = "in", device = cairo_pdf)



# detection

pdf(file = "detection_trendsT.pdf", width = 10, height = 6)

n = c("p[1]", "p[2]", "p[3]", "p[4]", "p[5]", "p[6]", "p[7]", "p[8]", "p[9]", "p[10]", "p[11]", "p[12]", "p[13]",  "p[14]", "p[15]", "p[16]",  "p[17]")

op <- par(mfrow=c(1,1), mar=c(4,4,1,1)+0.1, pch=4, lwd=2, las=1)

# Number of Occupied Sites
toplot <- jagsOut2$summary[n, c(1,3,7)]

plotsegraph <- function(loc, lower, upper, wiskwidth, color = "grey", linewidth = 2) {
  
  w <- wiskwidth/2
  segments(x0 = loc, x1 = loc, y0 = lower, y1 = upper, col = color, 
           lwd = linewidth)
  segments(x0 = loc - w, x1 = loc + w, y0 = upper, y1 = upper, 
           col = color, lwd = linewidth)  # upper whiskers
  segments(x0 = loc - w, x1 = loc + w, y0 = lower, y1 = lower, 
           col = color, lwd = linewidth)  # lower whiskers
}

op <- par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(newyears, toplot[,1], col = "black", pch = 21, bg = "black", cex = 1.5,
     xlim = c(year1, yearn), ylim = c(0, 1), ylab = "", xlab = "", axes = FALSE)
plotsegraph(newyears, toplot[,2], toplot[,3], 0.2, color = "black")
axis(1, year1:yearn)
axis(2) 

par(las = 0)
mtext("Year", side = 1, line = 2.5, cex = 1.5)
mtext("Detection Probability", side = 2, line = 3.7, cex = 1.5)

dev.off()
