rm(list = ls())
library(R2WinBUGS)
library(reshape2)
library(data.table)

#----------------
# Sampling data
#----------------

# read in data
df = read.csv("rawdat.csv", sep = ",")

ids = c(2,3,4,5,6,7,11,12,13,14,15,16,18,19,21,30,33,34,36,41,50,51,53,103,107,111,126,201,202)


# extract year from dates and merge to data
Ys = format(as.Date(df[,2], format = "%m/%d/%Y") - 180, "%Y")
df = cbind(df, Year = as.numeric(Ys))
df[is.na(df[,3]),3] <- 0

# omit NAs and extract years of interest
df = na.omit(df)
df = df[which(df$Year > 2008), ]

df = df[order(df[,1]),]
df = df[df$SITE %in% ids,]
df = df[with(df, order(SITE, Year)),]

# function to calculate within season sampling occasion and store in a separate column
j.list = vector()
j.list[1] = 1
count = 1

for (i in 2:length(df[, 1])) {
  if (df[i, 1] == df[i - 1, 1] && df[i, 4] == df[i - 1, 4]) {
    count = count + 1
    j.list[i] = count
  }
  else{
    count = 1
    j.list[i] = count
  }
}

df = cbind(df, j.list)
head(df)


# recast to have each sampling occasion as its own column
dd = dcast(df, Year + SITE ~  j.list, value.var = "AMBBIS")


# add rows for sites that were never sampled to balance dataframe
DT = as.data.table(dd)
setkey(DT, SITE, Year)
dbase = DT[CJ(unique(SITE), seq(min(Year), max(Year)))]

#--------------------
# Distance Matrix
#--------------------

# read in LatLon data and create a matrix with distances between ponds
LL = read.csv("OccLatLon.csv", header = T, sep = ",")
length(df[!duplicated(df$SITE), 1])
LL.dat = LL[LL$ID %in% df[!duplicated(df$SITE), 1], ]
LL.dat[, 1] = as.factor(LL.dat[, 1])
colnames(LL.dat) = c("Site", "lat", "lon")

ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}

dmat = GeoDistanceInMetresMatrix(LL.dat)
#-------------------------------------------------------------

# extract sites for which there is GPS data 
dbase = dbase[dbase$SITE %in% LL.dat$Site, ]

#----------------
# Occupancy Data
#----------------

# create 3D matrix for occupancy data
Ys = as.factor(dbase$Year)
Ss = as.factor(dbase$SITE)

y <- array(NA, dim = c(length(unique(Ss)), 15, length(unique(Ys))))	# sites, reps, years


# reformat occupancy data into the 3D matrix
Y1 = min(dbase$Year)
YN = max(dbase$Year)

for (i in 1:8) {
  sel.rows <- (dbase$Year - Y1) + 1 == i
  y[, , i] <- as.matrix(dbase[sel.rows, 3:17])
}

# Convert counts to binary presence/absence (0,1)
y[y > 1] = 1

# Look at the number of sites with detections for each year
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
apply(tmp, 2, sum, na.rm = TRUE)

y

#----------------------------
# Dynamic Occupancy Model Run
#----------------------------

# Bundle data
win.data <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3], dmat = dmat)

# Initial values
z = apply(y, c(1, 3), max, na.rm = TRUE)
z[z == "-Inf"] = 0 
inits <- function(){ list(z = z)}

# Parameters monitored
params <- c("PSI1", "muZ", "E", "g", "P", "n.occ", "Alpha")

# MCMC settings
ni <- 500000   # iterations
nt <- 100      # thinning
nb <- 10000    # burn in
nc <- 3        # chains

# Call WinBUGS from R
outAB <- bugs(win.data, inits, params, "SMDynOcc.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = "C:/Program Files (x86)/WinBUGS14", working.directory = getwd())
# save output
save(outAB, file = "Occupancy_Run2.RData")
