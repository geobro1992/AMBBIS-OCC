rm(list = ls())
library(R2WinBUGS)

#----------------
# Sampling data
#----------------
df = read.csv("OccDat.csv", sep = ",")
# extract year from dates and merge to data
Ys = format(as.Date(df$Date, format = "%d/%m/%Y"), "%Y")
df = cbind(df, Ys)

# omit NAs and extract most recent years
df = na.omit(df)
df = df[which(df$Year > 2007), ]

# function to create column for intra-annual sampling occasion
j.list = vector()
j.list[1] = 1
count = 1

for (i in 2:length(df[, 1])) {
  if(df[i,1] == df[i-1,1] && df[i,2] == df[i-1,2]){
    count = count + 1
    j.list[i] = count
  }
  else{
    count = 1
    j.list[i] = count
  }
}

df = cbind(df, j.list)


# recast to have each sampling occasion as its own column
library(reshape2)
dd = dcast(df, Year + Site ~  j.list, value.var = "Caps")

# add rows for sites that were never sampled to balance dataframe
library(data.table)
DT = as.data.table(dd)
setkey(DT, Site, Year)
dbase = DT[CJ(unique(Site), seq(min(Year), max(Year)))]

#--------------------
# Distance Matrix
#--------------------

# read in LatLon data and create a matrix with distances between ponds
LL = read.csv("OccLatLon.csv", header = T, sep = ",")
length(df[!duplicated(df$Site), 1])
LL.dat = LL[LL$ID %in% df[!duplicated(df$Site), 1], ]
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
dbase = dbase[dbase$Site %in% LL.dat$Site, ]


#----------------
# Occupancy Data
#----------------

# create 3D matrix for occupancy data
Ys = as.factor(dbase$Year)
Ss = as.factor(dbase$Site)

y <- array(NA, dim = c(length(unique(Ss)), 11, length(unique(Ys))))	# sites, reps, years


# reformat occupancy data into the 3D matrix
Y1 = min(dbase$Year)
YN = max(dbase$Year)

for (i in 1:9) {
  sel.rows <- (dbase$Year - Y1) + 1 == i
  y[, , i] <- as.matrix(dbase[sel.rows, 3:13])
}

y[y > 1] = 1

# Look at the number of sites with detections for each year
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
apply(tmp, 2, sum, na.rm = TRUE)



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
ni <- 2000000
nt <- 100
nb <- 10000
nc <- 3

# Call WinBUGS from R
outAB <- bugs(win.data, inits, params, "SMDynOcc.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = "C:/Program Files (x86)/WinBUGS14", working.directory = getwd())
save(outAB, file = "Occupancy_Run.RData")
