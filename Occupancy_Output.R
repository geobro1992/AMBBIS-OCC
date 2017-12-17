#------------------------------
# Model output
#------------------------------
library(R2WinBUGS)
library(reshape2)
library(data.table)

load(file = "Occupancy_Run2.RData")

ds = seq(1, 2000, length.out = 50)
A = median(outAB$sims.list$Alpha)
A.low = min(quantile(outAB$sims.list$Alpha, probs = c(0.05, 0.95)))
A.high = max(quantile(outAB$sims.list$Alpha, probs = c(0.05, 0.95)))
log(2)/A
log(2)/A.high
log(2)/A.low

con = exp(-A*ds)
con.low = exp(-A.low*ds)
con.high = exp(-A.high*ds)
plot(ds, con, xlab = "Distance (m)", ylab = "Connectivity", type = "l", lwd = 2)
lines(ds, con.low, lty = "dashed")
lines(ds, con.high, lty = "dashed")

colon = 1 - exp(-con)
col.high = 1 - exp(-con.high)
col.low = 1 - exp(-con.low)
plot(ds,colon, xlab = "Distance (m)", ylab = "Colonization Rate", type = "l", lwd = 2)
lines(ds, col.low, lty = "dashed")
lines(ds, col.high, lty = "dashed")


pdf("SMDynOcc_Posts.pdf", 50, 50)
plot(outAB, display.parallel = FALSE)
dev.off()


Pstar <- array(NA, dim = c(outAB$n.sims, 10))
x <- cbind(rep(1, 2700), rep(2, 2700), rep(3, 2700), rep(4, 2700), rep(5, 2700), rep(6, 2700), rep(7, 2700), rep(8, 2700), rep(9, 2700), rep(10, 2700)) 

for (i in 1:outAB$n.sims) {
  for (j in 1:10) {
    Pstar[i,j] <- 1 - (1 - mean(outAB$sims.list$P[i,]))^j
  } #j
} #i

boxplot(Pstar ~ x, col = "gray", las = 1, ylab = "Detection Probability", xlab = "Number of surveys", outline = FALSE)
abline(h = 0.95, lty = 2, lwd = 2)


plot(con, colon, xlab = "Connectivity", ylab = "Colonization Rate", type = "l", lwd = 2)
lines(con.high, col.low, lty = "dashed")
lines(con.low, col.high, lty = "dashed")
