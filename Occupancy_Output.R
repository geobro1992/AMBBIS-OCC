#------------------------------
# Model output
#------------------------------

outAB = load(file = "outAB.RData")


ds = seq(1, 3200, length.out = 50)
A = outAB$mean$Alpha
A.low = min(quantile(outAB$sims.list$Alpha, probs = c(0.05, 0.95)))
A.high = max(quantile(outAB$sims.list$Alpha, probs = c(0.05, 0.95)))

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
x <- cbind(rep(1, 13500), rep(2, 13500), rep(3, 13500), rep(4, 13500), rep(5, 13500), rep(6, 13500), rep(7, 13500), rep(8, 13500), rep(9, 13500), rep(10, 13500)) 

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