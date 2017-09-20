#--------------------------
# Dynamic Occupancy Model
#--------------------------

sink("SMDynOcc.txt")
cat("
    model {
    
    # Specify priors
    psi1 ~ dunif(-10, 10)
    logit(PSI1) <- psi1
    alpha ~ dunif(-10,10)
    Alpha <- exp(alpha)
    mu_p ~ dunif(-10,10)
    sd_p ~ dunif(0,5)
    
    for(k in 1:nyear-1){
    phi[k] ~ dunif(-4, -2)
    logit(E[k]) <- phi[k]
    }
    
    for(k in 1:nyear){
    p[k] ~ dnorm(mu_p, sd_p)
    logit(P[k]) <- p[k]
    }
    
    # Initial State and Site Connectivities    
    for(m in 1:nsite){
    z[m,1] ~ dbern(PSI1)
    
    for(n in 1:nsite){
    con[m,n] <- exp(-Alpha * dmat[m,n])
    }
    }
    
    # Ecological submodel: Define state conditional on parameters
    for(k in 2:nyear){
    for(i in 1:nsite){
    for(j in 1:nsite){
    disp[i,j,k-1] <- z[j,k-1] * con[i,j]
    }
    s[i,k-1] <- sum(disp[i,,k-1])
    g[i,k-1] <- 1 - exp(-s[i,k-1])
    
    #  Estar[i,k-1] <- E[k-1] *(1-g[i,k-1])
    muZ[i,k-1]<- max(min(z[i,k-1]*(1-E[k-1]) + (1-z[i,k-1])*g[i,k-1],0.99),0.001)
    z[i,k] ~ dbern(muZ[i,k-1])
    } #i
    } #k
    
    # Observation model
    for (i in 1:nsite){
    for (j in 1:nrep){
    for (k in 1:nyear){
    muy[i,j,k] <- z[i,k]*P[k]
    y[i,j,k] ~ dbern(muy[i,j,k])
    } #k
    } #j
    } #i
    
    # Derived parameters
    n.occ[1]<-sum(z[1:nsite,1])
    for (k in 2:nyear){
    n.occ[k] <- sum(z[1:nsite,k])
    }
    }
    ",fill = TRUE)
sink()