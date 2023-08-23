# ******************************************************
#
#   Analytical Paleobiology Workshop 2023
#
#   Module 2: Paleodiversity analyses
#   Day 4 | Thursday, August 24th
#
#   Emma Dunne (emma.dunne@fau.de)
#   Lisa Schnetz (lisa.schnetz@gmail.com)
# ______________________________________________________
#
#   3. Sampling standardisation
# 
# ******************************************************



# 0. Packages used in this script -----------------------------------------

library(tidyverse)



# 1. Prepare data ---------------------------------------------------------



dat <- Triassic_div

# Compute diversity indices manually
# Immediate values
S <- nrow(dat) # = number of species
N <- sum(Individuals) # N
cl <- length(levels(factor(Class))) # Number of classes
gen <- length(levels(factor(Genus))) # Number of genera
menh <- S/sqrt(N) # Menhinick Index
marg <- (S-1)/log(N) # Margalef Index




# Get one rarefied diversity value for 100 individuals
gensp <- paste(dat$Genus, dat$Species)
abu <- rep(gensp, dat$Individuals)

trial <- 1000 # subsampling trials
quota <- 100 # subsampling quota
div <- numeric()

for (i in 1:trial) {
  z <- sample(abu, quota) # subsampling without replacement
  div[i] <- length(levels(factor(z)))    
}  


# Expected number of species according to Good (1953)
E.Sm <- S - sum((1-psp)^quota) 

# Empirical rarefaction
# Expected number of species if x individuals are collected

trial <- 200 # subsampling trials
rardiv <- numeric()
erdiv <- numeric()
count <- 1 
sq <- seq(1, 651, by=10)
for (j in sq)  {    # loop for the quota
  
  div <- numeric(trial)
  
  for (i in 1:trial) { # loop for the trial
    z <- sample(abu, j) # subsampling without replacement
    div[i] <- length(levels(factor(z)))    
  }  
  
  # Rarefied diversity
  rardiv[count] <- mean(div)
  erdiv[count] <- sd(div)
  count <- count+1
}

plot(sq, rardiv, type="l")
segments(sq, rardiv-erdiv, sq, rardiv+erdiv)

op <- par(mfrow=c(2,1), lwd=2, mar=c(4,4,1,1))
plot(sq, rardiv, type="l", xlab="", ylab="S")
plot(sq, rardiv, log="xy", type="l", xlab="N", ylab="S")
x <- seq(600); y <- x^0.8
points(x,y, type="l", lty=2)
par(op)

