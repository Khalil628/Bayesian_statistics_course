### 2021-10-26 F Barraquand -- Bayesian Statistics Course for PhD Students

################# TD2 - linear model #######################

## From Chapter 5 in Kéry (2010)

###  1. The simplest model, for a single mean

pdf(width=10,height=6,file="Peregrine_male_mass.pdf")
### Data generation
# Generate two samples of body mass measurements of male peregrines
y10 <- rnorm(n = 10, mean = 600, sd = 30) # Sample of 10 birds
y1000 <- rnorm(n = 1000, mean = 600, sd = 30) # Sample of 1000 birds

# Plot data
xlim = c(450, 750)
par(mfrow = c(2,1))
hist(y10, col = 'grey ', xlim = xlim, main = 'Body mass (g) of 10 male peregrines')
hist(y1000, col = 'grey', xlim = xlim, main = ' Body mass (g) of 1000 male peregrines')
dev.off()

### Analysis using R
summary(lm(y1000 ~ 1))

###Analysis using JAGS
library(R2jags)		                # Load R2jags // we use R2jags while the last version of Kéry(2010) uses package jagsUI

# Bundle and summarize the data set passed to JAGS
str(jags.data <- list(mass = y1000, nobs = length(y1000)))

# Specify model in BUGS language
cat(file = "model.txt", "
model {

# Priors
population.mean ~ dunif(0,5000)         # Normal parameterized by precision
precision <- 1 / population.variance    # Precision = 1/variance
population.variance <- population.sd * population.sd
population.sd ~ dunif(0,100)

# Likelihood
for(i in 1:nobs){
  mass[i] ~ dnorm(population.mean, precision)
}

}
")

# Initial values
inits <- function(){list(population.mean = rnorm(1,600), population.sd = runif(1, 1, 30))}

# Parameters to be monitored (= to estimate)
params <- c("population.mean", "population.sd", "population.variance")

# MCMC settings
nc <- 3					# Number of chains
ni <- 1000				# Number of draws from posterior (for each chain)
nb <- 1					# Number of draws to discard as burn-in
nt <- 1					# Thinning rate

# Call JAGS from R, check convergence and summarize posteriors
out <- jags(jags.data, inits, parameters=params, model.file="model.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out, dig = 2)



### 2. Now we do the T-test // From Chapter 7 in Kéry (2010)

### T-test with equal variances

### Data generation
n1 <- 6				# Number of females
n2 <- 4				# Number of males
mu1 <- 105				# Population mean of females
mu2 <- 77.5				# Population mean of males
sigma <- 10*2.75				# Average population SD of both

n <- n1+n2				# Total sample size
y1 <- rnorm(n1, mu1, sigma)		# Data for females
y2 <- rnorm(n2, mu2, sigma)		# Date for males
y <- c(y1, y2)				# Aggregate both data sets
x <- rep(c(0,1), c(n1, n2))		# Indicator for male
boxplot(y ~ x, col = "grey", xlab = "Male", ylab = "Wingspan (cm)", las = 1)

n <- n1+n2				# Total sample size
alpha <- mu1				# Mean for females serves as the intercept
beta <- mu2-mu1				# Beta is the difference male-female
E.y <- alpha + beta*x			# Expectation
y.obs <- rnorm(n = n, mean = E.y, sd = sigma)	# Add random variation
boxplot(y.obs ~ x, col = "grey", xlab = "Male", ylab = "Wingspan (cm)", las = 1)

### Analysis using R
fit1 <- lm(y ~ x)			# Analysis of first data set
fit2 <- lm(y.obs ~ x)			# Analysis of second data set
summary(fit1)
summary(fit2)

anova(fit1)
anova(fit2)

#model.matrix(fit1)
#model.matrix(fit2)

# Bundle and summarize the data set passed to JAGS
str(jags.data <- list(x = x, y = y, n = n))  # passing the data

# [When I write ... you are meant to fill in the right code instead!]

# Specify model in BUGS language
cat(file = "ttest.txt", "
model {

# Priors
 mu1 ~ dnorm(0,0.001)			  # Beware dnorm(mean,precision)
 delta ~ dnorm(0,0.001)			# Large variance = Small precision
 tau <-pow(sigma,-2)
 sigma ~ dunif(0, 250) #dexp(1) #n'importe quoi

# Likelihood
 for (i in 1:n) {
    y[i] ~ dnorm(mu[i], tau) 
    mu[i] <- mu1 + delta * x[i]
    residual[i] <- y[i] - mu[i]		# Define residuals
 }

# Derived quantities: one of the greatest things about a Bayesian analysis
 mu2 <- mu1 + delta			# Difference in wingspan
}
")

# Inits function
inits <- function(){list(mu1=rnorm(1), delta=rnorm(1), sigma = rlnorm(1))}

# Parameters to estimate
params <- c("mu1","mu2", "delta", "sigma", "residual")

# MCMC settings
nc <- 3	;  ni <- 1000  ;  nb <- 1  ;  nt <- 1

# Call JAGS from R, check convergence and summarize posteriors
out <- jags(jags.data, inits, parameters=params, model.file="ttest.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out, dig = 2)

# Let's compare this to the classical function in R
# Does the 95% CI for delta resembles the 95% CRI ? 
# CRI 
quantile(out$BUGSoutput$sims.list$delta,c(0.025,0.975))
# CI
t.test(y1,y2,var.equal = TRUE)

## However with Bayes we can compute stuff like Pr(delta<0|data) which more handy than the p-values
mean(out$BUGSoutput$sims.list$delta < 0)

# Follow-up you want to get a better intuition behind tests -> decrease the sample size and increase the SD until the p-value is not significant. 



### 3. One-way analysis of variance

# (Fixed-effects ANOVA) // from Chapter 9 in Kéry (2010)
# We study how the snout-vent length (SVL), a measure of size in herpetology, varies across 5 different populations.
# We take 10 snakes from each population. 

# Data generation
ngroups <- 5				# Number of populations
nsample <- 10				# Number of snakes in each
pop.means <- c(50, 40, 45, 55, 60) 	# Population mean SVL
sigma <- 3				# Residual sd

n <- ngroups * nsample 			# Total number of data points
eps <- rnorm(n, 0, sigma)		# Residuals 
x <- rep(1:5, rep(nsample, ngroups)) 	# Indicator for population
means <- rep(pop.means, rep(nsample, ngroups))
X <- as.matrix(model.matrix(~ as.factor(x)-1)) # Create design matrix
X					# Inspect that
y <- as.numeric(X %*% as.matrix(pop.means) + eps) # assemble -- NOTE: as.numeric ESSENTIAL for WinBUGS
# %*% denotes matrix multiplication
boxplot(y~x, col="grey", xlab="Population", ylab="SVL", main="", las = 1)


## Maximum likelihood analysis using R
print(anova(lm(y~as.factor(x))))
cat("\n")
print(summary(lm(y~as.factor(x)))$coeff, dig = 3)
cat("Sigma:         ", summary(lm(y~as.factor(x)))$sigma, "\n")


## Bayesian analysis using JAGS

# Bundle and summarize the data set passed to JAGS
str(bdata <- list(y = y, x = x))

# Specify model in BUGS language
cat(file = "anova.txt", "
    model {
    
    # Priors
    for (i in 1:5){			# Implicitly define alpha as a vector
    alpha[i] ~ dnorm(50,0.01)
    }
    sigma ~ dunif(0,50)
    tau <- pow(sigma,-2)
    
    # Likelihood
    for (i in 1:50) {
    y[i] ~ dnorm(mu[i],tau)
    mu[i] <- alpha[x[i]]
    }
    
    # Derived quantities
    effect2 <- alpha[2] - alpha[1]
    effect3 <- alpha[3] - alpha[1]
    effect4 <- alpha[4] - alpha[1]
    effect5 <- alpha[5] - alpha[1]
    
    # Custom ``hypothesis'' test / Define your own contrasts
    test1 <- (effect2+effect3) + (effect4+effect5) # Equals zero when 2+3 = 4+5
    test2 <- effect5 - 2 * effect4 		# Equals zero when effect5 = 2*effect4
    }
    ")

# Inits function
inits <- function(){ list(alpha = rnorm(5, mean = mean(y)), sigma = rlnorm(1) )}

# Parameters to estimate
params <- c("alpha", "sigma", "effect2", "effect3", "effect4", "effect5", "test1", "test2")

# MCMC settings
 nc <- 3  ;  ni <- 1200  ;  nb <- 200  ;  nt <- 2

# Call JAGS, check convergence and summarize posteriors
out <- jags(bdata, inits, params, "anova.txt", n.thin = nt, n.chains = nc, 
n.burnin = nb, n.iter = ni)
print(out, dig = 3)
