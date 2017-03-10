#########################################################################
#                                                                       #
#     Fit the RDINA model to the Tatsuoka Fraction Subtraction data     #
#     Author: Dave Rackham                                              #
#     Created: 3/9/2017                                                 #
#                                                                       #
#########################################################################

#########################################################################
#                                                                       #
#     REMEMBER TO CHANGE YOUR WORKING DIRECTORY!                        #
#                                                                       #
#     Session > Set Working Directory > To Source File Location         #
#                                                                       #
#########################################################################

#------------------------------------------------------------------------
# Set up the environment
#------------------------------------------------------------------------

# Change working directory
setwd("~/dev/intro-to-dcms/Fraction Subtraction DINA")

# Source the helper files which do some of the heavy lifting
source("../Helpers/extractCodaVariables.R")
source("../Helpers/DCM Packages.R")

#------------------------------------------------------------------------
# Set up all of the simulation parameters
#------------------------------------------------------------------------

# Load the fraction subtraction Q-matrix
q <- as.matrix(fraction.subtraction.qmatrix)

# Specify number of examinees
I <- 536
J <- 20
K <- 8

# Load the fraction subtraction response data
responses <- as.matrix(fraction.subtraction.data)

# Create data object for JAGS
data <- list("I" = I, "J" = J, "K" = K, "resp" = responses)

# Specify which model to use
model <- 'Fraction-Subtraction-RDINA.jags'

# Specify which parameters to save
jags.params = c('fHat', 'dHat', 'alpha1', 'alpha2', 'alpha3', 'alpha4', 'alpha5', 'alpha6', 'alpha7', 'alpha8')

# Generate model
generateFractionSubtractionRDINA()

#------------------------------------------------------------------------
# Fit the model to the data
#------------------------------------------------------------------------
posterior <- rDINAJagsSim(data, jagsModel = model, jags.params = jags.params,
                          maxCores = 1, adaptSteps = 100, burnInSteps = 500,
                          numSavedSteps = 5000, thinSteps = 3)

#------------------------------------------------------------------------
# Prep the posterior for analysis
#------------------------------------------------------------------------

# Convert to coda object
codaSamples = as.mcmc.list(combine.mcmc(posterior)) # resulting codaSamples object has these indices: codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
codaList <- as.mcmc.list(posterior)

# Setup ggmcmc object
S <- ggs(as.mcmc.list(posterior))

#------------------------------------------------------------------------
# Get summaries of d and f
#------------------------------------------------------------------------
dHatSummary <- add.summary(posterior, vars = c("dHat"))
dHatSummary

fHatSummary <- add.summary(posterior, vars = c("fHat"))
fHatSummary

#------------------------------------------------------------------------
# Check d and f for convergence
#------------------------------------------------------------------------

# d
ggmcmc(S, file="dHat.pdf", plot=c("traceplot", "density", "autocorrelation"), family="dHat", param_page=4)

# f
ggmcmc(S, file="fHat.pdf", plot=c("traceplot", "density", "autocorrelation"), family="fHat", param_page=4)

# extract d and f means from posterior
dHat <- extractCodaVariables(x=codaSamples, params='dHat', exact=FALSE)
fHat <- extractCodaVariables(x=codaSamples, params='fHat', exact=FALSE)

alpha1 <- extractCodaVariables(x=codaSamples, params='alpha1', exact=FALSE)
alpha2 <- extractCodaVariables(x=codaSamples, params='alpha2', exact=FALSE)
alpha3 <- extractCodaVariables(x=codaSamples, params='alpha3', exact=FALSE)
alpha4 <- extractCodaVariables(x=codaSamples, params='alpha4', exact=FALSE)
alpha5 <- extractCodaVariables(x=codaSamples, params='alpha5', exact=FALSE)
alpha6 <- extractCodaVariables(x=codaSamples, params='alpha6', exact=FALSE)
alpha7 <- extractCodaVariables(x=codaSamples, params='alpha7', exact=FALSE)
alpha8 <- extractCodaVariables(x=codaSamples, params='alpha8', exact=FALSE)

g <- lapply(fHat[,1], function(x) exp(x) / (1 + exp(x)))

s <- list()
for (i in 1:nrow(fHat)){
  s[i] <- 1 - exp(fHat[i,1] + dHat[i, 1]) / (1 + exp(fHat[i,1] + dHat[i, 1]))
}

