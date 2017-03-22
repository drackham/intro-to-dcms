#########################################################################
#                                                                       #
#     Fit the RDINA model to the Tatsuoka Fraction Subtraction data     #
#     without items 3,7 and 9                                           #
#     Author: Dave Rackham                                              #
#     Created: 3/20/2017                                                #
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
J <- 17
K <- 8

# Load the fraction subtraction response data
responses <- as.matrix(fraction.subtraction.data)
responses.partial <- responses[,-c(3,7,9)]

# Create data object for JAGS
data <- list("I" = I, "J" = J, "K" = K, "resp" = responses.partial)

# Specify which model to use
model <- 'Partial-Fraction-Subtraction-RDINA.jags'

# Specify which parameters to save
jags.params = c('fHat', 'dHat', 'alpha1', 'alpha2', 'alpha3', 'alpha4', 'alpha5', 'alpha6', 'alpha7', 'alpha8')


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

#------------------------------------------------------------------------
# Extract means and and sds from posteriors
#------------------------------------------------------------------------

# extract d and f 
dHat <- extractCodaVariables(x=codaSamples, params='dHat', exact=FALSE)
fHat <- extractCodaVariables(x=codaSamples, params='fHat', exact=FALSE)

# extract alphas
alpha1 <- extractCodaVariables(x=codaSamples, params='alpha1', exact=FALSE)
alpha2 <- extractCodaVariables(x=codaSamples, params='alpha2', exact=FALSE)
alpha3 <- extractCodaVariables(x=codaSamples, params='alpha3', exact=FALSE)
alpha4 <- extractCodaVariables(x=codaSamples, params='alpha4', exact=FALSE)
alpha5 <- extractCodaVariables(x=codaSamples, params='alpha5', exact=FALSE)
alpha6 <- extractCodaVariables(x=codaSamples, params='alpha6', exact=FALSE)
alpha7 <- extractCodaVariables(x=codaSamples, params='alpha7', exact=FALSE)
alpha8 <- extractCodaVariables(x=codaSamples, params='alpha8', exact=FALSE)

# combine alphas into a mastery vector
alphas <- matrix(c(alpha1[,1], alpha2[,1], alpha3[,1], alpha4[,1], alpha5[,1], alpha6[,1], alpha7[,1], alpha8[,1]),
                 nrow = nrow(alpha1), ncol = 8)

#------------------------------------------------------------------------
# Convert d and f into g and s
#------------------------------------------------------------------------

g <- lapply(fHat[,1], function(x) exp(x) / (1 + exp(x)))

s <- list()
for (i in 1:nrow(fHat)){
  s[i] <- 1 - exp(fHat[i,1] + dHat[i, 1]) / (1 + exp(fHat[i,1] + dHat[i, 1]))
}

#------------------------------------------------------------------------
# Analyze specific items
#------------------------------------------------------------------------

# Item analysis
fHat1.all <- codaSamples[[1]][,1]
fHat3.all <- codaSamples[[1]][,3]
plot(density(fHat1.all))

gHat1.all <- exp(fHat1.all) / (1 + exp(fHat1.all))
gHat3.all <- exp(fHat3.all) / (1 + exp(fHat3.all))

plot(density(gHat3.all))
mean(gHat3.all)

#------------------------------------------------------------------------
# Analyze examinees who got most items wrong
#------------------------------------------------------------------------

# Get list of examinees who got most items wrong
total <- rowSums(responses)
responses.with.total <- data.frame(responses, total)
most.wrong <- which(responses.with.total$total < 2)

mastery.most.wrong <- alphas[most.wrong,]

# Correlation between Alphas
cor(alphas)
