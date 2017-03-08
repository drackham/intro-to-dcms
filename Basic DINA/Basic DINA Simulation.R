#########################################################################
#                                                                       #
#     Fit the DINA model to simulated data using the 'dina' package     #
#     Author: Dave Rackham                                              #
#     Created: 3/2/2017                                                 #
#                                                                       #
#########################################################################

#########################################################################
#                                                                       #
#           REMEMBER TO CHANGE YOUR WORKING DIRECTORY!                  #
#                                                                       #
#           Session > Set Working Directory > To Source File Location   #
#                                                                       #
#########################################################################
# Change working directory
setwd("~/dev/intro-to-dcms/Basic DINA")

# Source the helper file which does some of the heavy lifting
source('Basic DINA Simulation Setup.R')

#------------------------------------------------------------------------
# Inspect the dataset
#------------------------------------------------------------------------

# Table of response frequencies
table(data$resp)  # in R when you have a variable and you follow it with $
                  # the variable after is inside of a list in the first var
                  # so myList$varInList

# Mean and standard devition of guessing parameters
mean(data$g)
sd(data$g)

# Mean and standard devition of slipping parameters
mean(data$s)
sd(data$s)

#------------------------------------------------------------------------
# Fit the RDINA model using JAGS
#------------------------------------------------------------------------

# Note: this could take a few seconds!
posterior <- DINA_Gibbs(data$resp, attr.matrix, q, chain_length = chainLength)

#------------------------------------------------------------------------
# Check the slipping parameters for convergence
#------------------------------------------------------------------------

# Create a variable that is the posterior
slipping.posterior <- as.mcmc(t(posterior$SigS))

# Plot the posterior density of the slipping parameters
plot(slipping.posterior)

# Plot the autocorrelation of the slipping parameters
autocorr.plot(slipping.posterior)

#------------------------------------------------------------------------
# Check the guessing parameters for convergence
#------------------------------------------------------------------------

# Create a variable that is the posterior
guessing.posterior <- as.mcmc(t(posterior$GamS))

# Plot the posterior density of the slipping parameters
plot(guessing.posterior)

# Plot the autocorrelation of the slipping parameters
autocorr.plot(guessing.posterior)

#------------------------------------------------------------------------
# Analyze the posterior distribution
#------------------------------------------------------------------------

# Get the mean and standard devition of the guessing posteriors
guess.mean <- apply(posterior$GamS[,burnin:chainLength],1,mean)
guess.sd <- apply(posterior$GamS[,burnin:chainLength],1,sd)

# Get the mean and standard devition of the slipping posteriors
slip.mean <- 1-apply(posterior$SigS[,burnin:chainLength],1,mean)
slip.sd <- apply(posterior$SigS[,burnin:chainLength],1,sd)

# Summarize the posterior distriubtions for all parameters
output=cbind(guess.mean, guess.sd, slip.mean, slip.sd)
colnames(output) = c('g Est','g SE','1-s Est','1-s SE')
rownames(output) = colnames(data$resp)
print(output,digits=3)

#------------------------------------------------------------------------
# Calculate the root mean squared errors (RMSE)
#------------------------------------------------------------------------

# Guessing parameters
guess.rmse <- sqrt((output[,1] - data$g)^2)
mean(guess.rmse)


# Slipping parameters
slip.rmse <- sqrt((output[,3] - (1-data$s))^2)
mean(slip.rmse)
