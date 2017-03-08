#########################################################################
#                                                                       #
#     Fit the DINA model to simulated data using JAGS                   #
#     Author: Dave Rackham                                              #
#     Created: 3/6/2017                                                 #
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
setwd("~/dev/intro-to-dcms/JAGS DINA")

# Source the helper file which does some of the heavy lifting
source("../Helpers/extractCodaVariables.R")

#------------------------------------------------------------------------
# Load all of the needed packages
#------------------------------------------------------------------------

library('devtools')
# install_github("drackham/dcms", ref="e5931eb6f2262bf72ccb9ee973ed167764dc9e31")
library('dcms')
library('coda')
library('ggplot2')
library('ggmcmc')
library('parallel')
library('runjags')

#------------------------------------------------------------------------
# Set up all of the simulation parameters
#------------------------------------------------------------------------

model <- 'R-DINA-Non-Hierarchical.jags'
N <- 200
max_cores <- 1
iter <- 5000
chains <- 1

# Load and save the simulated data see: http://stackoverflow.com/questions/9083907/r-how-to-call-an-object-with-the-character-variable-of-the-same-name
data <- R_DINA_SimpleQ_LN.200

# Specify which model to use
generateRDINAJagsNonHierachical()

#------------------------------------------------------------------------
# Fit the model to the data
#------------------------------------------------------------------------

posterior <- rDINAJagsSim(data, jagsModel = model,
                    maxCores = max_cores, adaptSteps = 100, burnInSteps = 1000,
                    numSavedSteps = iter, thinSteps = 1)

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

# Make vector of wanted parameter names
wanted_pars <- c(paste0("dHat[", 1:data$J, "]"))

# Get estimated and generating values for wanted parameters
generating_values = c(data$d)
estimated_values <- dHat

# Assesmble a data frame to pass to ggplot()
sim_df <- data.frame(parameter = factor(wanted_pars, rev(wanted_pars)), row.names = NULL)
sim_df$middle <- estimated_values[, "mean"] - generating_values
sim_df$lower <- estimated_values[, "2.5%"] - generating_values
sim_df$upper <- estimated_values[, "97.5%"] - generating_values

# Plot the discrepancy
ggplot(sim_df) + aes(x = parameter, y = middle, ymin = lower, ymax = upper) +
  scale_x_discrete() + geom_abline(intercept = 0, slope = 0, color = "white") +
  geom_linerange() + geom_point(size = 2) + labs(y = "Discrepancy", x = NULL, title = "dHat Sim vs. Predicted Discrepency") +
  theme(panel.grid = element_blank()) + coord_flip()

# Make vector of wanted parameter names
wanted_pars <- c(paste0("fHat[", 1:data$J, "]"))

# Get estimated and generating values for wanted parameters
generating_values = c(data$f)
estimated_values <- fHat

# Assesmble a data frame to pass to ggplot()
sim_df <- data.frame(parameter = factor(wanted_pars, rev(wanted_pars)), row.names = NULL)
sim_df$middle <- estimated_values[, "mean"] - generating_values
sim_df$lower <- estimated_values[, "2.5%"] - generating_values
sim_df$upper <- estimated_values[, "97.5%"] - generating_values

# Plot the discrepancy
ggplot(sim_df) + aes(x = parameter, y = middle, ymin = lower, ymax = upper) +
  scale_x_discrete() + geom_abline(intercept = 0, slope = 0, color = "white") +
  geom_linerange() + geom_point(size = 2) + labs(y = "Discrepancy", x = NULL, title = "fHat Sim vs. Predicted Discrepency") +
  theme(panel.grid = element_blank()) + coord_flip()

# Check alpha recovery
plot(data$alphaIK[,1], alpha1[,1], main = "Alpha1 Recovery")
abline(a=0,b=1)

plot(data$alphaIK[,2], alpha2[,1], main = "Alpha2 Recovery")
abline(a=0,b=1)

# Check Prob(response) recovery
probHat <- boot::inv.logit(fHat + dHat)
prob <- boot::inv.logit(data$f + data$d)

diff <- (prob - probHat)
mean(diff[,1])

# Convert d and f to g and s
gSim <- data$g
sSim <- boot::inv.logit(data$f + data$d)

gHat <- boot::inv.logit(fHat[,1])
sHat <- boot::inv.logit(fHat[,1] + dHat[,1])

# Inspect simulated vs. predicted
plot(gHat, gSim, main = "Guess simulated vs. predicted", xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1)

plot(sHat, sSim, main = "Slip simulated vs. predicted", xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1)
