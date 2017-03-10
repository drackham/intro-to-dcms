#########################################################################
#                                                                       #
#     Set up the 'Basic DINA Simulation.R' environment                  #
#     Author: Dave Rackham                                              #
#     Created: 3/2/2017                                                 #
#                                                                       #
#########################################################################

# Load the simulated dataset and create the needed parameters
#------------------------------------------------------------------------
data <- get('R_DINA_SimpleQ_LN.1000') # This is a low-noise data set I simulated
                                      # See the dcmdata R package for details

# Load the Q-matrix
q <- simpleQ()

# Specify the number of attributes
K <- 2

# Creating matrix of possible attribute profiles
attr.matrix <- rep(0,K)
for(j in 1:K){
  temp = combn(1:K,m=j)
  temp.matrix = matrix(0,ncol(temp),K)
  for(j in 1:ncol(temp)) temp.matrix[j,temp[,j]] = 1
  attr.matrix = rbind(attr.matrix,temp.matrix)
}
attr.matrix <- as.matrix(attr.matrix)

# Specify the burnin and chain length for the MCMC analysis
burnin <- 500
chainLength <- 5000
