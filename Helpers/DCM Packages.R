#########################################################################
#                                                                       #
#     Set up the DCM packages environment                  #
#     Author: Dave Rackham                                              #
#     Created: 3/9/2017                                                 #
#                                                                       #
#########################################################################

# Install and load the needed libraries
#------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(devtools, dina, ggmcmc, runjags, coda, fraction.subtraction.data)

# Load the dcms and dcmdata packages written by Dave Rackham
install_github("drackham/dcms", ref="develop")
install_github("drackham/dcmdata", ref="develop")
library('dcms')
library('dcmdata')
