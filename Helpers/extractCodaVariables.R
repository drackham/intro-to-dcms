################ Extract coda variables #######################
# https://github.com/johnbaums/jagstools/blob/master/R/jagsresults.R

getMCMCSample <- function(mcmclist,i){
  chainext <- function(x,i) return(x[,i])
  return(as.mcmc.list(lapply(mcmclist, chainext, i = i)))
}


extractCodaVariables <- function(x, params, invert = FALSE, exact = TRUE, regex = FALSE,
                                 probs=c(0.025, 0.975), ...) {
  if(!regex) {
    params <- paste0(gsub('(?=\\.|\\[|\\])', '\\1\\\\', params, perl=TRUE),
                     '(\\[.*\\])?', collapse='|')
    if (exact) params <- paste("^",
                               gsub("\\|", "$|^", params),
                               "$", sep = "")
  } else {
    if(exact) warning('exact=TRUE ignored when regex=TRUE')
  }
  switch(is(x)[1],
         rjags={
           stop('Please convert object to type mcmc.list.')
#            For now only supporting mcmc.list types
#            nm <- dimnames(x$BUGSoutput$sims.array)[[3]]
#            i <- grep(params, nm, invert=invert, ...)
#            if(length(i) == 0) stop('No parameters match params', call.=FALSE)
#            samp <- x$BUGSoutput$sims.array[, , i, drop=FALSE]
#            rhat_neff <- x$BUGSoutput$summary[i, c('Rhat', 'n.eff'), drop=FALSE]
#            return(cbind(t(apply(
#              samp, 3, function(x)
#                c(mean=mean(x), sd=sd(x), quantile(x, probs=probs)))), rhat_neff))
         },
         mcmc.list={
           nm <- colnames(x[[1]])
           i <- grep(params, nm, invert=invert)

           if(length(i) == 0) stop('No parameters match params', call.=FALSE)
           # Use nm and i to extract the specific samples that we are interested in
           xExtracted <-as.mcmc.list(getMCMCSample(x, i[1]:tail(i, n=1)))

           # Collapse the extracted samples into one data frame
           t(apply(do.call(rbind, xExtracted), 2, function(z) { # 2 indicates apply over columns
           #t(apply(do.call(rbind, x), 2, function(z) { # 2 indicates apply over columns
             c(mean=mean(z), sd=sd(z), quantile(z, probs))
           }))
         },
         {
           stop('x must be an mcmc.list or rjags object.')
         }
  )
}