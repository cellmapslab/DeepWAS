library(dplyr)
library(readr)
library(stabs)

ncores <- 18

fit.stabs <- function(res, stabs.cutoff=0.7, stabs.PFER=0.5) {
  fits <- mclapply(res$features.uniq, function(feature) {
    print(feature)
    mat <- cbind(res$x[, unique(res$snp.eval.split[[feature]])],
                 res$covar)
    f <- try(stabsel(mat, res$y, cutoff = stabs.cutoff, PFER = stabs.PFER,
                 fitfun = glmnet.lasso, args.fitfun = list(alpha=1),
                 papply=lapply))
    if(class(f)=='try-error')return(NULL)
    f
  }, mc.cores = ncores)
  names(fits) <- fits$features.uniq
  fits
}

res <- readRDS('de1.Rds')

stabs.cutoffs <- c(0.7, 0.9)
stabs.PFERs <- c(0.1, 1.0)

result <- lapply(stabs.cutoffs, function(cutoff){
  lapply(stabs.PFERs, function(PFER){
    cat(sprintf('Running stabs with Pfer:%f and cutoff:%f...\n', PFER, cutoff))
    f <- fit.stabs(res, stabs.cutoff = cutoff, stabs.PFER = PFER)
    f
  })
})

saveRDS(result, 'de1-results.Rds', compress=F)