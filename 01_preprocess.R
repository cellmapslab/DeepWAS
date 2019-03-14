library(data.table)
library(svd)
library(Rtsne)


`%+%` <- paste0
vol <- './'
genfile <- vol %+% 'MS1_functional.gen'
covarfile <- vol %+% 'MS1_covariates_complete.txt'
phenfile <- vol %+% 'MS1.fam'

read.data <- function(genfile, phenfile, covarfile) {

  gen <- fread(genfile, header = F, sep=' ', data.table=F)
  phen <- fread(phenfile, data.table = F)
  covar <- fread(covarfile, data.table = F)

  covar$SEX   <- factor(covar$SEX, levels = c('male', 'female'))
  covar$PHENO <- factor(covar$PHENO, levels = c('ctrl', 'case'))

  nsnps <- nrow(gen)
  nsamples <- nrow(phen)

  dim(gen)
  dim(phen)
  dim(covar)

  stopifnot((ncol(gen)-6)/3 == nsamples)

  gen[1:10,1:10]

  dom.gen <- 2*as.matrix(gen[,seq(7, ncol(gen), 3)]) +
               as.matrix(gen[,seq(8, ncol(gen), 3)])

  dim(dom.gen) # 37251 x 2894

  # sort covariates based on FID
  covar <- covar[ match(phen$V1, covar$FID), ]

  stopifnot(all(covar$FID == phen$V1))
  stopifnot(all(as.numeric(covar$PHENO) == phen$V6))
  stopifnot(all(as.numeric(covar$SEX) == phen$V5))

  colnames(dom.gen) <- covar$FID
  rownames(dom.gen) <- gen$V3

  dom.gen <- as.data.frame(t(dom.gen))
  final.covar <- covar[, c(3, 5:12)]
  final.covar$SEX <- as.numeric(final.covar$SEX) - 1
  rownames(final.covar) <- rownames(dom.gen)
  final.pheno <- as.numeric(covar$PHENO) - 1

  # Load evalues ------------------------------------------------------------
  snp.eval <- fread(vol %+% 'evalueswithpos.tsv',
                    colClasses = c('snp', 'mineval.feature.uniq', 'mineval.feature', 'eval',
                                   'chr', 'pos', 'al1', 'al2'),
                    data.table=F)

  snp.eval$chrpos <- paste(snp.eval$chr, snp.eval$pos,  sep= ':')

  gen.pos.snp <- gen[, c('V3', 'V2', 'V4')]
  colnames(gen.pos.snp) <- c('gen.snp', 'gen.chr', 'gen.pos')
  gen.pos.snp$chrpos <- paste(gen.pos.snp$gen.chr, gen.pos.snp$gen.pos, sep=':')

  snp.eval <- merge(snp.eval, gen.pos.snp, by = 'chrpos')
  head(snp.eval)

  length(unique(snp.eval$snp))
  snp.eval <- snp.eval[snp.eval$gen.snp%in%colnames(dom.gen), ]
  length(unique(snp.eval$snp))
  snp.eval$snp <- snp.eval$gen.snp

  stopifnot(ncol(dom.gen) == length(unique(snp.eval$snp)))

  snp.eval.split <- split(snp.eval$snp, snp.eval$mineval.feature.uniq)
  snp.eval.split.len <- sapply(snp.eval.split, length)

  # Save names of DeepSEA features
  features.uniq <- unique(snp.eval$mineval.feature.uniq)
  feature.mapping <- unique(snp.eval[,c('mineval.feature.uniq', 'mineval.feature')])
  features <- feature.mapping[match(features.uniq, feature.mapping$mineval.feature.uniq),]$mineval.feature

  res <- list(x=dom.gen,
              covar=final.covar,
              y=final.pheno,
              snp.eval=snp.eval,
              snp.eval.split=snp.eval.split,
              features=features,
              features.uniq=features.uniq)
  res
}

saveRDS(res, 'de1.Rds', compress = F)
