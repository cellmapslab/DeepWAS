library(data.table)
library(dplyr)

exc <- function(x) setdiff(x, c('SEX', 'CHIP', paste0('C', 1:8)))
selected <- function(e){e<-e[!is.na(e)]; exc(unique(names(unlist(sapply(e, `[[`, 'selected')))))}
selected.single <- function(e){e<-e[!is.na(e)]; paste(exc(unique(names(e$selected))), collapse = ',')}

extract.df <- function(res, features, features.uniq) {
  
  cutoffs <- c(0.7, 0.9)
  pfers <- c(0.1, 1.0)
  
  do.call(rbind.data.frame, lapply(seq_along(res), function(i){
    
    cutoff <- cutoffs[i]
    
    l <- lapply(seq_along(res[[i]]), function(j){
      
      pfer <- pfers[j]
      
      feature.map <- data.frame(features=features,
                                features.uniq=features.uniq,
                                cell.line=sapply(strsplit(features, '\\|'), `[[`, 1),
                                tf=sapply(strsplit(features, '\\|'), `[[`, 2),
                                treatment=sapply(strsplit(features, '\\|'), `[[`, 3),
                                hits=gsub('`', '', sapply(res[[i]][[j]], selected.single)),
                                cutoff=cutoff,
                                pfer=pfer,
                                stringsAsFactors = F)
      
      feature.map <- transform(feature.map, hits=strsplit(hits, ',')) %>% unnest(hits)
      feature.map
    })
    do.call(rbind.data.frame, l)
  }))
}

snp.pos <- fread('MS1_functional_positions.gen', data.table=F, col.names = c('snp', 'chr', 'pos'))
snp.pos <- unique(snp.pos)


de1.res <- readRDS('de1-results.Rds')
de1 <- readRDS('de1.Rds')
de1.df <- extract.df(de1.res, de1$features, de1$features.uniq)
de1.df <- merge(de1.df, snp.pos, by.x = 'hits', by.y='snp', all.x = T)
write.table(de1.df, 'de1-results.tsv', row.names = F, quote = F, sep='\t')

