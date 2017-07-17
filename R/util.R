## phyndr samplr

phyndr_treedata <- function(x, n){
  phy <- x$phy
  dat <- x$dat
  new_phy <- phyndr::phyndr_sample_n(phy, n)
  td_set <- lapply(new_phy, function(y) {
    tmp <- filter(dat, binomial %in% y$tip.label)
    out <- group_by(tmp, binomial) %>% sample_n(1, replace=FALSE)
    rownames(out) <- out$binomial
    class(y) <- "phylo"
    geiger::treedata(y,out, warnings=FALSE)
  })
  td_set
}

## function for getting 95% HPD from mcmc.samples
hdr <- diversitree:::hdr

## function for determining what proportion of 95% HPD overlaps with 0
z.hdr <- function(x){
  hpd <- hdr(x)
  ## all samples in 95% HPD
  nf <- which(which(x >= hpd.lim[1]) %in% which(x <= hpd.lim[2]))
  ## lower than 0
  l <- length(which(x[nf] < 0))/length(x)
  ## greater than 0
  g <- length(which(x[nf] > 0))/length(x)
  list(lower.than.zero=l, greater.than.zero=g)
}

cache.mcmc <- function(x)
  diversitree:::get.cache(get.likelihood(x))


## function for summarizing table for xy 2 zw analysis
build.tab.xz <- function(cache){
  st <- cache$states
  
  mhom <- length(intersect(which(st[,1] == 0), which(st[,2] == 0)))
  mhet <- length(intersect(which(st[,1] == 1), which(st[,2] == 0)))
  fhom <- length(intersect(which(st[,1] == 0), which(st[,2] == 1)))
  fhet <- length(intersect(which(st[,1] == 1), which(st[,2] == 1)))
  uhom <- length(intersect(which(st[,1] == 0), which(is.na(st[,2]))))
  
  dat <- c(mhom, mhet, fhom, fhet, uhom)
  names(dat) <- c("male.hom", "male.het", "fem.hom", "fem.het", "unk.hom")
  dat
}

trans.multitrait <- function(x){
  cc <- cache.mcmc(x)
  lik <- make.musse.multitrait(cc$info$phy, data.frame(cc$states),
                               depth=c(0,0,1))
  r <- lapply(seq_len(nrow(x)), function(i) get.pars(x[i,c(2:11)], lik))
  r <- do.call(rbind,r)
  as.data.frame(r)
}

get.pars <- function(p, lik){
  p.full <- lik(as.numeric(p), pars.only=TRUE)
  p.full[grepl("^q", names(p.full))]
}

unfactor <- function(x)
  as.numeric(levels(x)[as.integer(x)])

state_summary_mk <- function(x, col.names){
  tmp <- lapply(x, function(y){
    cm <- cache.mcmc(y)
      a <- length(which(cm$states == 1))
      b <- length(which(cm$states == 2))
      data.frame(a=a, b=b)
  })
  tmp <- bind_rows(tmp)
  colnames(tmp) <- col.names
  summarise_all(tmp, funs(mean))
}







