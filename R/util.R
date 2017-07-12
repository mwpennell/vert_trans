## phyndr samplr

phyndr_samplr <- function(x, n){
  phy <- x$phy
  dat <- x$dat
  new_phy <- phyndr::phyndr_sample_n(phy, n)
  td_set <- lapply(new_phy, function(y) {
    tmp <- filter(dat, binomial %in% y$tip.label)
    out <- group_by(tmp, binomial) %>% sample_n(1, replace=FALSE)
    rownames(out) <- out$binomial
    geiger::treedata(y,out, warnings=FALSE)
  })
  td_set
}






