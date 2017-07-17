library(diversitree)

run_mkmulti_sample <- function(file_name, exp_pr=10, n_step=50000){
  td <- readRDS(file_name)
  out <- lapply(td, function(x) {
    phy <- x$phy
    dat <- x$dat
    ## multitrait version
    lik <- make.musse.multitrait(phy, dat, depth=c(0,0,1))
    p <- starting.point.musse.multitrait(phy, lik)
    prior.par <- make.prior.exponential(exp_pr)
    prior <- function(pars)
      prior.par(lik(pars,pars.only = TRUE))
    
    samp.tmp <- mcmc(lik, x.init = p, w=1, prior=prior,
                     nsteps=100, print.every = 0)
    w <- diff(sapply(samp.tmp[2:11], range))
    mcmc(lik, x.init = p, w=w, prior=prior,
                      nsteps=n_step, print.every=0)
  })
  stem <- strsplit(file_name, ".", fixed=TRUE)[[1]][1]
  saveRDS(out, paste0(stem, "_results", ".rds"))
}

## HET-HOM
run_mkmulti_sample("output/het-hom/fish_sampdata_10.rds")
run_mkmulti_sample("output/het-hom/squa_sampdata_10.rds")
