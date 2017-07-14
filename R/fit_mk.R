## Functions for running Mk models

library(diversitree)

run_mk_sample <- function(file_name, exp_pr=10, n_step=50000){
  td <- readRDS(file_name)
  out <- lapply(td, function(x) {
    phy <- x$phy
    dat <- x$dat
    lik <- make.mk2(phy, dat)
    p <- starting.point.bisse(phy, dat)
    pp <- p[c("lambda0", "lambda1")]
    names(pp) <- argnames(lik)
    prior <- make.prior.exponential(exp_pr)
    samp.tmp <- mcmc(lik, x.init = pp, w=1, prior=prior,
                     nsteps=100, print.every = 0)
    w <- diff(sapply(samp.tmp[2:3], range))
    mcmc(lik, x.init = pp, w=w, prior=prior,
                      nsteps=n_step, print.every=0)
  })
  stem <- strsplit(file_name, ".", fixed=TRUE)[[1]][1]
  saveRDS(out, paste0(stem, "_results", ".rds"))
}

run_mk_sample("output/esd-gsd/fish_sampdata_10.rds")
run_mk_sample("output/esd-gsd/squa_sampdata_10.rds")
