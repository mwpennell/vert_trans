## fit trap

library(diversitree)

run_mk3_sample <- function(file_name, exp_pr=10, n_step=20000){
  td <- readRDS(file_name)
  out <- lapply(td, function(x) {
    phy <- x$phy
    dat <- x$dat
    lik <- make.mkn(tree= phy, states = dat, k = 3)
    p <- starting.point.musse(phy, 3)
    pp <- p[argnames(lik)]
    names(pp) <- argnames(lik)
    prior <- make.prior.exponential(exp_pr)
    samp.tmp <- mcmc(lik, x.init = pp, w=1, prior=prior,
                     nsteps=100, print.every = 0)
    w <- diff(sapply(samp.tmp[2:7], range))
    mcmc(lik, x.init = pp, w=w, prior=prior,
         nsteps=n_step, print.every=0)
  })
  stem <- strsplit(file_name, ".", fixed=TRUE)[[1]][1]
  saveRDS(out, paste0(stem, "_results", ".rds"))
}

#run_mk3_sample("output/trap/fish_sampdata_2.rds")
run_mk3_sample("output/trap/squa_sampdata_2.rds")
