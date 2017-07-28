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
                      nsteps=n_step, print.every=500)
  })
  stem <- strsplit(file_name, ".", fixed=TRUE)[[1]][1]
  saveRDS(out, paste0(stem, "_results", ".rds"))
}

## Parallelize
run_mkmulti_single <- function(x, exp_pr=10, n_step=50000){
    phy <- x$phy
    dat <- x$dat
    ## multitrait version
    lik <- diversitree::make.musse.multitrait(phy, dat, depth=c(0,0,1))
    p <- diversitree::starting.point.musse.multitrait(phy, lik)
    prior.par <- diversitree::make.prior.exponential(exp_pr)
    prior <- function(pars)
      prior.par(lik(pars,pars.only = TRUE))
    
    samp.tmp <- diversitree::mcmc(lik, x.init = p, w=1, prior=prior,
                     nsteps=100, print.every = 0)
    w <- diff(sapply(samp.tmp[2:11], range))
    out <- diversitree::mcmc(lik, x.init = p, w=w, prior=prior,
         nsteps=n_step, print.every=1000)
    out
}

library(parallel)
n_cores <- detectCores() - 2
cl <- makeCluster(n_cores)
clusterExport(cl, "run_mkmulti_single")
ff <- readRDS("output/het-hom/fish_sampdata_10.rds")
out <- parLapply(cl, ff, function(x) run_mkmulti_single(x))
saveRDS(out, "output/het-hom/fish_sampdata_10_results.rds")

#aa <- readRDS("output/het-hom/amph_sampdata_10.rds")
#out <- parLapply(cl, ff, function(x) {
#  run_mkmulti_single(x)
#})
#saveRDS(out, "output/het-hom/amph_sampdata_10_results.rds")
stopCluster(cl)



#run_mkmulti_sample("output/het-hom/fish_sampdata_10.rds")
#run_mkmulti_sample("output/het-hom/amph_sampdata_10.rds")
#run_mkmulti_sample("output/het-hom/fish_sampdata_2.rds", n_step = 10000)
#run_mkmulti_sample("output/het-hom/amph_sampdata_2.rds", n_step = 10000)