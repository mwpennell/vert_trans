source("R/util.R")
source("R/plotting.R")

## ESD to GSD

esd.f <- readRDS("output/esd-gsd/fish_sampdata_10_results.rds")
esd.s <- readRDS("output/esd-gsd/squa_sampdata_10_results.rds")

## Remove burnin from all samples
esd.f <- lapply(esd.f, function(x) x[-seq_len(10000),])
esd.f <- dplyr::bind_rows(esd.f)
esd.s <- lapply(esd.s, function(x) x[-seq_len(10000),])
esd.s <- dplyr::bind_rows(esd.s)

esd.f$d <- esd.f$q10 - esd.f$q01

length(which(esd.f$d < 0))/nrow(esd.f)

esd.f$m <- esd.f$q10 / esd.f$q01
median(esd.f$m)

esd.s$d <- esd.s$q10 - esd.s$q01

length(which(esd.s$d < 0))/nrow(esd.s)

esd.s$m <- esd.s$q10 / esd.s$q01
median(esd.s$m)


p <- ggplot(esd.f, aes(x=d))
p <- p + geom_histogram(bins = 50, fill=cols[2], alpha=1) + xlim(-0.12,0.12)
p <- p + ylab("Posterior density") 
p <- p +xlab(bquote('Net rate of transitions '*my^-1*''))
p <- p + geom_vline(aes(xintercept=0), colour=cols[1], alpha=0.75)
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank())
p
ggsave(filename = "figs/esd-gsd-fish_v2.pdf")
q <- ggplot(esd.s, aes(x=d))
q <- q + geom_histogram(bins = 50, fill=cols[2], alpha=1) + xlim(-0.025, 0.025)
q <- q + ylab("Posterior density") 
q <- q + xlab(bquote('Net rate of transitions '*my^-1*''))
q <- q + geom_vline(aes(xintercept=0), colour=cols[1], alpha=0.75)
q <- q + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank(),
               axis.title.y=element_blank())
q
ggsave("figs/esd-gsd-squa_v2.pdf")


## With ambiguous species removed


esd.f <- readRDS("output/esd-gsd/fish_sampdata_10_results_alt.rds")
esd.s <- readRDS("output/esd-gsd/squa_sampdata_10_results_alt.rds")

## Remove burnin from all samples
esd.f <- lapply(esd.f, function(x) x[-seq_len(10000),])
esd.f <- dplyr::bind_rows(esd.f)
esd.s <- lapply(esd.s, function(x) x[-seq_len(10000),])
esd.s <- dplyr::bind_rows(esd.s)

esd.f$d <- esd.f$q10 - esd.f$q01

length(which(esd.f$d < 0))/nrow(esd.f)

esd.f$m <- esd.f$q10 / esd.f$q01
median(esd.f$m)

esd.s$d <- esd.s$q10 - esd.s$q01

length(which(esd.s$d < 0))/nrow(esd.s)

esd.s$m <- esd.s$q10 / esd.s$q01
median(esd.s$m)


p <- ggplot(esd.f, aes(x=d))
p <- p + geom_histogram(bins = 50, fill=cols[2], alpha=1) + xlim(-0.075,0.075)
p <- p + ylab("Posterior density") 
p <- p +xlab(bquote('Net rate of transitions '*my^-1*''))
p <- p + geom_vline(aes(xintercept=0), colour=cols[1], alpha=0.75)
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank())
p
ggsave(filename = "figs/esd-gsd-fish-alt.pdf")
q <- ggplot(esd.s, aes(x=d))
q <- q + geom_histogram(bins = 50, fill=cols[2], alpha=1) + xlim(-0.025, 0.025)
q <- q + ylab("Posterior density") 
q <- q + xlab(bquote('Net rate of transitions '*my^-1*''))
q <- q + geom_vline(aes(xintercept=0), colour=cols[1], alpha=0.75)
q <- q + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank(),
               axis.title.y=element_blank())
q
ggsave("figs/esd-gsd-squa-alt.pdf")
