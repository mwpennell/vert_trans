## analysis of xy vs zw

source("R/util.R")
source("R/plotting.R")

## XY to ZW

xy.f <- readRDS("output/xy-zw/fish_sampdata_10_results.rds")
xy.s <- readRDS("output/xy-zw/squa_sampdata_10_results.rds")
xy.a <- readRDS("output/xy-zw/amph_sampdata_10_results.rds")

## Remove burnin from all samples
xy.f <- lapply(xy.f, function(x) x[-seq_len(10000),])
xy.f <- dplyr::bind_rows(xy.f)

xy.f$d <- -(xy.f$q01 - xy.f$q10)

length(which(xy.f$d < 0))/nrow(xy.f)

p <- ggplot(xy.f, aes(x=d))
p <- p + geom_histogram(bins = 50, fill=cols[2], alpha=1) + xlim(c(-0.12, 0.12))
p <- p + ylab("Posterior density") 
p <- p +xlab(bquote('Net rate of transitions '*my^-1*''))
p <- p + geom_vline(aes(xintercept=0), colour=cols[1], alpha=0.75)
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank())
p
ggsave(filename = "figs/xy-zw-fish-v2.pdf")

## Remove burnin from all samples
xy.s <- lapply(xy.s, function(x) x[-seq_len(10000),])
xy.s <- dplyr::bind_rows(xy.s)

xy.s$d <- xy.s$q10 - xy.s$q01

length(which(xy.s$d < 0))/nrow(xy.s)

p <- ggplot(xy.s, aes(x=d))
p <- p + geom_histogram(bins = 50, fill=cols[2], alpha=1) + xlim(c(-0.007, 0.007))
p <- p +xlab(bquote('Net rate of transitions '*my^-1*''))
p <- p + geom_vline(aes(xintercept=0), colour=cols[1], alpha=0.75)
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank(),
               axis.title.y=element_blank())
p
ggsave(filename = "figs/xy-zw-squa-v2.pdf")

## Remove burnin from all samples
xy.a <- lapply(xy.a, function(x) x[-seq_len(10000),])
xy.a <- dplyr::bind_rows(xy.a)

xy.a$d <- xy.a$q10 - xy.a$q01

length(which(xy.a$d < 0))/nrow(xy.a)

p <- ggplot(xy.a, aes(x=d))
p <- p + geom_histogram(bins = 50, fill=cols[2], alpha=1) + xlim(c(-0.06, 0.06))
p <- p +xlab(bquote('Net rate of transitions '*my^-1*''))
p <- p + geom_vline(aes(xintercept=0), colour=cols[1], alpha=0.75)
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank(),
               axis.title.y=element_blank())
p
ggsave(filename = "figs/xy-zw-amph-v2.pdf")
