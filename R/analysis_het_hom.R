source("R/util.R")
source("R/plotting.R")

## Het-Hom

hh.f <- readRDS("output/het-hom/fish_sampdata_2_results.rds")
hh.a <- readRDS("output/het-hom/amph_sampdata_2_results.rds")

## Translate parameters
hh.f.mk <- lapply(hh.f, trans.multitrait)
hh.a.mk <- lapply(hh.a, trans.multitrait)

## Remove burnin from all samples, both translated and not
hh.f.mk <- lapply(hh.f.mk, function(x) x[-seq_len(2000),])
hh.f.mk <- dplyr::bind_rows(hh.f.mk)
hh.f <- lapply(hh.f, function(x) x[-seq_len(2000),])
hh.f <- dplyr::bind_rows(hh.f)

hh.a.mk <- lapply(hh.a.mk, function(x) x[-seq_len(2000),])
hh.a.mk <- dplyr::bind_rows(hh.a.mk)
hh.a <- lapply(hh.a, function(x) x[-seq_len(2000),])
hh.a <- dplyr::bind_rows(hh.a)


hh.f.mk$d <- hh.f.mk$q00.10 + hh.f.mk$q01.11 - hh.f.mk$q10.00 - hh.f.mk$q11.01

length(which(hh.f.mk$d < 0))/nrow(hh.f.mk)

p <- ggplot(hh.f.mk, aes(x=d))
p <- p + geom_histogram(bins = 50, fill=cols[3], alpha=1) 
p <- p + ylab("Posterior density") + xlab("Net rate of sex chromosome differentiation")
p <- p + geom_vline(aes(xintercept=0))
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank())
p 
ggsave(filename = "figs/het-hom-fish.pdf")

## Dependency on XY vs ZW systems
hh.f$g <- hh.f$qH01.G
p <- ggplot(hh.f, aes(x=g))
p <- p + geom_histogram(bins = 50, fill=cols[4], alpha=1)
p <- p + ggtitle("On differentiation")
p <- p + ylab("Posterior density") + xlab("Effect size of ZW vs. XY system")
p <- p + geom_vline(aes(xintercept=0))
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank())
p 
ggsave("figs/het-hom-fish-xy.pdf")

## Dependency on XY vs ZW systems
hh.f$g <- hh.f$qH10.G
p <- ggplot(hh.f, aes(x=g))
p <- p + geom_histogram(bins = 50, fill=cols[4], alpha=1)
p <- p + ggtitle("On turnover")
p <- p + ylab("Posterior density") + xlab("Effect size of ZW vs. XY system")
p <- p + geom_vline(aes(xintercept=0))
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank())
p 
ggsave("figs/hom-het-fish-xy.pdf")

hh.a.mk$d <- hh.a.mk$q00.10 + hh.a.mk$q01.11 - hh.a.mk$q10.00 - hh.a.mk$q11.01

length(which(hh.a.mk$d < 0))/nrow(hh.a.mk)

p <- ggplot(hh.a.mk, aes(x=d))
p <- p + geom_histogram(bins = 50, fill=cols[3], alpha=1) 
p <- p + ylab("Posterior density") + xlab("Net rate of sex chromosome differentiation")
p <- p + geom_vline(aes(xintercept=0))
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank())
p 
ggsave(filename = "figs/het-hom-amph.pdf")

## Dependency on XY vs ZW systems
hh.a$g <- hh.a$qH01.G
p <- ggplot(hh.a, aes(x=g))
p <- p + geom_histogram(bins = 50, fill=cols[4], alpha=1)
p <- p + ggtitle("On differentiation")
p <- p + ylab("Posterior density") + xlab("Effect size of ZW vs. XY system")
p <- p + geom_vline(aes(xintercept=0))
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank())
p 
ggsave("figs/het-hom-amph-xy.pdf")
length(which(hh.a$qH01.G > 0)) / nrow(hh.a)
## Dependency on XY vs ZW systems
hh.a$g <- hh.a$qH10.G
p <- ggplot(hh.a, aes(x=g))
p <- p + geom_histogram(bins = 50, fill=cols[4], alpha=1)
p <- p + ggtitle("On turnover")
p <- p + ylab("Posterior density") + xlab("Effect size of ZW vs. XY system")
p <- p + geom_vline(aes(xintercept=0))
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank())
p 
ggsave("figs/hom-het-amph-xy.pdf")
length(which(hh.a$qH10.G > 0)) / nrow(hh.a)


