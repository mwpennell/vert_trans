source("R/util.R")
source("R/plotting.R")

## Trap

hh.f <- readRDS("output/trap/fish_sampdata_10_results.rds")


hh.f <- lapply(hh.f, function(x) x[-seq_len(10000),])
hh.f <- dplyr::bind_rows(hh.f)

## Dependency of GSD -> ESD on heteromorphy
hh.f$df <- hh.f$q31 - hh.f$q21

p <- ggplot(hh.f, aes(x=df))
p <- p + geom_histogram(bins = 50, fill=cols[2], alpha=1) + xlim(-0.15, 0.15)
p <- p + ylab("Posterior density") 
p <- p +xlab(bquote('Net rate of transitions '*my^-1*''))
p <- p + geom_vline(aes(xintercept=0), color=cols[1], alpha=0.75)
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank())
p 
ggsave("figs/trap-fish-v2.pdf")

p <- ggplot(hh.f, aes(x=q31))
p <- p + geom_histogram(bins = 500, fill=cols[2], alpha=1) + xlim(0, 0.05)
p <- p + ylab("Posterior density") 
p <- p +xlab(bquote('Transition rate (heteromorphy to ESD) '*my^-1*''))
p <- p + geom_vline(aes(xintercept=0), color=cols[1], alpha=0.75)
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank())
p
ggsave("figs/rev/figureS2-heteromorhphy-esd.pdf")

p <- ggplot(hh.f, aes(x=q21))
p <- p + geom_histogram(bins = 500, fill=cols[2], alpha=1) + xlim(0, 0.05)
p <- p + ylab("Posterior density") 
p <- p +xlab(bquote('Transition rate (homomorphy to ESD) '*my^-1*''))
p <- p + geom_vline(aes(xintercept=0), color=cols[1], alpha=0.75)
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank())
p
ggsave("figs/rev/figureS3-homomorphy-esd.pdf")
