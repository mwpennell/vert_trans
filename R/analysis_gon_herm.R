source("R/util.R")
source("R/plotting.R")

## Gono to herm

gon.f <- readRDS("output/gon-herm/fish_sampdata_10_results.rds")

## Remove burnin from all samples
gon.f <- lapply(gon.f, function(x) x[-seq_len(10000),])
gon.f <- dplyr::bind_rows(gon.f)

gon.f$d <- gon.f$q01 - gon.f$q10

length(which(gon.f$d < 0))/nrow(gon.f)

gon.f$r <- gon.f$q01 / gon.f$q10
median(gon.f$r)

p <- ggplot(gon.f, aes(x=d))
p <- p + geom_histogram(bins = 50, fill=cols[2], alpha=1)
p <- p + ylab("Posterior density") + xlab("Net transition rate (Gonochorism to Hermaphroditism)")
p <- p + geom_vline(aes(xintercept=0), colour=cols[1], alpha=0.75)
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank())
p
ggsave(filename = "figs/gono-herm-fish.pdf")
