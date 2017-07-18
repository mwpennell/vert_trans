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

esd.f$d <- esd.f$q01 - esd.f$q10
profile.block(esd.f$d, col.line=cols["dg"], col.fill=cols["dg"],
              col.l=cols["esd"], col.g=cols["mhet"], cex.lab=1.25,
              cex.axis=1.25, xlim=c(-0.1,0.0025),
              xlab="Net transition rate from GSD to ESD", main="Fish")

length(which(esd.f$d < 0))/nrow(esd.f)

esd.f$m <- esd.f$q10 / esd.f$q01
median(esd.f$m)

esd.s$d <- esd.s$q01 - esd.s$q10
profile.block(esd.s$d, col.line=cols["dg"], col.fill=cols["dg"],
              col.l=cols["esd"], col.g=cols["mhet"], cex.lab=1.25,
              cex.axis=1.25, xlab="Net transition rate from GSD to ESD",
              xlim=c(-0.029,0.0025), main="Squamates")

length(which(esd.s$d < 0))/nrow(esd.s)

esd.s$m <- esd.s$q10 / esd.s$q01
median(esd.s$m)

pdf("figs/gsd-esd-draft.pdf")
par(mfrow=c(1,2))
profile.block(esd.f$d, col.line=cols["dg"], col.fill=cols["dg"],
              col.l=cols["esd"], col.g=cols["mhet"], 
              xlim=c(-0.1,0.0025),
              xlab="Net transition rate from GSD to ESD", main="Fish")
profile.block(esd.s$d, col.line=cols["dg"], col.fill=cols["dg"],
              col.l=cols["esd"], col.g=cols["mhet"], 
              xlab="Net transition rate from GSD to ESD",
              xlim=c(-0.029,0.0025), main="Squamates")
dev.off()

profile.block(esd.f$d, col.line=cols[1], col.fill=cols[2],
              col.l="white", col.g="white",  cex.lab=1,
              cex.axis=1, xlim=c(-0.1,0.0025), vline=0,
              xlab="Net transition rate from GSD to ESD", main="Fish")


p <- ggplot(esd.f, aes(x=d))
p <- p + geom_histogram(bins = 50, fill=cols[5], alpha=1) + xlim(-0.15,0.005)
p <- p + ylab("Posterior density") + xlab("Net transition rate (GSD to ESD)")
p <- p + geom_vline(aes(xintercept=0))
p <- p + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank())
p
ggsave(filename = "figs/esd-gsd-fish.pdf")
q <- ggplot(esd.s, aes(x=d))
q <- q + geom_histogram(bins = 50, fill=cols[5], alpha=1) + xlim(-0.029,0.0025)
q <- q + xlab("Net transition rate (GSD to ESD)")
q <- q + geom_vline(aes(xintercept=0))
q <- q + theme(panel.background=element_blank(), 
               axis.ticks.y = element_blank(),
               axis.text.y=element_blank(),
               axis.title.y=element_blank())
q
ggsave("figs/esd-gsd-squa.pdf")

