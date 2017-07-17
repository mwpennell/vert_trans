source("R/util.R")
source("R/plotting.R")

## ESD to GSD

gon.f <- readRDS("output/gon-herm/fish_sampdata_10_results.rds")

## Remove burnin from all samples
gon.f <- lapply(gon.f, function(x) x[-seq_len(10000),])
gon.f <- dplyr::bind_rows(gon.f)

gon.f$d <- gon.f$q01 - gon.f$q10
profile.block(gon.f$d, col.line=cols["dg"], col.fill=cols["dg"],
             col.l=cols["herm"], col.g=cols["mhet"], cex.lab=1.25,
             cex.axis=1.25,
             xlab="Net transition rate from gonochorism to hermaphroditism",
             main="Fish")

length(which(gon.f$d < 0))/nrow(gon.f)

gon.f$r <- gon.f$q01 / gon.f$q10
median(gon.f$r)

dev.off()
pdf("figs/gono-herm-draft.pdf")
profile.block(gon.f$d, col.line=cols["dg"], col.fill=cols["dg"],
              col.l=cols["herm"], col.g=cols["mhet"], cex.lab=1.25,
              cex.axis=1.25,
              xlab="Net transition rate from gonochorism to hermaphroditism",
              main="Fish")
dev.off()
