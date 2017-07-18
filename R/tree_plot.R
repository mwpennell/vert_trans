## tree plotting

cols <- c("#053061", "#78B7C5", "firebrick3", "#EE6363", 
          "#0B775E", "#D8B70A", "#999999")
names(cols) <- c("mhet", "mhom", "fhet", "fhom", "esd", "herm", "lg")
cols <- list(S=as.character(cols[1:7]))

## get taxonomy
all_dat <- read.csv("data/vert_db.csv")
tax <- all_dat[,1:5]
tax <- tax[!duplicated(tax$Genus),]

f <- readRDS("output/sd-all/fish_sampdata_10.rds")
f <- f[[1]]
phy.f <- ladderize(f$phy)
genus.f <- sub("_.+$", "", phy.f$tip.label)
fam.f <- as.character(sapply(genus.f, function(x)
{as.character(tax[which(tax$Genus == x), "Family"])}))
dat.f <- as.data.frame(f$dat)

diversitree::trait.plot(phy.f, dat=dat.f, cols=cols, 
                        class=fam.f, w=1/20, legend=FALSE, 
                        cex.lab=0.55, margin=1/3)

s <- readRDS("output/sd-all/squa_sampdata_10.rds")
s <- s[[1]]
phy.s <- ladderize(s$phy)
genus.s <- sub("_.+$", "", phy.s$tip.label)
fam.s <- as.character(sapply(genus.s, function(x)
{as.character(tax[which(tax$Genus == x), "Family"])}))
dat.s <- as.data.frame(s$dat)

diversitree::trait.plot(phy.s, dat=dat.s, cols=cols, 
                        class=fam.s, w=1/20, legend=FALSE, 
                        cex.lab=0.55, margin=1/3)

a <- readRDS("output/sd-all/amph_sampdata_10.rds")
a <- a[[1]]
phy.a <- ladderize(a$phy)
genus.a <- sub("_.+$", "", phy.a$tip.label)
fam.a <- as.character(sapply(genus.a, function(x)
{as.character(tax[which(tax$Genus == x), "Family"])}))
dat.a <- as.data.frame(a$dat)

diversitree::trait.plot(phy.a, dat=dat.a, cols=cols, 
                        class=fam.a, w=1/20, legend=FALSE, 
                        cex.lab=0.55, margin=1/3)

dev.off()
pdf(file="figs/tree-plot.pdf",height=4, width=12)
par(mfrow=c(1,3))
trait.plot(phy.f, dat=dat.f, cols=cols, class=fam.f, w=1/20, legend=FALSE, cex.lab=0.5, margin=1/2.5)
trait.plot(phy.s, dat=dat.s, cols=cols, class=fam.s, w=1/20, legend=FALSE, cex.lab=0.5, margin=1/2.5)
trait.plot(phy.a, dat=dat.a, cols=cols, class=fam.a, w=1/20, legend=FALSE, cex.lab=0.5, margin=1/2.5)
dev.off()



