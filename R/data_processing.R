## data cleaning and processing
library(ape)
library(dplyr)
library(phyndr)

all_dat <- read.csv("data/vert_db.csv")
all_dat <- mutate(all_dat, binomial=paste(Genus, species, sep="_"))
all_dat <- mutate(all_dat, binomial=gsub(" [A-Z]", replacement="", binomial))
all_dat <- mutate(all_dat, binomial=gsub("[-]", replacement="", binomial))
all_dat <- mutate(all_dat, binomial=gsub("*([0-9])", replacement="", binomial))

fish_dat <- filter(all_dat, Higher.taxonomic.group=="Fish")
squa_dat <- filter(all_dat, Order=="Squamata")
amph_dat <- filter(all_dat, Higher.taxonomic.group=="Amphibia")

fish_tre <- read.tree("data/fish.tre")
squa_tre <- read.tree("data/squa.tre")
amph_tre <- read.tree("data/amph.tre")

## make ultrametric
fish_tre <- phangorn::nnls.tree(cophenetic(fish_tre),fish_tre,rooted=TRUE)
squa_tre <- phangorn::nnls.tree(cophenetic(squa_tre),squa_tre,rooted=TRUE)
amph_tre <- phangorn::nnls.tree(cophenetic(amph_tre),amph_tre,rooted=TRUE)

fish_spp <- unique(fish_dat$binomial)
squa_spp <- unique(squa_dat$binomial)
amph_spp <- unique(amph_dat$binomial)

fish_match <- phyndr::phyndr_genus(fish_tre, fish_spp)
squa_match <- phyndr::phyndr_genus(squa_tre, squa_spp)
amph_match <- phyndr::phyndr_genus(amph_tre, amph_spp)

fish_out <- list(phy=fish_match, dat=fish_dat)
squa_out <- list(phy=squa_match, dat=squa_dat)
amph_out <- list(phy=amph_match, dat=amph_dat)

saveRDS(fish_out, "output/fish_proc.rds")
saveRDS(squa_out, "output/squa_proc.rds")
saveRDS(amph_out, "output/amph_proc.rds")

