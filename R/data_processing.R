## data cleaning and processing
library(ape)
library(dplyr)
library(phyndr)

all_dat <- read.csv("data/vert_db.csv")
all_dat <- mutate(all_dat, binomial=paste(Genus, species, sep="_"))
fish_dat <- filter(all_dat, Higher.taxonomic.group=="Fish")
squa_dat <- filter(all_dat, Order=="Squamata")
amph_dat <- filter(all_dat, Higher.taxonomic.group=="Amphibia")
rownames(fish_dat) <- fish_dat$binomial
rownames(squa_dat) <- squa_dat$binomial
rownames(amph_dat) <- amph_dat$binomial

fish_tre <- read.tree("data/fish.tre")
squa_tre <- read.tree("data/squa.tre")
amph_tre <- read.tree("data/amph.tre")

## make ultrametric
fish_tre <- phangorn::nnls.tree(cophenetic(fish_tre),fish_tre,rooted=TRUE)
squa_tre <- phangorn::nnls.tree(cophenetic(squa_tre),squa_tre,rooted=TRUE)
amph_tre <- phangorn::nnls.tree(cophenetic(amph_tre),amph_tre,rooted=TRUE)

fish_match <- phyndr::phyndr_genus(fish_tre, rownames(fish_dat))
squa_match <- phyndr::phyndr_genus(squa_tre, rownames(squa_dat))
amph_match <- phyndr::phyndr_genus(amph_tre, rownames(amph_dat))

saveRDS(fish_match, "output/fish_tree_proc.rds")
saveRDS(squa_match, "output/squa_tree_proc.rds")
saveRDS(amph_match, "output/amph_tree_proc.rds")

