build_het_hom <- function(file_name, n_samp){
  
  td <- readRDS(file_name)
  d  <- td$dat
  phy <- td$phy
  clade <- strsplit(file_name, split="_")[[1]][1]
  clade <- strsplit(clade, split="/")[[1]][2]
  
  ## pull out speices for which we have info on karyotype
  k <- grep("Karyotype", colnames(d))[1]
  h.names <- c("homomorphic", "XY", "ZW", "complex XY", "complex ZW")
  tmp <- which(d[,k] %in% h.names)
  dd <- d[tmp,]
  
  h.code <- function(x){
    if (x == "homomorphic"){
      y <- 0
    } else { y <- 1}
    y
  }
  
  h.states <- sapply(dd[,k], function(x) h.code(x))
  
  ## find genotypic states
  g <- grep("Genotypic", colnames(dd))[1]
  g.states <- factor(dd[,g])
  g.names <- c("male heterogametic", "female heterogametic")
  g.states <- factor(g.states, levels=g.names, labels=c(0,1))
  g.states <- as.numeric(levels(g.states)[as.integer(g.states)])
  
  ## build matrix for analysis
  hg <- cbind(h.states, g.states)
  rownames(hg) <- dd$binomial
  colnames(hg) <- c("H", "G")
  
  ## add binomial for phyndr
  hg <- as.data.frame(hg)
  hg$binomial <- rownames(hg)
  
  out <- list(phy=phy, data=hg)
  
  saveRDS(out, paste0("output/het-hom/", clade, "_fulldata", ".rds"))
  
  ## match to tree
  ptd <- phyndr_treedata(out, n_samp)
  
  ## get rid of binomial
  ptd <- lapply(ptd, function(x) {
    dat <- as.data.frame(x$data[,c("H", "G")])
    dat$G <- unfactor(dat$G)
    dat$H <- unfactor(dat$H)
    list(phy = x$phy, dat=dat)
  })
  
  saveRDS(ptd, paste0("output/het-hom/", clade, "_sampdata_", n_samp, ".rds"))
}

source("R/util.R")

## Build fish data
build_het_hom("output/fish_proc.rds", 10)
## Build amphibian data
build_het_hom("output/amph_proc.rds", 10)
