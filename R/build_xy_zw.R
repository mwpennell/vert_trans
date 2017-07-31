## XY to ZW analyses 
## minor appendum to build_all_tree

## collate sd info for big tree figure
source("R/util.R")

build_xy <- function(file_name, n_samp){
  
  td <- readRDS(file_name)
  d  <- td$dat
  phy <- td$phy
  clade <- strsplit(file_name, split="_")[[1]][1]
  clade <- strsplit(clade, split="/")[[1]][2]
  
  ## gonochorous systems
  s <- grep("Sexual.system", colnames(d))[1]
  s.names <- c("gonochorous", "hermaphrodite")
  tmp <- which(d[,s] %in% s.names)
  dd <- d[tmp,]
  
  ## get column for environment
  e <- grep("Environment", colnames(d))[1]
  
  ## get column for karyotype
  k <- grep("Karyotype", colnames(d))[1]
  
  ## get column for heterogamety
  g <- grep("Genotypic", colnames(d))[1]
  
  ## exclude polygenic
  if ("polygenic" %in% dd[,g])
    dd <- dd[-which(dd[,g] == "polygenic"),]
  
  states <- rep(NA, nrow(dd))
  names(states) <- dd$binomial
  
  ## get indices for each
  
  ## state 1: male heterogametic -- heteromorphic
  mhet <- which(dd[,g] == "male heterogametic")
  het <- which(dd[,k] != "homomorphic")
  mhet.het <- intersect(mhet, het)
  states[mhet.het] <- rep(0, length(mhet.het))
  
  ## state 2: male heterogametic -- homomorphic
  hom <- which(dd[,k] == "homomorphic")
  mhet.hom <- intersect(mhet, hom)
  states[mhet.hom] <- rep(0, length(mhet.hom))
  
  ## state 3: female heterogametic -- heteromorphic
  fhet <- which(dd[,g] == "female heterogametic")
  fhet.het <- intersect(fhet, het)
  states[fhet.het] <- rep(1, length(fhet.het))
  
  ## state 4: female heterogametic -- homomorphic
  fhet.hom <- intersect(fhet, hom)
  states[fhet.hom] <- rep(1, length(fhet.hom))
  
  ## build matrix with dummy variable for time being
  sd <- cbind.data.frame(states, dd$binomial)
  colnames(sd) <- c("xy", "binomial")
  
  ## NA omit
  sd <- na.omit(sd)
  
  out <- list(phy=phy, data=sd)
  
  saveRDS(out, paste0("output/xy-zw/", clade, "_fulldata", ".rds"))
  
  ## match to tree
  ptd <- phyndr_treedata(out, n_samp)
  ptd <- lapply(ptd, function(x) {
    dat <- as.numeric(x$data[,"xy"])
    names(dat) <- rownames(x$data)
    list(phy = x$phy, dat=dat)
  })
  
  saveRDS(ptd, paste0("output/xy-zw/", clade, "_sampdata_", n_samp, ".rds"))
}

## write sd data for all 3 clades
build_xy("output/fish_proc.rds", n_samp=10)
build_xy("output/squa_proc.rds", n_samp=10)
build_xy("output/amph_proc.rds", n_samp=10)
