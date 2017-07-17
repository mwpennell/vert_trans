## Hermaphroditism and Gonochorism

build_gon_herm <- function(file_name, n_samp){
  
  td <- readRDS(file_name)
  d  <- td$dat
  phy <- td$phy
  clade <- strsplit(file_name, split="_")[[1]][1]
  clade <- strsplit(clade, split="/")[[1]][2]
  
  ## only consider gonochorous systems
  s <- grep("Sexual.system", colnames(d))[1]
  s.names <- c("gonochorous", "hermaphrodite")
  tmp <- which(d[,s] %in% s.names)
  dd <- d[tmp,]
  
  g.states <- factor(dd[,s], levels=s.names, labels=c(0,1))
  g.states <- as.numeric(levels(g.states)[as.integer(g.states)])
  
  ## build matrix with dummy variable for time being
  gon <- cbind.data.frame(g.states, dd$binomial)
  colnames(gon) <- c("gon.herm", "binomial")
  
  ## remove NAs
  gon <- na.omit(gon)
  
  out <- list(phy=phy, data=gon)
  
  saveRDS(out, paste0("output/gon-herm/", clade, "_fulldata", ".rds"))
  
  ## match to tree
  ptd <- phyndr_treedata(out, n_samp)
  ptd <- lapply(ptd, function(x) {
    dat <- as.numeric(x$data[,"gon.herm"])
    names(dat) <- rownames(x$data)
    list(phy = x$phy, dat=dat)
  })
  
  saveRDS(ptd, paste0("output/gon-herm/", clade, "_sampdata_", n_samp, ".rds"))
}

source("R/util.R")

## Build fish data
build_gon_herm("output/fish_proc.rds", 10)