## Build data for ESD to GSD

build_esd_gsd <- function(file_name, n_samp){

  td <- readRDS(file_name)
  d  <- td$dat
  phy <- td$phy
  clade <- strsplit(file_name, split="_")[[1]][1]
  clade <- strsplit(clade, split="/")[[1]][2]
  
  ## only consider gonochorous systems
  s <- grep("Sexual.system", colnames(d))[1]
  s.names <- c("gonochorous")
  tmp <- which(d[,s] %in% s.names)
  dd <- d[tmp,]
  
  e <- grep("Environmental", colnames(d))[1]
  
  ## code blank as 0, everything else as 1
  e.code <- function(x){
    if (x == ""){
      y <- 0
    } else {
      y <- 1
    }
  }
  
  e.states <- sapply(dd[,e], function(x) e.code(x))
  
  ## build matrix with dummy variable for time being
  esd <- cbind.data.frame(e.states, dd$binomial)
  colnames(esd) <- c("esd.gsd", "binomial")
  
  ## remove NAs
  esd <- na.omit(esd)
  
  out <- list(phy=phy, data=esd)
  
  saveRDS(out, paste0("output/esd-gsd/", clade, "_fulldata", ".rds"))
  
  ## match to tree
  ptd <- phyndr_treedata(out, n_samp)
  ptd <- lapply(ptd, function(x) {
    dat <- as.numeric(x$data[,"esd.gsd"])
    names(dat) <- rownames(x$data)
    list(phy = x$phy, dat=dat)
  })

  saveRDS(ptd, paste0("output/esd-gsd/", clade, "_sampdata_", n_samp, ".rds"))
}


## Build fish data
build_esd_gsd("output/fish_proc.rds", 10)
## Build squamate data
build_esd_gsd("output/squa_proc.rds", 10)
