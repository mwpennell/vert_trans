build_trap <- function(file_name, n_samp){
  
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
  k <- grep("Karyotype", colnames(d))[1]
  h.names <- c("XY", "ZW", "complex XY", "complex ZW")
  
  ## code blank as 0, everything else as 1
  tr.code <- function(x, i){
    if (x[i,e] != "" & x[i,k] == ""){
      y <- 1
    } else if (x[i,e] != "" & x[i,k] != ""){ ## drop these
      y <- NA
    } else if (x[i,e] == "" & x[i,k] == "homomorphic"){
      y <- 2
    } else if (x[i,e] == "" & x[i,k] %in% h.names){
      y <- 3
    } else {
      y <- NA
    }
  }
  
  tr.states <- sapply(seq_len(nrow(dd)), function(x) tr.code(dd,x))
  
  ## build matrix for analysis
  hg <- cbind(tr.states, dd$binomial)
  rownames(hg) <- dd$binomial
  colnames(hg) <- c("esd.sc", "binomial")
  
  ## add binomial for phyndr
  hg <- as.data.frame(hg)
  hg$binomial <- rownames(hg)
  
  ## remove NAs
  hg <- na.omit(hg)
  
  out <- list(phy=phy, data=hg)
  
  saveRDS(out, paste0("output/trap/", clade, "_fulldata", ".rds"))
  
  ## match to tree
  ptd <- phyndr_treedata(out, n_samp)
  ptd <- lapply(ptd, function(x) {
    dat <- as.numeric(x$data[,"esd.sc"])
    names(dat) <- rownames(x$data)
    list(phy = x$phy, dat=dat)
  })
  
  saveRDS(ptd, paste0("output/trap/", clade, "_sampdata_", n_samp, ".rds"))
}

source("R/util.R")

## Build fish data
build_trap("output/fish_proc.rds", 10)
## Build amphibian data
build_trap("output/amph_proc.rds", 10)
