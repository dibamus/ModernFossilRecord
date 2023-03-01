###CHANGED FOR NODE DEPTH MEASUREMENT 
simulateDropoffs <- function(multiphy,subsetLevels,replicates,df){
  source("Scripts/Phyfunctions.R")
  require('parallel')
  require('doParallel')
  
  #generate 3d array (ntrees*subsetLevels*nreplicates)
  a <- array(dim = c(length(multiphy),
                     length(subsetLevels),
                     replicates),
             dimnames = list(names(multiphy),
                             names(subsetLevels),
                             paste0("r_",as.character(1:replicates))))
  #generate 4d array (ntrees*subsetLevels*maxnodes*replicates)
  d <- array(data = NA,
             dim = c(length(multiphy),
                     length(subsetLevels),#thisis the problematic line - where is "subsetList" coming from
                     24, #getting counts of nodes within bins now instead of measuring the depth of each node
                     replicates),
             dimnames = list(names(multiphy),
                             names(subsetLevels),
                             paste0("10^",seq(-3,2.75,by = 0.25),"_myr"),
                             paste0("r_",as.character(1:replicates))))
  
  cl <- makeCluster(parallel::detectCores() - 1)
  clusterExport(cl,
                varlist  = c("multiphy", "subsetLevels","replicates","df","a", "d"),
                envir = environment())
  clusterEvalQ(cl,library('tidyverse'))
  clusterEvalQ(cl,library('dplyr'))
  clusterEvalQ(cl,library('phytools'))
  clusterEvalQ(cl,library('ape'))
  clusterEvalQ(cl,library('phytools'))
  clusterEvalQ(cl,library('foreach'))
  clusterEvalQ(cl,library('doParallel'))
  clusterEvalQ(cl,source("Scripts/Phyfunctions.R"))
  
  registerDoParallel(cl)
  
  foreach(i = 1:replicates) %dopar% {
    s <- sample(df$fossilArea)
    f <- function(x,s){df$Tree_Name[which(s < x)]}
    subsetList <- lapply(subsetLevels, f, s)
    drop.tip.multiPhylo(multiphy, subsetList$km100)[[1]]$edge.length %>% sum()
  }
  
  results <- foreach(i = 1:replicates) %dopar% {
    s <- sample(df$fossilArea)
    f <- function(x,s){df$Tree_Name[which(s < x)]}
    subsetList <- lapply(subsetLevels, f, s)

    return(list(a = as.matrix(treedropoff(multiphy,subsetList)$lengths),
    d = treedropoff(multiphy,subsetList)$depths))
  }
  
  for (i in 1:replicates){
    a[,,i] <- results[[i]]$a
    d[,,,i]<- results[[i]]$d
  }  
  
  stopImplicitCluster()
  return(list(lengths = a, #return filled array of tree lengths
              depths = d)) #return filled array of node depth counts
  
}