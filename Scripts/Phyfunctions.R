###08-24-22 MODIFIED FUNCTIONS TO TRACK CHANGES IN NODE HEIGHT DISTRIBUTION
drop.tip.multiPhylo<-function(phy, tip){
  if(!inherits(phy,"multiPhylo"))
    stop("phy is not an object of class \"multiPhylo\".")
  else {
    trees<-lapply(phy,drop.tip,tip=tip)
    class(trees)<-"multiPhylo"
  }
  trees
}

tree.tl <- function(phy){ #measure total tree length
  sum(phy$edge.length)
}

tree.tl.multi <- function(multiphy) {# measure all tree lengths in multiphylo
  lapply(multiphy,tree.tl)
}

internal.nd <- function(phy) { #measure the depth of all internal nodes (tips have depth of 0) (Deprecated for hist version)
  nd <- node.depth.edgelength(phy)[-c(1:length(phy$tip.label))]
  nd <- max(nd) - nd
  nd
}

internal.nd.hist <- function(phy) { #get histogram counts of internal node depths (tips have depth of 0)
  nd <- node.depth.edgelength(phy)[-c(1:length(phy$tip.label))]
  nd <- max(nd) - nd
  nd[nd < 0.001]<- 0.0011 #replace small values for log range safety
  hist(nd %>% log10() , breaks = seq(-3,3, by = 0.25), plot = FALSE)$counts
}

###CHANGED FOR NODE DEPTH MEASUREMENT (see return line)
tree.tl.multi.sub <- function(multiphy,sub) {#measure total history loss when subsampling trees
  #multiphy = multiphylo object
  #sub = list of tips to drop from tree
  multiphy <- drop.tip.multiPhylo(multiphy, sub) #subsample trees
  tl <- tree.tl.multi(multiphy)
  nd <- t(sapply(multiphy, internal.nd.hist))  #needs to be transposed
  
  return(list(tl = tl, nd = nd))
}

###CHANGED FOR NODE DEPTH MEASUREMENT
#take a multiphylo and a list of subsets, measure tree length for each subsetting operation
treedropoff <- function(multiphy,subsetList) {
  require('parallel')
  m <- matrix(nrow = length(multiphy),
              ncol = length(subsetList))
  colnames(m) <- names(subsetList)
  
  #create a t X s X n matrix of node depths t = number of trees, s = # of subsets, n = 24
  # when using the new form, tree.tlmulti.sub returns counts of nodes in bins, so n = 24 bins
  # this means the array can be A LOT smaller - like not 20+ gb (as it was when each node depth was stored)
  d <- array(data = NA, #matrix of node depths t X s X n
             dim = c(length(multiphy),
                     length(subsetList),
                     24),
             dimnames = list(NULL,
                             names(subsetList),
                             paste0("10^",seq(-3,2.75,by = 0.25),"_myr")))
  
  for (i in 1:length(subsetList)){
    print(i)
    f <- tree.tl.multi.sub(multiphy, subsetList[[i]])
    m[,i] <- unlist(f$tl)
    #for the depths matrix, not the entire t X n slice is filled; 
    #only t-length(s[i]) X n data points are generated because length(s[i]) # of tips were removed

    d[c(1:length(multiphy)),i,c(1:dim(f$nd)[2])] <- f$nd
  }
  return(list(lengths = m, # returns a ntrees*nsubsets matrix of tree lengths
              depths = d))
}

#Finds how many species in a data frame have a fossil area larger than x
subsetFunction <-function(x,df){df$Tree_Name[which(df$fossilArea < x)]}

#
force.ultrametric.multi <- function(multiPhy){
  require('phytools')
  for (i in 1:length(multiPhy)){
    multiPhy[[i]] <- force.ultrametric(multiPhy[[i]], method = "extend")
  }
  multiPhy
}
