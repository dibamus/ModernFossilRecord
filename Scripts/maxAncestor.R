maxAncestor <- function(phy, tip.data){
  require('phytools')
  if(Ntip(phy) != length(tip.data)){
    paste0("Number of tips (",Ntip(phy),
           ") does not equal number of tip data points (",
           length(tip.data),")")
    return()
  }
  tip.data <- tip.data[phy$tip.label]
  
  Node.data<- setNames(c(tip.data,rep(NA, times = phy$Nnode)),
                       c(phy$tip.label, 1:phy$Nnode+Ntip(phy)))
  
  for(i in 1:phy$Nnode+Ntip(phy)){
      Node.data[i] <- max(Node.data[getDescendants(phy,i)], na.rm = TRUE)
  }
  return(list(tips = Node.data[1:Ntip(phy)],
            nodes = Node.data[1:phy$Nnode+Ntip(phy)]))
}