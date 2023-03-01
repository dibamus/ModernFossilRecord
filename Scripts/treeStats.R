#Get tree shape statistics for a subset of trees
treeStats <- function(df, phylo, ntrees = 10){
  rownames(df) <- df$Tree_Name
  for (i in 1:1000) {
    ab <- log10(df[phylo[[i]]$tip.label,"fossilArea"] + 0.1) + 1
    phylo[[i]]$tip.ab <- ab/max(ab)
  }
  
  trees <- phylo[sample(c(1:length(phylo)),ntrees)]
  
  sapply(1:10, FUN = function(x){maxlik.betasplit(trees[[x]])}) %>% 
    t %>% 
    data.frame
}
