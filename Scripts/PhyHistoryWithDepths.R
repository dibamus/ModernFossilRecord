#This function measures a multiphylo (phy) for a taxonomic group (taxon) to determine the distribution of phylogenetic diversity measures for the phylogenetic trees it represents. It then repeats that process for only those species whose $fossilArea measurement in (df) is >1 , >10, >100, and >1000, trimming trees to only include subsets that match that inclusion criterion.
# It then builds 100 simulated datasets by reshuffling all of the $fossilArea values between species and repeating the trimming and measurement process. It is currently set to simulate "rep" # of datasets
# Once the data have been generated, it uses a cohen's D test to determine whether the fossil record and the simulated fossil records represent different levels of phylogenetic diversity.
# Finally, it plots phylogenetic diversity from the real vs simulated datasets with a ribbon graph.
# The function returns a list of cohen's D values for the four different inclusion criteria.

PDloss <- function(df,phy,taxon,rep = 100){
  source("Scripts/Phyfunctions.R")
  source("Scripts/fossilAesthetics.R")
  source("Scripts/parallelSimulateDropoffs.R")
  require("effsize")
  require("ggplot2")
  require("scales")
  require("phytools")
  require("tidyverse")
  
  col <- fossAES[taxon][1]
  
  dir <- paste("Results/",taxon,"/",sep = "")
  
  #Measure Trees
  subsetLevels <- c(0,1,10,100,1000,10000,100000)
  names(subsetLevels) <- c("All","km1","km10","km100","km1000","km10,000","km100,000")
  
  subsetList <- lapply(subsetLevels, subsetFunction, df = df)
  
  #Measure inital lengthsof all trees (no spp dropped)
  intlengths <- tree.tl.multi(phy)
  #Measure length dropoff with real ranges (long step)
  truedata <- treedropoff(phy,subsetList)
  
  truelengths <- truedata$lengths
  saveRDS(truelengths, file = paste(dir,taxon,"HistoryObserved.rds",sep = ""))
  
  truedepths <- truedata$depths
  saveRDS(truedepths, file = paste(dir,taxon,"NodeDepthsObserved.rds",sep = ""))
  
  ####TREE LENGTH CODE####
  
  #Measure length dropoff with bootstrapped ranges (long step)
  simulatedlengths <- simulateDropoffs(phy,subsetLevels,rep,df) #generates rep simulated datasets
  
  saveRDS(simulatedlengths, file = paste(dir,taxon,"HistorySims.rds",sep = ""))
  
  for (i in 1:length(phy)) { # change history matrices to % completion (instead of length)
    truelengths[i,]<- truelengths[i,]/intlengths[[i]]
    simulatedlengths$lengths[i,,]<- simulatedlengths$lengths[i,,]/intlengths[[i]]
  }
  
  #Effect size of incomplete preservation
  
  effsize <- list(
    "0km"= cohen.d(truelengths[,1],simulatedlengths$lengths[,1,]),
    "1km" = cohen.d(truelengths[,2],simulatedlengths$lengths[,2,]), 
    "10km" = cohen.d(truelengths[,3],simulatedlengths$lengths[,3,]), 
    "100km" = cohen.d(truelengths[,4],simulatedlengths$lengths[,4,]), 
    "1000km" = cohen.d(truelengths[,5],simulatedlengths$lengths[,5,]),
    "10000km" = cohen.d(truelengths[,6],simulatedlengths$lengths[,6,]),
    "100000km" = cohen.d(truelengths[,7],simulatedlengths$lengths[,7,]))
  
  
  #Plot real vs random fossilization
  Tsd <- c(sd(truelengths[,1]),
           sd(truelengths[,2]),
           sd(truelengths[,3]),
           sd(truelengths[,4]),
           sd(truelengths[,5]),
           sd(truelengths[,6]),
           sd(truelengths[,7]))
  Tm <- c(mean(truelengths[,1]),
          mean(truelengths[,2]),
          mean(truelengths[,3]),
          mean(truelengths[,4]),
          mean(truelengths[,5]),
          mean(truelengths[,6]),
          mean(truelengths[,7]))
  Ssd <- c(sd(simulatedlengths$lengths[,1,]),
           sd(simulatedlengths$lengths[,2,]),
           sd(simulatedlengths$lengths[,3,]),
           sd(simulatedlengths$lengths[,4,]),
           sd(simulatedlengths$lengths[,5,]),
           sd(simulatedlengths$lengths[,6,]),
           sd(simulatedlengths$lengths[,7,]))
  Sm <- c(mean(simulatedlengths$lengths[,1,]),
          mean(simulatedlengths$lengths[,2,]),
          mean(simulatedlengths$lengths[,3,]),
          mean(simulatedlengths$lengths[,4,]),
          mean(simulatedlengths$lengths[,5,]),
          mean(simulatedlengths$lengths[,6,]),
          mean(simulatedlengths$lengths[,7,]))
  cutoff <- subsetLevels
  
  PD.df <- data.frame(Tm,Tsd,Sm,Ssd,cutoff)
  rownames(PD.df) <- names(subsetLevels)
  
  PDlossplot <- ggplot(PD.df, aes(cutoff, Tm,Sm))+
    geom_ribbon(aes(ymin=Sm - 2*Ssd, ymax=Sm + 2*Ssd, fill="simulated"), alpha =  0.5)+
    geom_line(aes(cutoff, Sm,color="simulated"), lwd=1, linetype = "dashed") +
    
    geom_ribbon(aes(ymin=Tm - 2*Tsd, ymax=Tm + 2*Tsd, fill="observed"), alpha =  0.5)+
    geom_line(aes(cutoff, Tm,color="observed"), lwd=1, linetype = "dashed") +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), #trans_breaks is retired - use log_breaks instead? or perhaps pretty_breaks?
                  labels = trans_format("log10", math_format(10^.x))) +
    
    ggtitle("PD loss from incomplete sampling of modern diversity") +
    scale_color_manual(values = c(simulated = "gray", observed =unlist(fossAES[taxon][1],use.names=FALSE)[2]))+
    scale_fill_manual(values = c(simulated = "gray", observed =unlist(fossAES[taxon][1],use.names=FALSE)[2]))+
    labs(color = "Mean SD + 95% CI", fill = "Mean SD + 95% CI", 
         x = bquote('Minimum fossil range size for inclusion '~(km^2)),
         y = bquote('Proportion of PD represented'))
  
  saveRDS(effsize, file = paste(dir,taxon,"FossilRecordEffect.R",sep = "")) 
  
  saveRDS(PDlossplot,file = paste(dir,taxon,"PDloss.",sep = ""))
  ggsave(filename = paste(dir,taxon,"PDloss.png",sep = ""),
        plot = PDlossplot,
        height = 7, width = 8)
  
  PDlossplot
  return(effsize)
}
