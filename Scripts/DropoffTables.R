# Dropoff/Exclusion Functions
GenerateDropoffs <- function(df,taxon){
  require('ggplot2')
  
  groupcol <- df$group
  groups <- levels(as.factor(groupcol))
  
  #getInclusionsTable - creates a table of numbers of species above a certain range size
  getInclusionsTable <- function(df,col="fossilArea",focal) {
    z <- length(focal)
    y <- length(which(df[focal,col] > 0))
    a <- length(which(df[focal,col] > 1)) # how many ranges above 1 km2?
    b <- length(which(df[focal,col] > 10)) # how many ranges above 10 km2?
    c <- length(which(df[focal,col] > 100)) # how many ranges above 100 km2?
    d <- length(which(df[focal,col] > 1000))# how many ranges above 1000 km2?
    e <- length(which(df[focal,col] > 10000)) # how many ranges above 10000 km2?
    f <- length(which(df[focal,col] > 100000)) # how many ranges above 100000 km2?
    c(z,y,a,b,c,d,e,f)
  }
  
  #functions for easy group charting
  inc <- function(group) {
    getInclusionsTable(df,"fossilArea",which(groupcol == group))
  }
  percinc <- function(group) {
    getInclusionsTable(df,"fossilArea",which(groupcol == group))/length(which(groupcol == group))
  }
  
  getInclusions <- function(df, col="fossilArea", focal, x) {
    length(which(df[focal,col] > x))
  }
  
  #### Dropoff charts ####
  
  fun <- function(x, group, perc = T) { # count how many ranges in the group are larger than x
    v <-getInclusions(df,"fossilArea",which(groupcol == group),x)
    if (perc) { #get proportion of total
      v <- v/length(which(groupcol == group))
    }
    v
  }
  
  # Sets up a vector of values and calls fun for each one to determine the ranges above that value in the group given
  dropoffinc <- function(group) { 
    group <- group
    samples <- c(1:9997)*10
    samples <- c(0,1,5,samples)
    sapply(samples,fun, group)
  }
  
  # Table of species diversity within groups for given cutoff values
  groupInclusion <- t(sapply(groups,inc))
  groupInclusion <- rbind(getInclusionsTable(df,focal = c(1:length(df$group))),groupInclusion)
  colnames(groupInclusion) <- c("Modern",">0km",">1km",">10km",">100km",">1000km",">10,000km",">100,000km")
  rownames(groupInclusion)[1] <-"All"
  
  #groupInclusionGen <- rbind(groupInclusion,All_Genera)
  
  
  # Table of record completeness within groups for given cutoff values
  groupCompleteness <- t(sapply(groups,percinc))
  groupCompleteness <- rbind(getInclusionsTable(df,focal = c(1:length(df$group)))/length(df$group),groupCompleteness)
  colnames(groupCompleteness) <- c("Modern",">0km",">1km",">10km",">100km",">1000km",">10,000km",">100,000km")
  rownames(groupCompleteness)[1] <-"All"
  
  # Create matrix for plotting dropoff in recovered completion with increasing cutoff value
  dropoff <- t(sapply(groups,dropoffinc))
  colnames(dropoff) <- as.character(c(0,1,5,c(1:9997)*10))
  
  #####Plot Dropoff ####
  library("reshape2")
  source("Scripts/fossilAesthetics.R")
  fossAES$survivingBird <- fossAES$Bird
  fossAES$survivingMammal <- fossAES$Mammal
  fossAES$survivingReptile <- fossAES$Reptile
  
  melted <- melt(dropoff, id.vars=colnames())
  colnames(melted)<- c("Group","Size","Percent")
  
  d <- ggplot(data=melted, aes(x=Size, y=Percent, group=Group, colour = Group, linetype = Group)) + 
    geom_line(size = 1) +
    scale_x_log10() +
    scale_color_manual(values = unlist(fossAES[taxon], use.names=FALSE)) +
    ylab("Percent Completeness") +
    xlab(bquote('Minimum range size for Inclusion'~(km^2))) +
    ggtitle("Effect of inclusion criteria on completness")
  
  #### Save Results ####
  dir <- paste("Results/",taxon,"/",sep = "")
  
  #Save Dropoff Plot
  png(paste(dir,taxon,"Dropoff.png",sep = ""))
  plot(d)
  dev.off()
  save(d, file = paste(dir,taxon,"Dropoffplot.R",sep = ""))
  print(paste("Dropoff plot saved to ",dir,taxon,"Dropoff.png",sep = ""))
  print(paste("Dropoff plot saved to ",dir,taxon,"Dropoff.R",sep = ""))
  
  #Save tables
  write.csv(groupInclusion,file = paste(dir,taxon,"GroupInclusion.csv",sep=""))
  print(paste("Group inclusion table saved to ",dir,taxon,"GroupInclusion.csv",sep = ""))
  
  write.csv(groupCompleteness,file = paste(dir,taxon,"GroupCompleteness.csv",sep=""))
  print(paste("Group completeness table saved to ",dir,taxon,"GroupCompleteness.csv",sep = ""))
}