fossBarPlot <- function(df,cdf,taxon){
  require(dplyr)
  require(ggplot2)
  require(reshape)
  require(tidyr)
  source("Scripts/fossilAesthetics.R")
  df<-df[-1,] #drop "All" row
  cdfAll<-cdf[1,]
  cdf<-cdf[-1,] #drop "All" row
  taxon <- taxon
  
  df2 <- melt(df, id='X')
  colnames(df2) <- list("Group","Inclusion","Diversity")
  
  cdf2 <- melt(cdf, id = 'X') %>%
    mutate(variable = factor(variable, levels = unique(variable)))
  colnames(cdf2) <- list("Group","Inclusion","Completeness")
  
  sumdiv <- sum(df$Modern)#get the total # of extant species
  
  p <- ggplot(df2, aes(x = Inclusion, y = Diversity)) +
    geom_col(aes(fill = Group), width = 0.7, alpha = 0.6) + # bar plot
    
    scale_fill_manual(values = unlist(fossAES[taxon], use.names=FALSE)) + 
    
    geom_line(data = cdf2, aes(x = Inclusion, y = sumdiv*Completeness, #plot completeness for individual categories
                               group = Group, color = Group, linetype = Group),
              size = 1) +
    geom_point(data = cdf2, aes(x = Inclusion, y = sumdiv*Completeness,color = Group, pch = Group),#points for each category completeness
               size = 2) +
    
    scale_color_manual(values = unlist(fossAES[taxon], use.names=FALSE)) +
    
    scale_y_continuous(sec.axis=sec_axis(
      ~.*(1/sumdiv),name="Percent of species diversity preserved in fossil record", labels=scales::percent))+
    
    xlab(bquote('Criteria for inclusion')) +
    ylab(bquote('Number of species')) +
    theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=1, hjust=1))+
    scale_x_discrete(labels= c("Extant",
                               bquote('FGR >0'~km^2),
                               bquote('FGR >1'~km^2),
                               bquote('FGR >10'~km^2),
                               bquote('FGR >100'~km^2),
                               bquote('FGR >1000'~km^2),
                               bquote('FGR >10,000'~km^2),
                               bquote('FGR >100,000'~km^2)
    ))
  
  
  dir <- paste("Results/",taxon,"/",sep = "")
  saveRDS(p, file = paste(dir,taxon,"RecordBarPlot.R",sep = ""))
  ggsave(filename = paste(dir,taxon,"RecordBarPlot.png",sep = ""),
         plot = p,
         height = 7, width = 8)
}
