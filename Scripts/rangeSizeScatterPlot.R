# Modern VS Fossil Ranges scatterplot

#Quick plot extant Range size vs Fossil Range Size
library("scales")
library("RColorBrewer")
library("ggplot2")
source("Scripts/fossilAesthetics.R")

p <- ggplot(data = df, aes(x = Area, y = fossilArea, colour = group , shape = group)) +
  { if (rain == TRUE)
    geom_segment(aes(x = Area, y = Area, 
                     xend = Area, yend = fossilArea, 
                     colour = group), 
                 alpha = 0.1)} +
  geom_point() +
  scale_color_manual(values = unlist(fossAES[taxon], use.names=FALSE)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(yintercept = 1,linetype = "longdash") +
  geom_hline(yintercept = 10,linetype = "twodash") +
  geom_hline(yintercept = 100,linetype = "dashed") +
  geom_hline(yintercept = 1000,linetype = "dotted") +
  coord_cartesian(xlim =c(0.5*min(df$Area),2*max(df$Area))) +
  ggtitle("Modern Range size vs Predicted Fossil Range size") +
  xlab(bquote('Modern Range Size'~(km^2))) +
  ylab(bquote('Predicted Fossil Range Size'~(km^2)))

#Save the plot
dir <- paste("Results/",taxon,"/",sep = "")

ggsave(filename = paste(dir,taxon,"RangeSizeScatterPlot.png",sep = ""),
       plot = p,
       height = 7, width = 10)

save(p, file = paste(dir,taxon,"RangeSizeScatterPlot.R",sep = ""))

print(paste("Scatter plot saved to ",dir,taxon,"RangeSizeScatterPlot.png",sep = ""))
print(paste("Scatter plot saved to ",dir,taxon,"RangeSizeScatterPlot.R",sep = ""))