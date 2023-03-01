require("ggplot2")
require("tidyr")
require("reshape2")
require("dplyr")
require("cowplot")
source("Scripts/fossilAesthetics.R")

taxon = "Amphibian"

#Load in Modern amphibian data
df <- read.csv(file = "Results/Amphibian/AmphibianGroupInclusion.csv", header = T)[2:4,]
df <- melt(df, id='X')
colnames(df) <- list("Group","Inclusion","Diversity")
df$Dataset <- "Ac"
df$Alpha = 1

#Load in Post-extinction amphibian data

xdf <- read.csv(file = "Results/SurvivingAmphibian/SurvivingAmphibianGroupInclusion.csv", header = T)[2:4,]
xdf <- melt(xdf, id='X')
colnames(xdf) <- list("Group","Inclusion","Diversity")
xdf$Dataset = "pAc"
xdf$Alpha = 0.5


combined <- rbind(df, xdf) %>% as.data.frame() 
map <- setNames(c("Extant","0","1","10","100","1000","10,000","100,000"),levels(combined$Inclusion))
combined$Inclusion <- factor(map[combined$Inclusion], levels = c("Extant","0","1","10","100","1000","10,000","100,000"))
combined <- filter(combined, Inclusion != "0")


p <- ggplot() +
  geom_col(data = combined,
           aes(x=Dataset, y=Diversity, fill = Group, alpha = Alpha), position = "stack") +
  
  xlab(bquote('Minimum FGR size for inclusion '(km^2))) +
  ylab(bquote('Number of species'))+
  
  theme_light() +
  theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=1, hjust=1),
        
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(0, "mm"),                       # remove spacing between facets
        strip.background = element_rect(size = 0.5),
        strip.text.x = element_text(size = 12,face = c("bold")),
        plot.margin = unit(c(0,0,0,0), "mm")) +
  scale_fill_manual(values = fossAES$Amphibian) +
  scale_alpha(range=c(0.5,1)) +
  facet_grid(~Inclusion, switch = "both") +
  guides(alpha = "none", fill = "none")

difference <- read.csv(file = "Results/Species-Genus_Diversity.csv", 
                       header = T)[1:6,] %>% 
  data.frame(row.names = 1) %>% t  %>% data.frame %>%
  mutate(species = Surviving_Amphibian_Species/Amphibian_species,
         genus = Surviving_Amphibian_Genera/Amphibian_Genera)

difference$Inclusion <- factor(row.names(difference), levels = row.names(difference))
difference <- melt(difference[,7:9], id.vars = c("Inclusion"))  %>% 
  filter(Inclusion != "X.0km")

 
exMag <- ggplot(data = difference, aes(x = Inclusion, y = 1 - value, group = variable, linetype = variable)) +
  geom_line(size = 1.3, color = fossAES$Amphibian[[2]]) +
  geom_point(size = 10, color = fossAES$Amphibian[[2]], pch = "*") +
  coord_cartesian(xlim = c(0.5,7.5),ylim = c(-0.02,0.45), expand = FALSE)+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5), color = "gray",size = 0.5)+
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L)) +
  geom_label(data = difference %>% filter(Inclusion == "Modern"),
            aes(x = 1, y = 1-value, label = variable),
            label.size = 0,
            hjust = 0,
            nudge_x = 0.1,
            nudge_y = 0.03) +
  ylab("Extinction Magnitude")+
  theme_light() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_line(color = gray, size = 0.5),
        plot.margin = unit(c(0.5,0.5,0,0.5), "mm")) +
  
  guides(linetype = "none")


combinedFig <- plot_grid(exMag,p, nrow = 2, align = "v", axis = "lrbt",rel_heights = c(0.5,1))

ggsave(filename = "Figures/Figure_5.png",
       plot = combinedFig,
       height = 6, width = 8)
