#revised version of Figure 4, showing absolute extiction magnitude and 
#proportional extinction of four tetrapod groups.

require("ggplot2")
require("tidyr")
require("reshape2")
require("dplyr")
require("stringr")
require("cowplot")
source("Scripts/fossilAesthetics.R")

difference <- read.csv(file = "Results/Species-Genus_Diversity.csv", 
                       header = T)

pAc <- which(grepl("Surviving",difference$X,fixed=T))
Ac <- which(!grepl("Surviving",difference$X,fixed=T))

prop <- 1-difference[pAc,-1]/difference[Ac,-1] %>% as.data.frame

prop$level <- rep(c("species","Genus","Family"),4)
prop$taxon <- str_replace(difference[Ac,1],"_species|_Genera|_Families","")
prop$group <- paste0(prop$taxon,prop$level)
prop$taxon <- factor(prop$taxon, levels = c("Bird","Reptile","Amphibian","Mammal"), ordered = TRUE)
prop <- prop[,-2] #remove the 0km level

#Table_3 - extinction magnitude

T3 <- data.frame(group = prop$group) %>% cbind(Diversity=difference$Modern[-pAc],prop[,1:7])

write.csv(T3, file = "Results/Table_3.csv")

# Figure 5

df <- melt(prop, id.vars = c("level","taxon","group"))

map <- setNames(c("Extant","1","10","100","1000","10,000","100,000"),unique(df$variable))
df$variable <- factor(map[df$variable], levels = map, ordered = TRUE)

exMag <- ggplot(data = df, aes(x = variable, y = value, group = group, linetype = level, color = taxon)) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = c(Mammal = fossAES$Mammal$Terrestrial,
                                Reptile = fossAES$Reptile$lizard,
                                Amphibian = fossAES$Amphibian$salamander,
                                Bird = fossAES$Bird$Nonpasserine))+
  scale_linetype_manual(values = c(species ="solid",Genus ="dashed", Family ="dotted"))+
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L)) +
  #coord_cartesian(xlim=c(1,7), expand = FALSE)+
  ylab("Extinction Magnitude")+
  xlab(bquote('Minimum FGR size for inclusion '(km^2)))+
  facet_wrap(taxon ~., )+
  theme_light()+
  theme(axis.text.x=element_text(color = "black", size=9, angle=45, vjust=1, hjust=1),
        panel.grid.minor.y = element_line(linetype = "dashed"),
        panel.grid.minor.x = element_blank(),
        strip.text = element_blank() , 
        strip.background = element_blank())+
    guides(linetype = "none",color = "none")


ExRateGraph <- readRDS("Figures/ExtinctionRates.rds")

Fig5 <- plot_grid(ExRateGraph + theme_light()+
                    xlab(bquote('Minimum FGR size for inclusion '(km^2)))+
                    theme(axis.text.x=element_text(color = "black", size=9, angle=45, vjust=1, hjust=1),
                          panel.grid.minor.y = element_line(linetype = "dashed"),
                          panel.grid.minor.x = element_blank()),
                  exMag,
                  
                  labels = c("A","B"),
                  axis = 'tb',
                  align = 'v',
                  rel_widths = c(3,5),
                  ncol =2)

ggsave(filename = "Figures/Figure_5_2col.pdf",
       plot = Fig5,
       height = 8, width = 15, units = "cm")
