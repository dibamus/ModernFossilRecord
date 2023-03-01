require(ggplot2)
source('Scripts/fossilAesthetics.R')

####Datasets for plotting ####

#Modern and fossil range sizes
densitydata <- rbind(rbind(select(Mammals,binomial,Area,class) %>% mutate(foss = "Extant"),
                     select(GARD,binomial,Area,class) %>% mutate(foss = "Extant"),
                     select(BDF,binomial,Area,class) %>% mutate(foss = "Extant"),
                     select(ADF,binomial,Area,class) %>% mutate(foss = "Extant")),
                     
                     rbind(select(Mammals,binomial,fossilArea,class) %>% mutate(foss = "Fossil"),
                           select(GARD,binomial,fossilArea,class) %>% mutate(foss = "Fossil"),
                           select(BDF,binomial,fossilArea,class) %>% mutate(foss = "Fossil"),
                           select(ADF,binomial,fossilArea,class) %>% mutate(foss = "Fossil")) %>%
                       `colnames<-`(c("binomial","Area","class","foss")))

densitydata$class <- recode(densitydata$class, #recoding this so the legend looks nice
                            "AVES" = "Bird",
                            "MAMMALIA" = "Mammal",
                            "REPTILE" = "Reptile",
                            "AMPHIBIA" = "Amphibian",) %>%
  factor(levels = c("Amphibian","Bird", "Mammal", "Reptile")) # stating the order that they should plot in

#this plots on a log scale, so replace small Areas with 0.1 to keep it tidy
densitydata$Area[which(densitydata$Area < 0.1)] <- 0.1


#### Density plots of Modern and fossil range sizes ####
density <- ggplot(data = densitydata, aes(group = class, fill = class, color = class, linetype = class))+
  geom_density(aes(x = Area), alpha = 0.5, size = 1.3) +
  
  xlab(bquote("Range Size "(km^2))) +
  ylab("Percent Species Diversity") +
  
  scale_fill_manual(values = c(Amphibian = fossAES$Amphibian$salamander,
                               Bird = fossAES$Bird$Nonpasserine,
                               Mammal = fossAES$Mammal$Terrestrial,
                               Reptile = fossAES$Reptile$lizard)) +
  scale_color_manual(values = c(Amphibian = fossAES$Amphibian$salamander,
                                Bird = fossAES$Bird$Nonpasserine,
                                Mammal = fossAES$Mammal$Terrestrial,
                                Reptile = fossAES$Reptile$lizard)) +
  
  scale_y_continuous(labels = function(x) scales::percent(abs(x))) +
  scale_x_log10(breaks = c(0.1,1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000),
                labels = c("0","1","10","100","1000","10,000","100,000","1,000,000","10,000,000","100,000,000","1,000,000,000"),
                expand = c(0,0,0,0)) +
  theme_light()+
  theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=1, hjust=1),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.9, 0.9),
        legend.title = element_blank()) +
  facet_wrap(~foss, nrow = 2)

#### save the file ####
ggsave(filename = paste("Figures/Figure_2.png",sep = ""),
       plot = density,
       height = 8, width = 6)