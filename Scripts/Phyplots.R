#Phyplots - combined phylogenetic history loss plots
library("tidyverse")
source("Scripts/fossilAesthetics.R")
require("scales")
library("ggpubr")
library("dplyr")

loadRData <- function(fileName){
  #loads an RData file, and returns it
  readRDS(fileName)
  get(ls()[ls() != "fileName"])
}

mam <- readRDS("Results/Mammal/MammalFossilRecordEffect.R")
amp <- readRDS("Results/Amphibian/AmphibianFossilRecordEffect.R")
rep <- readRDS("Results/Reptile/ReptileFossilRecordEffect.R")
bir <- readRDS("Results/Bird/BirdFossilRecordEffect.R")

getEstimateD <- function(df){
  sapply(c(1:length(df)), FUN = function(x){df[[x]]$estimate})
}


CohensD <- rbind(getEstimateD(rep),
                 getEstimateD(mam),
                 getEstimateD(amp),
                 getEstimateD(bir)) %>%
  `rownames<-` (c("Reptile","Mammal","Amphibian","Bird")) %>%
  `colnames<-` (c("0km","1km","10km","100km","1000km","10,000km","100,000km")) %>%
  as.data.frame()
# Calculated the final Cohen's D (100,000km) wrong in the script. The bug is fixed now
#however, I'm not sure why there's no difference between the km0 simulated and observed datasets. 
#Something to look into


#### Plots


Mamplot <- readRDS("Results/Mammal/MammalPDloss.R")
Mamplot$data$Taxon <- "Mammal"
Mamplot$data$D <- getEstimateD(mam)

Repplot <- readRDS("Results/Reptile/ReptilePDloss.R")
Repplot$data$Taxon <- "Reptile"
Repplot$data$D <- getEstimateD(rep)

Ampplot <- readRDS("Results/Amphibian/AmphibianPDloss.R")
Ampplot$data$Taxon <- "Amphibian"
Ampplot$data$D <- getEstimateD(amp)

Birplot <- readRDS("Results/Bird/BirdPDloss.R")
Birplot$data$Taxon <- "Bird"
Birplot$data$D <- getEstimateD(bir)

allData <- rbind(Repplot$data,
      Mamplot$data,
      Ampplot$data,
      Birplot$data) %>% filter(cutoff >0 )
  
allData$diff <- allData$Tm - allData$Sm
allData$large <- allData$D %>% abs > 0.8


#### Combine all PD Loss Plots ####

cols <- scale_color_manual(values = c(Mammal = fossAES$Mammal$Terrestrial,
                                      Reptile = fossAES$Reptile$lizard,
                                      Amphibian = fossAES$Amphibian$salamander,
                                      Bird = fossAES$Bird$Nonpasserine))
fils <- scale_fill_manual(values = c(Mammal = fossAES$Mammal$Terrestrial,
                                     Reptile = fossAES$Reptile$lizard,
                                     Amphibian = fossAES$Amphibian$salamander,
                                     Bird = fossAES$Bird$Nonpasserine))



combinedPD <- ggplot(allData, aes(cutoff, Tm,Sm)) +
  geom_ribbon(aes(ymin=Sm - 2*Ssd, ymax=Sm + 2*Ssd, group = Taxon, fill= Taxon), alpha =  0.3) +
  geom_line(aes(cutoff, Sm,color="gray", group = Taxon), lwd=1, linetype = "dashed") +

  geom_ribbon(aes(ymin=Tm - 2*Tsd, ymax=Tm + 2*Tsd, fill = Taxon), alpha =  0.6) +
  geom_line(aes(cutoff, Tm, color = Taxon), lwd=1) +
  scale_y_continuous(labels = scales::percent, breaks = 0.2*c(0:5)) +
  scale_x_log10(breaks = c(1,10,100,1000,10000,100000),
                labels = c("1","10","100","1000","10,000","100,000")) +
  cols +
  fils +
  labs(color = "Mean SD + 95% CI", fill = "Mean SD + 95% CI",
       x = bquote('Minimum FGR size for inclusion '(km^2)),
       y = bquote('Proportion of PD represented')) +
  guides(color = "none", fill = "none")+
  theme_light() + 
  theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=1, hjust=1),
        panel.border = element_blank(),
        panel.grid.minor.x = element_blank())



#### Combined PD difference Plot ####

differencePD <- ggplot(allData, aes(cutoff, y = diff, color = Taxon)) +
  geom_line(size = 1.3)+
  geom_point(data = filter(allData, large == TRUE), pch = "*", size = 9) +
  scale_x_log10(breaks = c(1,10,100,1000,10000,100000), 
              labels = c("1","10","100","1000","10,000","100,000")) +
  cols+
  scale_y_continuous(labels = scales::percent) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1.3) +
  #ggtitle("Difference between Observed and Simulated PD") +
  labs(x = bquote('Minimum FGR size for inclusion '(km^2)),
       y = bquote('Difference of Means (% of total PD)')) +
  guides(color = "none", size = "none", fill = "none") +
  theme_light() + 
  theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=1, hjust=1),
        panel.border = element_blank(),
        panel.grid.minor.x = element_blank())
  
combined <- ggarrange(combinedPD,differencePD) + guides(color = "none",fill = "none")

ggsave("Figures/Figure_4.png", plot = combined,
       width = 7,
       height = 5)
