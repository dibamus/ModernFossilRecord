### Calculating the extinction rate stuff

## E/MSY = extinctions per million species-years. background rates are estimated 
## from fossil extinctions that took place in million-year-or-more time bins. 
## For current rates, the proportion of species extinct in a comparatively very
## short time (one to a few centuries) is extrapolated to predict what the rate 
## would be over a million years. 

## Average fossil rate (about E/MSY≈1.8) for mammals
## MSY = 10,000 species per 100 years
## Multiply number of species by the years considered


## In our study, this is the proportion of species that would appear extinct 
## based on geographic range size over the Anthropocene 

## Extinction rates are very very high, but of course they are!! it's based on 
## geographic range size... IF we captured everything, fossil extinction rates
## SHOULD be much higher. The average for mammals is 1.8 
## The extinction rates in the past are based on the fossil record, so we 
## expect them to be small -- 


library(MASS) 
library(reshape2) 
library(reshape) 
library(tidyr)

ExRates <- function(){
  Species_Genus_Diversity <- read.csv("Results/Species-Genus_Diversity.csv")
  
  ## Cycle through a range of years 
  years = data.frame(1, 10, 100, 1000, 10000, 100000, 1000000)
  ##### Do based on IUCN data ####
  
  ## Because we are already doing a subset of the data here (i.e., how many species
  ## are most likely to go extinct), then I don't need to subtract from the modern. 
  ## I just use the raw numbers I get from the following calculations. 
  ## This will mean that the extinction rates DECREASE, bc the number of species does also.
  
  ## Amphibians ##
  
  amph_IUCN = Species_Genus_Diversity[1,2:9] - Species_Genus_Diversity[4,2:9]
  
  mod_amph_IUCN = amph_IUCN[1,1]
  mod_amph = Species_Genus_Diversity[1,2]
  
  ## Reptiles ##
  reptiles_IUCN =  GARD %>% filter(REDcat == "EN" | REDcat == "CR" | REDcat == "DD"
                                   | REDcat == "EX" | REDcat == "EW") 
  mod_rep = Species_Genus_Diversity[7,2]
  rep_IUCN = data.frame(NA, NA, NA, NA, NA, NA, NA, NA)
  colnames(rep_IUCN) <- c('Modern', '>0km', '1km', '>10km', '>100km', '>1000km', '>10000km', '>1000000km')
  mod_rep_IUCN = nrow(reptiles_IUCN)
  
  rep_IUCN[1,1] = nrow(reptiles_IUCN)
  rep_IUCN[1,2] = nrow(reptiles_IUCN %>% filter(fossilArea > 0))
  rep_IUCN[1,3] = nrow(reptiles_IUCN %>% filter(fossilArea > 1))
  rep_IUCN[1,4] = nrow(reptiles_IUCN %>% filter(fossilArea > 10))
  rep_IUCN[1,5] = nrow(reptiles_IUCN %>% filter(fossilArea > 100))
  rep_IUCN[1,6] = nrow(reptiles_IUCN %>% filter(fossilArea > 1000))
  rep_IUCN[1,7] = nrow(reptiles_IUCN %>% filter(fossilArea > 10000))
  rep_IUCN[1,8] = nrow(reptiles_IUCN %>% filter(fossilArea > 100000))
  
  ## Mammals ##
  mammals_IUCN =  Mammals %>% filter(category == "EN" | category == "CR" | category == "DD"
                                     | category == "EX" | category == "EW") 
  mod_mam = Species_Genus_Diversity[10,2]
  mam_IUCN = data.frame(NA, NA, NA, NA, NA, NA, NA, NA)
  colnames(mam_IUCN) <- c('Modern', '>0km', '1km', '>10km', '>100km', '>1000km', '>10000km', '>1000000km')
  mod_mam_IUCN = nrow(mammals_IUCN)
  
  mam_IUCN[1,1] = nrow(mammals_IUCN)
  mam_IUCN[1,2] = nrow(mammals_IUCN %>% filter(fossilArea > 0))
  mam_IUCN[1,3] = nrow(mammals_IUCN %>% filter(fossilArea > 1))
  mam_IUCN[1,4] = nrow(mammals_IUCN %>% filter(fossilArea > 10))
  mam_IUCN[1,5] = nrow(mammals_IUCN %>% filter(fossilArea > 100))
  mam_IUCN[1,6] = nrow(mammals_IUCN %>% filter(fossilArea > 1000))
  mam_IUCN[1,7] = nrow(mammals_IUCN %>% filter(fossilArea > 10000))
  mam_IUCN[1,8] = nrow(mammals_IUCN %>% filter(fossilArea > 100000))
  
  ## Birds ## 
  
  birds_IUCN =  BDF %>% filter(category == "EN" | category == "CR" | category == "DD"
                               | category == "EX" | category == "EW") 
  mod_bir = Species_Genus_Diversity[13,2]
  bir_IUCN = data.frame(NA, NA, NA, NA, NA, NA, NA, NA)
  colnames(bir_IUCN) <- c('Modern', '>0km', '1km', '>10km', '>100km', '>1000km', '>10000km', '>1000000km')
  mod_bir_IUCN = nrow(birds_IUCN)
  
  bir_IUCN[1,1] = nrow(birds_IUCN)
  bir_IUCN[1,2] = nrow(birds_IUCN %>% filter(fossilArea > 0))
  bir_IUCN[1,3] = nrow(birds_IUCN %>% filter(fossilArea > 1))
  bir_IUCN[1,4] = nrow(birds_IUCN %>% filter(fossilArea > 10))
  bir_IUCN[1,5] = nrow(birds_IUCN %>% filter(fossilArea > 100))
  bir_IUCN[1,6] = nrow(birds_IUCN %>% filter(fossilArea > 1000))
  bir_IUCN[1,7] = nrow(birds_IUCN %>% filter(fossilArea > 10000))
  bir_IUCN[1,8] = nrow(birds_IUCN %>% filter(fossilArea > 100000))
  
  
  ## For the first set of calculations, we are assuming that the difference 
  ## between the FGR bins are the species going extinct, 
  ## So the calculation should be (difference)/(modern*years/1mill)
  
  ## For the IUCN data, the number of species in each bin is what is going extinct
  ## So the calculation should be (# of species)/(modern*years/1mill)
  
  ### Calculate extinction rates ##
  years = data.frame(1, 10, 100, 1000, 10000, 100000, 1000000)
  
  
  for (i in 1:ncol(amph_IUCN)) {
    for (j in 2:ncol(years)){
      amph_IUCN[j,i] = (amph_IUCN[1, i])/(mod_amph*years[,j]/1000000)
    }
  }
  rownames(amph_IUCN) = c("1", "10", "100", "1000", "10000", "100000", "1000000")
  
  rep_IUCN_species_year = rep_IUCN[1,1]
  for (i in 1:ncol(rep_IUCN)) {
    for (j in 2:ncol(years)){
      rep_IUCN[j,i] = (rep_IUCN[1, i])/(mod_rep*years[,j]/1000000)
    }
  }
  
  mam_IUCN_species_year = mam_IUCN[1,1]
  for (i in 1:ncol(mam_IUCN)) {
    for (j in 2:ncol(years)){
      mam_IUCN[j,i] = (mam_IUCN[1, i])/(mod_mam*years[,j]/1000000)
    }
  }
  
  bir_IUCN_species_year = bir_IUCN[1,1]
  for (i in 1:ncol(bir_IUCN)) {
    for (j in 2:ncol(years)){
      bir_IUCN[j,i] = (bir_IUCN[1, i])/(mod_bir*years[,j]/1000000)
    }
  }
  
  
  ## Create data frame
  
  year = 3 ## This is 100 years based on the for loops
  
  Extinction_rates = data.frame(NA, NA, NA, NA, NA, NA, NA)
  colnames(Extinction_rates) <- c('modern','1', '10', '100', '1000', '10,000', '100,000')
  Extinction_rates[1, ] = amph_IUCN[year,c(1,3:8)]
  Extinction_rates[2, ] = rep_IUCN[year,c(1,3:8)]
  Extinction_rates[3, ] = mam_IUCN[year,c(1,3:8)]
  Extinction_rates[4, ] = bir_IUCN[year,c(1,3:8)]
  
  rownames(Extinction_rates) <- c('Amphibian','Reptile','Mammal','Bird')
  
  #Plot Extinction Rates
  source('Scripts/fossilAesthetics.R')
  
  Extinction_rates$rowname <- row.names(Extinction_rates)
  m <-melt(Extinction_rates, id.var = "rowname")
  m$variable <- as.factor(m$variable)
  m$rowname <- as.factor(m$rowname)
  
  cols <- scale_color_manual(values = c(Mammal = fossAES$Mammal$Terrestrial,
                                        Reptile = fossAES$Reptile$lizard,
                                        Amphibian = fossAES$Amphibian$salamander,
                                        Bird = fossAES$Bird$Nonpasserine))
  fils <- scale_fill_manual(values = c(Mammal = fossAES$Mammal$Terrestrial,
                                       Reptile = fossAES$Reptile$lizard,
                                       Amphibian = fossAES$Amphibian$salamander,
                                       Bird = fossAES$Bird$Nonpasserine))
  
  ExRateGraph <- ggplot(m, aes(group = rowname, color = rowname))+
    geom_line(data = m %>% filter(variable !='modern'), 
              aes(x = variable , y = value), 
              linewidth = 1.3)+
    geom_hline(data = m%>% filter(variable == 'modern'),
               aes(yintercept = value, color = rowname), 
               linetype = "dashed",
               linewidth = 1)+
    cols+
    labs(color = "Taxon",
         x = bquote('Minimum FGR size for inclusion '(km^2)),
         y = bquote('Extinction Rate (E/MSY)')) +
    guides(color = "none", fill = "none")+
    theme_light()+
    theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=1, hjust=1))
  return(list(graph = ExRateGraph,
              birdER = bir_IUCN,
              reptileER =rep_IUCN,
              mammalER = mam_IUCN,
              amphibianER = amph_IUCN))
}
Extinction_rates <- ExRates()
saveRDS(Extinction_rates, file = 'results/ExtinctionRates.rds')

saveRDS(ExRateGraph, file = 'Figures/ExtinctionRates.rds')
