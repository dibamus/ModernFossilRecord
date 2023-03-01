
setwd("proper working directory")

library('tidyverse')
library('dplyr')

#### HEAD SCRIPT####
#This script is meant to be the primary thread for running all R analyses for the paper. These analyses require:
# #The Nyberg & Howell sedimentary basins spatial dataset
# #the GARD reptile range dataset
# #IUCN spatial data for birds, mammals, and amphibians
# #IUCN API key
# #VertLife Phylogenies for mammals, birds, reptiles, & amphibians

####Read in Datasets####

## Amphibian
ADF <- readRDS("FossilAmphibians.rds")
ADF <- mutate(ADF, difference = fossilArea - Area) %>% #Calculate difference between fossil area and area
  mutate(proportion = fossilArea/Area) #calculate proportion of range in basins
amgroupnames <- c(ANURA = "frog", CAUDATA = "salamander",GYMNOPHIONA = "caecilian")
ADF$group <- amgroupnames[ADF$order_]
rm(amgroupnames)

## Mammals
Mammals <- as.data.frame(readRDS("FossilMammals.rds"))%>% 
  mutate(difference = fossilArea - Area) %>% #Calculate difference between fossil area and area
  mutate(proportion = fossilArea/Area) #calculate proportion of range in basins 
Mammals$group <- Mammals$marine
Mammals$group[which(Mammals$group == "true")] <- "Marine"
Mammals$group[which(Mammals$group == "false")] <- "Terrestrial"

## Reptiles
GARD <- readRDS("FossilReptiles.rds")

GARD <- mutate(GARD,difference = fossilArea - Area) %>% #Calculate difference between fossil area and area
  mutate(proportion = fossilArea/Area) #calculate proportion of range in basins
names(GARD)[1] <- "binomial"
GARD$class <- "REPTILE"
GARD$group[which(GARD$group == "worm lizard")] <- "lizard" #include amphisbaenids in lizards

## Birds

BDF <- readRDS("FossilBirds.rds")

BDF <- mutate(BDF, difference = fossilArea - Area) %>% #Calculate difference between fossil area and area
  mutate(proportion = fossilArea/Area) #calculate proportion of range in basins
BDF$group <- BDF$passerine



####General Stats####

frStats <- function(df){
  setNames(c(mean(df$Area),mean(df$fossilArea),mean(df$difference),mean(df$proportion)),c("Mean GR","Mean FGR","Difference FGR-GR","FGR/GR"))
}

genStats <- data.frame(rbind("Amphibian" = frStats(ADF),
                             "Bird" = frStats(BDF),
                             "Mammal" = frStats(Mammals),
                             "Reptile" = frStats(GARD)))

#### FIGURE 2 ####

source("Scripts/densityPlots.R") #generates histogram ("hist") density plot ("density"), and range size vs proportion of range in fossil record ("sizevsproportion") plots


####Generate Cutoff/Dropoff Data - Who is in the fossil record? ####
source("Scripts/DropoffTables.R", local = TRUE)

GenerateDropoffs(ADF, "Amphibian")
survivingAmphibians <- which(!ADF$category %in% c("EN","CR","EW","EX","DD"))
GenerateDropoffs(ADF[survivingAmphibians,], "survivingAmphibian")

GenerateDropoffs(Mammals, "Mammal")

GenerateDropoffs(GARD, "Reptile")

GenerateDropoffs(BDF, "Bird")


source("Scripts/BarPlot.R", local = TRUE)

fossBarPlot(df = read.csv(file = "Results/Reptile/ReptileGroupInclusion.csv", header = T),
                cdf = read.csv(file = "Results/Reptile/ReptileGroupCompleteness.csv", header = T),
                taxon = "Reptile")

fossBarPlot(df = read.csv(file = "Results/Amphibian/AmphibianGroupInclusion.csv", header = T),
              cdf = read.csv(file = "Results/Amphibian/AmphibianGroupCompleteness.csv", header = T),
              taxon = "Amphibian")

fossBarPlot(df = read.csv(file = "Results/Mammal/MammalGroupInclusion.csv", header = T),
              cdf = read.csv(file = "Results/Mammal/MammalGroupCompleteness.csv", header = T),
              taxon = "Mammal")

fossBarPlot(df = read.csv(file = "Results/Bird/BirdGroupInclusion.csv", header = T),
            cdf = read.csv(file = "Results/Bird/BirdGroupCompleteness.csv", header = T),
            taxon = "Bird")

#### FIGURE 3####
library("ggpubr")
loadPlot <- function(taxon){#Quickly loads PD loss figures for multiplotting - removes labels
  x<- readRDS(file = paste0("Results/",taxon,"/",taxon,"RecordBarPlot.R"))
  total <- sum(filter(x$data, Inclusion == "Modern")$Diversity)
  as.data.frame(cbind(x$data, Completeness = x$layers[[3]]$data$Completeness, Total = total))%>%
    filter(Inclusion != "X.0km")
}

Barplotdf <- rbind(loadPlot("Reptile") %>% mutate(order = "Reptile"),
                   loadPlot("Amphibian") %>% mutate(order = "Amphibian"),
                   loadPlot("Mammal") %>% mutate(order = "Mammal"),
                   loadPlot("Bird") %>% mutate(order = "Bird")) %>% 
  
  mutate(scaledpercent = Completeness * Total)
Barplotdf[Barplotdf == "croc"] <- "crocodylian"
Barplotdf[Barplotdf == "Rhynchocephalia"] <- "tuatara"

Barplotdf$order <- factor(Barplotdf$order, levels = c("Bird","Reptile","Amphibian","Mammal"))

diversity_barplot <- ggplot(data = Barplotdf %>%   
                              filter(Group != "tuatara")%>% #exclude rhynchocephalians (1 species) 
                              filter(Group != "crocodylian"), #exclude crocs (bar too small to show up)
                            aes(x = Inclusion, y = Diversity)) +
  geom_col(aes(fill = Group), width = 0.7) + # bar plot
  
  scale_fill_manual(breaks = names(c(fossAES[[1]][c(2,4,5)],fossAES[[2]],fossAES[[4]],fossAES[[5]])), 
                    values = unlist(c(fossAES[[1]][c(2,4,5)],fossAES[[2]],fossAES[[4]],fossAES[[5]]))) + 
  xlab(bquote('Minimum FGR size for inclusion '(km^2))) +
  ylab('Number of species') +
  theme_light() +
  geom_vline(xintercept = 1.5, color = "darkgray", size = 1.5) +
  theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=1, hjust=1),
        panel.grid.major.x = element_blank())+
  scale_x_discrete(labels= c("Extant",
                             "1",
                             "10",
                             "100",
                             "1000",
                             "10,000",
                             "100,000")) +
  guides(color = "none", fill = "none") +
  facet_wrap(~order,strip.position = "left")


completeness_line <- ggplot(data = Barplotdf %>% filter(Inclusion != "Modern"), 
                            aes(x = Inclusion, y = Completeness, 
                                group = Group, color = Group, label = Group)) +
  geom_line(size = 1.3) +
  # geom_text(data = Barplotdf %>% filter(Inclusion == "X.100.000km"), 
  #           aes(x = Inclusion, y = Completeness, group = Group), 
  #           color = "black",
  #           hjust = 0,
  #           nudge_x = 0.05,
  #           nudge_y = 0.02) +
  
  scale_color_manual(breaks = names(c(fossAES[[1]],fossAES[[2]],fossAES[[4]],fossAES[[5]])), 
                     values = unlist(c(fossAES[[1]],fossAES[[2]],fossAES[[4]],fossAES[[5]]))) +
  scale_y_continuous(labels = scales::percent) +
  xlab(bquote('Minimum FGR size for inclusion '(km^2))) +
  ylab('Percent species diversity preserved') +
  theme_light() + 
  theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=1, hjust=1))+
  scale_x_discrete(labels= c("1",
                             "10",
                             "100",
                             "1000",
                             "10,000",
                             "100,000"),
                   expand = c(0, 0)) +
  guides(color = "none") +
  facet_wrap(~order,strip.position = "left")

Figure_3 <- ggarrange(diversity_barplot,completeness_line)

#### TABLE 1 ####
require('rlang')

includedtaxa <- function(df,x,taxon){
  length(levels(as.factor(df[which(df$fossilArea>x),taxon])))
}

getinclist <- function(df,taxon = "binomial"){ # get list of how many species/genera/families are included in the fossil record at certain cutoff values
  cutofflevels <- c(0,1,10,100,1000,10000,100000)
  l <- c(length(unique(df[,taxon])),lapply(cutofflevels, includedtaxa, df=df, taxon = taxon))
  names(l) <- c("Modern",">0km",">1km",">10km",">100km",">1,000km",">10,000km",">100,000km")
  return(l)
}

incTaxa <- rbind(Amphibian_species = getinclist(ADF,"binomial"), #removing "family" since it is not in the GARD data
                 Amphibian_Genera = getinclist(ADF,"genus"),
                 Amphibian_Families = getinclist(ADF,"family"),
                 
                 Surviving_Amphibian_Species = getinclist(ADF[survivingAmphibians,],"binomial"),
                 Surviving_Amphibian_Genera = getinclist(ADF[survivingAmphibians,],"genus"),
                 Surviving_Amphibian_Families = getinclist(ADF[survivingAmphibians,],"family"),
                 
           Reptile_species = getinclist(GARD,"binomial"),
           Reptile_Genera = getinclist(GARD,"genus"),
           Reptile_Families = getinclist(GARD,"family"),
           
           Mammal_species = getinclist(Mammals,"binomial"),
           Mammal_Genera = getinclist(Mammals,"genus"),
           Mammal_Families = getinclist(Mammals,"family"),
           
           Bird_species = getinclist(BDF,"binomial"),
           Bird_Genera = getinclist(BDF,"genus"),
           Bird_Families = getinclist(BDF,"family")
           )

fun <- function(x){
  unlist(x)/x[[1]]
}
incTaxaComp <- t(apply(incTaxa,1,fun))
save(incTaxaComp,file = "Results/Species-Genus_Completeness.rda")
write.csv(signif(incTaxaComp,3), file = "Results/TABLE_1_Species-Genus_Completeness.csv")
write.csv(incTaxa, file = "Results/Species-Genus_Diversity.csv")


####PD Loss####
library("phytools")
library("apTreeshape")

PDgen <- function(df,phy,taxon,rep){
  df <- df
  taxon <- taxon
  phy <- phy
  rep <- rep
  source("Scripts/PhyHistoryWithDepths.R")
  PDloss(df,phy,taxon,rep)
}


##Squamates
sq.phy <- read.tree("Phylogenies/Squamates/squam_shl_new_Posterior_9755.1000.trees") #9755 species
sq.phy <- drop.tip.multiPhylo(sq.phy, setdiff(sq.phy$tree_9916$tip.label,GARD$Tree_Name))# 9755 species

sqaL <- phylosig(sq.phy[[1]], x = setNames(GARD$fossilArea,GARD$Tree_Name)[sq.phy[[1]]$tip.label],
                 method = "lambda")

#Phylogenetic signal lambda : 0.421916 
#logL(lambda) : -135232 

PDgen(GARD,sq.phy,"Reptile",rep=1000)
rm(sq.phy)

##Amphibians 
amp.phy <- read.tree("Phylogenies/Amphibians/amph_shl_new_Posterior_7238.1000.trees") #7238 species
amp.phy<- drop.tip.multiPhylo(amp.phy, setdiff(amp.phy[[1]]$tip.label,ADF$Tree_Name))#6636 species

ampL <- phylosig(amp.phy[[1]], x = setNames(ADF$fossilArea,ADF$Tree_Name)[amp.phy[[1]]$tip.label],
                 method = "lambda")

#Phylogenetic signal lambda : 0.405178 
#logL(lambda) : -88534.3 

PDgen(ADF,amp.phy,"Amphibian",rep=1000)
rm(amp.phy)

##Mammals
mam.phy <- read.nexus("Phylogenies/Mammal/output.nex")
mam.phy <- drop.tip.multiPhylo(mam.phy,setdiff(mam.phy[[1]]$tip.label, Mammals$Tree_Name))#5419/5911 species

mamL <- phylosig(mam.phy[[1]], x = setNames(Mammals$fossilArea,Mammals$Tree_Name)[mam.phy[[1]]$tip.label],
                 method = "lambda")
#Phylogenetic signal lambda : 0.899551 
#logL(lambda) : -88877.8 

PDgen(Mammals,mam.phy,"Mammal",rep=1000)
rm(mam.phy)

##Birds
bird.phy <- read.nexus("Phylogenies/Bird/output.nex") #9993 species
bird.phy <- drop.tip.multiPhylo(bird.phy,setdiff(bird.phy[[1]]$tip.label,BDF$Tree_Name))#9391/9993 species

#let's do a bit of phylogenetic signal testing

passers <- drop.tip.multiPhylo(bird.phy,BDF$Tree_Name[which(BDF$group == 'Nonpasserine')])
nonpass <- drop.tip.multiPhylo(bird.phy,BDF$Tree_Name[which(BDF$group == 'Passerine')])

birdL <- phylosig(bird.phy[[1]], x = setNames(BDF$fossilArea,BDF$Tree_Name)[bird.phy[[1]]$tip.label],
                  method = "lambda")

passL <- phylosig(passers[[1]], x = setNames(BDF$fossilArea,BDF$Tree_Name)[passers[[1]]$tip.label],
         method = "lambda")
#Phylogenetic signal lambda : 0.140828 
#logL(lambda) : -88329.8 

npasL <- phylosig(nonpass[[1]], x = setNames(BDF$fossilArea,BDF$Tree_Name)[nonpass[[1]]$tip.label],
                  method = "lambda")
#Phylogenetic signal lambda : 0.374 
#logL(lambda) : -63495.4 

PDgen(BDF,bird.phy,"Bird", rep = 1000)
rm(bird.phy)

#### FIGURE 4 ####
source("Scripts/Phyplots.R")
#### FIGURE 5
#### Figure S3 - Node depth ####
source("Scripts/nodePlots.R")

## Get Tree Statistics #### NOT CRITICAL ####

treestats.DF <- rbind(treeStats(ADF, amp.phy) %>% mutate(group = "Amphibian"),
                      treeStats(BDF, bird.phy) %>% mutate(group = "Bird"),
                      treeStats(Mammals, mam.phy) %>% mutate(group = "Mammal"),
                      treeStats(GARD, sq.phy) %>% mutate(group = "Squamate"))

  

####Generate Range Size Scatterplots #### NOT CRITICAL ####
rangeScatterPlot <- function(df,taxon, rain = TRUE){
  df <- df
  taxon <- taxon
  #print(df[1,])
  source("Scripts/rangeSizeScatterPlot.R", local = TRUE)
}

rangeScatterPlot(ADF, "Amphibian")
rangeScatterPlot(GARD, "Reptile")
rangeScatterPlot(Mammals,"Mammal")
rangeScatterPlot(BDF,"Bird")

#### Generate Supplementary Trees with Mapped Range Size cutoffs #### NOT CRITICAL####
library('phytools')
library('phangorn')
library('tidytree')
library('treeio')
library('ggtree')

taxon <- "Amphibian" 
df <- ADF

maxAncestorMap <- function(taxon,df){
  source('Scripts/fossilAesthetics.R')
  source('Scripts/maxAncestor.R')
  
  if(taxon == "Squamates"){
    fossAES <- list("Squamates" = list(lizard = "#3dab98",snake = "#3dab98"))
  }
  phy <- read.tree(paste0("Phylogenies/",taxon,"/",taxon,"_consensus.tre"))
  
  names <- df$Tree_Name[which(df$Tree_Name %in%  phy$tip.label)]
  df <- filter(df, Tree_Name %in% names)
  df <- df[which(!duplicated(df$Tree_Name)),] #remove any duplicated names
  df$roundFossArea <- log10(df$fossilArea)
  df$roundFossArea[which(df$roundFossArea < 0)] <- 0
  phy <- ape::drop.tip(phy, setdiff(phy$tip.label,df$Tree_Name))

  nodes <- maxAncestor(phy,setNames(df$roundFossArea,df$Tree_Name))
  nodestib <- tibble(node = 1:(Ntip(phy) + Nnode(phy)), `log FGR` = c(nodes$tips,nodes$nodes))
  # Map <- contMap(phy, nodes$tips,
  #                     method = "user", anc.states = nodes$nodes,
  #                     plot = FALSE)
  # colorMap <- setMap(Map,colors=c("#eeeeee",fossAES[[taxon]][[2]]))
  
  #svg(paste0("/Results/PhyloMaps/",taxon,".svg"))
  
  # Code of the plot
  tidytree <- full_join(phy, nodestib,by = "node")
  
  treeplot <- ggtree(tidytree, aes(color = `log FGR`), layout = "circular") +
    scale_color_continuous(low = "#eeeeee", high = fossAES[[taxon]][[2]]) +
    geom_tiplab(color = "black", 
                #hjust = -.5, 
                #align = TRUE,
                offset = -0.6,
                size = 0.5) +
    theme(legend.position="right")
  
  return(treeplot)
  
  # Close the graphics device
 # dev.off() 
}

maxAncestorMap("Squamates", GARD)
maxAncestorMap("Bird", BDF)
ampMap <- maxAncestorMap("Amphibian", ADF)
maxAncestorMap("Mammal", Mammals)


