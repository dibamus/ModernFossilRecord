#Cleaning mammal Data
setwd("C:/Users/Isaac/Box/Fossil Record")
library("sf")
library("rgdal")
library("tidyverse")
library("rgeos")
EDF <-read.csv("Geo_Intersects/Mammal_Range_Metadata.csv") #Extant Data Frame
EDF <- EDF[match(unique(EDF$binomial),EDF$binomial),]
EDF <- EDF[order(EDF$binomial),]
EDF$fossilArea <- 0

originalRanges <- read.csv("Geo_Intersects/MAMMALS_original_ranges_1.csv", header = T)
originalRanges <- originalRanges[order(originalRanges$binomial),]

MDF <- read.csv("Geo_Intersects/Mammals_Basins_Overlap.csv", header = T)
colnames(MDF)[4:5]<-c("Area","fossilArea")
MDF<-MDF[order(MDF$binomial),]


fossilized <- which(EDF$binomial%in%MDF$binomial)
EDF$fossilArea[fossilized] <- MDF$fossilArea

EDF$Area <- originalRanges$orig_area_km2

levels(as.factor(EDF$category)) # "CR" "DD" "EN" "EW" "EX" "LC" "NT" "VU"

IUCNlevels <- c("LC","NT","VU","EN","CR","EW","EX","DD")#order these from least to most extinct

EDF$category <- ordered(EDF$category, IUCNlevels)

EDF$Tree_Name <- sub(" ","_",EDF$binomial)

#get the synonyms
source("Scripts/getClearSynonyms.R")
library("phytools")
phy <- read.nexus("Phylogenies/Mammal/output.nex")
MammalTaxonomy <- get_clear_synonyms(EDF,phy)

saveRDS(MammalTaxonomy,file = "Results/Mammal/MammalTaxonomy.rda")

EDF$Tree_Name <- MammalTaxonomy$namematches$replacement

#fix species in MammalTaxonomy$unclear
EDF$Tree_Name[which(EDF$binomial == "Cheracebus purinus")] <- "Callicebus_purinus" 
EDF$Tree_Name[which(EDF$binomial == "Cheracebus torquatus")] <- "Callicebus_torquatus"
EDF$Tree_Name[which(EDF$binomial == "Hylopetes sagitta")] <- "Petinomys_sagitta" 
EDF$Tree_Name[which(EDF$binomial == "Hypsugo anthonyi")] <- "Pipistrellus_anthonyi" 
EDF$Tree_Name[which(EDF$binomial == "Hypsugo joffrei")] <- "Pipistrellus_joffrei" 
EDF$Tree_Name[which(EDF$binomial == "Lycalopex griseus")] <- "Pseudalopex_griseus" 
EDF$Tree_Name[which(EDF$binomial == "Lycalopex gymnocercus")] <- "Pseudalopex_gymnocercus" 
EDF$Tree_Name[which(EDF$binomial == "Plecturocebus caligatus")] <- "Callicebus_caligatus" 
EDF$Tree_Name[which(EDF$binomial == "Plecturocebus dubius")] <- "Callicebus_dubius"  

EDF <- EDF %>% mutate(inTree = sapply(Tree_Name, function(x){x %in% phy[[1]]$tip.label }), .after = Tree_Name)

#helper functions to look for species and genera in the tree
chphy <- function(x){ x %in% phy[[1]]$tip.label}
chgen <- function(x){phy[[1]]$tip.label[str_which(phy[[1]]$tip.label,x)]}

#look for potentially significant generic changes
commongenera <- table(filter(EDF, inTree == FALSE)$genus)

#Crocidura
chgen("Crocidura") %>% sort()
filter(EDF, EDF$genus == "Crocidura", inTree == FALSE)$binomial
EDF$Tree_Name[which(EDF$binomial == "Crocidura absconditus")] <- "Crocidura_abscondita" 

#Cryptotis
chgen("Cryptotis") %>% sort()
filter(EDF, EDF$genus == "Cryptotis", inTree == FALSE)$binomial
EDF$Tree_Name[which(EDF$binomial == "Cryptotis merus")] <- "Cryptotis_mera" 

#Dermanura
chgen("Dermanura") %>% sort()
filter(EDF, EDF$genus == "Dermanura", inTree == FALSE)$binomial
EDF$Tree_Name[which(EDF$binomial == "Dermanura azteca")] <- "Dermanura_aztecus" 
EDF$Tree_Name[which(EDF$binomial == "Dermanura cinerea")] <- "Dermanura_cinereus" 
EDF$Tree_Name[which(EDF$binomial == "Dermanura glauca")] <- "Dermanura_glaucus" 
EDF$Tree_Name[which(EDF$binomial == "Dermanura gnoma")] <- "Dermanura_gnomus" 
EDF$Tree_Name[which(EDF$binomial == "Dermanura rosenbergi")] <- "Dermanura_rosenbergii" 
EDF$Tree_Name[which(EDF$binomial == "Dermanura tolteca")] <- "Dermanura_toltecus"

#Heteromys
chgen("Heteromys") %>% sort()
chgen("Liomys") %>% sort()
filter(EDF, EDF$genus == "Heteromys", inTree == FALSE)$binomial
EDF[which(EDF$binomial %in% filter(EDF, EDF$genus == "Heteromys", inTree == FALSE)$binomial),]$Tree_Name <- chgen("Liomys") %>% sort()

# Hipposideros
chgen("Hipposideros") %>% sort()
chgen("Liomys") %>% sort()
chgen("Pipistrellus") %>% sort()
filter(EDF, EDF$genus == "Hipposideros", inTree == FALSE)$binomial
EDF$Tree_Name[which(EDF$binomial == "Hypsugo alaschanicus")] <- "Pipistrellus_alaschanicus"
EDF$Tree_Name[which(EDF$binomial == "Hypsugo anthonyi")] <- "Pipistrellus_anthonyi"
EDF$Tree_Name[which(EDF$binomial == "Hypsugo joffrei")] <- "Pipistrellus_joffrei"
EDF$Tree_Name[which(EDF$binomial == "Hypsugo savii")] <- "Pipistrellus_savii"
EDF$Tree_Name[which(EDF$binomial == "Hypsugo vordermanni")] <- "Pipistrellus_vordermanni"

#Change "Neotamias" Tree Names to "Tamias"
chgen("Neotamias")# not in tree
chgen("Tamias")
view(filter(EDF, genus == "Neotamias"))
EDF[str_which(EDF$Tree_Name, "Neotamias"),]$Tree_Name <- EDF[str_which(EDF$Tree_Name, "Neotamias"),]$Tree_Name %>% 
  str_replace("Neotamias_","Tamias_")

#Piliocolobus
chgen("Piliocolobus")
chgen("Procolobus") %>% sort()
chgen("Colobus") %>% sort()
chgen("Tropicolobus") %>% sort()
filter(EDF, EDF$genus == "Piliocolobus", inTree == FALSE)$binomial
EDF$Tree_Name[which(EDF$binomial == "Piliocolobus kirkii")] <- "Procolobus_kirkii" 
EDF$Tree_Name[which(EDF$binomial == "Piliocolobus kirkii")] <- "Procolobus_kirkii"

#Plecturocebus
chgen("Plecturocebus") %>% sort()
chgen("Callicebus") %>% sort()
chgen("Callithrix") %>% sort()
filter(EDF, EDF$genus == "Plecturocebus", inTree == FALSE)$binomial
EDF$Tree_Name[which(EDF$binomial == "Plecturocebus caligatus")] <- "Callicebus_caligatus"
EDF$Tree_Name[which(EDF$binomial == "Plecturocebus dubius")] <- "Callicebus_dubius"

#Tupaia
chgen("Tupaia") %>% sort()
chgen("Urogale") %>% sort()
filter(EDF, EDF$genus == "Tupaia", inTree == FALSE)$binomial
EDF$Tree_Name[which(EDF$binomial == "Tupaia everetti")] <- "Urogale_everetti"


EDF <- EDF %>% mutate(inTree = sapply(Tree_Name, function(x){x %in% phy[[1]]$tip.label }), .after = Tree_Name)

saveRDS(EDF, file = "FossilMammals.rds")
write.csv(EDF, file = "FossilMammals.csv")
