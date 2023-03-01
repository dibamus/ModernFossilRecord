#Amphibian Data cleaning
setwd("C:/Users/Isaac/Box/Fossil Record")

originalRanges <- read.csv("Results/Amphibian/AMPHIBIANS_Dissolve_TableToExcel.csv", header = T)

fossilRanges <- read.csv("Results/Amphibian/AMPHIBIANS_PairwiseIntersect_TableToExcel.csv", header = T)

taxonomy <- read.csv("Results/Amphibian/AMPHIBIANS_TableToExcel.csv", header = T)

#set up area and fossil area columns
colnames(originalRanges)[6] <- "Area"
originalRanges$fossilArea <- 0

#add fossilArea to originalranges df
originalRanges[match(fossilRanges$binomial,originalRanges$binomial),]$fossilArea <- fossilRanges$final_area_km2

#just get taxonomic data for species (trim subspecies, subpopulations, etc)
taxonomy <- taxonomy[match(unique(taxonomy$binomial),taxonomy$binomial),]

#get the dfs in the same order
originalRanges <- originalRanges[order(originalRanges$binomial),]
taxonomy <- taxonomy[order(taxonomy$binomial),]

ADF <- data.frame(binomial = originalRanges$binomial,
                  kingdom = 'ANIMALIA',
                  phylum = 'CHORDATA',
                  class = 'AMPHIBIA',
                  order_ = taxonomy$order_,
                  family = taxonomy$family,
                  genus = gsub(" .+","",originalRanges$binomial),
                  category = taxonomy$category,
                  Area = originalRanges$Area,
                  fossilArea = originalRanges$fossilArea,
                  Tree_Name = sub(" ","_",originalRanges$binomial))

IUCNlevels <- c("LC","NT","VU","EN","CR","EW","EX","DD")#order these from least to most extinct

ADF$category <- ordered(ADF$category, IUCNlevels)

#get the synonyms
source("Scripts/getClearSynonyms.R")
library("phytools")
phy <- read.tree("Phylogenies/Amphibians/amph_shl_new_Posterior_7238.1000.trees")
AmphibianTaxonomy <- get_clear_synonyms(ADF,phy)
saveRDS(AmphibianTaxonomy,file = "Results/Amphibian/AmphibianTaxonomy.rda")

ADF$Tree_Name <- AmphibianTaxonomy$namematches$replacement

View(AmphibianTaxonomy$unclear)
# match "unclear" taxa
#The Hylaranas
ADF$Tree_Name[which(ADF$binomial == "Amnirana albolabris")] <- "Hylarana_albolabris"
ADF$Tree_Name[which(ADF$binomial == "Amnirana longipes")] <- "Hylarana_longipes"
ADF$Tree_Name[which(ADF$binomial == "Chalcorana crassovis")] <- "Hylarana_crassovis"
ADF$Tree_Name[which(ADF$binomial == "Chalcorana campeni")] <- "Hylarana_campeni"
ADF$Tree_Name[which(ADF$binomial == "Dryophytes immaculatus")] <- "Hylarana_immaculata"
ADF$Tree_Name[which(ADF$binomial == "Dryophytes japonicus")] <- "Hylarana_japonica"
ADF$Tree_Name[which(ADF$binomial == "Dryophytes suweonensis")] <- "Hylarana_suweonensis"
ADF$Tree_Name[which(ADF$binomial == "Sylvirana nigrovittata")] <- "Hylarana_nigrovittata"

#Other
ADF$Tree_Name[which(ADF$binomial == "Espadarana audax")] <- "Centrolene_audax"
ADF$Tree_Name[which(ADF$binomial == "Megophrys carinense")] <- "Brachytarsophrys_carinense"
ADF$Tree_Name[which(ADF$binomial == "Nidirana adenopleura")] <- "Babina_adenopleura"
ADF$Tree_Name[which(ADF$binomial == "Rentapia everetti")] <- "Pedostibes_everetti"
#ADF$Tree_Name[which(ADF$binomial == "Sarcohyla labeculata")] <- "Plectrohyla calthula" #best guess
ADF$Tree_Name[which(ADF$binomial == "Sclerophrys cristigans")] <- "Amietophrynus_cristigans"
ADF$Tree_Name[which(ADF$binomial == "Sclerophrys togoensis")] <- "Amietophrynus_togoensis"
ADF$Tree_Name[which(ADF$binomial == "Uperodon anamaliensis")] <- "Ramanella_anamalaiensis"
ADF$Tree_Name[which(ADF$binomial == "Vitreorana ritae")] <- "Cochranella_ritae"

ADF <- ADF %>% mutate(inTree = sapply(Tree_Name, function(x){x %in% phy[[1]]$tip.label }), .after = Tree_Name)

saveRDS(ADF, file = "FossilAmphibians.rds")
write.csv(ADF, file = "FossilAmphibians.csv")
