setwd("C:/Users/User1/Box/Fossil Record") 
library('tidyverse')
#### RAW DATA PREP - Skip after 1st read-in ####
GARD <- read.csv(file ="Geo_Intersects/Reptile_Ranges.csv", header = T)
colnames(GARD)[2] <- "Area"
row.names(GARD) <- GARD$Binomial.C.80

reptileData <- read.csv("Geo_Intersects/Reptiles_Basins_Overlap.csv", header = T)
colnames(reptileData)[5] <- "Area"

GARD$fossilArea <- rep(0,length(GARD$Area))
GARD[reptileData$Binomial,"fossilArea"] <- reptileData$Area
GARD$genus <- sub(" .*", "", GARD$Binomial)
colnames(GARD) <- c("Binomial","Area","FID","TID","group","Value","fossilArea","genus")


load("IUCNcategories.rda")
IUCNcategories <- IUCNcategories[which(is.na(IUCNcategories$infra_rank)),] #pay no attention to below-species taxa
IUCNcategories <- IUCNcategories[which(is.na(IUCNcategories$population)),] #pay no attention to populations

GARD_IUCN <- IUCNcategories$scientific_name[which(IUCNcategories$scientific_name %in% GARD$Binomial)]

Leftout <- setdiff(IUCNcategories$scientific_name[which(IUCNcategories$order_name %in% c("SQUAMATA","TESTUDINES","CROCODYLIA","RHYNCHOCEPHALIA"))],GARD_IUCN)
length(Leftout) #719

GARD$REDcat <- rep(NA, length(GARD$Binomial))


GARD$REDcat[which(GARD$Binomial %in% GARD_IUCN)] <- IUCNcategories$category[which(IUCNcategories$scientific_name %in% GARD_IUCN)]

#check to make sure that the IUCN categories are correctly assigned:

names <- GARD$Binomial[which(!is.na(GARD$REDcat))]

all(GARD$REDcat[which(GARD$Binomial %in% names)] == IUCNcategories$category[which(IUCNcategories$scientific_name %in% names)])

GARD <- GARD[-10065,]

#get family data
library('taxize')

#the following genera break the taxoze workflow, remove them
exc <- c("Nebulifera","Python", "Tenuidactylus","Trachischium","Xenopeltis")
gens <-  unique(GARD$genus)[-(which(unique(GARD$genus) %in% exc))]

# find classification data for all these animals - TAKES A LONG TIME
#fams <- classification(gens, db = 'itis')

#save that so I never have to do it again.
saveRDS(fams,"ReptileFamilies.rds")

#now let's translate the taxize data into a usable form
famMap <- setNames(1:length(gens),gens)

for (i in names(famMap)){
  if(all(fams[i] %>% is.na)){
    famMap[i]<-NA
  }
  else{
    if(!('family' %in% fams[[i]]$rank)){
      famMap[i]<-NA
    }
    else{
      famMap[i] <- filter(fams[[i]], rank == 'family')$name
    }
  }
}

#add back in those genera that break taxize
famMap <- c(famMap,setNames(rep(NA,length(exc)),exc))

# manually enter data for genera taxize missed
famMap.na <- which(is.na(famMap))
famMap.na[1:36] <- c("Scincidae", #Afroablepharus
                     "Gekkonidae", #Altiphylax
                     "Diplodactylidae", #Amalosia
                     "Scincidae", #Androngo
                     "Scincidae", #Asymblepharus
                     "Scincidae", 
                     "Uropeltidae",
                     "Agamidae", 
                     "Sphaerodactylidae",
                     "Anguidae",
                     "Diplodactylidae",
                     "Scincidae",
                     "Diplodactylidae",
                     "Gekkonidae",
                     "Typhlopidae",
                     "Scincidae",
                     "Anguidae",
                     "Diplodactylidae",
                     "Scincidae",
                     "Diplodactylidae",
                     "Uropeltidae",
                     "Colubridae",
                     "Gymnophthalmidae",
                     "Gymnophthalmidae",
                     "Gekkonidae",
                     "Colubridae",
                     "Gekkonidae",
                     "Gekkonidae",
                     "Scincidae",
                     "Gymnophthalmidae",
                     "Colubridae",
                     "Diplodactylidae",
                     "Pythonidae",
                     "Gekkonidae",
                     "Colubridae",
                     "Xenopeltidae")

famMap[names(famMap.na)] <- famMap.na

#Map families onto Genera
GARD$family <- famMap[match(GARD$genus,names(famMap))]

GARD$family[which(GARD$family == "Trogonophiidae")] <- "Trogonophidae"

#Check for duplicates/misspellings
unique(GARD$family)%>% sort


GARD$Tree_Name <- gsub(" ","_",GARD$Binomial)

write.csv(GARD[,c("Binomial","genus","family","group","Tree_Name","Area","fossilArea","REDcat","FID","TID")], "FossilReptiles.csv")
saveRDS(GARD, file = "FossilReptiles.rds")
