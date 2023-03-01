library (tidyr)
library (dplyr)
library (readr)
library (ggplot2)

# #Load in Modern Fossil data
# df1 = read_csv(file = "FossilAmphibians.csv")
# df2 = read_csv(file = "FossilReptiles.csv")
# df3 = read_csv(file = "FossilMammals.csv")
# df4 = read_csv(file = "FossilBirds.csv")
# colnames(df2)[colnames(df2) == "Binomial"] ="binomial"
# 
# df1 = df1[c("binomial","genus","family","fossilArea","Area")] %>% mutate( taxon = "Amphibian")
# df2 = df2[c("binomial","genus","family","fossilArea","Area")] %>% mutate( taxon = "Reptile")
# df3 = df3[c("binomial","genus","family","fossilArea","Area")] %>% mutate( taxon = "Mammal")
# df4 = df4[c("binomial","genus","family","fossilArea","Area")] %>% mutate( taxon = "Bird")
# 
# all_modernfossil_verts <- rbind(df1, df2, df3, df4)
# 
# write_csv(all_modernfossil_verts, "all_modernfossil_verts.csv")
# 

#Load in all vertebrate potential fossil data
all_modernfossil_verts <- read_csv("all_modernfossil_verts.csv")

#Load in PBDB data
pbdb_verts = read_csv(file = "verts_PBDB_clean.csv")
#pbdb_verts <- pbdb_verts[which(pbdb_verts$is_extant == "extant"),]

# are species from our IUCN/GARD datasets in the PBDB dataset?
all_modernfossil_verts$PBDB_sp <- sapply(all_modernfossil_verts$binomial, 
                                         FUN = function(x){
                                           any(x %in% c(pbdb_verts$accepted_name,pbdb_verts$taxon_name))
                                           })

# are there extant species in the PBDB dataset that are not projected to fossilize?

missingsp <- which(
  !(filter(pbdb_verts, is_extant == "extant")["accepted_name"]%>% unlist %>% unique) %in% 
    all_modernfossil_verts$binomial)

missing <- pbdb_verts[missingsp,]

write_csv(missing, "missing_species.csv") # this is to be filled out and saved as an .xslx

#compare the two

same_sp = subset(all_modernfossil_verts, all_modernfossil_verts$binomial %in% pbdb_verts$accepted_name)

missing_sp <- subset(all_modernfossil_verts, !(pbdb_verts$accepted_name %in% all_modernfossil_verts$binomial ))

missing_gen <- unique(pbdb_verts$genus)[which(!(unique(pbdb_verts$genus) %in%  all_modernfossil_verts$genus))]

missing_fam <- unique(pbdb_verts$family)[which(!(unique(pbdb_verts$family) %>% toupper() %in%  
                                                   toupper(all_modernfossil_verts$family)))]
missingtaxa <- data.frame(taxon = c(missing_gen,missing_fam),
                          rank = c(rep("genus",length(missing_gen)),
                                   rep("family", length(missing_fam))))

write_csv(missingtaxa, "missing_genus_family.csv")

#percent of species in PBDB
length(which(all_modernfossil_verts$PBDB_sp))/dim(all_modernfossil_verts)[1]

length(which(all_modernfossil_verts$PBDB_sp))/length(unique(pbdb_verts$accepted_name))
