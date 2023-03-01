setwd("C:/Users/User1/Box/Fossil Record")
library("tidyverse")

#Load original ranges
originalRanges <- read.csv("Geo_Intersects/BIRDS_Dissolve_TableToExcel_1.csv", header = T)
#load fossil ranges
fossilRanges <- read.csv("Geo_Intersects/BIRDS_PairwiseIntersect_all_TableToExcel.csv", header = T)

#set up area and fossil area columns
colnames(originalRanges)[6] <- "Area"
originalRanges$fossilArea <- 0

#add fossilArea to originalranges df
originalRanges[match(fossilRanges$binomial,originalRanges$binomial),]$fossilArea <- fossilRanges$final_area_km2

#add taxonomic data from BOTW v5 (compatible with our bird geography data)

taxonomy <- read.csv("Results/Bird/HBW-BirdLife_List_of_Birds_v5.csv", header = T)
unmatched <- which(!(taxonomy$Scientific.name %in% originalRanges$binomial))
taxonomy <- taxonomy[-unmatched,]

#get the dfs in the same order
originalRanges <- originalRanges[order(originalRanges$binomial),]
taxonomy <- taxonomy[order(taxonomy$Scientific.name),]

#check if they have all the same species
all(originalRanges$binomial == taxonomy$Scientific.name) #looks like we're all set

BDF <- data.frame(binomial = originalRanges$binomial,
                  kingdom = 'ANIMALIA',
                  phylum = 'CHORDATA',
                  class = 'AVES',
                  order_ = taxonomy$Order,
                  family = taxonomy$Family.name,
                  genus = gsub(" .+","",originalRanges$binomial),
                  category = taxonomy$X2020.IUCN.Red.List.category,
                  Area = originalRanges$Area,
                  fossilArea = originalRanges$fossilArea,
                  Tree_Name = sub(" ","_",originalRanges$binomial),
                  passerine = 'Passerine')
BDF$passerine[which(BDF$order_ != 'PASSERIFORMES')] <- "Nonpasserine"

IUCNlevels <- c("LC","NT","VU","EN","CR","EW","EX","DD")#order these from least to most extinct

BDF$category <- ordered(BDF$category, IUCNlevels)

source("Scripts/getClearSynonyms.R")
library("phytools")
phy <- read.nexus("Phylogenies/Bird/output.nex")
BirdTaxonomy <- get_clear_synonyms(BDF,phy)
BirdTaxonomy <- readRDS(file = "Results/Bird/BirdTaxonomy.rda")
#saveRDS(BirdTaxonomy,file = "Results/Bird/BirdTaxonomy.rda")

BDF$Tree_Name <- BirdTaxonomy$namematches$replacement

#fix  1 in BirdTaxonomy$unclear
BDF$Tree_Name[which(BDF$binomial == "Setophaga pinus")] <- "Dendroica_pinus" #only 1 in BirdTaxonomy$unclear

BDF <- BDF %>% mutate(inTree = sapply(Tree_Name, function(x){x %in% phy[[1]]$tip.label }), .after = Tree_Name)

#helper functions to look for species and genera in the tree
chphy <- function(x){ x %in% phy[[1]]$tip.label}
chgen <- function(x){phy[[1]]$tip.label[str_which(phy[[1]]$tip.label,x)]}

#look for potentially significant generic changes
commongenera <- table(filter(BDF, inTree == FALSE)$genus)



chgen("Nesoclopeus")
BDF$Tree_Name[which(BDF$Tree_Name == "Hypotaenidia_woodfordi")] <- "Nesoclopeus_woodfordi"

#Hydrobates
chgen("Hydrobates")
chgen("Oceanodroma") %>% sort()
view(filter(BDF, genus == "Hydrobates"))

BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_castro")] <- "Oceanodroma_castro"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_furcatus")] <- "Oceanodroma_furcata"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_homochroa")] <- "Oceanodroma_homochroa"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_hornbyi")] <- "Oceanodroma_hornbyi"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_markhami")] <- "Oceanodroma_markhami"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_matsudairae")] <- "Oceanodroma_matsudairae"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_melania")] <- "Oceanodroma_melania"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_monorhis")] <- "Oceanodroma_monorhis"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_monteiroi")] <- "Oceanodroma_monteiroi"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_tethys")] <- "Oceanodroma_tethyse"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_tristrami")] <- "Oceanodroma_tristrami"



#Spatula
chgen("Spatula")
view(filter(BDF, genus == "Spatula"))
BDF$Tree_Name[which(BDF$genus == "Spatula")] <- str_replace(BDF$Tree_Name[which(BDF$genus == "Spatula")],
                                                            "Spatula_",
                                                            "Anas_")
# Aerodramus
chgen("Aerodramus") %>% sort()
chgen("Collocalia") %>% sort()
view(filter(BDF, genus == "Aerodramus", inTree == FALSE))
BDF$Tree_Name[which(BDF$Tree_Name == "Aerodramus_brevirostris")] <- "Collocalia_brevirostris"
BDF$Tree_Name[which(BDF$Tree_Name == "Aerodramus_fuciphagus")] <- "Collocalia_fuciphaga"
BDF$Tree_Name[which(BDF$Tree_Name == "Aerodramus_leucophaeus")] <- "Collocalia_leucophaea"
BDF$Tree_Name[which(BDF$Tree_Name == "Aerodramus_spodiopygius")] <- "Collocalia_spodiopygia"
BDF$Tree_Name[which(BDF$Tree_Name == "Aerodramus_terraereginae")] <- "Collocalia_terraereginae"
BDF$Tree_Name[which(BDF$Tree_Name == "Aerodramus_vanikorensis")] <- "Collocalia_vanikorensis"

# Aethopyga
chgen("Aethopyga") %>% sort()
view(filter(BDF, genus == "Aethopyga", inTree == FALSE))

#Alopeconas
chgen("Alopecoenas") %>% sort()
chgen("Gallicolumba") %>% sort()
view(filter(BDF, genus == "Alopecoenas", inTree == FALSE))
BDF$Tree_Name[which(BDF$Tree_Name == "Alopecoenas_beccarii")] <- "Gallicolumba_beccarii"
BDF$Tree_Name[which(BDF$Tree_Name == "Alopecoenas_canifrons")] <- "Gallicolumba_canifrons"
BDF$Tree_Name[which(BDF$Tree_Name == "Alopecoenas_erythropterus")] <- "Gallicolumba_erythroptera"
BDF$Tree_Name[which(BDF$Tree_Name == "Alopecoenas_hoedtii")] <- "Gallicolumba_hoedtii"
BDF$Tree_Name[which(BDF$Tree_Name == "Alopecoenas_jobiensis")] <- "Gallicolumba_jobiensis"
BDF$Tree_Name[which(BDF$Tree_Name == "Alopecoenas_kubaryi")] <- "Gallicolumba_kubaryi"
BDF$Tree_Name[which(BDF$Tree_Name == "Alopecoenas_rubescens")] <- "Gallicolumba_rubescens"
BDF$Tree_Name[which(BDF$Tree_Name == "Alopecoenas_sanctaecrucis")] <- "Gallicolumba_sanctaecrucis"
BDF$Tree_Name[which(BDF$Tree_Name == "Alopecoenas_stairi")] <- "Gallicolumba_stairi"
BDF$Tree_Name[which(BDF$Tree_Name == "Alopecoenas_xanthonurus")] <- "Gallicolumba_xanthonura"

#Amazilia
chgen("Amazilia") %>% sort()
chgen("Hylocharis") %>% sort()
chgen("Chrysuronia") %>% sort()
chgen("Lepidopyga") %>% sort()
chgen("Damophila") %>% sort()
filter(BDF, genus == "Amazilia", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Amazilia_chrysura")] <- "Hylocharis_chrysura"
BDF$Tree_Name[which(BDF$Tree_Name == "Amazilia_cyanus")] <- "Hylocharis_cyanus"
BDF$Tree_Name[which(BDF$Tree_Name == "Amazilia_eliciae")] <- "Hylocharis_eliciae"
BDF$Tree_Name[which(BDF$Tree_Name == "Amazilia_grayi")] <- "Hylocharis_grayi"
BDF$Tree_Name[which(BDF$Tree_Name == "Amazilia_humboldtii")] <- "Hylocharis_humboldtii"
BDF$Tree_Name[which(BDF$Tree_Name == "Amazilia_sapphirina")] <- "Hylocharis_sapphirina"
BDF$Tree_Name[which(BDF$Tree_Name == "Amazilia_chrysura")] <- "Hylocharis_chrysura"
BDF$Tree_Name[which(BDF$Tree_Name == "Amazilia_oenone")] <- "Chrysuronia_oenone"
BDF$Tree_Name[which(BDF$Tree_Name == "Amazilia_coeruleogularis")] <- "Lepidopyga_coeruleogularis"
BDF$Tree_Name[which(BDF$Tree_Name == "Amazilia_goudoti")] <- "Lepidopyga_goudoti"
BDF$Tree_Name[which(BDF$Tree_Name == "Amazilia_lilliae")] <- "Lepidopyga_lilliae"
BDF$Tree_Name[which(BDF$Tree_Name == "Amazilia_julie")] <- "Damophila_julie"

#Antigone
chgen("Grus") %>% sort()
filter(BDF, genus == "Antigone", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Grus_antigone")] <- "Antigone_antigone"
BDF$Tree_Name[which(BDF$Tree_Name == "Grus_canadensis")] <- "Antigone_canadensis"
BDF$Tree_Name[which(BDF$Tree_Name == "Grus_rubicunda")] <- "Antigone_rubicunda"
BDF$Tree_Name[which(BDF$Tree_Name == "Grus_vipio")] <- "Antigone_vipio"

#Ardenna
chgen("Puffinus") %>% sort()
filter(BDF, genus == "Ardenna", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Ardenna_bulleri")] <- "Puffinus_bulleri"
BDF$Tree_Name[which(BDF$Tree_Name == "Ardenna_carneipes")] <- "Puffinus_carneipes"
BDF$Tree_Name[which(BDF$Tree_Name == "Ardenna_creatopus")] <- "Puffinus_creatopus"
BDF$Tree_Name[which(BDF$Tree_Name == "Ardenna_gravis")] <- "Puffinus_gravis"
BDF$Tree_Name[which(BDF$Tree_Name == "Ardenna_grisea")] <- "Puffinus_griseus"
BDF$Tree_Name[which(BDF$Tree_Name == "Ardenna_pacifica")] <- "Puffinus_pacificus"
BDF$Tree_Name[which(BDF$Tree_Name == "Ardenna_bulleri")] <- "Puffinus_bulleri"
BDF$Tree_Name[which(BDF$Tree_Name == "Ardenna_tenuirostris")] <- "Puffinus_tenuirostris"

#Argya
chgen("Turdoides") %>% sort()
filter(BDF, genus == "Argya", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Argya_altirostris")] <- "Turdoides_altirostris"
BDF$Tree_Name[which(BDF$Tree_Name == "Argya_aylmeri")] <- "Turdoides_aylmeri"
BDF$Tree_Name[which(BDF$Tree_Name == "Argya_caudata")] <- "Turdoides_caudata"
BDF$Tree_Name[which(BDF$Tree_Name == "Argya_earlei")] <- "Turdoides_earlei"
BDF$Tree_Name[which(BDF$Tree_Name == "Argya_fulva")] <- "Turdoides_fulva"
BDF$Tree_Name[which(BDF$Tree_Name == "Argya_malcolmi")] <- "Turdoides_malcolmi"
BDF$Tree_Name[which(BDF$Tree_Name == "Argya_rubiginosa")] <- "Turdoides_rubiginosa"
BDF$Tree_Name[which(BDF$Tree_Name == "Argya_squamiceps")] <- "Turdoides_squamiceps"
BDF$Tree_Name[which(BDF$Tree_Name == "Argya_subrufa")] <- "Turdoides_subrufa"

#Asthenes
chgen("Schizoeaca") %>% sort()
filter(BDF, genus == "Asthenes", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Asthenes_coryi")] <- "Schizoeaca_coryi"
BDF$Tree_Name[which(BDF$Tree_Name == "Asthenes_fuliginosa")] <- "Schizoeaca_fuliginosa"
BDF$Tree_Name[which(BDF$Tree_Name == "Asthenes_griseomurina")] <- "Schizoeaca_griseomurina"
BDF$Tree_Name[which(BDF$Tree_Name == "Asthenes_palpebralis")] <- "Schizoeaca_palpebralis"
BDF$Tree_Name[which(BDF$Tree_Name == "Asthenes_harteri")] <- "Schizoeaca_harteri"
BDF$Tree_Name[which(BDF$Tree_Name == "Asthenes_moreirae")] <- "Schizoeaca_moreirae"
BDF$Tree_Name[which(BDF$Tree_Name == "Asthenes_perijana")] <- "Schizoeaca_perijana"
BDF$Tree_Name[which(BDF$Tree_Name == "Asthenes_harteri")] <- "Schizoeaca_harteri"
BDF$Tree_Name[which(BDF$Tree_Name == "Asthenes_usheri")] <- "Schizoeaca_usheri"

#Calidris
chgen("Limicola") %>% sort()
chgen("Philomachus") %>% sort()
chgen("Eurynorhynchus") %>% sort()
chgen("Aphriza") %>% sort()
chgen("Tryngites") %>% sort()
filter(BDF, genus == "Calidris", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Calidris_falcinellus")] <- "Limicola_falcinellus"
BDF$Tree_Name[which(BDF$Tree_Name == "Calidris_pugnax")] <- "Philomachus_pugnax"
BDF$Tree_Name[which(BDF$Tree_Name == "Calidris_pygmeus")] <- "Eurynorhynchus_pygmeus"
BDF$Tree_Name[which(BDF$Tree_Name == "Calidris_virgata")] <- "Aphriza_virgata"
BDF$Tree_Name[which(BDF$Tree_Name == "Calidris_subruficollis")] <- "Tryngites_subruficollis"

#Calliope
chgen("Luscinia") %>% sort()
filter(BDF, genus == "Calliope", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Calliope_calliope")] <- "Luscinia_calliope"
BDF$Tree_Name[which(BDF$Tree_Name == "Calliope_obscura")] <- "Luscinia_obscura"
BDF$Tree_Name[which(BDF$Tree_Name == "Calliope_pectardens")] <- "Luscinia_pectardens"
BDF$Tree_Name[which(BDF$Tree_Name == "Calliope_pectoralis")] <- "Luscinia_pectoralis"

#Carpodacus
chgen("Carpodacus") %>% sort()
filter(BDF, genus == "Carpodacus", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Carpodacus_sillemi")] <- "Leucosticte_sillemi"

#Ceblepyris
chgen("Coracina") %>% sort()
filter(BDF, genus == "Ceblepyris", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Ceblepyris_caesius")] <- "Coracina_caesia"
BDF$Tree_Name[which(BDF$Tree_Name == "Ceblepyris_cinereus")] <- "Coracina_cinerea"
BDF$Tree_Name[which(BDF$Tree_Name == "Ceblepyris_graueri")] <- "Coracina_graueri"
BDF$Tree_Name[which(BDF$Tree_Name == "Ceblepyris_pectoralis")] <- "Coracina_pectoralis"

#Ceratopipra 
chgen("Pipra") %>% sort()
filter(BDF, genus == "Ceratopipra", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Ceratopipra_chloromeros")] <- "Pipra_chloromeros"
BDF$Tree_Name[which(BDF$Tree_Name == "Ceratopipra_cornuta")] <- "Pipra_cornuta"
BDF$Tree_Name[which(BDF$Tree_Name == "Ceratopipra_erythrocephala")] <- "Pipra_erythrocephala"
BDF$Tree_Name[which(BDF$Tree_Name == "Ceratopipra_mentalis")] <- "Pipra_mentalis"
BDF$Tree_Name[which(BDF$Tree_Name == "Ceratopipra_rubrocapilla")] <- "Pipra_rubrocapilla"

#Cercomacroides
chgen("Cercomacra") %>% sort()
filter(BDF, genus == "Cercomacroides", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Cercomacroides_laeta" )] <- "Cercomacra_laeta"
BDF$Tree_Name[which(BDF$Tree_Name == "Cercomacroides_parkeri" )] <- "Cercomacra_parkeri"
BDF$Tree_Name[which(BDF$Tree_Name == "Cercomacroides_serva" )] <- "Cercomacra_serva"
BDF$Tree_Name[which(BDF$Tree_Name == "Cercomacroides_tyrannina" )] <- "Cercomacra_tyrannina"

#Ceyx
chgen("Alcedo") %>% sort()
filter(BDF, genus == "Ceyx", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Ceyx_azureus" )] <- "Alcedo_azurea"

#Chalcites
chgen("Chrysococcyx") %>% sort()
filter(BDF, genus == "Chalcites", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Chalcites_basalis" )] <- "Chrysococcyx_basalis"
BDF$Tree_Name[which(BDF$Tree_Name == "Chalcites_crassirostris" )] <- "Chrysococcyx_crassirostris"
BDF$Tree_Name[which(BDF$Tree_Name == "Chalcites_lucidus" )] <- "Chrysococcyx_lucidus"
BDF$Tree_Name[which(BDF$Tree_Name == "Chalcites_megarhynchus")] <- "Rhamphomantis_megarhynchus"
BDF$Tree_Name[which(BDF$Tree_Name == "Chalcites_meyerii" )] <- "Chrysococcyx_meyeri"
BDF$Tree_Name[which(BDF$Tree_Name == "Chalcites_minutillus" )] <- "Chrysococcyx_minutillus"
BDF$Tree_Name[which(BDF$Tree_Name == "Chalcites_osculans" )] <- "Chrysococcyx_osculans"
BDF$Tree_Name[which(BDF$Tree_Name == "Chalcites_ruficollis" )] <- "Chrysococcyx_ruficollis"

#Cinnyris
chgen("Nectarinia") %>% sort()
filter(BDF, genus == "Cinnyris", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Cinnyris_chalcomelas" )] <- "Nectarinia_chalcomelas"
BDF$Tree_Name[which(BDF$Tree_Name == "Cinnyris_coccinigastrus" )] <- "Nectarinia_coccinigaster"
BDF$Tree_Name[which(BDF$Tree_Name == "Cinnyris_fuelleborni" )] <- "Nectarinia_fuelleborni"
BDF$Tree_Name[which(BDF$Tree_Name == "Cinnyris_neergaardi" )] <- "Nectarinia_neergardi"
BDF$Tree_Name[which(BDF$Tree_Name == "Cinnyris_prigoginei" )] <- "Nectarinia_chalcomelas"
BDF$Tree_Name[which(BDF$Tree_Name == "Cinnyris_reichenowi" )] <- "Nectarinia_reichenowi"
BDF$Tree_Name[which(BDF$Tree_Name == "Cinnyris_tsavoensis" )] <- "Nectarinia_tsavoensis" 

#Cyanoderma
chgen("Stachyris") %>% sort()
filter(BDF, genus == "Cyanoderma", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Cyanoderma_chrysaeum" )] <- "Stachyris_chrysaea" 
BDF$Tree_Name[which(BDF$Tree_Name == "Cyanoderma_melanothorax" )] <- "Stachyris_melanothorax" 
BDF$Tree_Name[which(BDF$Tree_Name == "Cyanoderma_pyrrhops" )] <- "Stachyris_pyrrhops"
BDF$Tree_Name[which(BDF$Tree_Name == "Cyanoderma_ruficeps" )] <- "Stachyris_ruficeps"
BDF$Tree_Name[which(BDF$Tree_Name == "Cyanoderma_rufifrons" )] <- "Stachyris_rufifrons"

#Cyornis
chgen("Rhinomyias") %>% sort()
filter(BDF, genus == "Cyornis", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Cyornis_colonus" )] <- "Rhinomyias_colonus"
BDF$Tree_Name[which(BDF$Tree_Name == "Cyornis_olivaceus" )] <- "Rhinomyias_olivaceus"
BDF$Tree_Name[which(BDF$Tree_Name == "Cyornis_pallidipes" )] <- "Cyornis_pallipes"
BDF$Tree_Name[which(BDF$Tree_Name == "Cyornis_ruficauda" )] <- "Cyornis_ruficauda"
BDF$Tree_Name[which(BDF$Tree_Name == "Cyornis_umbratilis" )] <- "Cyornis_umbratilis"

#Dessonornis
chgen("Cossypha") %>% sort()
filter(BDF, genus == "Dessonornis", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Dessonornis_archeri" )] <- "Cossypha_archeri"
BDF$Tree_Name[which(BDF$Tree_Name == "Dessonornis_caffer" )] <- "Cossypha_caffra"
BDF$Tree_Name[which(BDF$Tree_Name == "Dessonornis_humeralis" )] <- "Cossypha_humeralis"

#Dyaphorophyia
chgen("Platysteira") %>% sort()
filter(BDF, genus == "Dyaphorophyia", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Dyaphorophyia_blissetti" )] <- "Platysteira_blissetti"
BDF$Tree_Name[which(BDF$Tree_Name == "Dyaphorophyia_castanea" )] <- "Platysteira_castanea"
BDF$Tree_Name[which(BDF$Tree_Name == "Dyaphorophyia_chalybea" )] <- "Platysteira_chalybea"
BDF$Tree_Name[which(BDF$Tree_Name == "Dyaphorophyia_concreta" )] <- "Platysteira_concreta"
BDF$Tree_Name[which(BDF$Tree_Name == "Dyaphorophyia_jamesoni" )] <- "Platysteira_jamesoni"
BDF$Tree_Name[which(BDF$Tree_Name == "Dyaphorophyia_tonsa" )] <- "Platysteira_tonsa"

#Edolisoma
chgen("Edolisoma")
chgen("Lalage")%>% sort()
chgen("Coracina") %>% sort()
view(filter(BDF, genus == "Edolisoma"))

BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_anale")] <- "Coracina_analis"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_ceramense")] <- "Coracina_ceramensis"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_coerulescens")] <- "Coracina_coerulescens"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_dispar")] <- "Coracina_dispar"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_dohertyi")] <- "Coracina_dohertyi"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_holopolium")] <- "Coracina_holopolia"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_incertum")] <- "Coracina_incerta"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_melas")] <- "Coracina_melas"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_mindanense")] <- "Coracina_mindanensis"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_montanum")] <- "Coracina_montana"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_morio")] <- "Coracina_morio"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_ostentum")] <- "Coracina_ostenta"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_parvulum")] <- "Coracina_parvula"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_schisticeps")] <- "Coracina_schisticeps"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_sula")] <- "Coracina_sula"
BDF$Tree_Name[which(BDF$Tree_Name == "Edolisoma_tenuirostre")] <- "Coracina_tenuirostris"


#Erythrogenys
chgen("Pomatorhinus") %>% sort()
filter(BDF, genus == "Erythrogenys", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Erythrogenys_erythrocnemis" )] <- "Pomatorhinus_erythrocnemis"
BDF$Tree_Name[which(BDF$Tree_Name == "Erythrogenys_erythrogenys" )] <- "Pomatorhinus_erythrogenys"
BDF$Tree_Name[which(BDF$Tree_Name == "Erythrogenys_gravivox" )] <- "Pomatorhinus_gravivox"
BDF$Tree_Name[which(BDF$Tree_Name == "Erythrogenys_hypoleucos" )] <- "Pomatorhinus_hypoleucos"
BDF$Tree_Name[which(BDF$Tree_Name == "Erythrogenys_mcclellandi" )] <- "Pomatorhinus_mcclellandi"
BDF$Tree_Name[which(BDF$Tree_Name == "Erythrogenys_swinhoei" )] <- "Pomatorhinus_swinhoei"

#Erythropitta
chgen("Pitta") %>% sort()
filter(BDF, genus == "Erythropitta", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Erythropitta_arquata" )] <- "Pitta_arcuata"
BDF$Tree_Name[which(BDF$Tree_Name == "Erythropitta_dohertyi" )] <- "Pitta_dohertyi"
BDF$Tree_Name[which(BDF$Tree_Name == "Erythropitta_granatina" )] <- "Pitta_granatina"
BDF$Tree_Name[which(BDF$Tree_Name == "Erythropitta_kochi" )] <- "Pitta_kochi"
BDF$Tree_Name[which(BDF$Tree_Name == "Erythropitta_kochi" )] <- "Pitta_kochi"
BDF$Tree_Name[which(BDF$Tree_Name == "Erythropitta_venusta" )] <- "Pitta_venusta"

#Eupsittula
chgen("Eupsittula") %>% sort()
filter(BDF, genus == "Eupsittula", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Eupsittula_aurea" )] <- "Aratinga_aurea"
BDF$Tree_Name[which(BDF$Tree_Name == "Eupsittula_cactorum" )] <- "Aratinga_cactorum"
BDF$Tree_Name[which(BDF$Tree_Name == "Eupsittula_canicularis" )] <- "Aratinga_canicularis"
BDF$Tree_Name[which(BDF$Tree_Name == "Eupsittula_pertinax" )] <- "Aratinga_pertinax"

#Ficedula
chgen("Ficedula") %>% sort()
filter(BDF, genus == "Ficedula", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Ficedula_hodgsoni" )] <- "Muscicapella_hodgsoni"
BDF$Tree_Name[which(BDF$Tree_Name == "Ficedula_erithacus" )] <- "Ficedula_hodgsonii"

#Fraseria
chgen("Muscicapa") %>% sort()
chgen("Myioparus") %>% sort()
filter(BDF, genus == "Fraseria", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Fraseria_caerulescens" )] <- "Muscicapa_caerulescens"
BDF$Tree_Name[which(BDF$Tree_Name == "Fraseria_griseigularis" )] <- "Myioparus_griseigularis"
BDF$Tree_Name[which(BDF$Tree_Name == "Fraseria_lendu" )] <- "Muscicapa_lendu"
BDF$Tree_Name[which(BDF$Tree_Name == "Fraseria_olivascens" )] <- "Muscicapa_olivascens"
BDF$Tree_Name[which(BDF$Tree_Name == "Fraseria_plumbea" )] <- "Myioparus_plumbeus"
BDF$Tree_Name[which(BDF$Tree_Name == "Fraseria_tessmanni" )] <- "Muscicapa_tessmanni"

#Garrulax
chgen("Garrulax") %>% sort()
chgen("Pterorhinus") %>% sort()
chgen("Ianthocincla") %>% sort()
chgen("Babax") %>% sort()
chgen("Kaznakowia") %>% sort()
filter(BDF, genus == "Garrulax", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Garrulax_koslowi" )] <- "Babax_koslowi"
BDF$Tree_Name[which(BDF$Tree_Name == "Garrulax_lanceolatus" )] <- "Babax_lanceolatus"
BDF$Tree_Name[which(BDF$Tree_Name == "Garrulax_waddelli" )] <- "Babax_waddelli"

#Heleia
chgen("Heleia") %>% sort()
chgen("Lophozosterops") %>% sort()
chgen("Oculocincta") %>% sort()
filter(BDF, genus == "Heleia", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Heleia_dohertyi" )] <- "Lophozosterops_dohertyi"
BDF$Tree_Name[which(BDF$Tree_Name == "Heleia_goodfellowi" )] <- "Lophozosterops_goodfellowi"
BDF$Tree_Name[which(BDF$Tree_Name == "Heleia_javanica" )] <- "Lophozosterops_javanicus"
BDF$Tree_Name[which(BDF$Tree_Name == "Heleia_pinaiae" )] <- "Lophozosterops_pinaiae"
BDF$Tree_Name[which(BDF$Tree_Name == "Heleia_squamiceps" )] <- "Lophozosterops_squamiceps"
BDF$Tree_Name[which(BDF$Tree_Name == "Heleia_squamifrons" )] <- "Oculocincta_squamifrons"
BDF$Tree_Name[which(BDF$Tree_Name == "Heleia_superciliaris" )] <- "Lophozosterops_superciliaris"
BDF$Tree_Name[which(BDF$Tree_Name == "Heleia_wallacei" )] <- "Zosterops_wallacei"

#Hierococcyx
chgen("Cuculus") %>% sort()
filter(BDF, genus == "Hierococcyx", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Hierococcyx_vagans" )] <- "Cuculus_vagans"
BDF$Tree_Name[which(BDF$Tree_Name == "Hierococcyx_varius" )] <- "Cuculus_varius"

#Hydrobates
chgen("Hydrobates")
chgen("Oceanodroma") %>% sort()
filter(BDF, genus == "Hydrobates", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_castro")] <- "Oceanodroma_castro"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_furcatus")] <- "Oceanodroma_furcata"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_homochroa")] <- "Oceanodroma_homochroa"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_hornbyi")] <- "Oceanodroma_hornbyi"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_markhami")] <- "Oceanodroma_markhami"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_matsudairae")] <- "Oceanodroma_matsudairae"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_melania")] <- "Oceanodroma_melania"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_monorhis")] <- "Oceanodroma_monorhis"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_monteiroi")] <- "Oceanodroma_monteiroi"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_tethys")] <- "Oceanodroma_tethyse"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrobates_tristrami")] <- "Oceanodroma_tristrami"

#Hydrornis
chgen("Pitta") %>% sort()
filter(BDF, genus == "Hydrornis", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrornis_baudii")] <- "Pitta_baudii"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrornis_caeruleus")] <- "Pitta_caerulea"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrornis_cyaneus")] <- "Pitta_cyanea"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrornis_elliotii")] <- "Pitta_elliotii"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrornis_guajanus")] <- "Pitta_guajana"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrornis_gurneyi")] <- "Pitta_gurneyi"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrornis_nipalensis")] <- "Pitta_nipalensis"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrornis_oatesi")] <- "Pitta_oatesi"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrornis_phayrei")] <- "Pitta_phayrei"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrornis_schneideri")] <- "Pitta_schneideri"
BDF$Tree_Name[which(BDF$Tree_Name == "Hydrornis_soror")] <- "Pitta_soror"

#Hypotaenidia
chgen("Hypotaenidia")
chgen("Gallirallus") %>% sort()
filter(BDF, genus == "Hypotaenidia", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Hypotaenidia_insignis")] <- "Gallirallus_insignis"
BDF$Tree_Name[which(BDF$Tree_Name == "Hypotaenidia_okinawae")] <- "Gallirallus_okinawae"
BDF$Tree_Name[which(BDF$Tree_Name == "Hypotaenidia_owstoni")] <- "Gallirallus_owstoni"
BDF$Tree_Name[which(BDF$Tree_Name == "Hypotaenidia_philippensis")] <- "Gallirallus_philippensis"
BDF$Tree_Name[which(BDF$Tree_Name == "Hypotaenidia_rovianae")] <- "Gallirallus_rovianae"
BDF$Tree_Name[which(BDF$Tree_Name == "Hypotaenidia_sylvestris")] <- "Gallirallus_sylvestris"
BDF$Tree_Name[which(BDF$Tree_Name == "Hypotaenidia_torquata")] <- "Gallirallus_torquatus"

#Icterus
chgen("Icterus") %>% sort()
filter(BDF, genus == "Icterus", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Icterus_bullockiorum")] <- "Icterus_bullockii"

#Kittacincla
chgen("Kittacincla") %>% sort()
filter(BDF, genus == "Kittacincla", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Kittacincla_cebuensis")] <- "Copsychus_cebuensis"
BDF$Tree_Name[which(BDF$Tree_Name == "Kittacincla_nigra")] <- "Copsychus_niger"

#Lalage
chgen("Lalage") %>% sort()
chgen("Coracina") %>% sort()
filter(BDF, genus == "Lalage", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Lalage_fimbriata")] <- "Coracina_fimbriata"
BDF$Tree_Name[which(BDF$Tree_Name == "Lalage_melanoptera")] <- "Coracina_melanoptera"
BDF$Tree_Name[which(BDF$Tree_Name == "Lalage_melaschistos")] <- "Coracina_melaschistos"
BDF$Tree_Name[which(BDF$Tree_Name == "Lalage_newtoni")] <- "Coracina_newtoni"
BDF$Tree_Name[which(BDF$Tree_Name == "Lalage_polioptera")] <- "Coracina_polioptera"
BDF$Tree_Name[which(BDF$Tree_Name == "Lalage_typica")] <- "Coracina_typica"

#Larvivora
chgen("Larvivora") %>% sort()
chgen("Luscinia") %>% sort()
filter(BDF, genus == "Larvivora", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Larvivora_brunnea")] <- "Luscinia_brunnea"
BDF$Tree_Name[which(BDF$Tree_Name == "Larvivora_cyane")] <- "Luscinia_cyane"
BDF$Tree_Name[which(BDF$Tree_Name == "Larvivora_ruficeps")] <- "Luscinia_ruficeps"
BDF$Tree_Name[which(BDF$Tree_Name == "Larvivora_sibilans")] <- "Luscinia_sibilans"

#Laterallus
chgen("Laterallus") %>% sort()
chgen("Hapalocrex") %>% sort()
chgen("Atlantisia") %>% sort()
chgen("Porzana") %>% sort()
filter(BDF, genus == "Laterallus", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Laterallus_flaviventer")] <- "Porzana_flaviventer"
BDF$Tree_Name[which(BDF$Tree_Name == "Laterallus_rogersi")] <- "Atlantisia_rogersi"
BDF$Tree_Name[which(BDF$Tree_Name == "Laterallus_spilonota")] <- "Laterallus_spilonotus"
BDF$Tree_Name[which(BDF$Tree_Name == "Laterallus_spilopterus")] <- "Porzana_spiloptera"

#Leistes
chgen("Leistes") %>% sort()
chgen("Sturnella") %>% sort()
filter(BDF, genus == "Leistes", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Leistes_bellicosus")] <- "Sturnella_bellicosa"
BDF$Tree_Name[which(BDF$Tree_Name == "Leistes_defilippii")] <- "Sturnella_defilippii"
BDF$Tree_Name[which(BDF$Tree_Name == "Leistes_loyca")] <- "Sturnella_loyca"
BDF$Tree_Name[which(BDF$Tree_Name == "Leistes_militaris")] <- "Sturnella_militaris"
BDF$Tree_Name[which(BDF$Tree_Name == "Leistes_superciliaris")] <- "Sturnella_superciliaris"

#Lophorina
chgen("Lophorina") %>% sort()
chgen("Ptiloris") %>% sort()
filter(BDF, genus == "Lophorina", inTree == FALSE)$Tree_Name

BDF$Tree_Name[which(BDF$Tree_Name == "Lophorina_intercedens")] <- "Ptiloris_intercedens"
BDF$Tree_Name[which(BDF$Tree_Name == "Lophorina_magnifica")] <- "Ptiloris_magnificus"
BDF$Tree_Name[which(BDF$Tree_Name == "Lophorina_paradisea")] <- "Ptiloris_paradiseus"
BDF$Tree_Name[which(BDF$Tree_Name == "Lophorina_victoriae")] <- "Ptiloris_victoriae"

#Melaenornis
chgen("Melaenornis") %>% sort()
chgen("Dioptrornis") %>% sort()
chgen("Namibornis") %>% sort()
chgen("Sigelus") %>% sort()
filter(BDF, genus == "Melaenornis", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Melaenornis_brunneus")] <- "Dioptrornis_brunneus"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaenornis_chocolatinus")] <- "Dioptrornis_chocolatinus"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaenornis_fischeri")] <- "Dioptrornis_fischeri"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaenornis_herero")] <- "Namibornis_herero"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaenornis_semipartitus")] <- "Empidornis_semipartitus"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaenornis_silens")] <- "Sigelus_silens"

#Melaniparus
chgen("Melaniparus") %>% sort()
chgen("Parus") %>% sort()
filter(BDF, genus == "Melaniparus", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Melaniparus_afer")] <- "Parus_afer"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaniparus_albiventris")] <- "Parus_albiventris"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaniparus_amabilis")] <- "Parus_amabilis"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaniparus_cinerascens")] <- "Parus_cinerascen"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaniparus_fasciiventer")] <- "Parus_fasciiventer"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaniparus_fringillinus")] <- "Parus_fringillinus"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaniparus_funereus")] <- "Parus_funereus"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaniparus_griseiventris")] <- "Parus_griseiventris"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaniparus_guineensis")] <- "Parus_guineensis"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaniparus_leucomelas")] <- "Parus_leucomelas"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaniparus_leuconotus")] <- "Parus_leuconotus"
BDF$Tree_Name[which(BDF$Tree_Name == "Melaniparus_thruppi")] <- "Parus_thruppi"

#Mixornis
chgen("Mixornis") %>% sort()
chgen("Macronous") %>% sort()
filter(BDF, genus == "Mixornis", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Mixornis_bornensis")] <- "Macronous_bornensis"
BDF$Tree_Name[which(BDF$Tree_Name == "Mixornis_gularis")] <- "Macronouss_gularis"
BDF$Tree_Name[which(BDF$Tree_Name == "Mixornis_kelleyi")] <- "Macronous_kelleyi"

#Myiothlypis
chgen("Myiothlypis") %>% sort()
chgen("Basileuterus") %>% sort()
chgen("Phaeothlypis") %>% sort()
filter(BDF, genus == "Myiothlypis", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Basileuterus_basilicus")] <- "Myiothlypis_basilica"
BDF$Tree_Name[which(BDF$Tree_Name == "Myiothlypis_chlorophrys")] <- "Basileuterus_chlorophrys"
BDF$Tree_Name[which(BDF$Tree_Name == "Myiothlypis_roraimae")] <- "Basileuterus_bivittatus"

#Myrmelastes
chgen("Myrmelastes") %>% sort()
chgen("Schistocichla") %>% sort()
chgen("Myrmeciza") %>% sort()
filter(BDF, genus == "Myrmelastes", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Myrmelastes_brunneiceps")] <- "Schistocichla_brunneiceps"
BDF$Tree_Name[which(BDF$Tree_Name == "Myrmelastes_caurensis")] <- "Schistocichla_caurensis"
BDF$Tree_Name[which(BDF$Tree_Name == "Myrmelastes_humaythae")] <- "Schistocichla_humaythae"
BDF$Tree_Name[which(BDF$Tree_Name == "Myrmelastes_hyperythrus")] <- "Myrmeciza_hyperythra"
BDF$Tree_Name[which(BDF$Tree_Name == "Myrmelastes_leucostigma")] <- "Schistocichla_leucostigma"
BDF$Tree_Name[which(BDF$Tree_Name == "Myrmelastes_rufifacies")] <- "Schistocichla_rufifacies"
BDF$Tree_Name[which(BDF$Tree_Name == "Myrmelastes_saturatus")] <- "Schistocichla_saturata"
BDF$Tree_Name[which(BDF$Tree_Name == "Myrmelastes_schistaceus")] <- "Schistocichla_schistacea"

#Myrmoderus
chgen("Myrmoderus") %>% sort()
chgen("Myrmeciza") %>% sort()
filter(BDF, genus == "Myrmoderus", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Myrmoderus_ferrugineus")] <- "Myrmeciza_ferruginea"
BDF$Tree_Name[which(BDF$Tree_Name == "Myrmoderus_loricatus")] <- "Myrmeciza_loricata"
BDF$Tree_Name[which(BDF$Tree_Name == "Myrmoderus_ruficauda")] <- "Myrmeciza_ruficauda"
BDF$Tree_Name[which(BDF$Tree_Name == "Myrmoderus_squamosus")] <- "Myrmeciza_squamosa"

#Pachycephala
chgen("Pachycephala") %>% sort()
chgen("Colluricincla") %>% sort()
filter(BDF, genus == "Pachycephala", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Pachycephala_tenebrosa")] <- "Colluricincla_tenebrosa"

#Pachysylvia
chgen("Pachysylvia") %>% sort()
chgen("Hylophilus") %>% sort()
chgen("Vireo") %>% sort()
filter(BDF, genus == "Pachysylvia", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Pachysylvia_aurantiifrons")] <- "Hylophilus_aurantiifrons"
BDF$Tree_Name[which(BDF$Tree_Name == "Pachysylvia_decurtata")] <- "Hylophilus_decurtatus"
BDF$Tree_Name[which(BDF$Tree_Name == "Pachysylvia_hypochrysea")] <- "Vireo_hypochryseus"
BDF$Tree_Name[which(BDF$Tree_Name == "Pachysylvia_hypoxantha")] <- "Hylophilus_hypoxanthus"
BDF$Tree_Name[which(BDF$Tree_Name == "Pachysylvia_muscicapina")] <- "Hylophilus_muscicapinus"
BDF$Tree_Name[which(BDF$Tree_Name == "Pachysylvia_aurantiifrons")] <- "Hylophilus_aurantiifrons"
BDF$Tree_Name[which(BDF$Tree_Name == "Pachysylvia_semibrunnea")] <- "Hylophilus_semibrunneus"

#Phylloscopus
chgen("Phylloscopus") %>% sort()
chgen("Seicercus") %>% sort()
filter(BDF, genus == "Phylloscopus", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Phylloscopus_burkii")] <- "Seicercus_burkii"
BDF$Tree_Name[which(BDF$Tree_Name == "Phylloscopus_castaniceps")] <- "Seicercus_castaniceps"
BDF$Tree_Name[which(BDF$Tree_Name == "Phylloscopus_grammiceps")] <- "Seicercus_grammiceps"
BDF$Tree_Name[which(BDF$Tree_Name == "Phylloscopus_intermedius")] <- "Seicercus_affinis"
BDF$Tree_Name[which(BDF$Tree_Name == "Phylloscopus_montis")] <- "Seicercus_montis"
BDF$Tree_Name[which(BDF$Tree_Name == "Phylloscopus_omeiensis")] <- "Seicercus_omeiensis"
BDF$Tree_Name[which(BDF$Tree_Name == "Phylloscopus_poliogenys")] <- "Seicercus_poliogenys"
BDF$Tree_Name[which(BDF$Tree_Name == "Phylloscopus_soror")] <- "Seicercus_soror"
BDF$Tree_Name[which(BDF$Tree_Name == "Phylloscopus_tephrocephalus")] <- "Seicercus_tephrocephalus"
BDF$Tree_Name[which(BDF$Tree_Name == "Phylloscopus_valentini")] <- "Seicercus_valentini"
BDF$Tree_Name[which(BDF$Tree_Name == "Phylloscopus_whistleri")] <- "Seicercus_whistleri"

#Picoides
chgen("Picoides") %>% sort()
chgen("Dendrocopos") %>% sort()
filter(BDF, genus == "Picoides", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Picoides_canicapillus")] <- "Dendrocopos_canicapillus"
BDF$Tree_Name[which(BDF$Tree_Name == "Picoides_kizuki")] <- "Dendrocopos_kizuki"
BDF$Tree_Name[which(BDF$Tree_Name == "Picoides_maculatus")] <- "Dendrocopos_maculatus"
BDF$Tree_Name[which(BDF$Tree_Name == "Picoides_moluccensis")] <- "Dendrocopos_moluccensis"
BDF$Tree_Name[which(BDF$Tree_Name == "Picoides_nanus")] <- "Dendrocopos_nanus"
BDF$Tree_Name[which(BDF$Tree_Name == "Picoides_ramsayi")] <- "Dendrocopos_ramsayi"
BDF$Tree_Name[which(BDF$Tree_Name == "Picoides_temminckii")] <- "Dendrocopos_temminckii"

#Poecile
chgen("Poecile") %>% sort()
chgen("Parus") %>% sort()
filter(BDF, genus == "Poecile", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Poecile_davidi")] <- "Parus_davidi"
BDF$Tree_Name[which(BDF$Tree_Name == "Poecile_hypermelaenus")] <- "Parus_hypermelaenus"
BDF$Tree_Name[which(BDF$Tree_Name == "Poecile_lugubris")] <- "Parus_lugubris"
BDF$Tree_Name[which(BDF$Tree_Name == "Poecile_palustris")] <- "Parus_palustris"
BDF$Tree_Name[which(BDF$Tree_Name == "Poecile_superciliosus")] <- "Parus_superciliosus"

#Pogonornis
chgen("Pogonornis") %>% sort()
chgen("Lybius") %>% sort()
filter(BDF, genus == "Pogonornis", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Pogonornis_bidentatus")] <- "Lybius_bidentatus"
BDF$Tree_Name[which(BDF$Tree_Name == "Pogonornis_dubius")] <- "Lybius_dubius"
BDF$Tree_Name[which(BDF$Tree_Name == "Pogonornis_melanopterus")] <- "Lybius_melanopterus"
BDF$Tree_Name[which(BDF$Tree_Name == "Pogonornis_minor")] <- "Lybius_minor"
BDF$Tree_Name[which(BDF$Tree_Name == "Pogonornis_rolleti")] <- "Lybius_rolleti"

#Pogonotriccus
chgen("Pogonotriccus") %>% sort()
chgen("Phylloscartes") %>% sort()
filter(BDF, genus == "Pogonotriccus", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Pogonotriccus_chapmani")] <- "Phylloscartes_chapmani"
BDF$Tree_Name[which(BDF$Tree_Name == "Pogonotriccus_eximius")] <- "Phylloscartes_eximius"
BDF$Tree_Name[which(BDF$Tree_Name == "Pogonotriccus_lanyoni")] <- "Phylloscartes_lanyoni"
BDF$Tree_Name[which(BDF$Tree_Name == "Pogonotriccus_ophthalmicus")] <- "Phylloscartes_ophthalmicus"
BDF$Tree_Name[which(BDF$Tree_Name == "Pogonotriccus_orbitalis")] <- "Phylloscartes_orbitalis"
BDF$Tree_Name[which(BDF$Tree_Name == "Pogonotriccus_poecilotis")] <- "Phylloscartes_poecilotis"
BDF$Tree_Name[which(BDF$Tree_Name == "Pogonotriccus_venezuelanus")] <- "Phylloscartes_venezuelanus"

#Poodytes
chgen("Poodytes") %>% sort()
chgen("Megalurus") %>% sort()
chgen("Eremiornis") %>% sort()
filter(BDF, genus == "Poodytes", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Poodytes_albolimbatus")] <- "Megalurus_chapmani"
BDF$Tree_Name[which(BDF$Tree_Name == "Poodytes_carteri")] <- "Eremiornis_carteri"
BDF$Tree_Name[which(BDF$Tree_Name == "Poodytes_gramineus")] <- "Megalurus_gramineus"

#Psittacara
chgen("Psittacara") %>% sort()
chgen("Heliolais") %>% sort()
chgen("Spiloptila") %>% sort()
filter(BDF, genus == "Psittacara", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Prinia_erythroptera")] <- "Heliolais_erythropterus"
BDF$Tree_Name[which(BDF$Tree_Name == "Prinia_rufifrons")] <- "Spiloptila_rufifrons"

#Psilopogon
chgen("Psilopogon") %>% sort()
chgen("Aratinga") %>% sort()
filter(BDF, genus == "Psilopogon", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Psittacara_acuticaudatus")] <- "Aratinga_acuticaudatus"
BDF$Tree_Name[which(BDF$Tree_Name == "Psittacara_chloropterus")] <- "Aratinga_chloroptera"
BDF$Tree_Name[which(BDF$Tree_Name == "Psittacara_erythrogenys")] <- "Aratinga_erythrogenys"
BDF$Tree_Name[which(BDF$Tree_Name == "Psittacara_euops")] <- "Aratinga_euops"
BDF$Tree_Name[which(BDF$Tree_Name == "Psittacara_finschi")] <- "Aratinga_finschi"
BDF$Tree_Name[which(BDF$Tree_Name == "Psittacara_frontatus")] <- "Aratinga_frontatus"
BDF$Tree_Name[which(BDF$Tree_Name == "Psittacara_holochlorus")] <- "Aratinga_holochlora"
BDF$Tree_Name[which(BDF$Tree_Name == "Psittacara_leucophthalmus")] <- "Aratinga_leucophthalma"
BDF$Tree_Name[which(BDF$Tree_Name == "Psittacara_mitrata")] <- "Aratinga_mitrata"
BDF$Tree_Name[which(BDF$Tree_Name == "Psittacara_rubritorquis")] <- "Aratinga_rubritorquis"
BDF$Tree_Name[which(BDF$Tree_Name == "Psittacara_wagleri")] <- "Aratinga_wagleri"

#Pteroglossus
chgen("Pteroglossus") %>% sort()
chgen("Aratinga") %>% sort()
filter(BDF, genus == "Pteroglossus", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Pteroglossus_beauharnaisii")] <- "Pteroglossus_beauharnaesii"

#Ramphiculus 
chgen("Ramphiculus") %>% sort()
chgen("Ptilinopus") %>% sort()
filter(BDF, genus == "Ramphiculus", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Ramphiculus_jambu")] <- "Ptilinopus_jambu"
BDF$Tree_Name[which(BDF$Tree_Name == "Ramphiculus_marchei")] <- "Ptilinopus_marchei"
BDF$Tree_Name[which(BDF$Tree_Name == "Ramphiculus_leclancheri")] <- "Ptilinopus_leclancheri"
BDF$Tree_Name[which(BDF$Tree_Name == "Ramphiculus_merrilli")] <- "Ptilinopus_merrilli"
BDF$Tree_Name[which(BDF$Tree_Name == "Ramphiculus_occipitalis")] <- "Ptilinopus_occipitalis"

#Schoeniparus
chgen("Schoeniparus") %>% sort()
chgen("Alcippe") %>% sort()
filter(BDF, genus == "Schoeniparus", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Schoeniparus_brunneus")] <- "Alcippe_brunnea"
BDF$Tree_Name[which(BDF$Tree_Name == "Schoeniparus_castaneceps")] <- "Alcippe_castaneceps"
BDF$Tree_Name[which(BDF$Tree_Name == "Schoeniparus_cinereus")] <- "Alcippe_cinerea"
BDF$Tree_Name[which(BDF$Tree_Name == "Schoeniparus_dubius")] <- "Alcippe_dubia"
BDF$Tree_Name[which(BDF$Tree_Name == "Schoeniparus_klossi")] <- "Alcippe_klossi"
BDF$Tree_Name[which(BDF$Tree_Name == "Schoeniparus_rufogularis")] <- "Alcippe_rufogularis"
BDF$Tree_Name[which(BDF$Tree_Name == "Schoeniparus_variegaticeps")] <- "Alcippe_variegaticeps"

#Scleroptila
chgen("Scleroptila") %>% sort()
chgen("Francolinus") %>% sort()
filter(BDF, genus == "Scleroptila", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Scleroptila_afra")] <- "Francolinus_africanus"
BDF$Tree_Name[which(BDF$Tree_Name == "Scleroptila_streptophora")] <- "Francolinus_streptophora"

#Sittiparus
chgen("Sittiparus") %>% sort()
chgen("Parus") %>% sort()
filter(BDF, genus == "Sittiparus", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Sittiparus_semilarvatus")] <- "Parus_semilarvatus"

#Sporophila
chgen("Sporophila") %>% sort()
chgen("Oryzoborus") %>% sort()
filter(BDF, genus == "Sporophila", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Sporophila_atrirostris")] <- "Oryzoborus_atrirostris"
BDF$Tree_Name[which(BDF$Tree_Name == "Sporophila_funereus")] <- "Oryzoborus_funerea"

#Sylvia
chgen("Sylvia") %>% sort()
chgen("Horizorhinus") %>% sort()
chgen("Lioptilus") %>% sort()
filter(BDF, genus == "Sylvia", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Sylvia_dohrni")] <- "Horizorhinus_dohrni"
BDF$Tree_Name[which(BDF$Tree_Name == "Sylvia_nigricapillus")] <- "Lioptilus_nigricapillus"
BDF$Tree_Name[which(BDF$Tree_Name == "Sylvia_subcoerulea")] <- "Sylvia_subcaerulea"

#Tangara
chgen("Tangara") %>% sort()
chgen("Thraupis") %>% sort()
filter(BDF, genus == "Tangara", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Tangara_argentea")] <- "Thraupis_cyanoptera"

#Trochalopteron
chgen("Trochalopteron") %>% sort()
chgen("Garrulax") %>% sort()
filter(BDF, genus == "Trochalopteron", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Trochalopteron_chrysopterum")] <- "Garrulax_chrysopterus"
BDF$Tree_Name[which(BDF$Tree_Name == "Trochalopteron_imbricatum")] <- "Garrulax_imbricatus"
BDF$Tree_Name[which(BDF$Tree_Name == "Trochalopteron_melanostigma")] <- "Garrulax_melanostigma"
BDF$Tree_Name[which(BDF$Tree_Name == "Trochalopteron_ngoclinhense")] <- "Garrulax_ngoclinhensis"
BDF$Tree_Name[which(BDF$Tree_Name == "Trochalopteron_peninsulae")] <- "Garrulax_peninsulae"
BDF$Tree_Name[which(BDF$Tree_Name == "Trochalopteron_variegatum")] <- "Garrulax_variegatus"

#Turdinus
chgen("Turdinus") %>% sort()
chgen("Gypsophila") %>% sort()
chgen("Napothera") %>% sort()
filter(BDF, genus == "Turdinus", inTree == FALSE)$Tree_Name
BDF$Tree_Name[which(BDF$Tree_Name == "Turdinus_brevicaudatus")] <- "Napothera_brevicaudata"
BDF$Tree_Name[which(BDF$Tree_Name == "Turdinus_brevicaudatus")] <- "Thraupis_cyanoptera"
BDF$Tree_Name[which(BDF$Tree_Name == "Turdinus_crassus")] <- "Napothera_crassa"
BDF$Tree_Name[which(BDF$Tree_Name == "Turdinus_marmoratus")] <- "Turdinus_marmorata"

#Zosterops
chgen("Zosterops") %>% sort()
view(filter(BDF, genus == "Zosterops", inTree == FALSE))
BDF$Tree_Name[which(BDF$Tree_Name == "Zosterops_chloronothos")] <- "Zosterops_chloronothus"

#Zapornia
chgen("Zapornia")
chgen("Amaurornis") %>% sort()
chgen("Porzana") %>% sort()
view(filter(BDF, genus == "Zapornia", inTree == FALSE))
BDF$Tree_Name[which(BDF$Tree_Name == "Zapornia_akool")] <- "Amaurornis_akool"
BDF$Tree_Name[which(BDF$Tree_Name == "Zapornia_atra")] <- "Porzana_atra"
BDF$Tree_Name[which(BDF$Tree_Name == "Zapornia_bicolor")] <- "Amaurornis_bicolor"
BDF$Tree_Name[which(BDF$Tree_Name == "Zapornia_flavirostra")] <- "Amaurornis_flavirostra"
BDF$Tree_Name[which(BDF$Tree_Name == "Zapornia_fusca")] <- "Porzana fusca"
BDF$Tree_Name[which(BDF$Tree_Name == "Zapornia_olivieri")] <- "Amaurornis_olivieri"
BDF$Tree_Name[which(BDF$Tree_Name == "Zapornia_parva")] <- "Porzana_parva"
BDF$Tree_Name[which(BDF$Tree_Name == "Zapornia_paykullii")] <- "Porzana_paykullii"
BDF$Tree_Name[which(BDF$Tree_Name == "Zapornia_pusilla")] <- "Porzana_pusilla"
BDF$Tree_Name[which(BDF$Tree_Name == "Zapornia_tabuensis")] <- "Porzana_tabuensis"
BDF$Tree_Name[which(BDF$Tree_Name == "Zapornia_tabuensis")] <- "Porzana_tabuensis"


#add "inTree" column
BDF <- BDF %>% mutate(inTree = sapply(Tree_Name, function(x){x %in% phy[[1]]$tip.label }), .after = Tree_Name)

saveRDS(BDF, file = "FossilBirds.rds")
write.csv(BDF, file = "FossilBirds.csv")
