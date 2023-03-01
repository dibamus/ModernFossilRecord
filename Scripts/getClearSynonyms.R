#This code uses the taxize package and the ITIS database to check for synonyms
# of the species names found in a data frame and cross-reference them with 
# the species names in a tree
get_clear_synonyms <- function(df,phy){
  require("phytools")
  require("taxize")
  require("tidyverse")

  df_unique <- df$Tree_Name[!(df$Tree_Name %in% phy[[1]]$tip.label)]
  phy_unique <- phy[[1]]$tip.label[!(phy[[1]]$tip.label %in% df$Tree_Name)][-1]
  
  get_syn <- function(x){
    t <- synonyms(x, db = "itis")
    if(is.na(t)){
      return("not found")
    }
    
    if(nrow(t[[1]])==0){
      return(x)
    }
    else{
      acc <- NULL
      if("acc_name" %in% colnames(t[[1]])){
        acc<- t[[1]]$acc_name %>% unique() %>% str_replace(" ","_")
      }
      return(c(acc,t[[1]]$syn_name %>% unlist() %>% str_replace(" ","_") %>% unique()))
    }
    
  }
  
  df_unique_syn <- sapply(df_unique, get_syn)  #select the name that matches the input and is valid
  
  
  syn_found <- function(x, source, ref){
    source[x][[1]][which(source[x][[1]] %in% ref)] 
  }
  df_unique_syn_found <- sapply(names(df_unique_syn), syn_found, 
                                source = df_unique_syn, ref = phy_unique)
  
  df_replacements <- df_unique_syn_found[which(sapply(df_unique_syn_found, length) == 1)]
  
  all_df_names <- df$Tree_Name
  
  match_syn <- function(x) {
    all_df_names[match(x,all_df_names)] <- df_replacements[x][[1]]
  }
  
  match_syn <- function(x) {
    if(!(x %in% names(df_replacements))){
      return(x)
    }
    else{
      return(df_replacements[x][[1]])
    }
  }
  
  namematches <- data.frame("original" = df$Tree_Name,"replacement" = sapply(df$Tree_Name, FUN = match_syn))
  
  unclear <- df_unique_syn_found[which(sapply(df_unique_syn_found, length) > 1)]
  
  nosyn <- df_unique_syn_found[which(sapply(df_unique_syn_found, length) == 0)]
  
  return(list(namematches = namematches,
         unclear = unclear,
         nosyn = nosyn))
}
