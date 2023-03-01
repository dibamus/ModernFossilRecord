matchIUCNtoFossil <- function(fileName){

  loadRData <- function(rda.file){ #got this code from stackoverflow https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
    #load up the .rda fossil data so that it can be renamed in endangered.fossil
    load(rda.file)
    get(ls()[ls() != "rda.file"])
  }
  load("IUCNcategories.rda")
  #List which species polygons have fossil record
  ls <-apply(loadRData(fileName), 1, any)
  
  #get a list of species names
  spnames <- sub("^(\\S*\\s+\\S+).*", "\\1", names(ls))
  
  #make df
  end.foss <- as.data.frame(ls)
  end.foss$binom <- spnames
  
  #find unique binomials
  binom <- levels(as.factor(spnames))
  fossil <- sapply(binom,
                   function(spname){
                   any(end.foss$ls[which(end.foss$binom == spname)])
                   }
  )
  REDcat <- allnames$cat[which(allnames$allnames %in% binom)]
  names(REDcat) <- allnames$allnames[which(allnames$allnames %in% binom)]
  
  sp.foss <- data.frame(
    fossil = fossil,
    REDcat = rep("DD",length(binom)), #initialize all as "DD" - some do not match spatial data (taxonomic changes)
    row.names = binom
  )
  
  sp.foss[names(REDcat),2] <- REDcat
  
  #return dataframe
  sp.foss
}
