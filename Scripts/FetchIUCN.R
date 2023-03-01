setwd("~/Fossil Record")
library(rredlist)

#From rredlist example code - download entire database
out <- rl_sp(all = TRUE, key = key)
length(out)
vapply(out, "[[", 1, "count")
IUCNcategories <- do.call(rbind, lapply(out, "[[", "result"))

save(IUCNcategories, file = "IUCNcategories.rda")
