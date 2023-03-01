#phylogenetic signal

library('phytools')
bird <- read.tree("Phylogenies/Bird/Bird_consensus.tre")
reptile <- read.tree("Phylogenies/Squamates/Squamates_consensus.tre")
amphibian <- read.tree("Phylogenies/Amphibians/Amphibians_consensus.tre")
mammal <- read.tree("Phylogenies/Mammal/Mammal_consensus.tre")


WilksLambda <- function(phy,taxon = NA,df){
  if (!is.na(taxon)){
    df <- filter(df, group == taxon)
  }
  
  f <- phylosig(phy, x = setNames(df$fossilArea,df$Tree_Name)[phy$tip.label],
                method = "lambda")
  e <- phylosig(phy, x = setNames(df$Area,df$Tree_Name)[phy$tip.label],
                method = "lambda") 
  return(setNames(c(f$lambda,f$logL,e$lambda,e$logL),
                  nm =  c("FGR Lambda", "FGR Log Likelihood",
                          "Geographic Range Lambda","Geographic Range Log Likelihood")))
}

a <- WilksLambda(bird,df=BDF)
b <- WilksLambda(bird,taxon = "Passerine", df=BDF)
c <- WilksLambda(bird,taxon = "Nonpasserine", df=BDF)
      
d <- WilksLambda(reptile,df=GARD)
g <- WilksLambda(reptile,taxon = "lizard",df=GARD)
h <- WilksLambda(reptile,taxon = "snake",df=GARD)
      
i <- WilksLambda(amphibian,df=ADF)
j <- WilksLambda(amphibian,taxon = "frog",df=ADF)
k <- WilksLambda(amphibian,taxon = "salamander",df=ADF)
l <- WilksLambda(amphibian,taxon = "caecilian",df=ADF)
      
m <- WilksLambda(mammal,df=Mammals)
n <- WilksLambda(mammal,taxon = "Terrestrial",df=Mammals)
o <- WilksLambda(mammal,taxon = "Marine", df=Mammals)


Lambdas <- rbind(a,b,c,d,g,h,i,j,k,l,m,n,o)

rownames(Lambdas) <- c("Bird","Passerine","Non-Passerine","Squamate","Lizard","Snake",
                      "Amphibian","Frog","Salamander","Caecilian","Mammal","Terrestrial Mammal",
                      "Marine Mammal")
write.csv(Lambdas, file ="Results/PhylogeneticSignalTable.csv")                 
