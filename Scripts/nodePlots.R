#Plots nodedepth data for a diven dataset (ds)
nodeDepthDF <- function(ds){
  require('reshape2')
  require('tidyverse')
  require('ggplot2')
  require('forcats')
  obsdata <- readRDS(paste0("Results/",ds,"/",ds,"NodeDepthsObserved.rds"))
  simdata <- readRDS(paste0("Results/",ds,"/",ds,"HistorySims.rds"))$depths
  
  obsdata_means <- colMeans(obsdata)
  simdata_means <- rowMeans(simdata[,,,],dims = 3) %>% colMeans()
  
  
  df<- rbind(melt(obsdata_means, as.is = TRUE) %>% mutate(dataset = "observed"),
             melt(simdata_means) %>% mutate(dataset = "simulated"),
             melt(obsdata_means) %>% mutate(value = value - melt(simdata_means)$value, dataset = "difference")) %>%
    mutate(node_age = fct_rev(X2))
   
  return(df)
}
#get DF of all node depth measurements
df <- rbind(nodeDepthDF("Amphibian") %>% mutate(taxon = "Amphibian"),
         nodeDepthDF("Reptile") %>% mutate(taxon = "Squamate"),
         nodeDepthDF("Mammal") %>% mutate(taxon = "Mammal"),
         nodeDepthDF("Bird") %>% mutate(taxon = "Bird"))

df$taxon <- factor(df$taxon, levels = c("Bird","Squamate","Amphibian","Mammal"))

levels(df$node_age) <- 10^seq(from =2.75, to =-3, by = -0.25)
levels(df$X1) <- c("All","1","10","100","1000","10,000","100,000")
colnames(df)[1] <- "Minimum FGR size (km^2)"

df$yearslost <- df$value*as.numeric(df$node_age)


df$node_age <- as.numeric(as.character(df$node_age))

geointervals <- data.frame("start" = c(0.01,2.6,23,66,201.3),
                           "end" = c(2.6,23,66,201.3,251.9),
                           "Period" = factor(c("Quaternary","Neogene","Paleogene","Cretaceous","Jurassic"),
                                           levels = c("Quaternary","Neogene","Paleogene","Cretaceous","Jurassic")))

helpfultext <- data.frame(x = 0.5, y = c(50, -50),
                          taxon = "Bird",
                          label = c("better node recovery", "worse node recovery"),
                          label2 = c("better PD recovery","wose PD recovery"))
geoLabels <- geointervals %>%
  mutate(taxon = "Bird") %>%
  filter( Period != "Jurassic")
helpfultext$taxon <- factor(helpfultext$taxon, levels = c("Bird","Squamate","Amphibian","Mammal"))
geoLabels$taxon <- factor(geoLabels$taxon, levels = c("Bird","Squamate","Amphibian","Mammal"))


library('scales')  
depthdifferences <- function(yvar) {
  yvar <- sym(yvar)
  ggplot() +
  geom_hline(yintercept = 0) +
    geom_rect(data = geointervals,
              aes(xmin = start,
                  xmax = end,
                  ymin = 0,
                  ymax = Inf,
                  fill = Period),
              alpha = 0.45)  +
    geom_rect(data = geointervals,
              aes(xmin = start,
                  xmax = end,
                  ymin = -Inf,
                  ymax = 0,
                  fill = Period),
              alpha = 0.25) +
  

  scale_fill_manual(values = c("#fbee5a", "#f6cc5e", "#f1a962", "#ec8766", "#e7656a")) +
  
  geom_line(data = df %>% filter(dataset == "difference", 
                                 `Minimum FGR size (km^2)` != "All",
                                 node_age >= 10^-2), 
            aes(x = node_age, 
                y = !!yvar, #change back to "value" to get node number
                color = `Minimum FGR size (km^2)`, 
                group = `Minimum FGR size (km^2)`),
            size = 1.5) +
  # demarcate Geologic Periods
  geom_vline(xintercept = geointervals$end, linetype = "dashed", color="#999999", size=1) +
  
  scale_x_continuous(trans  = scales::compose_trans("log10", "reverse"),
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)),
                     limits = c(251.9,0.01),
                     expand = c(0, 0))+
  xlab("node age (Ma)") +
  facet_wrap(.~taxon)
}

#generate node number plot
depthdifferences("value") +
  
  #Text annotations
  geom_text(data = helpfultext, size =9, aes(x,y,label=label),hjust = "left") +
  geom_text(data = geoLabels, size=9,
            aes(y = -100, x = end - end*0.2, label = Period),
            angle = 90,
            hjust = "left") +
  
  guides(fill = "none") +
  ylab("Difference of nodes recovered by geographically structured VS random preservation") +
  
  ggtitle(paste("Geographic structuring of phylogenetic node recovery potential"))

#generate history lost plot
depthdifferences("yearslost") +
  
  geom_text(data = helpfultext, size =9, aes(x,y*10,label=label2),hjust = "left") +
  geom_text(data = geoLabels, size =9,
            aes(y = -1000, x = end - end*0.2, label = Period),
            angle = 90,
            hjust = "left") +
  guides(fill = "none") +
  ylab("Difference in phylogenetic diversity (Ma) recovered by geographically structured VS random preservation") +
  
  ggtitle(paste("Geographic structuring of phylogenetic diversity recovery potential"))

#for publication as figure 5
fig5 <- depthdifferences("yearslost") +  
  geom_text(data = helpfultext, size =9 * (5/14), aes(x,y*10,label=label2),hjust = "left") +
  geom_text(data = geoLabels, size =9 * (5/14),
            nudge_x = 0.1,
            aes(y = -1000, x = end - end*0.2, label = Period),
            angle = 90,
            hjust = "left") +
  guides(fill = "none")+
  ylab(expression(paste("\u394"," phylogenetic diversity (Ma) recovered"))) +
  theme(strip.text.x.top = element_blank(),
        
        legend.title = element_blank())

ggsave("Figures/Figure_4_1col.pdf", plot=fig5, width = 15, height = 10, units="cm")
