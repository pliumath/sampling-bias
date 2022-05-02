library(ape)
library(diversitree)
library(phangorn)
library(phyloTop)
library(ggplot2)
library(reshape2)
library(scales)
library(egg)
library(gridBase)
library(grid)
library(latex2exp)

##########################################
##########################################
#load true trees

load("truetree01.Rdata")

h <- history.from.sim.discrete(phy, 0:1)

plot(h, phy,show.node.label = F,show.tip.label = F, show.tip.state = F, label.offset = 0.1, cex=1 ,cols = c('#EDB120','#0072BD'))

n0 <- sum(phy$tip.state == 0)
n1 <- sum(phy$tip.state == 1)

N0 <- sum(phy$node.state == 0)
N1 <- sum(phy$node.state == 1)


##########################################
##########################################
#estimates with subsampled tips

times <- 50 #downsampling times
k <- seq(0.05,0.95,0.1) #the location-A proportion
l <- length(k)

ntips <- 100 #number of tips in the downsapled trees

n0 <- sum(phy$tip.state == 0)
n1 <- sum(phy$tip.state == 1)

mrate <- matrix(0,nrow = times,ncol = l)

for (i in 1:l) {
  
  print(i)
  
  for (j in 1:times) {
    
    print(j)
    
    to.remove <- rep(FALSE,phy$Nnode+1) #preallocate tips to remove
    
    sn0 <- floor(k[i]*ntips) #number of selected location-A tips with a given proportion
    sn1 <- ntips - sn0
    
    rand.remove0 <- sample(which(phy$tip.state == 0),n0-sn0) #sample tips to remove
    rand.remove1 <- sample(which(phy$tip.state == 1),n1-sn1)
    to.remove[rand.remove0] <- TRUE
    to.remove[rand.remove1] <- TRUE
    
    downsampled <- prune(phy,to.remove)
    
    sh <- history.from.sim.discrete(downsampled, 0:1)
    dans <- ace(downsampled$tip.state, downsampled, type = "d") #reconstruct ancestral states on the downsampled tree
    
    mrate[j,i] <- dans$rates
    
  }
  
}


##########################################
##########################################
##visualization (Fig 3)


mr4plot <- as.data.frame(mrate)

mr4plot  <- melt(mr4plot)
mr4plot$variable <- sapply(rep(k, each=times),as.character)

p <- ggplot(mr4plot, aes(x=variable, y=value)) + 
  geom_jitter(shape=16, position=position_jitter(0.1),size=0.75) + 
  geom_vline(aes(xintercept = n0/(phy$Nnode+1)*(l+1)), color = '#000000', show.legend = NA) + 
  geom_hline(aes(yintercept = 0.1), color = '#000000', show.legend = NA) + 
  geom_violin(alpha =0.5) + 
  labs(x = "Proportion of location-A (amber) tips in a downsampled tree", y = "Migration rate") +
  scale_y_continuous(trans='log10',limits = c(0.01,100)) +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle(TeX("Migration rate $\\alpha = 0.1$")) + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,width = unit(16, "cm"),height = unit(9, "cm"))

grid.draw(tp)
