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

##########################################
##########################################
#genarate a true tree

flag <- 0
while (flag == 0) {
  
  numcands <- 5
  
  ntips <- 100 #number of tips in the downsapled trees
  
  lmd <- 4 #birth rate
  mu <- 1 #death rate
  alph <- 0.9 #migration rate low-0.1,0.3, high-0.7,0.9
  
  pars <- c(lmd, lmd, mu, mu, alph, alph) #pars: birth rate 0, birth rate 1, death rate 0, death rate 1, migration rate 0->1, migration rate 1->0
  
  phys <- trees(pars, type="bisse", n=numcands, max.taxa= 4*lmd*ntips, max.t=Inf, include.extinct=TRUE) #generate a truth tree candidate
  
  for (i in 1:length(phys)) {
    
    cand <- phys[[i]]
    
    if(!is.null(cand) && cand$Nnode > ntips) {
      
      to.drop <- startsWith(cand$tip.label,"sp") #drop the extant tips
      cand <- prune(cand,to.drop)
      
      rootstate <- cand$node.state[1] #swap node states if the root is blue
      
      n0 <- sum(cand$tip.state == 0) #number of location-A tips
      n1 <- sum(cand$tip.state == 1) #number of location-B tips
      
      if (n0 >= ntips & n1 >= ntips & n0>=n1 & rootstate == 0) {
        
        phy <- cand
        flag <- 1
        break
        
      }
      
    }
    
  }
  
}


h <- history.from.sim.discrete(phy, 0:1)

plot(h, phy,show.node.label = F,show.tip.label = F, show.tip.state = F, label.offset = 0.1, cex=1 ,cols = c('#EDB120','#0072BD'))

sum(phy$tip.state == 0)
sum(phy$tip.state == 1)

sum(phy$tip.state == 0) / (sum(phy$tip.state == 0) + sum(phy$tip.state == 1))

sum(phy$node.state == 0)
sum(phy$node.state == 1)

##########################################
##########################################
#estimates with downsampled trees


times <- 50 #number of downsampled trees
k <- seq(0.05,0.95,0.1) #proportion of location-A (amber) tips
l <- length(k)


rtnd <- matrix("",nrow = times,ncol = l)
rtstacc <- matrix(0,nrow = times,ncol = l)


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
    
    rtnd[j,i] <- downsampled$node.label[1]
    rtstacc[j,i] <- dans$lik.anc[1,1]
    
  }
  
}


##########################################
##########################################
##visualization


rootdata <- as.data.frame(rtstacc)
rootdata  <- melt(rootdata)
rootdata$variable <- sapply(rep(k, each=times),as.character)

rootnode <- as.data.frame(rtnd)
rootnode <- melt(rootnode,0)

rootdata[,3] <- rootnode$value
names(rootdata) <- c("variable","value","rootnode")

p <- ggplot(rootdata, aes(x=variable, y=value)) + 
  geom_jitter(shape=16, position=position_jitter(0.1),size=2.5, alpha = 0.5,aes(color=rootnode)) + geom_vline(aes(xintercept = n0/(phy$Nnode+1)*(l+1)), color = '#000000', show.legend = NA) + 
  scale_color_manual(values=c('#000000', '#D95319', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'))+
  labs(x = "Proportion of location-A (amber) tips in a downsampled tree", y = "Accuracy", color = "Root Node") +
  ylim(0,1) +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Reconstructed root node and absolute accuracy of the root in downsampled trees") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12)) + theme(legend.position = "bottom")

tp <- set_panel_size(p,width  = unit(16, "cm"),height = unit(9, "cm"))

grid.draw(tp)


