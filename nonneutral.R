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
#functions

nullmodaccu <- function(downsampled,dans){
  
  N <- downsampled$Nnode 
  p = sum(downsampled$tip.state == 1)/(N+1)
  
  results <- c()
  
  accuracy <-  cbind(abs(dans$lik.anc[,2] - rep(p,N)),1-abs(dans$lik.anc[,2] - rep(p,N)))
  
  return(mean(accuracy[,2])) 
  
}

##########################################
##########################################
#generate true trees

flag <- 0
while (flag == 0) {

  numcands <- 5

  ntips <- 100 #number of tips in the downsapled trees

  lmd <- 4 #birth rate
  mu <- 1 #death rate
  alph <- 0.3 #migration rate

  pars <- c(lmd+4, lmd, mu, mu, alph, alph) #pars: birth rate 0, birth rate 1, death rate 0, death rate 1, migration rate 0->1, migration rate 1->0

  phys <- trees(pars, type="bisse", n=numcands, max.taxa= 8*lmd*ntips, max.t=Inf, include.extinct=TRUE) #generate a truth tree candidate

  for (i in 1:length(phys)) {

    cand <- phys[[i]]

    if(!is.null(cand) && cand$Nnode > ntips) {

      to.drop <- startsWith(cand$tip.label,"sp") #drop the extant tips
      cand <- prune(cand,to.drop)

      rootstate <- cand$node.state[1] #swap node states if the root is blue

      n0 <- sum(cand$tip.state == 0) #number of location-A tips
      n1 <- sum(cand$tip.state == 1) #number of location-B tips

      if (n0 >= ntips & n1 >= ntips & rootstate == 0) {

        phy <- cand
        flag <- 1
        break

      }

    }

  }

}

h <- history.from.sim.discrete(phy, 0:1)

plot(h, phy,show.node.label = F,show.tip.label = F, show.tip.state = F, label.offset = 0.1, cex=1 ,cols = c('#EDB120','#0072BD'))

n0 <- sum(phy$tip.state == 0)
n1 <- sum(phy$tip.state == 1)

N0 <- sum(phy$node.state == 0)
N1 <- sum(phy$node.state == 1)


##########################################
##########################################
#estimates with subsampled tips

times <- 50 #number of downsampled trees
k <- seq(0.05,0.95,0.1) #proportion of location-A (amber) tips
l <- length(k)

ntips <- 100 #number of tips in the downsapled trees

meanaccu <- matrix(0,nrow = times,ncol = l) #preallocating
adjaccu <- matrix(0,nrow = times,ncol = l)

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
    
    
    accuracy <-  cbind(abs(dans$lik.anc[,2] - downsampled$node.state),1-abs(dans$lik.anc[,2] - downsampled$node.state))
    meanaccu[j,i] <- mean(accuracy[,2]) #absolute accuracy
    
    adjaccu[j,i] <-  mean(accuracy[,2]) - nullmodaccu(downsampled)  #relative accuracy
    
  }
  
}



##########################################
##########################################
##visualization


ma4p <- as.data.frame(meanaccu)
ma4p <- melt(ma4p)
ma4p$variable <- sapply(rep(k, each=times),as.character)
ad4p <- as.data.frame(adjaccu)
ad4p <- melt(ad4p)
ad4p$variable <- sapply(rep(k, each=times),as.character)
ma4p <- rbind(ma4p,ad4p)
ma4p[,3] <- c(rep("mean",times*l), rep("rel",times*l))
names(ma4p) <- c("prop","accu","type")

type.labs <- c("Absolute accuracy","Relative accuracy")
names(type.labs) <- c("mean","rel")

p <- ggplot(ma4p, aes(x=prop, y=accu)) + 
  geom_jitter(shape=16, position=position_jitter(0.1),size=0.75) + 
  geom_vline(aes(xintercept = n0/(phy$Nnode+1)*(l+1)), color = '#000000', show.legend = NA) + 
  geom_violin(alpha =0.5) + 
  facet_grid(type~., labeller = labeller(type = type.labs)) +
  labs(x = "Proportion of location-A (amber) tips in a downsampled tree", y = "Accuracy") +
  ylim(0,1) +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Absolute and relative accuracy of downsampled trees \nfor a true tree with non-neutral speciation rates (location-A lineages branch faster)") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,width = unit(16, "cm"),height = unit(9, "cm"))

grid.draw(tp)




p <- ggplot(ma4p, aes(x=prop, y=accu)) + 
  geom_jitter(shape=16, position=position_jitter(0.1),size=0.75) + 
  geom_vline(aes(xintercept = n0/(phy$Nnode+1)*(l+1)), color = '#000000', show.legend = NA) + 
  geom_hline(aes(yintercept = 0), color = '#000000', show.legend = NA) + 
  geom_violin(alpha =0.5) + 
  facet_grid(type~., labeller = labeller(type = type.labs)) +
  labs(x = "Proportion of location-A (amber) tips in a downsampled tree", y = "Accuracy") +
  ylim(-0.5,0.5) +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Absolute and relative accuracy of downsampled trees \nfor a true tree with non-neutral speciation rates (location-A lineages branch faster)") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,width = unit(16, "cm"),height = unit(9, "cm"))

grid.draw(tp)


