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
#generate true trees

KK <- 25 #number of true trees

trees <- list()

for (kk in 1:KK) {
  
  print(kk)
  
  flag <- 0
  while (flag == 0) {
    
    numcands <- 25
    
    ntips <- 100 #number of tips in the downsapled trees
    
    lmd <- 4 #birth rate
    mu <- 1 #death rate
    alph <- 0.9 #migration rate
    
    pars <- c(lmd, lmd, mu, mu, alph, alph) #pars: birth rate 0, birth rate 1, death rate 0, death rate 1, migration rate 0->1, migration rate 1->0
    
    phys <- trees(pars, type="bisse", n=numcands, max.taxa= 3*lmd*ntips, max.t=Inf, include.extinct=TRUE) #generate a truth tree candidate
    
    for (i in 1:length(phys)) {
      
      cand <- phys[[i]]
      
      if(!is.null(cand) && cand$Nnode > ntips) {
        
        to.drop <- startsWith(cand$tip.label,"sp") #drop the extant tips
        cand <- prune(cand,to.drop)
        
        rootstate <- cand$node.state[1] #swap node states if the root is location-B (blue)
        
        n0 <- sum(cand$tip.state == 0) #number of location-A tips
        n1 <- sum(cand$tip.state == 1) #number of location-B tips
        
        if (n0 >= ntips & n1 >= ntips & n0>=n1 & rootstate == 0) {
          
          trees[[kk]] <- cand
          flag <- 1
          break
          
        }
        
      }
      
    }
    
  }
  
}


##########################################
##########################################
#estimates with downsampled trees

times <- 50 #number of downsampled trees
k <- seq(0.05,0.95,0.1) #proportion of location-A (amber) tips
l <- length(k)

ntips <- 100  #number of tips in the downsapled trees

Accu <- matrix(0,nrow = KK,ncol = l)
Node <- matrix(0,nrow = KK,ncol = 4)

for (kk in 1:KK) {
  
  print(kk)
  
  phy <- trees[[kk]]
  
  h <- history.from.sim.discrete(phy, 0:1)
  
  plot(h, phy,show.node.label = F,show.tip.label = F, show.tip.state = F, label.offset = 0.1, cex=1 ,cols = c('#EDB120','#0072BD'))
  
  n0 <- sum(phy$tip.state == 0)
  n1 <- sum(phy$tip.state == 1)
  
  N0 <- sum(phy$node.state == 0)
  N1 <- sum(phy$node.state == 1)
  
  Node[kk,] <- c(n0,n1,N0,N1)
  
  meanaccu <- matrix(0,nrow = times,ncol = l)
  
  for (i in 1:l) {
    
    print(i)
    
    for (j in 1:times) {
      
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
      
    }
    
  }
  
  Accu[kk,] <- colMeans(meanaccu)
  
}


##########################################
##########################################
##visualization (Fig S2)

data4j <- as.data.frame(Accu)
data4j <- melt(data4j)
data4j$variable <- rep(k, each=KK)

pr <- Node[,1]/(Node[,1]+Node[,2])
tst <- t.test(pr)
tst <- c(tst$conf.int[1],tst$estimate,tst$conf.int[2])
pr4p <- as.data.frame(pr)
pr4p[,2] <- rep(0,length(pr))
names(pr4p) <- c("x","y")

p <- ggplot(data4j,aes(x = variable, y= value)) + 
  geom_vline(aes(xintercept = tst[2]), color = '#000000', show.legend = NA) + 
  geom_vline(aes(xintercept = tst[1]), color = '#000000', show.legend = NA, linetype="dashed", alpha = 0.75) + 
  geom_vline(aes(xintercept = tst[3]), color = '#000000', show.legend = NA, linetype="dashed", alpha = 0.75) + 
  geom_point(data = pr4p, aes(x= x, y= y), shape= 3) +
  geom_jitter(shape=16, position=position_jitter(0.005),size=0.75) +
  geom_smooth(color = '#000000') +
  scale_x_continuous(breaks=k) +
  labs(x = "Proportion of location-A (amber) tips in a downsampled tree", y = "Accuracy") +
  ylim(0,1) + 
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Average of absolute accuracy of 50 downsampled trees of 25 true trees with migration rate 0.9") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,width = unit(16, "cm"),height = unit(9, "cm"))

grid.draw(tp)
