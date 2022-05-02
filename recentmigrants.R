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

nullmodaccu <- function(downsampled){
  
  N <- downsampled$Nnode 
  
  n0 <- round(sum(downsampled$tip.state == 0)*(N)/(N+1)) #get yellow proportion
  n1 <- round(sum(downsampled$tip.state == 1)*(N)/(N+1))
  
  results <- c()
  
  for (i in 1:50) {
    
    rand.state <- sample(1:N,n1) #assign random states to internal nodes
    states <- rep(0,N)
    states[rand.state] <- 1
    nms <- names(downsampled$node.state)
    names(states) <- nms
    downsampled$node.state <- states
    
    sh <- history.from.sim.discrete(downsampled, 0:1)
    
    dans <- ace(downsampled$tip.state, downsampled, type = "d") #reconstruct ancestral states on the downsampled tree
    
    accuracy <-  cbind(abs(dans$lik.anc[,2] - downsampled$node.state),1-abs(dans$lik.anc[,2] - downsampled$node.state))
    results[i] <- mean(accuracy[,2])
    
  }
  
  return(mean(results)) 
  
}


##########################################
##########################################
#load true trees

load("truetree09.Rdata")

h <- history.from.sim.discrete(phy, 0:1)

plot(h, phy,show.node.label = F,show.tip.label = F, show.tip.state = F, label.offset = 0.1, cex=1 ,cols = c('#EDB120','#0072BD'))

n0 <- sum(phy$tip.state == 0)
n1 <- sum(phy$tip.state == 1)

N0 <- sum(phy$node.state == 0)
N1 <- sum(phy$node.state == 1)


mns <- unique(phy$hist[startsWith(phy$hist[,8],"ex"),8]) # recent migrants information
tips <- as.data.frame(cbind(phy$tip.state == 1,rep(0,length(phy$tip.state))))
tips[mns,2] <- 1
tips[,3] <- phy$tip.label
names(tips) <- c("state","ismns","tips")

##########################################
##########################################
#estimates with downsampled trees


times <- 50 #number of downsampled trees
k <- seq(0.05,0.95,0.1) #recent migrants proportion
l <- length(k)

ntips <- 80 #number of tips in the downsapled trees

meanaccu <- matrix(0,nrow = times,ncol = l) #preallocating
adjaccu <- matrix(0,nrow = times,ncol = l)

dstrees <- list()

for (i in 1:l) {
  
  print(i)
  
  for (j in 1:times) {
    
    print(j)
    
    to.remove <- rep(TRUE,phy$Nnode+1) #preallocate tips to remove
    
    sn0 <- (1-k[i])*ntips/2 #number of selected non-recent-migrant tips with a given proportion
    sn1 <- k[i]*ntips/2  #number of selected recent-migrant tips with a given proportion
    
    rand.sel00 <- sample(which(tips[,1] == 0 & tips[,2] == 0),min(1*sn0,sum(tips[,1] == 0 & tips[,2] == 0))) #sample tips
    rand.sel01 <- sample(which(tips[,1] == 0 & tips[,2] == 1),min(1*sn1,sum(tips[,1] == 0 & tips[,2] == 1)))
    
    rand.sel10 <- sample(which(tips[,1] == 1 & tips[,2] == 0),min(1*sn0,sum(tips[,1] == 1 & tips[,2] == 0))) #sample tips
    rand.sel11 <- sample(which(tips[,1] == 1 & tips[,2] == 1),min(1*sn1,sum(tips[,1] == 1 & tips[,2] == 1)))
    
    to.remove[rand.sel00] <- FALSE
    to.remove[rand.sel01] <- FALSE
    
    to.remove[rand.sel10] <- FALSE
    to.remove[rand.sel11] <- FALSE
    
    
    downsampled <- prune(phy,to.remove)
    
    dstrees[[(i-1)*times+j]]  <- downsampled
    
    sh <- history.from.sim.discrete(downsampled, 0:1)
    dans <- ace(downsampled$tip.state, downsampled, type = "d") #reconstruct ancestral states on the downsampled tree
    
    accuracy <-  cbind(abs(dans$lik.anc[,2] - downsampled$node.state),1-abs(dans$lik.anc[,2] - downsampled$node.state))
    meanaccu[j,i] <- mean(accuracy[,2])
    
    adjaccu[j,i] <-  mean(accuracy[,2]) - nullmodaccu(downsampled)
    
  }
  
}

##########################################
##########################################
##visualization (Fig 2, S6)

meanaccu4plot <- as.data.frame(meanaccu)

meanaccu4plot <- melt(meanaccu4plot)
meanaccu4plot$variable <- sapply(rep(k, each=times),as.character)

p <- ggplot(meanaccu4plot, aes(x=variable, y=value)) + 
  geom_jitter(shape=16, position=position_jitter(0.1),size=0.75) + 
  geom_vline(aes(xintercept = (sum(tips[,2] == 1)/dim(tips)[[1]])*(l+1)), color = '#000000', show.legend = NA) + 
  geom_violin(alpha =0.5) + 
  labs(x = "Proportion of location-A (amber) tips in a downsampled tree", y = "Accuracy") +
  ylim(0,1) +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Absolute accuracy of downsampled trees with 50% location-A tips") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,width = unit(16, "cm"),height = unit(9, "cm"))

grid.draw(tp)




adjaccu4plot <- as.data.frame(adjaccu)

adjaccu4plot <- melt(adjaccu4plot)
adjaccu4plot$variable <- sapply(rep(k, each=times),as.character)

p <- ggplot(adjaccu4plot, aes(x=variable, y=value)) + 
  geom_jitter(shape=16, position=position_jitter(0.1),size=0.75) + 
  geom_vline(aes(xintercept = (sum(tips[,2] == 1)/dim(tips)[[1]])*(l+1)), color = '#000000', show.legend = NA) + 
  geom_hline(aes(yintercept = 0), color = '#000000', show.legend = NA) + 
  geom_violin(alpha =0.5) + 
  labs(x = "Proportion of recent migrants in a downsampled tree", y = "Accuracy") +
  ylim(-0.5,0.5) +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Relative accuracy of downsampled trees with 50% location-A tips") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,width = unit(16, "cm"),height = unit(9, "cm"))

grid.draw(tp)



ptree <- dstrees[[14]] 
ph <- history.from.sim.discrete(ptree, 0:1)

px <- ptree$tip.state
pans <- ace(px, ptree, type = "d")

paccu <-  cbind(abs(pans$lik.anc[,2] - ptree$node.state),1-abs(pans$lik.anc[,2] - ptree$node.state))

par(mfrow=c(1,2))

plot(ph, ptree,show.node.label = F,show.tip.label = T, show.tip.state = FALSE,label.offset = 0.1, cex=1, cols = c('#EDB120','#0072BD'))

plot(ptree, type = "p", label.offset = 0.1, cex = 1, show.node.label = F, show.tip.label = F)
co <- c('#D95319', '#77AC30')
nodelabels(pie = paccu, piecol = co, cex = 0.5)

ptree$tip.label[ptree$tip.label %in% mns]

#894x1142

ptree <- dstrees[[460]]
ph <- history.from.sim.discrete(ptree, 0:1)

px <- ptree$tip.state
pans <- ace(px, ptree, type = "d")

paccu <-  cbind(abs(pans$lik.anc[,2] - ptree$node.state),1-abs(pans$lik.anc[,2] - ptree$node.state))

par(mfrow=c(1,2))

plot(ph, ptree,show.node.label = F,show.tip.label = T, show.tip.state = FALSE,label.offset = 0.1, cex=1, cols = c('#EDB120','#0072BD'))

plot(ptree, type = "p", label.offset = 0.1, cex = 1, show.node.label = F, show.tip.label = F)
co <- c('#D95319', '#77AC30')
nodelabels(pie = paccu, piecol = co, cex = 0.5)

ptree$tip.label[!(ptree$tip.label %in% mns)]

#894x1142



