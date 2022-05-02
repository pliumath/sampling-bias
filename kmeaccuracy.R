library(ape)
library(diversitree)
library(phangorn)
library(phyloTop)
library(adephylo)
library(ggplot2)
library(reshape2)
library(scales)
library(egg)
library(gridBase)
library(grid)


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
#find the kmes
hist <- h$history

cands <- names(which(sapply(hist,function(x){m = dim(x)[1]; return(m != 1 & x[m,2]==1) != 1}) & startsWith(names(hist),"nd")==TRUE))
ind <- sapply(cands,function(x){which(phy$node.label==x)})

des <- listTips(phy)

data <- as.data.frame(cands)
data[,2] <- sapply(ind, function(x){sum(phy$tip.state[des[[x]]]==1)})

names(data) <- c("cands","numbluetips")

kmes <- data[data[,2]>=15,1]

##########################################
##########################################
#estimates with subsampled tips

times <- 50 #number of downsampled trees
k <- seq(0.05,0.95,0.1) #proportion of location-A (amber) tips
l <- length(k)

ntips <- 100 #number of tips in the downsapled trees

res <- matrix(0,nrow = times,ncol = l)
res0 <- matrix(0,nrow = times,ncol = l)
unobsc <- matrix(0,nrow = times,ncol = l)

estmig <- matrix(0,nrow = times,ncol = l)

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
    
    estmig[j,i] <- dans$rates
    
    unob <- unlist(sapply(kmes, function(x){ind <- which(downsampled$node.label == x); return(dans$lik.anc[ind,2])}))
    
    pkmes <- downsampled$node.label[unlist(sapply(kmes, function(x){ind <- which(downsampled$node.label == x); return( Ancestors(downsampled,ind+downsampled$Nnode+1,'parent'))})) - downsampled$Nnode - 1]
    punob <- unlist(sapply(pkmes, function(x){ind <- which(downsampled$node.label == x); return(dans$lik.anc[ind,1])}))
    
    unobsc[j,i] <- length(unob)
    
    if (length(unob) != 0) {
      # res[j,i] <- sum(unob*punob)/length(kmes)
      # res0[j,i] <- sum(unob)/length(kmes)
      
      res[j,i] <- sum(unob)/length(kmes)
      res0[j,i] <- sum(punob)/length(kmes)
    }
    
  }
  
}



##########################################
##########################################
##visualization (Fig S8)

res4plot <- as.data.frame(res)
res4plot <- melt(res4plot)


unob4plot <- as.data.frame(unobsc)
unob4plot <- sapply(melt(unob4plot),as.character)

res4plot[,3] <- unob4plot[,2]
res4plot$variable <- sapply(rep(k, each=times),as.character)



res4plot0 <- as.data.frame(res0)
res4plot0 <- melt(res4plot0)

res4plot0 <- as.data.frame(res0)
res4plot0 <- melt(res4plot0)

unob4plot0 <- as.data.frame(unobsc)
unob4plot0 <- sapply(melt(unob4plot0),as.character)

res4plot0[,3] <- unob4plot0[,2]
res4plot0$variable <- sapply(rep(k, each=times),as.character)

res4plot <- rbind(res4plot,res4plot0)

res4plot[,4] <- c(rep("Child node (location B)",500),rep("Parent node (location A)",500))

names(res4plot) <- c("prop","accu","numunob","pc")

res4plot$pc <- factor(res4plot$pc,levels = c("Parent node (location A)", "Child node (location B)"))



p <- ggplot(res4plot, aes(x=prop, y=accu)) + 
  geom_point(aes(color = numunob),position=position_jitter(0.1),shape=16, size=2.5, alpha = 0.5) + geom_vline(aes(xintercept = n0/(phy$Nnode+1)*(l+1)), color = '#000000', show.legend = NA) + 
  facet_grid(pc ~.) +
  scale_color_manual(values=c( '#000000','#4DBEEE' ,'#7E2F8E','#77AC30', '#D95319'))+
  labs(x = "Proportion of location-A (amber) tips in a downsampled tree", y = "Accuracy", color = "Number of \nobserved \nKMEs") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Average absolute accuracy of the parent node and the child node \nover all KMEs in the true tree with migration rate 0.1") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,width  = unit(16, "cm"),height = unit(9, "cm"))

grid.draw(tp)
