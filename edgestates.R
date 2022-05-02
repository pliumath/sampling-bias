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
#load true trees

load("truetree09.Rdata")

h <- history.from.sim.discrete(phy, 0:1)

plot(h, phy,show.node.label = F,show.tip.label = F, show.tip.state = F, label.offset = 0.1, cex=1 ,cols = c('#EDB120','#0072BD'))

n0 <- sum(phy$tip.state == 0)
n1 <- sum(phy$tip.state == 1)

N0 <- sum(phy$node.state == 0)
N1 <- sum(phy$node.state == 1)



##########################################
##########################################
#estimates with subsampled tips

nodesinfo <- function(downsampled,dans){
  
  #information on nodes
  dinfo <- cbind(1:(downsampled$Nnode*2+1),c(downsampled$tip.state,downsampled$node.state),c(downsampled$tip.state,dans$lik.anc[,2]))
  
  #true parent state
  tpstates <- dinfo[Ancestors(downsampled,1:(downsampled$Nnode*2+1),'parent'),2]
  tpstates <- as.matrix(tpstates,length(tpstates),1)
  
  #reconstructed parent state
  rpstates <- sapply(Ancestors(downsampled,1:(downsampled$Nnode*2+1),'parent'), function(x){
    
    if (x == 0) {
      y = -1
    }
    if (x != 0) {
      y = dinfo[dinfo[,1] == x,3]
    }
    
    return(y)
    
  })
  rpstates <- as.matrix(rpstates,length(rpstates),1)
  
  dinfo <- cbind(1:(downsampled$Nnode*2+1),c(downsampled$tip.state,downsampled$node.state),c(downsampled$tip.state,dans$lik.anc[,2]),Ancestors(downsampled,1:(downsampled$Nnode*2+1),'parent'),rpstates)
  dinfo <- dinfo[dinfo[,4] != 0,]
  dinfo <- as.data.frame(dinfo)
  dinfo[,6] <- dinfo[,5]
  dinfo[,5] <- tpstates
  
  names(dinfo) <- c("index","truestate","recstate","parnode","trueparstate","recparstate")
  
  return(dinfo)
  
}


times <- 50 #downsampling times
k <- seq(0.05,0.95,0.1) #the location-A proportion
l <- length(k)

ntips <- 100 #number of tips in the downsampled trees

n0 <- sum(phy$tip.state == 0)
n1 <- sum(phy$tip.state == 1)

tesyy <- matrix(0,nrow = times,ncol = l)
tesyb <- matrix(0,nrow = times,ncol = l)
tesby <- matrix(0,nrow = times,ncol = l)
tesbb <- matrix(0,nrow = times,ncol = l)

resyy <- matrix(0,nrow = times,ncol = l)
resyb <- matrix(0,nrow = times,ncol = l)
resby <- matrix(0,nrow = times,ncol = l)
resbb <- matrix(0,nrow = times,ncol = l)


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
    
    dans <- ace(downsampled$tip.state, downsampled, type = "d") #reconstruct ancestral states on the downsampled tree
    
    dinfo <- nodesinfo(downsampled,dans)
    
    tesyy[j,i] <- sum(dinfo[,5] == 0 & dinfo[,2] == 0)
    tesyb[j,i] <- sum(dinfo[,5] == 0 & dinfo[,2] == 1)
    tesby[j,i] <- sum(dinfo[,5] == 1 & dinfo[,2] == 0)
    tesbb[j,i] <- sum(dinfo[,5] == 1 & dinfo[,2] == 1)
    
    resyy[j,i] <- sum((1-dinfo[,6])*(1-dinfo[,3]))
    resyb[j,i] <- sum((1-dinfo[,6])*(dinfo[,3]))
    resby[j,i] <- sum((dinfo[,6])*(1-dinfo[,3]))
    resbb[j,i] <- sum((dinfo[,6])*(dinfo[,3]))
    
  }
  
}



##########################################
##########################################
##visualization (Fig S3 S4)

tesyy4p <- as.data.frame(tesyy)
tesyy4p <- melt(tesyy4p)
tesyy4p[,1] <- sapply(rep(k, each=times),as.character)
tesyy4p[,3] <- rep("Location-A parent node",l*times)
tesyy4p[,4] <- rep("Location-A child node",l*times)

tesyb4p <- as.data.frame(tesyb)
tesyb4p <- melt(tesyb4p)
tesyb4p[,1] <- sapply(rep(k, each=times),as.character)
tesyb4p[,3] <- rep("Location-A parent node",l*times)
tesyb4p[,4] <- rep("Location-B child node",l*times)

tesby4p <- as.data.frame(tesby)
tesby4p <- melt(tesby4p)
tesby4p[,1] <- sapply(rep(k, each=times),as.character)
tesby4p[,3] <- rep("Location-B parent node",l*times)
tesby4p[,4] <- rep("Location-A child node",l*times)

tesbb4p <- as.data.frame(tesbb)
tesbb4p <- melt(tesbb4p)
tesbb4p[,1] <- sapply(rep(k, each=times),as.character)
tesbb4p[,3] <- rep("Location-B parent node",l*times)
tesbb4p[,4] <- rep("Location-B child node",l*times)

tes4p <- as.data.frame(rbind(tesyy4p,tesyb4p,tesby4p,tesbb4p))
names(tes4p) <- c("prop","num","pare","chil")

tes4p$pare <- factor(tes4p$pare,levels = c("Location-A parent node", "Location-B parent node"))
tes4p$chil <- factor(tes4p$chil,levels = c("Location-A child node", "Location-B child node"))

p <- ggplot(tes4p, aes(x=prop, y=num)) + 
  geom_jitter(shape=16, position=position_jitter(0.1),size=0.75) + 
  geom_vline(aes(xintercept = n0/(phy$Nnode+1)*(l+1)), color = '#000000', show.legend = NA) + 
  geom_violin(alpha =0.5) + 
  facet_grid(chil ~ pare) +
  labs(x = "Proportion of location-A (amber) tips in a downsampled tree", y = "Number of edges") +
  ylim(-1,201) +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle(TeX("Number of edges (with true locations) in each state for the true tree ($\\alpha = 0.9$)")) + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,width = unit(8, "cm"),height = unit(6, "cm"))

grid.draw(tp)





resyy4p <- as.data.frame(resyy)
resyy4p <- melt(resyy4p)
resyy4p[,1] <- sapply(rep(k, each=times),as.character)
resyy4p[,3] <- rep("Location-A parent node",l*times)
resyy4p[,4] <- rep("Location-A child node",l*times)

resyb4p <- as.data.frame(resyb)
resyb4p <- melt(resyb4p)
resyb4p[,1] <- sapply(rep(k, each=times),as.character)
resyb4p[,3] <- rep("Location-A parent node",l*times)
resyb4p[,4] <- rep("Location-B child node",l*times)

resby4p <- as.data.frame(resby)
resby4p <- melt(resby4p)
resby4p[,1] <- sapply(rep(k, each=times),as.character)
resby4p[,3] <- rep("Location-B parent node",l*times)
resby4p[,4] <- rep("Location-A child node",l*times)

resbb4p <- as.data.frame(resbb)
resbb4p <- melt(resbb4p)
resbb4p[,1] <- sapply(rep(k, each=times),as.character)
resbb4p[,3] <- rep("Location-B parent node",l*times)
resbb4p[,4] <- rep("Location-B child node",l*times)

res4p <- as.data.frame(rbind(resyy4p,resyb4p,resby4p,resbb4p))
names(res4p) <- c("prop","num","pare","chil")

res4p$pare <- factor(res4p$pare,levels = c("Location-A parent node", "Location-B parent node"))
res4p$chil <- factor(res4p$chil,levels = c("Location-A child node", "Location-B child node"))


ces4p <- res4p
ces4p[,2] <- ces4p[,2] - tes4p[,2]
names(ces4p) <- c("prop","num","pare","chil")

ces4p$pare <- factor(ces4p$pare,levels = c("Location-A parent node", "Location-B parent node"))
ces4p$chil <- factor(ces4p$chil,levels = c("Location-A child node", "Location-B child node"))

p <- ggplot(ces4p, aes(x=prop, y=num)) + 
  geom_jitter(shape=16, position=position_jitter(0.1),size=0.75) + 
  geom_vline(aes(xintercept = n0/(phy$Nnode+1)*(l+1)), color = '#000000', show.legend = NA) + 
  geom_hline(aes(yintercept = 0), color = '#000000', show.legend = NA) + 
  geom_violin(alpha =0.5) + 
  facet_grid(chil ~ pare) +
  labs(x = "Proportion of location-A (amber) tips in a downsampled tree", y = "Difference") +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle(TeX("Difference in number of edges for the true tree ($\\alpha = 0.9$)")) + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12))

tp <- set_panel_size(p,width = unit(8, "cm"),height = unit(6, "cm"))

grid.draw(tp)





