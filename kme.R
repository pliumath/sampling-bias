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
#functions

corres <- function(phy){
  
  cor <- as.data.frame(1:(2*phy$Nnode+1))
  cor[,2] <- c(phy$tip.label,phy$node.label)
  cor[,3] <- c(phy$tip.state,phy$node.state)
  names(cor) <- c("ind","names","states")
  return(cor)
  
}

tran <- function(cor,nodes){
  
  if (is.numeric(nodes)) {
    return(sapply(nodes,function(x){cor[x,2]}))
  }
  
  if (is.character(nodes) ) {
    return(sapply(nodes,function(x){cor[cor[,2]==x,1]}))
  }
  
}

pnodes <- function(phy,cor,nodes){
  
  ind <- tran(cor,nodes)
  pind <- Ancestors(phy,ind,'parent')
  pnodes <- tran(cor,pind)
  return(pnodes)
  
}

cnodes <- function(phy,cor,nodes){
  
  ind <- tran(cor,nodes)
  cind <- Descendants(phy,ind,'children')
  cnodes <- tran(cor,cind)
  return(cnodes)
  
}

ftruekmes <- function(phy,cor,cond,D){
  
  h <- history.from.sim.discrete(phy, 0:1)
  hist <- h$history
  
  #child nodes with edges from amber to blue
  cands <- names(which(sapply(hist,function(x){m = dim(x)[1]; return(m != 1 & x[1,2]==0 & x[m,2]==1)}) & startsWith(names(hist),"nd")==TRUE))
  ind <- sapply(cands,function(x){which(phy$node.label==x)})
  
  #all tips of internal nodes
  des <- listTips(phy)
  
  res<- as.data.frame(cands)
  res[,2] <- sapply(ind, function(x){sum(phy$tip.state[des[[x]]]==1)})
  res <- res[res[,2]>=cond,]
  
  res[,3] <- pnodes(phy,cor,res[,1])
  
  kmes <- as.data.frame(cbind(res[,3],res[,1],res[,2]))
  kmes[,4] <- diag(D[tran(cor,res[,3]),tran(cor,res[,1])])
  names(kmes) <- c("pnodes","cnodes","numbltips","dist")
  
  return(kmes)
  
}


fdeskmes <- function(phy,cor,nodes){
  
  inds <- tran(cor,nodes)
  dinds <- Descendants(phy,inds,"all")
  deskmes <- lapply(dinds,function(x){tran(cor,x)})
  
  return(deskmes)
  
}


appearance <- function(downsampled,dcor,node,truekmes,deskmes,D){
  
  dcnodes <- cnodes(downsampled,dcor,node)
  
  if (is.na(dcnodes[1])) {
    
    return(0)
    
  }else {
    
    if (truekmes[truekmes[,1] == node,2] %in% dcnodes) {
      
      return(truekmes[truekmes[,1] == node,4])
      
    } else {
      
      cn <- truekmes[truekmes[,1] == node,2]
      cind <- tran(cor,dcnodes[sapply(dcnodes, function(x){x[1] %in% deskmes[[cn]]})])
      pind <- tran(cor,node)
      return(D[pind,cind])
      
    }
    
  }
  
  
}



##########################################
##########################################
#load true trees

flag <- 0
while (flag == 0) {
  
  numcands <- 25
  
  ntips <- 100 #number of tips in the downsapled trees
  
  lmd <- 4 #birth rate
  mu <- 1 #death rate
  alph <- 0.1 #migration rate low-0.1,0.3, high-0.7,0.9
  
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

plot(h, phy,show.node.label = T,show.tip.label = F, show.tip.state = F, label.offset = 0.1, cex=1 ,cols = c('#EDB120','#0072BD'))
axisPhylo(backward = F)

n0 <- sum(phy$tip.state == 0)
n1 <- sum(phy$tip.state == 1)

N0 <- sum(phy$node.state == 0)
N1 <- sum(phy$node.state == 1)

##########################################
##########################################
#find true kmes
cor <- corres(phy)

D <- dist.nodes(phy)

truekmes <- ftruekmes(phy,cor,15,D)

deskmes <- fdeskmes(phy,cor,truekmes[,2])

bluetipskmes <- lapply(deskmes, function(y){y[sapply(y, function(x){cor[cor[,2]==x,3]}) == 1 & startsWith(y,"ex")]})


##########################################
##########################################
#estimates with subsampled tips

times <- 100 #number of downsampled trees
k <- seq(0.05,0.95,0.1) #proportion of location-A (amber) tips
l <- length(k)

ntips <- 100  #number of tips in the downsapled trees

res <- list()
rtp <- list()
for (rr in truekmes[,1]) {
  
  res[[rr]] <- matrix(0,nrow = times,ncol = l)
  rtp[[rr]] <- matrix(0,nrow = times,ncol = l)
  
}


for (i in 1:l) {
  
  print(i)
  
  for (j in 1:times) {
    
    print(j)
    
    to.remove <- rep(FALSE,phy$Nnode+1) #preallocate tips to remove
    
    sn0 <- floor(k[i]*ntips) #number of selected amber tips with a given proportion
    sn1 <- ntips - sn0
    
    rand.remove0 <- sample(which(phy$tip.state == 0),n0-sn0) #sample tips to remove
    rand.remove1 <- sample(which(phy$tip.state == 1),n1-sn1)
    to.remove[rand.remove0] <- TRUE
    to.remove[rand.remove1] <- TRUE
    
    downsampled <- prune(phy,to.remove)
    
    sh <- history.from.sim.discrete(downsampled, 0:1)
    
    dcor <- corres(downsampled)
    
    for (rr in truekmes[,1]) {
      
      rtp[[rr]][j,i] <- sum(downsampled$tip.label %in% bluetipskmes[[truekmes[truekmes[,1] == rr,2]]])
      res[[rr]][j,i] <- appearance(downsampled,dcor,rr,truekmes,deskmes,D)
      
    }
    
    
  }
  
}


##########################################
##########################################
##visualization (Fig 4, S7)


unob4p <- data.frame()
data4p <- data.frame()
for (nms in names(res)) {
  
  temp <- res[[nms]]
  temp2 <- rtp[[nms]]
  da <- data.frame()
  
  for (ll in 1:l) {
    
    unob4p <- rbind(unob4p,sum(temp2[,ll] >=5 ))
    
    tempres <- as.data.frame(table(temp[,ll]))
    tempres[,3] <- rep(k[ll],length(table(temp[,ll])))
    tempres[,1] <- as.numeric(names(table(temp[,ll])))
    
    
    da <- rbind(da,tempres)
    
  }
  
  da[,4] <- rep(nms,dim(da)[1])
  
  data4p <- rbind(data4p,da)
  
}

data4p <- data4p[data4p[,1] != "0",]
names(data4p) <- c("v4","v1","v2","v3")
data4p$v3 <- factor(data4p$v3,levels = rev(names(res)))
for (nms in names(res)) {
  
  base <- truekmes[truekmes[,1] == nms,4]
  data4p[data4p[,4] == nms,1] = data4p[data4p[,4] == nms,1] - base
  
}
data4p[data4p[,1] < 10^-15,1]<- 0


unob4p[,2] <- rep(k,dim(truekmes)[1])
unob4p[,3] <- rep(names(res),each = 10)
names(unob4p) <- c("v1","v2","v3")
unob4p$v3 <- factor(unob4p$v3,levels = rev(names(res)))

hdr4p <- unob4p
hdr4p[,1] <- 100

v3.labs <- c("KME 1","KME 2","KME 3","KME 4")
names(v3.labs) <- rev(names(res))

p <- ggplot() +
  geom_bar(data = hdr4p, aes(x = v2+0.015, y = v1), stat = "identity", width = 0.03, alpha = 0.25) + 
  geom_bar(data = unob4p, aes(x = v2+0.015, y = v1) , stat = "identity", width = 0.03) + 
  geom_bar(data = hdr4p, aes(x = v2-0.015, y = v1), stat = "identity", width = 0.03, alpha = 0.25) + 
  geom_bar(data = data4p, aes(x = v2-0.015, y = v1, fill = v4), stat = "identity", width = 0.03) +
  facet_grid(v3 ~., labeller = labeller(v3 = v3.labs)) +
  scale_fill_gradient(low = '#77AC30',high = '#D95319', limits = c(0,1.75) ,breaks = c(0, 0.5, 1, 1.5)) +
  scale_x_continuous(breaks=k) +
  labs(x = "Proportion of location-A (amber) tips in a downsampled tree", y = "Count", fill = "Observed\nvia tips \nObscured \nDistance") +
  ylim(0,times) +
  theme(
    text = element_text(size = 12,family = "serif"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) + ggtitle("Number of downsampled trees in which a KME is \nobserved, erred or obscured for the true tree with migration rate 0.1") + theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(size=12)) + theme(legend.position = "right")

tp <- set_panel_size(p,width  = unit(16, "cm"),height = unit(9/2, "cm"))

grid.draw(tp)

