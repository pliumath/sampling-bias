library("strap")
library("stringr")
library("ips")
library("ggtree")
require(treeio)
require(ggplot2)
require(ggtree)
library("reshape2")
library("ggstance")
library("diversitree")
library("phytools")
library("ggimage")
library("ggpubr")
library("dplyr")
library("patchwork")
setwd("C:/Users/yexua/OneDrive/desktop/ggtree/123/ggtree/Ebola_outbreak_analysis")

###############################################
#read the tree file 
#insert group names (states) into the tree
#preparation before applying ace function
###############################################

#set seed that shows figures in the paper
set.seed(10403)

#read full tree
treefull <- read.nexus("meannode.tree")

#set negative branch to 0
for(i in 1:522){
  if(treefull$edge.length[i] < 0){
    treefull$edge.length[i] = 0
  }
}

#assign tip labels to the location using name file

groups <- read.table("tip_location.txt", sep='',
                     col.names = c('ID','group'),
                     header = FALSE, stringsAsFactors = FALSE)
gdata <- as.data.frame(groups)
sl <- grep("SierraLeone",gdata$group)
gdata$group[sl] = "Sierra Leone"

#transform full tree in order to apply ace function

tree1full <- full_join(as_tibble(treefull),gdata,by = c('label'='ID'))
tree2full <- as.treedata(tree1full)
tree3full <- as.phylo(tree1full)

#try to make node names

try_make_node <- tree3full
try_tree <- makeNodeLabel(tree3full,method = "number")

###############################################
#apply ace function 
###############################################

t <- tree3full$tip.label
t <- tree1full$group[1:262]
ansfull <- ace(t,treefull, type = 'd')

###############################################
#use lik.anc data to find the internal states
###############################################

jfull <- apply(ansfull$lik.anc, 1, which.max)
for(i in 263:523){
  if(jfull[i-262]==1){
    tree1full$group[i] = "Guinea"
    try_tree$node.label[i-262] = "Guinea"
  }
  else if(jfull[i-262]==2){
    tree1full$group[i] = "Liberia"
    try_tree$node.label[i-262] = "Liberia"
  }
  else if(jfull[i-262]==3){
    tree1full$group[i] = "Mali"
    try_tree$node.label[i-262] = "Mali"
  }
  else if(jfull[i-262]==4){
    tree1full$group[i] = "Sierra Leone"
    try_tree$node.label[i-262] = "Sierra Leone"
  }
  else{
    tree1full$group[i] = "Unknown"
    try_tree$node.label[i-262] = "Unknown"
  }
}

###############################################
#plot the full tree
###############################################

plot_treefull <- as.treedata(tree1full)

true_tree_circular<-ggtree(plot_treefull,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE),branch.length='none', layout='circular')+
  geom_point(aes(color=group))+theme_tree2()+ggtitle("topology: full tree \nASR: full tree") + 
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12))+labs(colour="location")

true_tree_treeshape<-ggtree(plot_treefull,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE))+
  geom_point(aes(color=group))+theme_tree2()+ggtitle("a)\ntopology: full tree \nASR: full tree") +
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12))+labs(colour="location", x = "Time")

#insert likelihood pie chart to the tree

ancstats<-as.data.frame(ansfull$lik.anc)
ancstats$node <- 1:try_tree$Nnode+Ntip(try_tree)
true_tree_pies <- nodepie(ancstats,cols = 1:5)
lik_true <- inset(true_tree_treeshape, true_tree_pies,width = 0.1,height = 0.1)
plot(lik_true)

################################################
#use true_tree to determine which tip where there was a recent state change (parent or parent's parent)
#extract data from the true tree
###############################################

dt <- true_tree_circular$data
tip <- c()

################################################
#consider three cases, state changes in 1 or 2 or 3 generations
################################################

#parent (1 generation)
#for(i in 1:262){
#  state <- dt$group[i]
#  if(dt$group[dt$parent[i]]!=state){
#    tip <- append(tip,i)
#  }
#}

#parent & parent's parent (2 generations)
for(i in 1:262){
  state <- dt$group[i]
  if(dt$group[dt$parent[i]]!=state|dt$group[dt$parent[dt$parent[i]]]!=state){
    tip <- append(tip,i)
  }
}

#parent & parent's parent & parent's parents' parent (3 generations)

#for(i in 1:262){
#  state <- dt$group[i]
#  if(dt$group[dt$parent[i]]!=state|dt$group[dt$parent[dt$parent[i]]]!=state|dt$group[dt$parent[dt$parent[dt$parent[i]]]]!=state){
#    tip <- append(tip,i)
#  }
#}

#find tips wish to drop

#checking
#for(i in 1:length(tip)){
#  print(dt$group[tip[i]])
#}

d_tip <- c(1:262)
d_tip <- d_tip[-tip]
#randomly downsample 80% of the tips

n <- floor(length(d_tip)*0.8)
ran_tra <- sample.int(length(d_tip),n)
tip <- sort(ran_tra)
tree1 <- full_join(as_tibble(treefull),gdata,by = c('label'='ID'))
tree2 <- as.phylo(tree1)
drop_traverlers <- function(tree1,tree2,tip){
  d <- c()
  rm_chr <- c()
  for(i in 1:length(tip)){
    d[i] <- tree2$tip.label[tip[i]]
    rm_chr[i] <- tree1$group[tip[i]]
  }
  name <- tip
  tip <- gdata$group
  tip <- tip[-name]
  
  tree_ds <- drop.tip(tree2, d,trim.internal = TRUE)
  
  for(i in 1:length(tree_ds$edge.length)){
    if(tree_ds$edge.length[i] < 0){
      tree_ds$edge.length[i] = 0
    }
  }
  results <- list("tree"=tree_ds, "tip"=tip,"drop"=d)
  return(results)
}
results <- drop_traverlers(tree1,tree2,tip)
tree_ds <- results$tree
tip <- results$tip
drop <- results$drop
#apply ace function to downsampled tree

ans <- ace(tip,tree_ds, type = 'd')
rfp <- ans$rates
ID <- tree_ds$tip.label
group <- tip
df <- data.frame(ID,group)
tree_ds_tib <- full_join(as_tibble(tree_ds),df,by = c('label'='ID'))
j <- apply(ans$lik.anc, 1, which.max)
len <- length(j)
num_tip <- len+1
for(i in num_tip+1:(num_tip+len)){
  if(i-num_tip < num_tip){
    #print(i-num_tip)
    if(j[i-num_tip]==1){
      tree_ds_tib$group[i] = "Guinea"
    }
    else if(j[i-num_tip]==2){
      tree_ds_tib$group[i] = "Liberia"
    }
    else if(j[i-num_tip]==3){
      tree_ds_tib$group[i] = "Mali"
    }
    else if(j[i-num_tip]==4){
      tree_ds_tib$group[i] = "Sierra Leone"
    }
    else{
      tree_ds_tib$group[i] = "Unknown"
    }
  }
}

###############################################
#plot the downsampled tree
###############################################

plot_ds <- as.treedata(tree_ds_tib)
down_travelers_circular<-ggtree(plot_ds,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE),branch.length='none', layout='circular')+
  geom_point(aes(color=group))+theme_tree2()+
  ggtitle("topology: downsampled tree \nASR: downsampled tree,\ndrop 80% non-travellers")+ 
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12))+labs(colour="location")

down_travelers_treeshape<-ggtree(plot_ds,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE))+
  geom_point(aes(color=group))+theme_tree2() +
  ggtitle("e)\ntopology: downsampled tree \nASR: downsampled tree,\ndrop 80% non-travellers")+ 
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12))+labs(colour="location",x="Time")

#insert likelihood pie chart to the tree

ancstatsdown<-as.data.frame(ans$lik.anc)
ancstatsdown$node <- 1:tree_ds$Nnode+Ntip(tree_ds)
down_travelers1pies <- nodepie(ancstatsdown,cols = 1:5)
lik_ds <- inset(down_travelers_treeshape, down_travelers1pies,width = 0.1,height = 0.1)
plot(lik_ds)

###############################################
#construct downsampled true tree likelihood
###############################################

treephylofull <- as.phylo(tree1full)
treedropfull <- drop.tip(treephylofull, drop,trim.internal=TRUE,subtree=TRUE)
try_drop <- drop.tip(try_tree,drop,trim.internal=TRUE)
t_drop <-as_tibble(try_drop)
drop_data <- gdata
for(i in 1:length(drop)){
  drop_data <- drop_data[-c((grep(drop[i],drop_data$ID,fixed = TRUE))),]
}
ID <- t_drop$label
group <- c(tip,try_drop$node.label)
df_drop <- data.frame(ID,group)
t_drop <- full_join(as_tibble(try_drop),df_drop,by=c('label'='ID'))

###############################################
#plot the downsampled true tree
###############################################

plot_ds_true <- as.treedata(try_drop)
true_down_traverlers_circular<-ggtree(plot_ds_true,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE),branch.length='none', layout='circular')+
  geom_point(aes(color=group))+theme_tree2()+
  ggtitle("topology: downsampled tree \nASR: full tree")+ 
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12.0))+labs(colour="location")

true_down_traverlers_treeshape<-ggtree(plot_ds_true,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE))+
  geom_point(aes(color=group))+theme_tree2() +
  ggtitle("d)\ntopology: downsampled tree \nASR: full tree")+ 
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12.0))+labs(colour="location",x="Time")

#+geom_text2(aes(subset=!isTip, label=node), hjust=-.3)
###############################################
#use lik.anc data to find the internal states for downsampled full tree
###############################################

tree1full_lik <- tree1full
try_tree_lik <- try_tree

for(i in 263:523){
  tree1full_lik$group[i] = toString(ansfull$lik.anc[i-262,])  
  try_tree_lik$node.label[i-262] = toString(ansfull$lik.anc[i-262,])  
}
treephylofull_lik <- as.phylo(tree1full_lik)
treedropfull_lik <- drop.tip(treephylofull_lik, drop,trim.internal=TRUE,subtree=TRUE)
try_drop_lik <- drop.tip(try_tree_lik,drop,trim.internal=TRUE)
t_drop_lik <-as_tibble(try_drop_lik)
drop_ancstats <- data.frame(Guinea=numeric(),Liberia=numeric(),Mali=numeric(),SierraLeone=numeric(), Unknown=numeric())
for(i in 1:length(try_drop_lik$node.label)){
  drop_ancstats[i,] <- as.numeric(strsplit(try_drop_lik$node.label[i],",")[[1]])
}
drop_ancstats$node <- 1:try_drop_lik$Nnode+Ntip(try_drop_lik)

true_tree1pies_lik <- nodepie(drop_ancstats,cols = 1:5)
lik_ds_true <- inset(true_down_traverlers_treeshape, true_tree1pies_lik,width = 0.1,height = 0.1)
plot14 <- plot(lik_ds_true)

###############################################
#arrange and save the plots
###############################################

p = ggarrange(true_tree_circular,true_down_traverlers_circular,down_travelers_circular,nrow=1,common.legend = TRUE,legend = "right")
#p
#ggsave("ebola_pie.pdf",width = 15,height = 20)
p1 = ggarrange(true_tree_treeshape,down_travelers_treeshape,true_down_traverlers_treeshape,nrow=1,common.legend = TRUE,legend = "right")
#p1
#ggsave("ebola_treeshape.pdf",width = 15,height = 20)
#p2 = ggarrange(plot(lik_true),plot(lik_ds),plot(lik_ds_true),nrow=1,common.legend = TRUE,legend = "right")
#p2
#ggsave("ebola_treeshape.pdf",width = 15,height = 20)