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

tree1full_lik <- tree1full
try_tree_lik <- try_tree

for(i in 263:523){
  tree1full_lik$group[i] = toString(ansfull$lik.anc[i-262,])  
  try_tree_lik$node.label[i-262] = toString(ansfull$lik.anc[i-262,])  
}


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
  geom_point(aes(color=group))+theme_tree2()+ggtitle("a)\ntopology: full tree \nASR: full tree")+
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12))+labs(colour="location", x = "Time")

#insert likelihood pie chart to the tree

ancstats<-as.data.frame(ansfull$lik.anc)
ancstats$node <- 1:try_tree$Nnode+Ntip(try_tree)
true_tree_pies <- nodepie(ancstats,cols = 1:5)
lik_true <- inset(true_tree_treeshape, true_tree_pies,width = 0.1,height = 0.1)
plot(lik_true)

###############################################
#Downsample Guinea or Sierra Leone
###############################################

gui_id <- grep("Guinea",groups$group)
sie_id <- grep("Sierra Leone",groups$group)

tree1 <- full_join(as_tibble(tree),gdata,by = c('label'='ID'))
tree3 <- as.phylo(tree1)


ran_gen <- function(gui_per,sie_per,tree1,tree3){
  
  n_gui <- floor(length(gui_id)*gui_per)
  n_sie <- floor(length(sie_id)*sie_per)
  
  ran_gui <- sample.int(length(gui_id),n_gui)
  ran_sie <- sample.int(length(sie_id),n_sie)
  ran_gui_sort <- sort(ran_gui)
  ran_sie_sort <- sort(ran_sie)
  
  drop_id <- c()
  if(length(ran_gui_sort)>0){
    for(i in 1:length(ran_gui_sort)){
      drop_id[i] = gui_id[ran_gui_sort[i]]
    }
  }
  
  if(length(ran_sie_sort)>0){
    for(i in 1:length(ran_sie_sort)){
      drop_id[i] = sie_id[ran_sie_sort[i]]
    }
  }
  tip <- gdata$group
  
  if(length(drop_id)>=1){
    d <- c()
    rm_chr <- c()
    for(i in 1:length(drop_id)){
      d[i] <- tree3$tip.label[drop_id[i]]
      rm_chr[i] <- tree1$group[drop_id[i]]
    }
    tip <- tip[-drop_id]}
  
  tree5 <- drop.tip(tree3, d,trim.internal = TRUE)
  
  for(i in 1:length(tree5$edge.length)){
    if(tree5$edge.length[i] < 0){
      tree5$edge.length[i] = 0
    }
  }
  results <- list("tree"=tree5, "tip"=tip,"drop"=d)
  return(results)
}
###############################################

#change here for Sierra Leone or Guinea
gui_per <- 0.8
sie_per <- 0.0

###############################################
results <- ran_gen(gui_per=gui_per,sie_per=sie_per,tree1=tree1,tree3=tree3)
tree5 <- results$tree
tip <- results$tip
drop <- results$drop

ans <- ace(tip,tree5, type = 'd')
ansgui <- ans
anssie <- ans

rfp <- ans$rates

ID <- tree5$tip.label
group <- tip
df <- data.frame(ID,group)

tree10 <- full_join(as_tibble(tree5),df,by = c('label'='ID'))

j <- apply(ans$lik.anc, 1, which.max)
len <- length(j)
num_tip <- len+1
for(i in num_tip+1:(num_tip+len)){
  if(i-num_tip < num_tip){
    #print(i-num_tip)
    if(j[i-num_tip]==1){
      tree10$group[i] = "Guinea"
    }
    else if(j[i-num_tip]==2){
      tree10$group[i] = "Liberia"
    }
    else if(j[i-num_tip]==3){
      tree10$group[i] = "Mali"
    }
    else if(j[i-num_tip]==4){
      tree10$group[i] = "Sierra Leone"
    }
    else{
      tree10$group[i] = "Unknown"
    }
  }
}

###############################################
#plot the downsampled tree depends on Guinea or Sierra Leone
###############################################

plot_ds <- as.treedata(tree10)
down_random_gui<-ggtree(plot_ds,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE))+
  geom_point(aes(color=group))+theme_tree2()+ggtitle("c)\nTopology: downsampled tree \nASR:downsample 80% Guinea")+
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12))+labs(colour="location")

ancstatsdown_gui<-as.data.frame(ans$lik.anc)
ancstatsdown_gui$node <- 1:tree5$Nnode+Ntip(tree5)
down_travelers1pies_gui <- nodepie(ancstatsdown_gui,cols = 1:5)
lik_gui <- inset(down_random_gui, down_travelers1pies_gui,width = 0.05,height = 0.05)
plot(lik_gui)

#plot_ds <- as.treedata(tree10)
#down_random_gui<-ggtree(plot_ds,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE))+
#  geom_point(aes(color=group))+theme_tree2()+ggtitle("c)\nTopology: downsampled tree \nASR:downsample 80% Sierra Leone")+ geom_text2(aes(subset=!isTip, label=node), hjust=-.1)+
#  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12))+labs(colour="location")

#ancstatsdown_sie<-as.data.frame(ans$lik.anc)
#ancstatsdown_sie$node <- 1:tree5$Nnode+Ntip(tree5)
#down_travelers1pies_sie <- nodepie(ancstatsdown_sie,cols = 1:5)
#lik_sie <- inset(down_random_sie, down_travelers1pies_gui,width = 0.05,height = 0.05)
#plot(lik_sie)

#############################################################

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
plot_ds_true <- as.treedata(try_drop)

true_down_circular<-ggtree(plot_ds_true,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE),branch.length='none', layout='circular')+
  geom_point(aes(color=group))+theme_tree2()+
  ggtitle("topology: downsampled tree \nASR: full tree")+ 
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12.0))+labs(colour="location")

true_down_treeshape<-ggtree(plot_ds_true,right=TRUE,mrsd = "2015-01-31",options(ignore.negative.edge=TRUE))+
  geom_point(aes(color=group))+theme_tree2()+ geom_text2(aes(subset=!isTip, label=node), hjust=-.1)+
  ggtitle("b)\nTopology: downsampled tree \nASR: full tree")+ 
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12.0))+labs(colour="location",x="Time")
#+geom_text2(aes(subset=!isTip, label=node), hjust=-.1)

####################################
treephylofull_lik <- as.phylo(tree1full_lik)
treedropfull_lik <- drop.tip(treephylofull_lik, drop,trim.internal=TRUE,subtree=TRUE)
try_drop_lik <- drop.tip(try_tree_lik,drop,trim.internal=TRUE)
t_drop_lik <-as_tibble(try_drop_lik)
drop_ancstats <- data.frame(Guinea=numeric(),Liberia=numeric(),Mali=numeric(),SierraLeone=numeric(), Unknown=numeric())
for(i in 1:length(try_drop_lik$node.label)){
  drop_ancstats[i,] <- as.numeric(strsplit(try_drop_lik$node.label[i],",")[[1]])
}

drop_ancstats$node <- 1:try_drop_lik$Nnode+Ntip(try_drop_lik)

p = ggarrange(true_tree_treeshape,true_down_treeshape,down_random_gui,nrow=1,common.legend = TRUE,legend = "right")

