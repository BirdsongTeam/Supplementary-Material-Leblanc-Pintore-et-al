######################################################################################################################
##############Foot adaptation to climbing in oven birds and woodcreepers (Furnariida)#################################
######################################################################################################################







#Loading packages ########
library(ape)
library(geomorph)
library(Morpho)
library(readr)
library(geomorph)
library(Morpho)
library(ggplot2)
library(ggalt)
library(ggrepel)
library(phytools)


############ set the working directory to SM8d_Data_file_loading



classification <- read_delim("classification.csv", 
                             delim = ";", escape_double = FALSE, col_types = cols(`Climbing index` = col_number(), `Tmt shape` = col_character(),
                                                                                  `BodyMass-Value` = col_number(), 
                                                                                  LogBodyMass = col_number(), `ForStrat-ground` = col_number(), 
                                                                                  `ForStrat-understory` = col_number(), 
                                                                                  `ForStrat-midhigh` = col_number(), 
                                                                                  `ForStrat-canopy` = col_number()), 
                             locale = locale(decimal_mark = ","), 
                             trim_ws = TRUE)
View(classification)

###########Loading the GPA alignment for each Digit X phalanx X 
load("D1P1")
load("D1P2")
load("D2P1")
load("D2P2")
load("D2P3")
load("D3P1")
load("D3P2")
load("D3P3")
load("D3P4")
load("D4P1")
load("D4P2")
load("D4P3")
load("D4P4")
load("D4P5")
load("tmt")




##########Generalized Procruste alignment 
GPA_D1P1 <- procSym(D1P1$dataslide)
GPA_D1P2 <- procSym(D1P2$dataslide)
GPA_D2P1 <- procSym(D2P1$dataslide)
GPA_D2P2 <- procSym(D2P2$dataslide)
GPA_D2P3 <- procSym(D2P3$dataslide)
GPA_D3P1 <- procSym(D3P1$dataslide)
GPA_D3P2 <- procSym(D3P2$dataslide)
GPA_D3P3 <- procSym(D3P3$dataslide)
GPA_D3P4 <- procSym(D3P4$dataslide)
GPA_D4P1 <- procSym(D4P1$dataslide)
GPA_D4P2 <- procSym(D4P2$dataslide)
GPA_D4P3 <- procSym(D4P3$dataslide)
GPA_D4P4 <- procSym(D4P4$dataslide)
GPA_D4P5 <- procSym(D4P5$dataslide)
GPA_tmt <- procSym(ptsarray)



#####classification correction for each phalanx as some bones are missing/damaged for some specimens #####

classification_D2P3 <- classification[-c(12,21),]############# delete specimens for digit 2 phalanx 3
cbind(classification_D2P3$ID,row.names(GPA_D2P3$PCscores))####Double check


classification_D3P3_D4P4<- classification[-c(29,55),]############# delete specimens for phalanges D3P3, D4P4
cbind(classification_D3P3_D4P4$ID,row.names(GPA_D3P3$PCscores))####Double check
cbind(classification_D3P3_D4P4$ID,row.names(GPA_D4P4$PCscores))####Double check


classification_D3P4 <- classification[-c(12,29),]############# delete specimens for digit 3 phalanx 4
cbind(classification_D3P4$ID,row.names(GPA_D3P4$PCscores))####Double check

classification_D4P3 <- classification[-c(55,23,40),] ########### delete specimens for digit 4 phalanx 3 
cbind(classification_D4P3$ID,row.names(GPA_D4P3$PCscores))####Double check



classification_inter_terminal <- classification[-c(55),] ###########Others intermediary phalanges D2P3,D3P2,D4P2,D4P3 and ungual phalanx digit 4 
cbind(classification_inter_terminal$ID,row.names(GPA_D4P5$PCscores))####Double check


##############loading label tree
label_tree <- read_csv("label_tree.csv")
label_tree_D2P3  <- label_tree[-c(12,21),]
label_tree_D3P4  <- label_tree[-c(12,29),]
label_D2P2_D3P2_D4P2_D4P5 <- label_tree[-c(55),]
label_D4P3<- label_tree[-c(55,23,40),]
label_tree_D3P3_D4P4<- label_tree[-c(55,29),]

##############loading  tree
tr2<-read.tree(file="tree_furna")######## All proximal phalanges + Digit 1 Phalanx 2
tr2_D2P3<-read.tree(file="tree_furna_D2P3.txt")#################Digit 2 Phalanx 3
tr2_D3P4<-read.tree(file="tree_furna_D3P4.txt")#################Digit 3 Phalanx 4
tr2_D2P2_D3P2_D4P2_D4P5<-read.tree(file="tree_furna_D2P2_D3P2_D4P2_D4P5.txt")##################Digit 2 phalanx 2, Digit 3 Phalanx 2, Digit 4 Phalanx 2, Digit 4 Phalanx 5
tr2_D4P3<- read.tree(file="tree_furna_D4P3.txt")##########Digit  4 phalanx 3
tr2_D3D3_D4D4<-read.tree(file="tree_furna_D3D3_D4D4.txt")###########Digit 3 Phalanx 3, Digit 4 Phalanx 4





############### formalize tree 

tr2=compute.brlen(tr2, method=10)
plot.phylo(tr2, type = "cladogram", use.edge.length = FALSE)
tres<-multi2di(tr2)

tr2_D2P3=compute.brlen(tr2_D2P3, method=10)
plot.phylo(tr2_D2P3, type = "cladogram", use.edge.length = FALSE)
tres_D2P3<-multi2di(tr2_D2P3)

tr2_D3P4=compute.brlen(tr2_D3P4, method=10)
plot.phylo(tr2_D3P4, type = "cladogram", use.edge.length = FALSE)
tres_D3P4<-multi2di(tr2_D3P4)

tr2_D2P2_D3P2_D4P2_D4P5=compute.brlen(tr2_D2P2_D3P2_D4P2_D4P5, method=10)
plot.phylo(tr2_D2P2_D3P2_D4P2_D4P5, type = "cladogram", use.edge.length = FALSE)
tres_D2P2_D3P2_D4P2_D4P5<-multi2di(tr2_D2P2_D3P2_D4P2_D4P5)

tr2_D4P3=compute.brlen(tr2_D4P3, method=10)
plot.phylo(tr2_D4P3, type = "cladogram", use.edge.length = FALSE)
tres_D4P3<-multi2di(tr2_D4P3)

tr2_D3P3_D4P4=compute.brlen(tr2_D3D3_D4D4, method=10)
plot.phylo(tr2_D3D3_D3D4, type = "cladogram", use.edge.length = FALSE)
tres_D3P3_D4P4<-multi2di(tr2_D3D3_D4D4)



########### Variance extraction ################

########tmt

GPA_tmt$Variance


#########Digit I
GPA_D1P1$Variance
GPA_D1P2$Variance

###########Digit II

GPA_D2P1$Variance
GPA_D2P2$Variance
GPA_D2D3$Variance


#########Digit III
GPA_D3P1$Variance
GPA_D3P2$Variance
GPA_D3P3$Variance
GPA_D3P4$Variance



###########Digit IV
GPA_D4P1$Variance
GPA_D4P2$Variance
GPA_D4P3$Variance
GPA_D4P4$Variance
GPA_D4P5$Variance





################################# Quick generated plot ################################


PCA_D1P1 <- {plot(GPA_D1P1$PCscores[,1],GPA_D1P1$PCscores[,2])
text(GPA_D1P1$PCscores[,1],GPA_D1P1$PCscores[,2],labels=classification$sp2)}

PCA_D1P2 <- {plot(GPA_D1P2$PCscores[,1],GPA_D1P2$PCscores[,2])
  text(GPA_D1P2$PCscores[,1],GPA_D1P2$PCscores[,2],labels=classification$sp2)}

PCA_D2P1 <- {plot(GPA_D2P1$PCscores[,1],GPA_D2P1$PCscores[,2])
  text(GPA_D2P1$PCscores[,1],GPA_D2P1$PCscores[,2],labels=classification$sp2)}

PCA_D2P2 <- {plot(GPA_D2P2$PCscores[,1],GPA_D2P2$PCscores[,2])
  text(GPA_D2P2$PCscores[,1],GPA_D2P2$PCscores[,2],labels=Classification_inter_terminal$sp2)}

PCA_D2P3 <- {plot(GPA_D2P3$PCscores[,1],GPA_D2P3$PCscores[,2])
  text(GPA_D2P3$PCscores[,1],GPA_D2P3$PCscores[,2],labels=classification_D2P3$sp2)}

PCA_D3P1 <- {plot(GPA_D3P1$PCscores[,1],GPA_D3P1$PCscores[,2])
  text(GPA_D3P1$PCscores[,1],GPA_D3P1$PCscores[,2],labels=classification$sp2)}

PCA_D3P2 <- {plot(GPA_D3P2$PCscores[,1],GPA_D3P2$PCscores[,2])
  text(GPA_D3P2$PCscores[,1],GPA_D3P2$PCscores[,2],labels=Classification_inter_terminal$sp2)}

PCA_D3P3 <- {plot(GPA_D3P3$PCscores[,1],GPA_D3P3$PCscores[,2])
  text(GPA_D3P3$PCscores[,1],GPA_D3P3$PCscores[,2],labels=classification_D3P3_D4P4$sp2)}

PCA_D3P4 <- {plot(GPA_D3P4$PCscores[,1],GPA_D3P4$PCscores[,2])
  text(GPA_D3P4$PCscores[,1],GPA_D3P4$PCscores[,2],labels=classification_D3P4$sp2)}

PCA_D4P1 <- {plot(GPA_D4P1$PCscores[,1],GPA_D4P1$PCscores[,2])
  text(GPA_D4P1$PCscores[,1],GPA_D4P1$PCscores[,2],labels=classification$sp2)}

PCA_D4P2 <- {plot(GPA_D4P2$PCscores[,1],GPA_D4P2$PCscores[,2])
  text(GPA_D4P2$PCscores[,1],GPA_D4P2$PCscores[,2],labels=Classification_inter_terminal$sp2)}

PCA_D4P3 <- {plot(GPA_D4P3$PCscores[,1],GPA_D4P3$PCscores[,2])
  text(GPA_D4P3$PCscores[,1],GPA_D4P3$PCscores[,2],labels=Classification_inter_terminal$sp2)}

PCA_D4P4 <- {plot(GPA_D4P4$PCscores[,1],GPA_D4P4$PCscores[,2])
  text(GPA_D4P4$PCscores[,1],GPA_D4P4$PCscores[,2],labels=classification_D3P3_D4P4)}

PCA_D4P5 <- {plot(GPA_D4P5$PCscores[,1],GPA_D4P5$PCscores[,2])
  text(GPA_D4P5$PCscores[,1],GPA_D4P5$PCscores[,2],labels=Classification_inter_terminal$sp2)}


PCA_tmt <- {plot(GPA_tmt$PCscores[,1],GPA_tmt$PCscores[,2])
  text(GPA_tmt$PCscores[,1],GPA_tmt$PCscores[,2],labels=classification$sp2)}




##################gg plot###################

library(ggplot2)
library(ggalt)
library(ggrepel)

dat_tmt<- data.frame(GPA_tmt$PCscores[,1],GPA_tmt$PCscores[,2])

dat_d1p1 <- data.frame(GPA_D1P1$PCscores[,1],GPA_D1P1$PCscores[,2])
dat_d1p2 <- data.frame(GPA_D1P2$PCscores[,1],GPA_D1P2$PCscores[,2])

dat_d2p1 <- data.frame(GPA_D2P1$PCscores[,1],GPA_D2P1$PCscores[,2])
dat_d2p2 <- data.frame(GPA_D2P2$PCscores[,1],GPA_D2P2$PCscores[,2])
dat_d2p3 <- data.frame(GPA_D2P3$PCscores[,1],GPA_D2P3$PCscores[,2])

dat_d3p1 <- data.frame(GPA_D3P1$PCscores[,1],GPA_D3P1$PCscores[,2])
dat_d3p2 <- data.frame(GPA_D3P2$PCscores[,1],GPA_D3P2$PCscores[,2])
dat_d3p3 <- data.frame(GPA_D3P3$PCscores[,1],GPA_D3P3$PCscores[,2])
dat_d3p4 <- data.frame(GPA_D3P4$PCscores[,1],GPA_D3P4$PCscores[,2])

dat_d4p1 <- data.frame(GPA_D4P1$PCscores[,1],GPA_D4P1$PCscores[,2])
dat_d4p2 <- data.frame(GPA_D4P2$PCscores[,1],GPA_D4P2$PCscores[,2])
dat_d4p3 <- data.frame(GPA_D4P3$PCscores[,1],GPA_D4P3$PCscores[,2])
dat_d4p4 <- data.frame(GPA_D4P4$PCscores[,1],GPA_D4P4$PCscores[,2])
dat_d4p5 <- data.frame(GPA_D4P5$PCscores[,1],GPA_D4P5$PCscores[,2])


#################### ggplot plot


#############Preparing the each plot

plot_tmt <- ggplot(dat_tmt , aes(dat_tmt[,1], dat_tmt[,2],color=classification$`Perching performance` ,label = classification$sp2,shape=classification$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification$`Perching performance`,fill=classification$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)

plot_d1p1 <- ggplot(dat_d1p1 , aes(dat_d1p1[,1], dat_d1p1[,2],color=classification$`Perching performance` ,label = classification$sp2,shape=classification$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification$`Perching performance`,fill=classification$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)
plot_d1p2 <- ggplot(dat_d1p2 , aes(dat_d1p2[,1], dat_d1p2[,2],color=classification$`Perching performance` ,label = classification$sp2,shape=classification$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification$`Perching performance`,fill=classification$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)

plot_d2p1 <- ggplot(dat_d2p1 , aes(dat_d2p1[,1], dat_d2p1[,2],color=classification$`Perching performance` ,label = classification$sp2,shape=classification$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification$`Perching performance`,fill=classification$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)
plot_d2p2 <- ggplot(dat_d2p2 , aes(dat_d2p2[,1], dat_d2p2[,2],color=classification_inter_terminal$`Perching performance` ,label = classification_inter_terminal$sp2,shape=classification_inter_terminal$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification_inter_terminal$`Perching performance`,fill=classification_inter_terminal$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)
plot_d2p3 <- ggplot(dat_d2p3 , aes(dat_d2p3[,1], dat_d2p3[,2],color=classification_D2P3$`Perching performance` ,label = classification_D2P3$sp2,shape=classification_D2P3$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification_D2P3$`Perching performance`,fill=classification_D2P3$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)

plot_d3p1 <- ggplot(dat_d3p1 , aes(dat_d3p1[,1], dat_d3p1[,2],color=classification$`Perching performance` ,label = classification$sp2,shape=classification$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification$`Perching performance`,fill=classification$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)
plot_d3p2 <- ggplot(dat_d3p2 , aes(dat_d3p2[,1], dat_d3p2[,2],color=classification_inter_terminal$`Perching performance` ,label = classification_inter_terminal$sp2,shape=classification_inter_terminal$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification_inter_terminal$`Perching performance`,fill=classification_inter_terminal$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)
plot_d3p3 <- ggplot(dat_d3p3 , aes(dat_d3p3[,1], dat_d3p3[,2],color=classification_D3P3_D4P4$`Perching performance` ,label = classification_D3P3_D4P4$sp2,shape=classification_D3P3_D4P4$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification_D3P3_D4P4$`Perching performance`,fill=classification_D3P3_D4P4$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)
plot_d3p4 <- ggplot(dat_d3p4 , aes(dat_d3p4[,1], dat_d3p4[,2],color=classification_D3P4$`Perching performance` ,label = classification_D3P4$sp2,shape=classification_D3P4$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification_D3P4$`Perching performance`,fill=classification_D3P4$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)

plot_d4p1 <- ggplot(dat_d4p1 , aes(dat_d4p1[,1], dat_d4p1[,2],color=classification$`Perching performance` ,label = classification$sp2,shape=classification$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification$`Perching performance`,fill=classification$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)
plot_d4p2 <- ggplot(dat_d4p2 , aes(dat_d4p2[,1], dat_d4p2[,2],color=classification_inter_terminal$`Perching performance` ,label = classification_inter_terminal$sp2,shape=classification_inter_terminal$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification_inter_terminal$`Perching performance`,fill=classification_inter_terminal$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)

plot_d4p3 <- ggplot(dat_d4p3 , aes(dat_d4p3[,1], dat_d4p3[,2],color=classification_D4P3$`Perching performance` ,label = classification_D4P3$sp2,shape=classification_D4P3$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification_D4P3$`Perching performance`,fill=classification_D4P3$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)

plot_d4p4 <- ggplot(dat_d4p4 , aes(dat_d4p4[,1], dat_d4p4[,2],color=classification_D3P3_D4P4$`Perching performance` ,label = classification_D3P3_D4P4$sp2,shape=classification_D3P3_D4P4$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification_D3P3_D4P4$`Perching performance`,fill=classification_D3P3_D4P4$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)

plot_d4p5 <- ggplot(dat_d4p5 , aes(dat_d4p5[,1], dat_d4p5[,2],color=classification_inter_terminal$`Perching performance` ,label = classification_inter_terminal$sp2,shape=classification_inter_terminal$Locomotion)) +
  geom_point(size=7) + geom_text_repel(color="black") + labs(title = "") + geom_encircle(aes(group=classification_inter_terminal$`Perching performance`,fill=classification_inter_terminal$`Perching performance`),colour="darkred", s_shape = 1, expand = 0,
                                                                                         alpha = 0.20, show.legend = FALSE)+xlab("PC1")+ylab("PC2")+scale_shape_manual(values=c(8, 19,17))+scale_color_manual(values=c("#193EE2", "#009E73", "#F0E442"))+scale_fill_manual(values=c("#56B4E9", "#009E73", "#F0E442"))+ theme_classic()+ geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.8,alpha=0.1) +geom_vline(xintercept = 0, linetype="dashed",color = "black", size=0.8,alpha=0.1)


##########################Plotting

####tarsometatarsus
plot_tmt 

####digit 1
plot_d1p1
plot_d1p2 
####digit 2
plot_d2p1
plot_d2p2
plot_d2p3
####digit 3
plot_d3p1
plot_d3p2
plot_d3p3
plot_d3p4
####digit 4
plot_d4p1
plot_d4p2
plot_d4p3
plot_d4p4
plot_d4p5







############################Phylomorphospace

PC_phy <-GPA_tmt$PCscores[,1:2]
row.names(PC_phy)<- label_tree$Sp

plotTree(tres,node.numbers=TRUE)
tres<-paintSubTree(tres,56,state="black")
tres<-paintSubTree(painted,62,state="blue")
tres<-paintSubTree(painted,77,state="red")
cols<-colnames(tres$mapped.edge)
names(cols)<-cols
plotSimmap(tres,cols,pts=FALSE)

phylomorphospace(tres,PC_phy,node.size=c(0,1.2),colors=cols,label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")



PC_phy <-GPA_D1P1$PCscores[,1:2]
row.names(PC_phy)<- label_tree$Sp
phylomorphospace(tres,PC_phy,node.size=c(0,1.2),label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")


PC_phy <-GPA_D1P2$PCscores[,1:2]
cbind(classification$sp2,label_tree$Sp)
row.names(PC_phy)<- label_tree$Sp
phylomorphospace(tres,PC_phy,node.size=c(0,1.2),label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")


PC_phy <-GPA_D2P1$PCscores[,1:2]
cbind(classification$sp2,label_tree$Sp)
row.names(PC_phy)<- label_tree$Sp
phylomorphospace(tres,PC_phy,node.size=c(0,1.2),label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")


PC_phy <-GPA_D2P2$PCscores[,1:2]
cbind(Classification_inter_terminal$sp2,label_D2P2_D3P2_D4P2_D4P5$Sp)
row.names(PC_phy)<- label_D2P2_D3P2_D4P2_D4P5$Sp
phylomorphospace(tres_D2P2_D3P2_D4P2_D4P5,PC_phy,node.size=c(0,1.2),label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")


PC_phy <-GPA_D2P3$PCscores[,1:2]
cbind(classification_D2P3$sp2,label_tree_D2P3$Sp)
row.names(PC_phy)<- label_tree_D2P3$Sp
phylomorphospace(tres_D2P3,PC_phy,node.size=c(0,1.2),label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")


PC_phy <-GPA_D3P1$PCscores[,1:2]
row.names(PC_phy)<- label_tree$Sp
phylomorphospace(tres,colors=cols,PC_phy,node.size=c(0,1.2),label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")


PC_phy <-GPA_D3P2$PCscores[,1:2]
cbind(Classification_inter_terminal$sp2,label_D2P2_D3P2_D4P2_D4P5$Sp)
row.names(PC_phy)<- label_D2P2_D3P2_D4P2_D4P5$Sp
phylomorphospace(tres_D2P2_D3P2_D4P2_D4P5,PC_phy,node.size=c(0,1.2),label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")


PC_phy <-GPA_D3P3$PCscores[,1:2]
cbind(classification_D3D3_D4D4$sp2,label_tree_D3P3_D4P4)
row.names(PC_phy)<- label_tree_D3P3_D4P4$Sp
phylomorphospace(tres_D3P3_D4P4,PC_phy,node.size=c(0,1.2),label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")


PC_phy <-GPA_D3P4$PCscores[,1:2]
cbind(classification_D3D4$sp2,label_tree_D3P4$Sp)
row.names(PC_phy)<- label_tree_D3P4$Sp
phylomorphospace(tres_D3P4,PC_phy,node.size=c(0,1.2),label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")


PC_phy <-GPA_D4P1$PCscores[,1:2]
cbind(classification$sp2,label_tree$Sp)
row.names(PC_phy)<- label_tree$Sp
phylomorphospace(tres,PC_phy,node.size=c(0,1.2),label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")


PC_phy <-GPA_D4P2$PCscores[,1:2]
cbind(Classification_inter_terminal$sp2,label_D2P2_D3P2_D4P2_D4P5$Sp)
row.names(PC_phy)<- label_D2P2_D3P2_D4P2_D4P5$Sp
phylomorphospace(tres_D2P2_D3P2_D4P2_D4P5,PC_phy,node.size=c(0,1.2),label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")


PC_phy <-GPA_D4P3$PCscores[,1:2]
cbind(classification_D4P3$sp2,label_D4P3$Sp)
row.names(PC_phy)<- label_D4P3$Sp
phylomorphospace(tres_D4P3,PC_phy,node.size=c(0,1.2),label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")


PC_phy <-GPA_D4P4$PCscores[,1:2]
cbind(classification_D3P3_D4P4$sp2,label_tree_D3D3_D4D4$Sp)
row.names(PC_phy)<- label_tree_D3P3_D4P4$Sp
phylomorphospace(tres_D3P3_D4P4,PC_phy,node.size=c(0,1.2),label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")


PC_phy <-GPA_D4P5$PCscores[,1:2]
cbind(Classification_inter_terminal$sp2,label_D2P2_D3P2_D4P2_D4P5$Sp)
row.names(PC_phy)<- label_D2P2_D3P2_D4P2_D4P5$Sp
phylomorphospace(tres_D2P2_D3P2_D4P2_D4P5,PC_phy,node.size=c(0,1.2),label="off",lwd=2,)
abline(h=0, v=0, col="grey75", lty="dashed")

 

