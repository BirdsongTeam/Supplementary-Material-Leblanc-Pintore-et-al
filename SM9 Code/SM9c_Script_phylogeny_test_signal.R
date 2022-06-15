library(ape)
library(readr)
library(geomorph)


label_tree <- read_csv("label_tree.csv")
label_tree_D2P3  <- label_tree[-c(12,21),]
label_tree_D3P4  <- label_tree[-c(12,29),]
label_D2P2_D3P2_D4P2_D4P5 <- label_tree[-c(55),]
label_D4P3<- label_tree[-c(55,23,40),]
label_tree_D3P3_D4P4<- label_tree[-c(55,29),]






tr2<-read.tree(file="tree_furna")######## All proximal phalanges + Digit 1 Phalanx 2
tr2_D2P3<-read.tree(file="tree_furna_D2P3.txt")#################Digit 2 Phalanx 3
tr2_D3P4<-read.tree(file="tree_furna_D3P4.txt")#################Digit 3 Phalanx 4
tr2_D2P2_D3P2_D4P2_D4P5<-read.tree(file="tree_furna_D2P2_D3P2_D4P2_D4P5.txt")##################Digit 2 phalanx 2, Digit 3 Phalanx 2, Digit 4 Phalanx 2, Digit 4 Phalanx 5
tr2_D4P3<- read.tree(file="tree_furna_D4P3.txt")##########Digit  4 phalanx 3
tr2_D3P3_D4P4<-read.tree(file="tree_furna_D3D3_D4D4.txt")###########Digit 3 Phalanx 3, Digit 4 Phalanx 4





############### formalize tree 

tr2=compute.brlen(tr2, method=10)
plot.phylo(tr2, type = "cladogram", use.edge.length = FALSE)
tr2[["edge.length"]]
tr2$tip.label %in% label_tree$Sp
tres<-multi2di(tr2)

tr2_D2P3=compute.brlen(tr2_D2P3, method=10)
plot.phylo(tr2_D2P3, type = "cladogram", use.edge.length = FALSE)
tr2_D2P3[["edge.length"]]
tr2_D2P3$tip.label %in% label_tree_claw2$Sp
tres_D2P3<-multi2di(tr2_D2P3)

tr2_D3P4=compute.brlen(tr2_D3P4, method=10)
plot.phylo(tr2_D3P4, type = "cladogram", use.edge.length = FALSE)
tr2_D3P4[["edge.length"]]
tr2_D3P4$tip.label %in% label_tree_D3P4$Sp
tres_D3P4<-multi2di(tr2_D3P4)

tr2_D2P2_D3P2_D4P2_D4P5=compute.brlen(tr2_D2P2_D3P2_D4P2_D4P5, method=10)
plot.phylo(tr2_D2P2_D3P2_D4P2_D4P5, type = "cladogram", use.edge.length = FALSE)
tr2_D2P2_D3P2_D4P2_D4P5[["edge.length"]]
tr2_D2P2_D3P2_D4P2_D4P5$tip.label %in% label_D2P2_D3P2_D4P2_D4P5$Sp
tres_D2P2_D3P2_D4P2_D4P5<-multi2di(tr2_D2P2_D3P2_D4P2_D4P5)

tr2_D4P3=compute.brlen(tr2_D4P3, method=10)
plot.phylo(tr2_D4P3, type = "cladogram", use.edge.length = FALSE)
tr2_D4P3[["edge.length"]]
tr2_D4P3$tip.label %in% label_D2P2_D3P2_D4P2_D4P5$Sp
tres_D4P3<-multi2di(tr2_D4P3)

tr2_D3P3_D4P4=compute.brlen(tr2_D3P3_D4P4, method=10)
plot.phylo(tr2_D3D3_D4D4, type = "cladogram", use.edge.length = FALSE)
tr2_D3P3_D4P4[["edge.length"]]
tres_D3P3_D4P4<-multi2di(tr2_D3P3_D4P4)





#############tmt 

#####global test

PS.shape<-physignal(PC_score_tmt, tres, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

tmt_GPA <- GPA_tmt$rotated
dimnames(tmt_GPA)[[3]] <- label_tree$Sp
tmt_Csize <- as.matrix(GPA_tmt$size)
row.names(tmt_Csize) <- label_tree$Sp
compgls<-procD.pgls(tmt_GPA~GPA_tmt$size, tres ,logsz = TRUE, iter=999)
compgls$aov.table

######PC test

cor.test(GPA_tmt$PCscores[,1],GPA_tmt$size)
cor.test(GPA_tmt$PCscores[,2],GPA_tmt$size)

PS.shape<-physignal(PC_score_tmt[,1], tres, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_tmt[,2], tres, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)


############# D1D1 Pcscore


PC_score_D1D1 <-GPA_D1D1$PCscores
row.names(PC_score_D1D1 ) <- label_tree$Sp


#####global test

PS.shape<-physignal(PC_score_D1D1, tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue


D1D1_GPA <- GPA_D1D1$rotated
dimnames(D1D1_GPA)[[3]] <- label_tree$Sp
D1D1_Csize <- as.matrix(GPA_D1D1$size)
row.names(D1D1_Csize) <- label_tree$Sp
compgls<-procD.pgls(D1D1_GPA~log(GPA_D1D1$size), tres , iter=999)
compgls$aov.table

#####PC test

cor.test(GPA_D1D1$PCscores[,1],GPA_D1D1$size)
cor.test(GPA_D1D1$PCscores[,2],GPA_D1D1$size)

PS.shape<-physignal(PC_score_D1D1[,1], tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue

PS.shape<-physignal(PC_score_D1D1[,2], tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue

############# D1D2 Pcscore

PC_score_D1D2 <-GPA_D1D2$PCscores
row.names(PC_score_D1D2) <- label_tree$Sp

#####global test

PS.shape<-physignal(PC_score_D1D2, tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue

D1D2_GPA <- GPA_D1D2$rotated
dimnames(D1D2_GPA)[[3]] <- label_tree$Sp
D1D2_Csize <- as.matrix(GPA_D1D2$size)
row.names(D1D2_Csize) <- label_tree$Sp
compgls<-procD.pgls(D1D2_GPA~GPA_D1D2$size, tres ,logsz = TRUE, iter=999)
compgls$aov.table


#####PC test


PS.shape<-physignal(PC_score_D1D2[,1], tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue

PS.shape<-physignal(PC_score_D1D2[,2], tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue



cor.test(GPA_D1D2$PCscores[,1],GPA_D1D2$size)
cor.test(GPA_D1D2$PCscores[,2],GPA_D1D2$size)


##################D2D1###############
PC_score_D2D1 <-GPA_D2D1$PCscores
row.names(PC_score_D2D1) <- label_tree$Sp

#####global test

D2D1_GPA <- GPA_D2D1$rotated
dimnames(D2D1_GPA)[[3]] <- label_tree$Sp
D2D1_Csize <- log(as.matrix(GPA_D2D1$size))
row.names(D2D1_Csize) <- label_tree$Sp

PS.shape<-physignal(PC_score_D2D1, tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue



compgls<-procD.pgls(D2D1_GPA~GPA_D2D1$size, tres ,logsz = TRUE, iter=999)
compgls$aov.table

#####PC test

PS.shape<-physignal(PC_score_D2D1[,1], tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue

PS.shape<-physignal(PC_score_D2D1[,2], tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue


cor.test(GPA_D2D1$PCscores[,1],GPA_D2D1$size)
cor.test(GPA_D2D1$PCscores[,2],GPA_D2D1$size)

################D2D2
PC_score_D2D2 <-GPA_D2D2$PCscores
row.names(PC_score_D2D2) <- label_D2P2_D3P2_D4P2_D4P5$Sp

#####global test

D2D2_GPA <- GPA_D2D2$rotated
dimnames(D2D2_GPA)[[3]] <- label_D2P2_D3P2_D4P2_D4P5$Sp
D2D2_Csize <- log(as.matrix(GPA_D2D2$size))
row.names(D2D2_Csize) <- label_D2P2_D3P2_D4P2_D4P5$Sp
compgls<-procD.pgls(D2D2_GPA~GPA_D2D2$size, phy=tres_D2P2_D3P2_D4P2_D4P5 ,logsz = TRUE, iter=999)
compgls$aov.table
summary(compgls)


PS.shape<-physignal(PC_score_D2D2, tres_D2P2_D3P2_D4P2_D4P5, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)


#####PC test

cor.test(GPA_D2D2$PCscores[,1],D2D2_Csize)
cor.test(GPA_D2D2$PCscores[,2],D2D2_Csize)

PS.shape<-physignal(PC_score_D2D2[,1],tres_D2P2_D3P2_D4P2_D4P5,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)


PS.shape<-physignal(PC_score_D2D2[,2], tres_D2P2_D3P2_D4P2_D4P5, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

################D2D3
PC_score_D2D3 <-GPA_D2D3$PCscores
row.names(PC_score_D2D3) <- label_tree_tree_D2P3$Sp

#####global test

PS.shape<-physignal(PC_score_D2D3,tres_D2P3,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

D2D3_GPA <- GPA_D2D3$rotated
dimnames(D2D3_GPA)[[3]] <- label_tree_tree_D2P3$Sp
D2D3_Csize <- log(as.matrix(GPA_D2D3$size))
row.names(D2D3_Csize) <- label_tree_claw2$Sp
compgls<-procD.pgls(D2D3_GPA~GPA_D2D3$size, tres_D2P3 ,logsz = TRUE, iter=999)
compgls$aov.table


#####PC test


cor.test(GPA_D2D3$PCscores[,1],D2D3_Csize)
cor.test(GPA_D2D3$PCscores[,2],D2D3_Csize)

PS.shape<-physignal(PC_score_D2D3[,1], tres_D2P3, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D2D3[,2], tres_D2P3, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

################D3D1
PC_score_D3D1 <-GPA_D3D1$PCscores
row.names(PC_score_D3D1) <- label_tree$Sp


#####global test
D3D1_GPA <- GPA_D3D1$rotated
dimnames(D3D1_GPA)[[3]] <- label_tree$Sp
D3D1_Csize <- log(as.matrix(GPA_D3D1$size))
row.names(D3D1_Csize) <- label_tree$Sp
compgls<-procD.pgls(D3D1_GPA~GPA_D3D1$size, tres ,logsz = TRUE, iter=999)
compgls$aov.table


PS.shape<-physignal(PC_score_D3D1,tres,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

#####PC test

cor.test(GPA_D3D1$PCscores[,1],D3D1_Csize)
cor.test(GPA_D3D1$PCscores[,2],D3D1_Csize)

PS.shape<-physignal(PC_score_D3D1[,1], tres, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D3D1[,2], tres, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

################D3D2
PC_score_D3D2 <-GPA_D3D2$PCscores
row.names(PC_score_D3D2) <- label_D2P2_D3P2_D4P2_D4P5$Sp


#####global test
D3D2_GPA <- GPA_D3D2$rotated
dimnames(D3D2_GPA)[[3]] <- label_D2P2_D3P2_D4P2_D4P5$Sp
D3D2_Csize <- log(as.matrix(GPA_D3D2$size))
row.names(D3D2_Csize) <- label_D2P2_D3P2_D4P2_D4P5$Sp
compgls<-procD.pgls(D3D2_GPA~GPA_D3D2$size, tres_D2P2_D3P2_D4P2_D4P5 ,logsz = TRUE, iter=999)
compgls$aov.table


PS.shape<-physignal(PC_score_D3D2,tres,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)


#####PC test
PS.shape<-physignal(PC_score_D3D2[,1], tres_D2P2_D3P2_D4P2_D4P5 , iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D3D2[,2], tres_D2P2_D3P2_D4P2_D4P5 , iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

cor.test(GPA_D3D2$PCscores[,1],D3D2_Csize)
cor.test(GPA_D3D2$PCscores[,2],D3D2_Csize)


################D3D3
PC_score_D3D3 <-GPA_D3D3$PCscores
row.names(PC_score_D3D3) <- label_tree_D3D3_D4D4$Sp

tres_D3D3_D3D4$tip.label

#####global test
D3D3_GPA <- GPA_D3D3$rotated
dimnames(D3D3_GPA)[[3]] <- label_tree_D3D3_D4D4$Sp
D3D3_Csize <- log(as.matrix(GPA_D3D3$size))
row.names(D3D3_Csize) <- label_tree_D3D3_D4D4$Sp
compgls<-procD.pgls(D3D3_GPA~GPA_D3D3$size, tres_D3D3_D4D4 ,logsz = TRUE, iter=999)
compgls$aov.table



PS.shape<-physignal(PC_score_D3D3,tres_D3D3_D4D4,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

#####PC test
PS.shape<-physignal(PC_score_D3D3[,1], tres_D3D3_D4D4 , iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D3D3[,2], tres_D3D3_D4D4, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

cor.test(GPA_D3D3$PCscores[,1],D3D3_Csize)
cor.test(GPA_D3D3$PCscores[,2],D3D3_Csize)


################D3D4
PC_score_D3D4 <-GPA_D3D4$PCscores
row.names(PC_score_D3D4) <- label_tree_D3P4$Sp


#####global test
D3D4_GPA <- GPA_D3D4$rotated
dimnames(D3D4_GPA)[[3]] <- label_tree_D3P4$Sp
D3D4_Csize <- log(as.matrix(GPA_D3D4$size))
row.names(D3D4_Csize) <- label_tree_claw2$Sp
compgls<-procD.pgls(D3D4_GPA~GPA_D3D4$size, tr2_D3P4 ,logsz = TRUE, iter=999)
compgls$aov.table


PS.shape<-physignal(PC_score_D3D4, tr2_D3P4, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)


#####global test
PS.shape<-physignal(PC_score_D3D4[,1], tr2_D3P4, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D3D4[,2], tr2_D3P4, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)


cor.test(GPA_D3D4$PCscores[,1],D3D4_Csize)
cor.test(GPA_D3D4$PCscores[,2],D3D4_Csize)





################D4D1
PC_score_D4D1 <-GPA_D4D1$PCscores
row.names(PC_score_D4D1) <- label_tree$Sp


#####global test
D4D1_GPA <- GPA_D4D1$rotated
dimnames(D4D1_GPA)[[3]] <- label_tree$Sp
D4D1_Csize <- log(as.matrix(GPA_D4D1$size))
row.names(D4D1_Csize) <- label_tree$Sp
compgls<-procD.pgls(D4D1_GPA~GPA_D4D1$size, tres ,logsz = TRUE, iter=999)
compgls$aov.table



PS.shape<-physignal(PC_score_D4D1,tres,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

#####PC test
PS.shape<-physignal(PC_score_D4D1[,1], tres, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D4D1[,2], tres, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

cor.test(GPA_D4D1$PCscores[,1],D4D1_Csize)
cor.test(GPA_D4D1$PCscores[,2],D4D1_Csize)
################D4D2
PC_score_D4D2 <-GPA_D4D2$PCscores
row.names(PC_score_D4D2) <- label_D2P2_D3P2_D4P2_D4P5$Sp


#####global test
D4D2_GPA <- GPA_D4D2$rotated
dimnames(D4D2_GPA)[[3]] <- label_D2P2_D3P2_D4P2_D4P5$Sp
D4D2_Csize <- log(as.matrix(GPA_D4D2$size))
row.names(D4D2_Csize) <- label_D2P2_D3P2_D4P2_D4P5$Sp
compgls<-procD.pgls(D4D2_GPA~GPA_D4D2$size, tres_D2P2_D3P2_D4P2_D4P5 ,logsz = TRUE, iter=999)
compgls$aov.table



PS.shape<-physignal(PC_score_D4D2,tres_D2P2_D3P2_D4P2_D4P5,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

#####gPC test
PS.shape<-physignal(PC_score_D4D2[,1], tres_D2P2_D3P2_D4P2_D4P5, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D4D2[,2], tres_D2P2_D3P2_D4P2_D4P5, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

cor.test(GPA_D4D2$PCscores[,1],D4D2_Csize)
cor.test(GPA_D4D2$PCscores[,2],D4D2_Csize)


################D4D3
PC_score_D4D3 <-GPA_D4D3$PCscores
row.names(PC_score_D4D3) <- label_D4P3$Sp

#####global test
PS.shape<-physignal(PC_score_D4D3,tres_D4P3,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

D4D3_GPA <- GPA_D4D3$rotated
dimnames(D4D3_GPA)[[3]] <- label_D4P3$Sp
D4D3_Csize <- log(as.matrix(GPA_D4D3$size))
row.names(D4D3_Csize) <- label_D4P3$Sp
compgls<-procD.pgls(D4D3_GPA~GPA_D4D3$size, tres_D4P3 ,logsz = TRUE, iter=999)
compgls$aov.table


#####PC test

cor.test(GPA_D4D3$PCscores[,1],D4D3_Csize)
cor.test(GPA_D4D3$PCscores[,2],D4D3_Csize)

PS.shape<-physignal(PC_score_D4D3[,1], tres_D4P3, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D4D3[,2], tres_D4P3, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)




################D4D4
PC_score_D4D4 <-GPA_D4D4$PCscores
row.names(PC_score_D4D4) <- label_tree_D3D3_D4D4$Sp



#####global test
D4D4_GPA <- GPA_D4D4$rotated
dimnames(D4D4_GPA)[[3]] <- label_tree_D3D3_D4D4$Sp
D4D4_Csize <- log(as.matrix(GPA_D4D4$size))
row.names(D4D4_Csize) <- label_tree_D3D3_D4D4$Sp
compgls<-procD.pgls(D4D4_GPA~GPA_D4D4$size, tres_D3D3_D4D4 ,logsz = TRUE, iter=999)
compgls$aov.table



PS.shape<-physignal(PC_score_D4D4, tres_D3D3_D4D4, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

#####PC test
PS.shape<-physignal(PC_score_D4D4[,1], tres_D3D3_D3D4, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D4D4[,2], tres_D3D3_D3D4, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

cor.test(GPA_D4D4$PCscores[,1],D4D4_Csize)
cor.test(GPA_D4D4$PCscores[,2],D4D4_Csize)


################D4D5
PC_score_D4D5 <-GPA_D4D5$PCscores
row.names(PC_score_D4D5) <- label_D2P2_D3P2_D4P2_D4P5$Sp

#####global test
D4D5_GPA <- GPA_D4D5$rotated
dimnames(D4D5_GPA)[[3]] <- label_D2P2_D3P2_D4P2_D4P5$Sp
D4D5_Csize <- log(as.matrix(GPA_D4D5$size))
row.names(D4D5_Csize) <- label_D2P2_D3P2_D4P2_D4P5$Sp
compgls<-procD.pgls(D4D5_GPA~GPA_D4D5$size, tres_D2P2_D3P2_D4P2_D4P5,logsz = TRUE, iter=999)
compgls$aov.table



PS.shape<-physignal(PC_score_D4D5, tres_D2P2_D3P2_D4P2_D4P5, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)



#####PC test
PS.shape<-physignal(PC_score_D4D5[,1], tres_D2P2_D3P2_D4P2_D4P5, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D4D5[,2], tres_D2P2_D3P2_D4P2_D4P5, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

cor.test(GPA_D4D5$PCscores[,1],D4D5_Csize)
cor.test(GPA_D4D5$PCscores[,2],D4D5_Csize)

#
