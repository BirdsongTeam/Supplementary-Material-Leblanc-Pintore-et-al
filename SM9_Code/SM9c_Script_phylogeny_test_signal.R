library(ape)
library(readr)
library(geomorph)


#############Loading label

label_tree <- read_csv("label_tree.csv")
label_tree_D2P3  <- label_tree[-c(12,21),]
label_tree_D3P4  <- label_tree[-c(12,29),]
label_D2P2_D3P2_D4P2_D4P5 <- label_tree[-c(55),]
label_D4P3<- label_tree[-c(55,23,40),]
label_tree_D3P3_D4P4<- label_tree[-c(55,29),]


############Loading tree

tr2<-read.tree(file="tree_furna")######## All proximal phalanges + Digit 1 Phalanx 2
tr2_D2P3<-read.tree(file="tree_furna_D2P3.txt")#################Digit 2 Phalanx 3
tr2_D3P4<-read.tree(file="tree_furna_D3P4.txt")#################Digit 3 Phalanx 4
tr2_D2P2_D3P2_D4P2_D4P5<-read.tree(file="tree_furna_D2P2_D3P2_D4P2_D4P5.txt")##################Digit 2 phalanx 2, Digit 3 Phalanx 2, Digit 4 Phalanx 2, Digit 4 Phalanx 5
tr2_D4P3<- read.tree(file="tree_furna_D4P3.txt")##########Digit  4 phalanx 3
tr2_D3P3_D4P4<-read.tree(file="tree_furna_D3P3_D4P4.txt")###########Digit 3 Phalanx 3, Digit 4 Phalanx 4





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
plot.phylo(tr2_D3P3_D4P4, type = "cladogram", use.edge.length = FALSE)
tr2_D3P3_D4P4[["edge.length"]]
tres_D3P3_D4P4<-multi2di(tr2_D3P3_D4P4)





#############tmt 

#####global test
PC_score_tmt <-GPA_tmt$PCscores
row.names(PC_score_tmt ) <- label_tree$Sp


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


compgls<-procD.pgls(PC_score_tmt[,1]~log(GPA_tmt$size), tres, iter=999)
compgls$aov.table

compgls<-procD.pgls(PC_score_tmt[,2]~log(GPA_tmt$size), tres, iter=999)
compgls$aov.table



PS.shape<-physignal(PC_score_tmt[,1], tres, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_tmt[,2], tres, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)


############# D1D1 Pcscore


PC_score_D1P1 <-GPA_D1P1$PCscores
row.names(PC_score_D1P1 ) <- label_tree$Sp


#####global test

PS.shape<-physignal(PC_score_D1P1, tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue


D1P1_GPA <- GPA_D1P1$rotated
dimnames(D1P1_GPA)[[3]] <- label_tree$Sp
D1P1_Csize <- as.matrix(GPA_D1P1$size)
row.names(D1P1_Csize) <- label_tree$Sp
compgls<-procD.pgls(D1P1_GPA~log(GPA_D1P1$size), tres , iter=999)
compgls$aov.table

#####PC test

compgls<-procD.pgls(PC_score_D1P1[,1]~log(GPA_D1P1$size), tres, iter=999)
compgls$aov.table

compgls<-procD.pgls(PC_score_D1P1[,2]~log(GPA_D1P1$size), tres, iter=999)
compgls$aov.table



cor.test(GPA_D1P1$PCscores[,1],log(GPA_D1P1$size))
cor.test(GPA_D1P1$PCscores[,2],log(GPA_D1P1$size))

PS.shape<-physignal(PC_score_D1P1[,1], tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue

PS.shape<-physignal(PC_score_D1P1[,2], tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue

############# D1D2 Pcscore

PC_score_D1P2 <-GPA_D1P2$PCscores
row.names(PC_score_D1P2) <- label_tree$Sp

#####global test

PS.shape<-physignal(PC_score_D1P2, tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue

D1P2_GPA <- GPA_D1P2$rotated
dimnames(D1P2_GPA)[[3]] <- label_tree$Sp
D1P2_Csize <- as.matrix(GPA_D1P2$size)
row.names(D1P2_Csize) <- label_tree$Sp
compgls<-procD.pgls(D1P2_GPA~GPA_D1P2$size, tres ,logsz = TRUE, iter=999)
compgls$aov.table


#####PC test


PS.shape<-physignal(PC_score_D1P2[,1], tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue

PS.shape<-physignal(PC_score_D1P2[,2], tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue


compgls<-procD.pgls(PC_score_D1P2[,1]~log(GPA_D1P2$size), tres, iter=999)
compgls$aov.table

compgls<-procD.pgls(PC_score_D1P2[,2]~log(GPA_D1P2$size), tres, iter=999)
compgls$aov.table



##################D2D1###############
PC_score_D2P1 <-GPA_D2P1$PCscores
row.names(PC_score_D2P1) <- label_tree$Sp

#####global test

D2P1_GPA <- GPA_D2P1$rotated
dimnames(D2P1_GPA)[[3]] <- label_tree$Sp
D2P1_Csize <- log(as.matrix(GPA_D2P1$size))
row.names(D2P1_Csize) <- label_tree$Sp

PS.shape<-physignal(PC_score_D2P1, tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue



compgls<-procD.pgls(D2P1_GPA~log(GPA_D2P1$size), tres ,logsz = TRUE, iter=999)
compgls$aov.table

#####PC test

PS.shape<-physignal(PC_score_D2P1[,1], tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue

PS.shape<-physignal(PC_score_D2P1[,2], tres, iter=1000)
PS.shape$phy.signal
PS.shape$pvalue


compgls<-procD.pgls(PC_score_D2P1[,1]~log(GPA_D1P2$size), tres, iter=999)
compgls$aov.table

compgls<-procD.pgls(PC_score_D2P1[,2]~log(GPA_D1P2$size), tres, iter=999)
compgls$aov.table





################D2D2
PC_score_D2P2 <-GPA_D2P2$PCscores
row.names(PC_score_D2P2) <- label_D2P2_D3P2_D4P2_D4P5$Sp

#####global test

D2P2_GPA <- GPA_D2P2$rotated
dimnames(D2P2_GPA)[[3]] <- label_D2P2_D3P2_D4P2_D4P5$Sp
D2P2_Csize <- log(as.matrix(GPA_D2P2$size))
row.names(D2P2_Csize) <- label_D2P2_D3P2_D4P2_D4P5$Sp
compgls<-procD.pgls(D2P2_GPA~GPA_D2P2$size, phy=tres_D2P2_D3P2_D4P2_D4P5 ,logsz = TRUE, iter=999)
compgls$aov.table
summary(compgls)


PS.shape<-physignal(PC_score_D2D2, tres_D2P2_D3P2_D4P2_D4P5, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)


#####PC test

compgls<-procD.pgls(PC_score_D2P2[,1]~log(GPA_D2P2$size), tres_D2P2_D3P2_D4P2_D4P5, iter=999)
compgls$aov.table

compgls<-procD.pgls(PC_score_D2P2[,2]~log(GPA_D2P2$size), tres_D2P2_D3P2_D4P2_D4P5, iter=999)
compgls$aov.table


PS.shape<-physignal(PC_score_D2D2[,1],tres_D2P2_D3P2_D4P2_D4P5,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)


PS.shape<-physignal(PC_score_D2D2[,2], tres_D2P2_D3P2_D4P2_D4P5, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

################D2P3
PC_score_D2P3 <-GPA_D2P3$PCscores
row.names(PC_score_D2P3) <- label_tree_D2P3$Sp

#####global test

PS.shape<-physignal(PC_score_D2P3,tres_D2P3,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

D2P3_GPA <- GPA_D2P3$rotated
dimnames(D2P3_GPA)[[3]] <- label_tree_D2P3$Sp
D2P3_Csize <- log(as.matrix(GPA_D2P3$size))
row.names(D2P3_Csize) <- label_tree_D2P3$Sp
compgls<-procD.pgls(D2P3_GPA~GPA_D2P3$size, tres_D2P3 ,logsz = TRUE, iter=999)
compgls$aov.table


#####PC test

compgls<-procD.pgls(PC_score_D2P3[,1]~log(GPA_D2P3$size), tres_D2P3, iter=999)
compgls$aov.table

compgls<-procD.pgls(PC_score_D2P3[,2]~log(GPA_D2P3$size), tres_D2P3, iter=999)
compgls$aov.table

PS.shape<-physignal(PC_score_D2P3[,1], tres_D2P3, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D2P3[,2], tres_D2P3, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

################D3D1
PC_score_D3P1 <-GPA_D3P1$PCscores
row.names(PC_score_D3P1) <- label_tree$Sp


#####global test
D3P1_GPA <- GPA_D3P1$rotated
dimnames(D3P1_GPA)[[3]] <- label_tree$Sp
D3P1_Csize <- log(as.matrix(GPA_D3P1$size))
row.names(D3P1_Csize) <- label_tree$Sp
compgls<-procD.pgls(D3P1_GPA~GPA_D3P1$size, tres ,logsz = TRUE, iter=999)
compgls$aov.table


PS.shape<-physignal(PC_score_D3P1,tres,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

#####PC test


compgls<-procD.pgls(PC_score_D3P1[,1]~log(GPA_D3P1$size), tres, iter=999)
compgls$aov.table

compgls<-procD.pgls(PC_score_D3P1[,2]~log(GPA_D3P1$size), tres, iter=999)
compgls$aov.table



PS.shape<-physignal(PC_score_D3P1[,1], tres, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D3P1[,2], tres, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

################D3P2
PC_score_D3P2 <-GPA_D3P2$PCscores
row.names(PC_score_D3P2) <- label_D2P2_D3P2_D4P2_D4P5$Sp


#####global test
D3P2_GPA <- GPA_D3P2$rotated
dimnames(D3P2_GPA)[[3]] <- label_D2P2_D3P2_D4P2_D4P5$Sp
D3P2_Csize <- log(as.matrix(GPA_D3P2$size))
row.names(D3P2_Csize) <- label_D2P2_D3P2_D4P2_D4P5$Sp
compgls<-procD.pgls(D3P2_GPA~GPA_D3P2$size, tres_D2P2_D3P2_D4P2_D4P5 ,logsz = TRUE, iter=999)
compgls$aov.table


PS.shape<-physignal(PC_score_D3P2,tres_D2P2_D3P2_D4P2_D4P5,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)


#####PC test
PS.shape<-physignal(PC_score_D3P2[,1], tres_D2P2_D3P2_D4P2_D4P5 , iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D3P2[,2], tres_D2P2_D3P2_D4P2_D4P5 , iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

compgls<-procD.pgls(PC_score_D3P2[,1]~log(GPA_D3P2$size),tres_D2P2_D3P2_D4P2_D4P5, iter=999)
compgls$aov.table

compgls<-procD.pgls(PC_score_D3P2[,2]~log(GPA_D3P2$size),tres_D2P2_D3P2_D4P2_D4P5, iter=999)
compgls$aov.table


################D3P3
PC_score_D3P3 <-GPA_D3P3$PCscores
row.names(PC_score_D3P3) <- label_tree_D3P3_D4P4$Sp

tres_D3P3_D3P4$tip.label

#####global test
D3P3_GPA <- GPA_D3P3$rotated
dimnames(D3P3_GPA)[[3]] <- label_tree_D3P3_D4P4$Sp
D3P3_Csize <- log(as.matrix(GPA_D3P3$size))
row.names(D3P3_Csize) <- label_tree_D3P3_D4P4$Sp
compgls<-procD.pgls(D3P3_GPA~GPA_D3P3$size, tres_D3P3_D4P4 ,logsz = TRUE, iter=999)
compgls$aov.table



PS.shape<-physignal(PC_score_D3P3,tres_D3P3_D4P4,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

#####PC test
PS.shape<-physignal(PC_score_D3P3[,1], tres_D3P3_D4P4 , iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D3P3[,2], tres_D3P3_D4P4, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

compgls<-procD.pgls(PC_score_D3P3[,1]~log(GPA_D3P3$size),tres_D2P2_D3P2_D4P2_D4P5, iter=999)
compgls$aov.table

compgls<-procD.pgls(PC_score_D3P3[,2]~log(GPA_D3P3$size),tres_D2P2_D3P2_D4P2_D4P5, iter=999)
compgls$aov.table



################D3P4
PC_score_D3P4 <-GPA_D3P4$PCscores
row.names(PC_score_D3P4) <- label_tree_D3P4$Sp


#####global test
D3P4_GPA <- GPA_D3P4$rotated
dimnames(D3P4_GPA)[[3]] <- label_tree_D3P4$Sp
D3P4_Csize <- log(as.matrix(GPA_D3P4$size))
row.names(D3P4_Csize) <- label_tree_claw2$Sp
compgls<-procD.pgls(D3P4_GPA~GPA_D3P4$size, tr2_D3P4 ,logsz = TRUE, iter=999)
compgls$aov.table


PS.shape<-physignal(PC_score_D3P4, tr2_D3P4, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)


#####global test
PS.shape<-physignal(PC_score_D3P4[,1], tr2_D3P4, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D3P4[,2], tr2_D3P4, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)


compgls<-procD.pgls(PC_score_D3P4[,1]~log(GPA_D3P4$size),tr2_D3P4, iter=999)
compgls$aov.table

compgls<-procD.pgls(PC_score_D3P4[,2]~log(GPA_D3P4$size),tr2_D3P4, iter=999)
compgls$aov.table






################D4D1
PC_score_D4P1 <-GPA_D4P1$PCscores
row.names(PC_score_D4P1) <- label_tree$Sp


#####global test
D4P1_GPA <- GPA_D4P1$rotated
dimnames(D4P1_GPA)[[3]] <- label_tree$Sp
D4P1_Csize <- log(as.matrix(GPA_D4P1$size))
row.names(D4P1_Csize) <- label_tree$Sp
compgls<-procD.pgls(D4P1_GPA~GPA_D4P1$size, tres ,logsz = TRUE, iter=999)
compgls$aov.table



PS.shape<-physignal(PC_score_D4P1,tres,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

#####PC test
PS.shape<-physignal(PC_score_D4P1[,1], tres, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D4P1[,2], tres, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)


compgls<-procD.pgls(PC_score_D4P1[,1]~log(GPA_D4P1$size), tres, iter=999)
compgls$aov.table


compgls<-procD.pgls(PC_score_D4P1[,2]~log(GPA_D4P1$size), tres, iter=999)
compgls$aov.table



################D4P2
PC_score_D4P2 <-GPA_D4P2$PCscores
row.names(PC_score_D4P2) <- label_D2P2_D3P2_D4P2_D4P5$Sp


#####global test
D4P2_GPA <- GPA_D4P2$rotated
dimnames(D4P2_GPA)[[3]] <- label_D2P2_D3P2_D4P2_D4P5$Sp
D4P2_Csize <- log(as.matrix(GPA_D4P2$size))
row.names(D4P2_Csize) <- label_D2P2_D3P2_D4P2_D4P5$Sp
compgls<-procD.pgls(D4P2_GPA~GPA_D4P2$size, tres_D2P2_D3P2_D4P2_D4P5 ,logsz = TRUE, iter=999)
compgls$aov.table



PS.shape<-physignal(PC_score_D4P2,tres_D2P2_D3P2_D4P2_D4P5,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

#####gPC test
PS.shape<-physignal(PC_score_D4P2[,1], tres_D2P2_D3P2_D4P2_D4P5, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D4P2[,2], tres_D2P2_D3P2_D4P2_D4P5, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

compgls<-procD.pgls(PC_score_D4P2[,1]~log(GPA_D4P2$size), tres_D2P2_D3P2_D4P2_D4P5, iter=999)
compgls$aov.table

compgls<-procD.pgls(PC_score_D4P2[,2]~log(GPA_D4P2$size), tres_D2P2_D3P2_D4P2_D4P5, iter=999)
compgls$aov.table

################D4P3
PC_score_D4P3 <-GPA_D4P3$PCscores
row.names(PC_score_D4P3) <- label_D4P3$Sp

#####global test
PS.shape<-physignal(PC_score_D4P3,tres_D4P3,iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

D4P3_GPA <- GPA_D4P3$rotated
dimnames(D4P3_GPA)[[3]] <- label_D4P3$Sp
D4P3_Csize <- log(as.matrix(GPA_D4P3$size))
row.names(D4P3_Csize) <- label_D4P3$Sp
compgls<-procD.pgls(D4P3_GPA~GPA_D4P3$size, tres_D4P3 ,logsz = TRUE, iter=999)
compgls$aov.table


#####PC test

compgls<-procD.pgls(PC_score_D4P3[,1]~log(GPA_D4P3$size), tres_D4P3, iter=999)
compgls$aov.table

compgls<-procD.pgls(PC_score_D4P3[,2]~log(GPA_D4P3$size), tres_D4P3, iter=999)
compgls$aov.table

PS.shape<-physignal(PC_score_D4P3[,1], tres_D4P3, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D4P3[,2], tres_D4P3, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)




################D4P4
PC_score_D4P4 <-GPA_D4P4$PCscores
row.names(PC_score_D4P4) <- label_tree_D3P3_D4P4$Sp



#####global test
D4P4_GPA <- GPA_D4P4$rotated
dimnames(D4P4_GPA)[[3]] <- label_tree_D3P3_D4P4$Sp
D4P4_Csize <- log(as.matrix(GPA_D4P4$size))
row.names(D4P4_Csize) <- label_tree_D3P3_D4P4$Sp
compgls<-procD.pgls(D4P4_GPA~GPA_D4P4$size, tres_D3P3_D4P4 ,logsz = TRUE, iter=999)
compgls$aov.table



PS.shape<-physignal(PC_score_D4P4, tres_D3P3_D4P4, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

#####PC test
PS.shape<-physignal(PC_score_D4P4[,1], tres_D3P3_D3P4, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D4P4[,2], tres_D3P3_D3P4, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

compgls<-procD.pgls(PC_score_D4P4[,1]~log(GPA_D4P4$size),  tres_D3P3_D4P4, iter=999)
compgls$aov.table

compgls<-procD.pgls(PC_score_D4P4[,2]~log(GPA_D4P4$size),  tres_D3P3_D4P4, iter=999)
compgls$aov.table

################D4P5
PC_score_D4P5 <-GPA_D4P5$PCscores
row.names(PC_score_D4P5) <- label_D2P2_D3P2_D4P2_D4P5$Sp

#####global test
D4P5_GPA <- GPA_D4P5$rotated
dimnames(D4P5_GPA)[[3]] <- label_D2P2_D3P2_D4P2_D4P5$Sp
D4P5_Csize <- log(as.matrix(GPA_D4P5$size))
row.names(D4P5_Csize) <- label_D2P2_D3P2_D4P2_D4P5$Sp
compgls<-procD.pgls(D4P5_GPA~GPA_D4P5$size, tres_D2P2_D3P2_D4P2_D4P5,logsz = TRUE, iter=999)
compgls$aov.table



PS.shape<-physignal(PC_score_D4P5, tres_D2P2_D3P2_D4P2_D4P5, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)



#####PC test
PS.shape<-physignal(PC_score_D4P5[,1], tres_D2P2_D3P2_D4P2_D4P5, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

PS.shape<-physignal(PC_score_D4P5[,2], tres_D2P2_D3P2_D4P2_D4P5, iter=1000)
print(PS.shape$phy.signal)
print(PS.shape$pvalue)

compgls<-procD.pgls(PC_score_D4P5[,1]~log(GPA_D4P5$size),tres_D2P2_D3P2_D4P2_D4P5, iter=999)
compgls$aov.table

compgls<-procD.pgls(PC_score_D4P5[,2]~log(GPA_D4P5$size),tres_D2P2_D3P2_D4P2_D4P5, iter=999)
compgls$aov.table

#