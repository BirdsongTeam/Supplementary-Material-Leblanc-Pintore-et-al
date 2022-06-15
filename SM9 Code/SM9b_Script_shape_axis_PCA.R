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

##################################### LOADING 

###########Loading the GPA alignment for each Digit X phalanx X 
load("D1P1")
D1P1 <- D1D1
load("D1P2")
D1P2 <- D1D2
load("D2P1")
D2P1 <- D2D1
load("D2P2")
D2P2 <- D2D2
load("D2P3")
D2P3 <- D2D3
load("D3P1")
D3P1 <- D3D1
load("D3P2")
D3P2 <- D3D2
load("D3P3")
D3P3 <- D3D3
load("D3P4")
D3P4 <- D3D4
load("D4P1")
D4P1 <- D4D1
load("D4P2")
D4P2 <- D4D2
load("D4P3")
D4P3 <- D4D3
load("D4P4")
D4P4 <- D4D4
load("D4P5")
D4P5 <- D4D5
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

###################################### put the correct directory pathway
template.tmt <- ply2mesh(filename = "./A8409_pd_tmt.ply")
template.D1P1 <- ply2mesh(filename = "./A8409_pd_d1_1.ply")
template.D2P1 <- ply2mesh(filename = "./A8409_pd_d2_1.ply")
template.D3P1 <- ply2mesh(filename = "./A8409_pd_d3_1.ply")
template.D4P1<- ply2mesh(filename = "./A8409_pd_d4_1.ply")

template.D1P2 <- ply2mesh(filename = "./A8409_pd_d1_2.ply")
template.D2P3 <- ply2mesh(filename = "./A8409_pd_d2_3.ply")
template.D3P4 <- ply2mesh(filename = "./A8409_pd_d3_4.ply")
template.D4P5 <- ply2mesh(filename = "./A8409_pd_d4_5.ply")


template.D2P2 <- ply2mesh(filename = "./A8409_pd_d2_2.ply")

template.D3P2 <- ply2mesh(filename = "./A8409_pd_d3_2.ply")
template.D3P3 <- ply2mesh(filename = "./A8409_pd_d3_3.ply")
template.D4P2 <- ply2mesh(filename = "./A8409_pd_d4_2.ply")
template.D4P3 <- ply2mesh(filename = "./A8409_pd_d4_3.ply")
template.D4P4 <- ply2mesh(filename = "./A8409_pd_d4_4.ply")





#############tmt

meanshape<-GPA_tmt$mshape
meshape_tmt<-tps3d(template.tmt,ptsarray[,,1], as.matrix(GPA_tmt$mshape))
plot3d(meshape_tmt, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_tmt<-GPA_tmt$PCscores[,1]
PC2_tmt<-GPA_tmt$PCscores[,2]
pPC1_tmt<-shape.predictor(GPA_tmt$rotated, x=PC1_tmt, intercept=FALSE, min=min(PC1_tmt), max=max(PC1_tmt))
pPC2_tmt<-shape.predictor(GPA_tmt$rotated, x=PC2_tmt, intercept=FALSE, min=min(PC2_tmt), max=max(PC2_tmt))
##PC1
PC1minshape_tmt<-pPC1_tmt$min
PC1maxshape_tmt<-pPC1_tmt$max
minshape_tmt<-tps3d(meshape_tmt, as.matrix(GPA_tmt$mshape), as.matrix(PC1minshape_tmt))
maxshape_tmt<-tps3d(meshape_tmt, as.matrix(GPA_tmt$mshape), as.matrix(PC1maxshape_tmt))


##########plot th shape deformation along the PC1
open3d()
############Minimum
plot3d(minshape_tmt, col="cornsilk", axes=FALSE, box=FALSE, specular=1,, emission=1,add=T)
############Maximum
plot3d(maxshape_tmt, col="cornsilk", axes=FALSE, box=FALSE, specular=1,alpha=0.25, emission=1 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC1maxshape_tmt,PC1minshape_tmt, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici


PC2minshape_tmt<-pPC2_tmt$min
PC2maxshape_tmt<-pPC2_tmt$max
minshape2_tmt<-tps3d(meshape_tmt, as.matrix(GPA_tmt$mshape), as.matrix(PC2minshape_tmt))
maxshape2_tmt<-tps3d(meshape_tmt, as.matrix(GPA_tmt$mshape), as.matrix(PC2maxshape_tmt))
##PC2

open3d()

plot3d(minshape2_tmt, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.2, add=T)
plot3d(maxshape2_tmt, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,,add=T)
deformGrid3d(PC2maxshape_tmt,PC2minshape_tmt, col1="brown1", col2="dodgerblue4", type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici


#############D1P1

meanshape<-GPA_D1P1$mshape
meshape_D1P1<-tps3d(template.D1P1,D1P1$dataslide[,,1], as.matrix(GPA_D1P1$mshape))
plot3d(meshape_D1P1, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_D1P1<-GPA_D1P1$PCscores[,1]
PC2_D1P1<-GPA_D1P1$PCscores[,2]
pPC1_D1P1<-shape.predictor(GPA_D1P1$rotated, x=PC1_D1P1, intercept=FALSE, min=min(PC1_D1P1), max=max(PC1_D1P1))
pPC2_D1P1<-shape.predictor(GPA_D1P1$rotated, x=PC2_D1P1, intercept=FALSE, min=min(PC2_D1P1), max=max(PC2_D1P1))
##PC1
PC1minshape_D1P1<-pPC1_D1P1$min
PC1maxshape_D1P1<-pPC1_D1P1$max
minshape_D1P1<-tps3d(meshape_D1P1, as.matrix(GPA_D1P1$mshape), as.matrix(PC1minshape_D1P1))
maxshape_D1P1<-tps3d(meshape_D1P1, as.matrix(GPA_D1P1$mshape), as.matrix(PC1maxshape_D1P1))

open3d()
plot3d(minshape_D1P1, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.25,add=T)
plot3d(maxshape_D1P1, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC1maxshape_D1P1,PC1minshape_D1P1, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici


PC2minshape_D1P1<-pPC2_D1P1$min
PC2maxshape_D1P1<-pPC2_D1P1$max
minshape2_D1P1<-tps3d(meshape_D1P1, as.matrix(GPA_D1P1$mshape), as.matrix(PC2minshape_D1P1))
maxshape2_D1P1<-tps3d(meshape_D1P1, as.matrix(GPA_D1P1$mshape), as.matrix(PC2maxshape_D1P1))
##PC2

open3d()

plot3d(minshape2_D1P1, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.2, add=T)
plot3d(maxshape2_D1P1, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,add=T)
deformGrid3d(PC2maxshape_D1P1,PC2minshape_D1P1, col1="brown1", col2="dodgerblue4", type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici
filename <- paste("pi2", formatC(i, digits = 1, flag = "0"), ".png", sep = "")


#############D1P2

meanshape<-GPA_D1P2$mshape
meshape_D1P2<-tps3d(template.D1P2,D1P2$dataslide[,,1], as.matrix(GPA_D1P2$mshape))
plot3d(meshape_D1P2, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_D1P2<-GPA_D1P2$PCscores[,1]
PC2_D1P2<-GPA_D1P2$PCscores[,2]
pPC1_D1P2<-shape.predictor(GPA_D1P2$rotated, x=PC1_D1P2, intercept=FALSE, min=min(PC1_D1P2), max=max(PC1_D1P2))
pPC2_D1P2<-shape.predictor(GPA_D1P2$rotated, x=PC2_D1P2, intercept=FALSE, min=min(PC2_D1P2), max=max(PC2_D1P2))
##PC1
PC1minshape_D1P2<-pPC1_D1P2$min
PC1maxshape_D1P2<-pPC1_D1P2$max
minshape_D1P2<-tps3d(meshape_D1P2, as.matrix(GPA_D1P2$mshape), as.matrix(PC1minshape_D1P2))
maxshape_D1P2<-tps3d(meshape_D1P2, as.matrix(GPA_D1P2$mshape), as.matrix(PC1maxshape_D1P2))

open3d()
bg3d("black")
plot3d(minshape_D1P2, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,add=T)
plot3d(maxshape_D1P2, col="cornsilk", axes=FALSE, box=FALSE, specular=1,,alpha=0.25, emission=1 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC1maxshape_D1P2,PC1minshape_D1P2, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici


PC2minshape_D1P2<-pPC2_D1P2$min
PC2maxshape_D1P2<-pPC2_D1P2$max
minshape2_D1P2<-tps3d(meshape_D1P2, as.matrix(GPA_D1P2$mshape), as.matrix(PC2minshape_D1P2))
maxshape2_D1P2<-tps3d(meshape_D1P2, as.matrix(GPA_D1P2$mshape), as.matrix(PC2maxshape_D1P2))
##PC2

open3d()
bg3d("black")
plot3d(minshape2_D1P2, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1, add=T)
plot3d(maxshape2_D1P2, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.2,add=T)
deformGrid3d(PC2maxshape_D1P2,PC2minshape_D1P2, col1="brown1", col2="dodgerblue4", type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici
filename <- paste("pi2", formatC(i, digits = 1, flag = "0"), ".png", sep = "")


#############D2P1

meanshape<-GPA_D2P1$mshape
meshape_D2P1<-tps3d(template.D2P1,as.matrix(D2P1$dataslide[,,1]), as.matrix(GPA_D2P1$mshape))
plot3d(meshape_D2P1, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_D2P1<-GPA_D2P1$PCscores[,1]
PC2_D2P1<-GPA_D2P1$PCscores[,2]
pPC1_D2P1<-shape.predictor(GPA_D2P1$rotated, x=PC1_D2P1, intercept=FALSE, min=min(PC1_D2P1), max=max(PC1_D2P1))
pPC2_D2P1<-shape.predictor(GPA_D2P1$rotated, x=PC2_D2P1, intercept=FALSE, min=min(PC2_D2P1), max=max(PC2_D2P1))
##PC1
PC1minshape_D2P1<-pPC1_D2P1$min
PC1maxshape_D2P1<-pPC1_D2P1$max
minshape_D2P1<-tps3d(meshape_D2P1, as.matrix(GPA_D2P1$mshape), as.matrix(PC1minshape_D2P1))
maxshape_D2P1<-tps3d(meshape_D2P1, as.matrix(GPA_D2P1$mshape), as.matrix(PC1maxshape_D2P1))

open3d()
plot3d(minshape_D2P1, col="cornsilk", axes=FALSE, box=FALSE,mag=2,alpha=1,add=T)
plot3d(maxshape_D2P1, col="cornsilk", axes=FALSE, box=FALSE,mag=2 ,alpha=0.25 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC2maxshape_D2P1,PC2minshape_D2P1, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.001, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici

plotRefToTarget(GPA_D2P1$rotated[,,38],GPA_D2P1$rotated[,,40],method = "vector",gridPars=gridPar(pt.bg = "green", pt.size = 1),mag=1.5)


PC2minshape_D2P1<-pPC2_D2P1$min
PC2maxshape_D2P1<-pPC2_D2P1$max
minshape2_D2P1<-tps3d(meshape_D2P1, as.matrix(GPA_D2P1$mshape), as.matrix(PC2minshape_D2P1))
maxshape2_D2P1<-tps3d(meshape_D2P1, as.matrix(GPA_D2P1$mshape), as.matrix(PC2maxshape_D2P1))
##PC2

open3d()

plot3d(minshape2_D2P1, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.75, add=T)
plot3d(maxshape2_D2P1, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.75,add=T)

####################D2P2

meanshape<-GPA_D2P2$mshape
meshape_D2P2<-tps3d(template.D2P2,as.matrix(D2P2$dataslide[,,1]), as.matrix(GPA_D2P2$mshape))
plot3d(meshape_D2P2, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_D2P2<-GPA_D2P2$PCscores[,1]
PC2_D2P2<-GPA_D2P2$PCscores[,2]
pPC1_D2P2<-shape.predictor(GPA_D2P2$rotated, x=PC1_D2P2, intercept=FALSE, min=min(PC1_D2P2), max=max(PC1_D2P2))
pPC2_D2P2<-shape.predictor(GPA_D2P2$rotated, x=PC2_D2P2, intercept=FALSE, min=min(PC2_D2P2), max=max(PC2_D2P2))
##PC1
PC1minshape_D2P2<-pPC1_D2P2$min
PC1maxshape_D2P2<-pPC1_D2P2$max
minshape_D2P2<-tps3d(meshape_D2P2, as.matrix(GPA_D2P2$mshape), as.matrix(PC1minshape_D2P2))
maxshape_D2P2<-tps3d(meshape_D2P2, as.matrix(GPA_D2P2$mshape), as.matrix(PC1maxshape_D2P2))

open3d()
plot3d(minshape_D2P2, col="cornsilk", axes=FALSE, box=FALSE,mag=2,alpha=0.25,add=T)
plot3d(maxshape_D2P2, col="cornsilk", axes=FALSE, box=FALSE,mag=2 ,alpha=1 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC2maxshape_D2P2,PC2minshape_D2P2, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.001, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici

plotRefToTarget(GPA_D2P2$rotated[,,38],GPA_D2P2$rotated[,,40],method = "vector",gridPars=gridPar(pt.bg = "green", pt.size = 1),mag=1.5)


PC2minshape_D2P2<-pPC2_D2P2$min
PC2maxshape_D2P2<-pPC2_D2P2$max
minshape2_D2P2<-tps3d(meshape_D2P2, as.matrix(GPA_D2P2$mshape), as.matrix(PC2minshape_D2P2))
maxshape2_D2P2<-tps3d(meshape_D2P2, as.matrix(GPA_D2P2$mshape), as.matrix(PC2maxshape_D2P2))
##PC2

open3d()

plot3d(minshape2_D2P2, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.75, add=T)
plot3d(maxshape2_D2P2, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.75,add=T)
deformGrid3d(PC2maxshape_D2P2,PC2minshape_D2P2, col1="brown1", col2="dodgerblue4", type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici
filename <- paste("pi2", formatC(i, digits = 1, flag = "0"), ".png", sep = "")


#############D2P3

meanshape<-GPA_D2P3$mshape
meshape_D2P3<-tps3d(template.D2P3,D2P3$dataslide[,,1], as.matrix(GPA_D2P3$mshape))
plot3d(meshape_D2P3, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_D2P3<-GPA_D2P3$PCscores[,1]
PC2_D2P3<-GPA_D2P3$PCscores[,2]
pPC1_D2P3<-shape.predictor(GPA_D2P3$rotated, x=PC1_D2P3, intercept=FALSE, min=min(PC1_D2P3), max=max(PC1_D2P3))
pPC2_D2P3<-shape.predictor(GPA_D2P3$rotated, x=PC2_D2P3, intercept=FALSE, min=min(PC2_D2P3), max=max(PC2_D2P3))
##PC1
PC1minshape_D2P3<-pPC1_D2P3$min
PC1maxshape_D2P3<-pPC1_D2P3$max
minshape_D2P3<-tps3d(meshape_D2P3, as.matrix(GPA_D2P3$mshape), as.matrix(PC1minshape_D2P3))
maxshape_D2P3<-tps3d(meshape_D2P3, as.matrix(GPA_D2P3$mshape), as.matrix(PC1maxshape_D2P3))

open3d()
plot3d(minshape_D2P3, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,add=T)
plot3d(maxshape_D2P3, col="cornsilk", axes=FALSE, box=FALSE, specular=1,,alpha=0.25, emission=1 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC1maxshape_D2P3,PC1minshape_D2P3, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici


PC2minshape_D2P3<-pPC2_D2P3$min
PC2maxshape_D2P3<-pPC2_D2P3$max
minshape2_D2P3<-tps3d(meshape_D2P3, as.matrix(GPA_D2P3$mshape), as.matrix(PC2minshape_D2P3))
maxshape2_D2P3<-tps3d(meshape_D2P3, as.matrix(GPA_D2P3$mshape), as.matrix(PC2maxshape_D2P3))
##PC2

open3d()

plot3d(minshape2_D2P3, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.2, add=T)
plot3d(maxshape2_D2P3, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,add=T)
deformGrid3d(PC2maxshape_D2P3,PC2minshape_D2P3, col1="brown1", col2="dodgerblue4", type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici


#############D3P1

meanshape<-GPA_D3P1$mshape
meshape_D3P1<-tps3d(template.D3P1,D3P1$dataslide[,,1], as.matrix(GPA_D3P1$mshape))
plot3d(meshape_D3P1, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_D3P1<-GPA_D3P1$PCscores[,1]
PC2_D3P1<-GPA_D3P1$PCscores[,2]
pPC1_D3P1<-shape.predictor(GPA_D3P1$rotated, x=PC1_D3P1, intercept=FALSE, min=min(PC1_D3P1), max=max(PC1_D3P1))
pPC2_D3P1<-shape.predictor(GPA_D3P1$rotated, x=PC2_D3P1, intercept=FALSE, min=min(PC2_D3P1), max=max(PC2_D3P1))
##PC1
PC1minshape_D3P1<-pPC1_D3P1$min
PC1maxshape_D3P1<-pPC1_D3P1$max
minshape_D3P1<-tps3d(meshape_D3P1, as.matrix(GPA_D3P1$mshape), as.matrix(PC1minshape_D3P1))
maxshape_D3P1<-tps3d(meshape_D3P1, as.matrix(GPA_D3P1$mshape), as.matrix(PC1maxshape_D3P1))

open3d()
plot3d(minshape_D3P1, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,add=T)
plot3d(maxshape_D3P1, col="cornsilk", axes=FALSE, box=FALSE, specular=1,,alpha=1 emission=1 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC1maxshape_D3P1,PC1minshape_D3P1, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici


PC2minshape_D3P1<-pPC2_D3P1$min
PC2maxshape_D3P1<-pPC2_D3P1$max
minshape2_D3P1<-tps3d(meshape_D3P1, as.matrix(GPA_D3P1$mshape), as.matrix(PC2minshape_D3P1))
maxshape2_D3P1<-tps3d(meshape_D3P1, as.matrix(GPA_D3P1$mshape), as.matrix(PC2maxshape_D3P1))
##PC2

open3d()
plot3d(minshape2_D3P1, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,, add=T)
plot3d(maxshape2_D3P1, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=1,add=T)
deformGrid3d(PC2maxshape_D3P1,PC2minshape_D3P1, col1="brown1", col2="dodgerblue4", type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici


##################D3P2
meanshape<-GPA_D3P2$mshape
meshape_D3P2<-tps3d(template.D3P2,as.matrix(D3P2$dataslide[,,1]), as.matrix(GPA_D3P2$mshape))
plot3d(meshape_D3P2, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_D3P2<-GPA_D3P2$PCscores[,1]
PC2_D3P2<-GPA_D3P2$PCscores[,2]
pPC1_D3P2<-shape.predictor(GPA_D3P2$rotated, x=PC1_D3P2, intercept=FALSE, min=min(PC1_D3P2), max=max(PC1_D3P2))
pPC2_D3P2<-shape.predictor(GPA_D3P2$rotated, x=PC2_D3P2, intercept=FALSE, min=min(PC2_D3P2), max=max(PC2_D3P2))
##PC1
PC1minshape_D3P2<-pPC1_D3P2$min
PC1maxshape_D3P2<-pPC1_D3P2$max
minshape_D3P2<-tps3d(meshape_D3P2, as.matrix(GPA_D3P2$mshape), as.matrix(PC1minshape_D3P2))
maxshape_D3P2<-tps3d(meshape_D3P2, as.matrix(GPA_D3P2$mshape), as.matrix(PC1maxshape_D3P2))

open3d()
plot3d(minshape_D3P2, col="cornsilk", axes=FALSE, box=FALSE,mag=2,alpha=0.25,add=T)
plot3d(maxshape_D3P2, col="cornsilk", axes=FALSE, box=FALSE,mag=2 ,alpha=1 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC2maxshape_D3P2,PC2minshape_D3P2, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.001, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici

plotRefToTarget(GPA_D3P2$rotated[,,38],GPA_D3P2$rotated[,,40],method = "vector",gridPars=gridPar(pt.bg = "green", pt.size = 1),mag=1.5)


PC2minshape_D3P2<-pPC2_D3P2$min
PC2maxshape_D3P2<-pPC2_D3P2$max
minshape2_D3P2<-tps3d(meshape_D3P2, as.matrix(GPA_D3P2$mshape), as.matrix(PC2minshape_D3P2))
maxshape2_D3P2<-tps3d(meshape_D3P2, as.matrix(GPA_D3P2$mshape), as.matrix(PC2maxshape_D3P2))
##PC2

open3d()

plot3d(minshape2_D3P2, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.75, add=T)
plot3d(maxshape2_D3P2, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.75,add=T)
deformGrid3d(PC2maxshape_D3P2,PC2minshape_D3P2, col1="brown1", col2="dodgerblue4", type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici


#################D3P3


meanshape<-GPA_D3P3$mshape
meshape_D3P3<-tps3d(template.D3P3,as.matrix(D3P3$dataslide[,,1]), as.matrix(GPA_D3P3$mshape))
plot3d(meshape_D3P3, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_D3P3<-GPA_D3P3$PCscores[,1]
PC2_D3P3<-GPA_D3P3$PCscores[,2]
pPC1_D3P3<-shape.predictor(GPA_D3P3$rotated, x=PC1_D3P3, intercept=FALSE, min=min(PC1_D3P3), max=max(PC1_D3P3))
pPC2_D3P3<-shape.predictor(GPA_D3P3$rotated, x=PC2_D3P3, intercept=FALSE, min=min(PC2_D3P3), max=max(PC2_D3P3))
##PC1
PC1minshape_D3P3<-pPC1_D3P3$min
PC1maxshape_D3P3<-pPC1_D3P3$max
minshape_D3P3<-tps3d(meshape_D3P3, as.matrix(GPA_D3P3$mshape), as.matrix(PC1minshape_D3P3))
maxshape_D3P3<-tps3d(meshape_D3P3, as.matrix(GPA_D3P3$mshape), as.matrix(PC1maxshape_D3P3))

open3d()
plot3d(minshape_D3P3, col="cornsilk", axes=FALSE, box=FALSE,mag=2,alpha=0.25,add=T)
plot3d(maxshape_D3P3, col="cornsilk", axes=FALSE, box=FALSE,mag=2 ,alpha=1 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC2maxshape_D3P3,PC2minshape_D3P3, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.001, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici

plotRefToTarget(GPA_D3P3$rotated[,,38],GPA_D3P3$rotated[,,40],method = "vector",gridPars=gridPar(pt.bg = "green", pt.size = 1),mag=1.5)


PC2minshape_D3P3<-pPC2_D3P3$min
PC2maxshape_D3P3<-pPC2_D3P3$max
minshape2_D3P3<-tps3d(meshape_D3P3, as.matrix(GPA_D3P3$mshape), as.matrix(PC2minshape_D3P3))
maxshape2_D3P3<-tps3d(meshape_D3P3, as.matrix(GPA_D3P3$mshape), as.matrix(PC2maxshape_D3P3))
##PC2

open3d()

plot3d(minshape2_D3P3, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.75, add=T)
plot3d(maxshape2_D3P3, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.75,add=T)
deformGrid3d(PC2maxshape_D3P3,PC2minshape_D3P3, col1="brown1", col2="dodgerblue4", type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici




#############D3P4

meanshape<-GPA_D3P4$mshape
meshape_D3P4<-tps3d(template.D3P4,D3P4$dataslide[,,1], as.matrix(GPA_D3P4$mshape))
plot3d(meshape_D3P4, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_D3P4<-GPA_D3P4$PCscores[,1]
PC2_D3P4<-GPA_D3P4$PCscores[,2]
pPC1_D3P4<-shape.predictor(GPA_D3P4$rotated, x=PC1_D3P4, intercept=FALSE, min=min(PC1_D3P4), max=max(PC1_D3P4))
pPC2_D3P4<-shape.predictor(GPA_D3P4$rotated, x=PC2_D3P4, intercept=FALSE, min=min(PC2_D3P4), max=max(PC2_D3P4))
##PC1
PC1minshape_D3P4<-pPC1_D3P4$min
PC1maxshape_D3P4<-pPC1_D3P4$max
minshape_D3P4<-tps3d(meshape_D3P4, as.matrix(GPA_D3P4$mshape), as.matrix(PC1minshape_D3P4))
maxshape_D3P4<-tps3d(meshape_D3P4, as.matrix(GPA_D3P4$mshape), as.matrix(PC1maxshape_D3P4))

open3d()
plot3d(minshape_D3P4, col="cornsilk", axes=FALSE, box=FALSE, specular=1, alpha=0.25,emission=1,add=T)
plot3d(maxshape_D3P4, col="cornsilk", axes=FALSE, box=FALSE, specular=1,alpha=0.75, emission=1 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC1maxshape_D3P4,PC1minshape_D3P4, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici



PC2minshape_D3P4<-pPC2_D3P4$min
PC2maxshape_D3P4<-pPC2_D3P4$max
minshape2_D3P4<-tps3d(meshape_D3P4, as.matrix(GPA_D3P4$mshape), as.matrix(PC2minshape_D3P4))
maxshape2_D3P4<-tps3d(meshape_D3P4, as.matrix(GPA_D3P4$mshape), as.matrix(PC2maxshape_D3P4))
##PC2

open3d()

plot3d(minshape2_D3P4, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.2, add=T)
plot3d(maxshape2_D3P4, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,add=T)
deformGrid3d(PC2maxshape_D3P4,PC2minshape_D3P4, col1="brown1", col2="dodgerblue4", type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici


#############D4P1

meanshape<-GPA_D4P1$mshape
meshape_D4P1<-tps3d(template.D4P1,D4P1$dataslide[,,1], as.matrix(GPA_D4P1$mshape))
plot3d(meshape_D4P1, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_D4P1<-GPA_D4P1$PCscores[,1]
PC2_D4P1<-GPA_D4P1$PCscores[,2]
pPC1_D4P1<-shape.predictor(GPA_D4P1$rotated, x=PC1_D4P1, intercept=FALSE, min=min(PC1_D4P1), max=max(PC1_D4P1))
pPC2_D4P1<-shape.predictor(GPA_D4P1$rotated, x=PC2_D4P1, intercept=FALSE, min=min(PC2_D4P1), max=max(PC2_D4P1))
##PC1
PC1minshape_D4P1<-pPC1_D4P1$min
PC1maxshape_D4P1<-pPC1_D4P1$max
minshape_D4P1<-tps3d(meshape_D4P1, as.matrix(GPA_D4P1$mshape), as.matrix(PC1minshape_D4P1))
maxshape_D4P1<-tps3d(meshape_D4P1, as.matrix(GPA_D4P1$mshape), as.matrix(PC1maxshape_D4P1))

open3d()
plot3d(minshape_D4P1, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,add=T)
plot3d(maxshape_D4P1, col="cornsilk", axes=FALSE, box=FALSE, specular=1,,alpha=0.25, emission=1 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC1maxshape_D4P1,PC1minshape_D4P1, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici


PC2minshape_D4P1<-pPC2_D4P1$min
PC2maxshape_D4P1<-pPC2_D4P1$max
minshape2_D4P1<-tps3d(meshape_D4P1, as.matrix(GPA_D4P1$mshape), as.matrix(PC2minshape_D4P1))
maxshape2_D4P1<-tps3d(meshape_D4P1, as.matrix(GPA_D4P1$mshape), as.matrix(PC2maxshape_D4P1))
##PC2

open3d()
plot3d(minshape2_D4P1, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.2, add=T)
plot3d(maxshape2_D4P1, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,add=T)
deformGrid3d(PC2maxshape_D4P1,PC2minshape_D4P1, col1="brown1", col2="dodgerblue4", type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici

############D4P2

meanshape<-GPA_D4P2$mshape
meshape_D4P2<-tps3d(template.D4P2,as.matrix(D4P2$dataslide[,,1]), as.matrix(GPA_D4P2$mshape))
plot3d(meshape_D4P2, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_D4P2<-GPA_D4P2$PCscores[,1]
PC2_D4P2<-GPA_D4P2$PCscores[,2]
pPC1_D4P2<-shape.predictor(GPA_D4P2$rotated, x=PC1_D4P2, intercept=FALSE, min=min(PC1_D4P2), max=max(PC1_D4P2))
pPC2_D4P2<-shape.predictor(GPA_D4P2$rotated, x=PC2_D4P2, intercept=FALSE, min=min(PC2_D4P2), max=max(PC2_D4P2))
##PC1
PC1minshape_D4P2<-pPC1_D4P2$min
PC1maxshape_D4P2<-pPC1_D4P2$max
minshape_D4P2<-tps3d(meshape_D4P2, as.matrix(GPA_D4P2$mshape), as.matrix(PC1minshape_D4P2))
maxshape_D4P2<-tps3d(meshape_D4P2, as.matrix(GPA_D4P2$mshape), as.matrix(PC1maxshape_D4P2))

open3d()
plot3d(minshape_D4P2, col="cornsilk", axes=FALSE, box=FALSE,mag=2,alpha=0.25,add=T)
plot3d(maxshape_D4P2, col="cornsilk", axes=FALSE, box=FALSE,mag=2 ,alpha=1 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC2maxshape_D4P2,PC2minshape_D4P2, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.001, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici

plotRefToTarget(GPA_D4P2$rotated[,,38],GPA_D4P2$rotated[,,40],method = "vector",gridPars=gridPar(pt.bg = "green", pt.size = 1),mag=1.5)


PC2minshape_D4P2<-pPC2_D4P2$min
PC2maxshape_D4P2<-pPC2_D4P2$max
minshape2_D4P2<-tps3d(meshape_D4P2, as.matrix(GPA_D4P2$mshape), as.matrix(PC2minshape_D4P2))
maxshape2_D4P2<-tps3d(meshape_D4P2, as.matrix(GPA_D4P2$mshape), as.matrix(PC2maxshape_D4P2))
##PC2

open3d()

plot3d(minshape2_D4P2, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.75, add=T)
plot3d(maxshape2_D4P2, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.75,add=T)
deformGrid3d(PC2maxshape_D4P2,PC2minshape_D4P2, col1="brown1", col2="dodgerblue4", type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici
filename <- paste("pi2", formatC(i, digits = 1, flag = "0"), ".png", sep = "")

#############D4P3


meanshape<-GPA_D4P3$mshape
meshape_D4P3<-tps3d(template.D4P3,as.matrix(D4P3$dataslide[,,1]), as.matrix(GPA_D4P3$mshape))
plot3d(meshape_D4P3, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_D4P3<-GPA_D4P3$PCscores[,1]
PC2_D4P3<-GPA_D4P3$PCscores[,2]
pPC1_D4P3<-shape.predictor(GPA_D4P3$rotated, x=PC1_D4P3, intercept=FALSE, min=min(PC1_D4P3), max=max(PC1_D4P3))
pPC2_D4P3<-shape.predictor(GPA_D4P3$rotated, x=PC2_D4P3, intercept=FALSE, min=min(PC2_D4P3), max=max(PC2_D4P3))
##PC1
PC1minshape_D4P3<-pPC1_D4P3$min
PC1maxshape_D4P3<-pPC1_D4P3$max
minshape_D4P3<-tps3d(meshape_D4P3, as.matrix(GPA_D4P3$mshape), as.matrix(PC1minshape_D4P3))
maxshape_D4P3<-tps3d(meshape_D4P3, as.matrix(GPA_D4P3$mshape), as.matrix(PC1maxshape_D4P3))

open3d()
plot3d(minshape_D4P3, col="cornsilk", axes=FALSE, box=FALSE,mag=2,alpha=0.25,add=T)
plot3d(maxshape_D4P3, col="cornsilk", axes=FALSE, box=FALSE,mag=2 ,alpha=1 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC2maxshape_D4P3,PC2minshape_D4P3, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.001, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici

plotRefToTarget(GPA_D4P3$rotated[,,38],GPA_D4P3$rotated[,,40],method = "vector",gridPars=gridPar(pt.bg = "green", pt.size = 1),mag=1.5)


PC2minshape_D4P3<-pPC2_D4P3$min
PC2maxshape_D4P3<-pPC2_D4P3$max
minshape2_D4P3<-tps3d(meshape_D4P3, as.matrix(GPA_D4P3$mshape), as.matrix(PC2minshape_D4P3))
maxshape2_D4P3<-tps3d(meshape_D4P3, as.matrix(GPA_D4P3$mshape), as.matrix(PC2maxshape_D4P3))
##PC2

open3d()

plot3d(minshape2_D4P3, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.75, add=T)
plot3d(maxshape2_D4P3, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.75,add=T)
deformGrid3d(PC2maxshape_D4P3,PC2minshape_D4P3, col1="brown1", col2="dodgerblue4", type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici

###################D4P4

meanshape<-GPA_D4P4$mshape
meshape_D4P4<-tps3d(template.D4P4,as.matrix(D4P4$dataslide[,,1]), as.matrix(GPA_D4P4$mshape))
plot3d(meshape_D4P4, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_D4P4<-GPA_D4P4$PCscores[,1]
PC2_D4P4<-GPA_D4P4$PCscores[,2]
pPC1_D4P4<-shape.predictor(GPA_D4P4$rotated, x=PC1_D4P4, intercept=FALSE, min=min(PC1_D4P4), max=max(PC1_D4P4))
pPC2_D4P4<-shape.predictor(GPA_D4P4$rotated, x=PC2_D4P4, intercept=FALSE, min=min(PC2_D4P4), max=max(PC2_D4P4))
##PC1
PC1minshape_D4P4<-pPC1_D4P4$min
PC1maxshape_D4P4<-pPC1_D4P4$max
minshape_D4P4<-tps3d(meshape_D4P4, as.matrix(GPA_D4P4$mshape), as.matrix(PC1minshape_D4P4))
maxshape_D4P4<-tps3d(meshape_D4P4, as.matrix(GPA_D4P4$mshape), as.matrix(PC1maxshape_D4P4))

open3d()
plot3d(minshape_D4P4, col="cornsilk", axes=FALSE, box=FALSE,mag=2,alpha=0.25,add=T)
plot3d(maxshape_D4P4, col="cornsilk", axes=FALSE, box=FALSE,mag=2 ,alpha=1 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC2maxshape_D4P4,PC2minshape_D4P4, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.001, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici

plotRefToTarget(GPA_D4P4$rotated[,,38],GPA_D4P4$rotated[,,40],method = "vector",gridPars=gridPar(pt.bg = "green", pt.size = 1),mag=1.5)


PC2minshape_D4P4<-pPC2_D4P4$min
PC2maxshape_D4P4<-pPC2_D4P4$max
minshape2_D4P4<-tps3d(meshape_D4P4, as.matrix(GPA_D4P4$mshape), as.matrix(PC2minshape_D4P4))
maxshape2_D4P4<-tps3d(meshape_D4P4, as.matrix(GPA_D4P4$mshape), as.matrix(PC2maxshape_D4P4))
##PC2

open3d()

plot3d(minshape2_D4P4, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.75, add=T)
plot3d(maxshape2_D4P4, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.75,add=T)
deformGrid3d(PC2maxshape_D4P4,PC2minshape_D4P4, col1="brown1", col2="dodgerblue4", type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici

#############D4P5

meanshape<-GPA_D4P5$mshape
meshape_D4P5<-tps3d(template.D4P5,D4P5$dataslide[,,1], as.matrix(GPA_D4P5$mshape))
plot3d(meshape_D4P5, col="ivory", axes=FALSE, box=FALSE, alpha=0.8 ,add=T)

#MIN et MAX shapes
PC1_D4P5<-GPA_D4P5$PCscores[,1]
PC2_D4P5<-GPA_D4P5$PCscores[,2]
pPC1_D4P5<-shape.predictor(GPA_D4P5$rotated, x=PC1_D4P5, intercept=FALSE, min=min(PC1_D4P5), max=max(PC1_D4P5))
pPC2_D4P5<-shape.predictor(GPA_D4P5$rotated, x=PC2_D4P5, intercept=FALSE, min=min(PC2_D4P5), max=max(PC2_D4P5))
##PC1
PC1minshape_D4P5<-pPC1_D4P5$min
PC1maxshape_D4P5<-pPC1_D4P5$max
minshape_D4P5<-tps3d(meshape_D4P5, as.matrix(GPA_D4P5$mshape), as.matrix(PC1minshape_D4P5))
maxshape_D4P5<-tps3d(meshape_D4P5, as.matrix(GPA_D4P5$mshape), as.matrix(PC1maxshape_D4P5))

open3d()
plot3d(minshape_D4P5, col="cornsilk", axes=FALSE, box=FALSE, specular=1,alpha=0.25, emission=1,add=T)
plot3d(maxshape_D4P5, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1 ,add=T)
#ADD LANDMARKS
deformGrid3d(PC1maxshape_D4P5,PC1minshape_D4P5, col1="brown1", col2="dodgerblue4",lines = T, type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici


PC2minshape_D4P5<-pPC2_D4P5$min
PC2maxshape_D4P5<-pPC2_D4P5$max
minshape2_D4P5<-tps3d(meshape_D4P5, as.matrix(GPA_D4P5$mshape), as.matrix(PC2minshape_D4P5))
maxshape2_D4P5<-tps3d(meshape_D4P5, as.matrix(GPA_D4P5$mshape), as.matrix(PC2maxshape_D4P5))
##PC2

open3d()

plot3d(minshape2_D4P5, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,alpha=0.2, add=T)
plot3d(maxshape2_D4P5, col="cornsilk", axes=FALSE, box=FALSE, specular=1, emission=1,add=T)
deformGrid3d(PC2maxshape_D4P5,PC2minshape_D4P5, col1="brown1", col2="dodgerblue4", type="s", size=0.0005, align=T, createMesh=T, add=T) ### tu devras surement adapter le size ici

