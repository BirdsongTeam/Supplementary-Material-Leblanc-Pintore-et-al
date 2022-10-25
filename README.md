This Readme.txt file was generated on 10/25/2022 by Killian Leblanc and Pauline Provini
-------------------
GENERAL INFORMATION
-------------------

Title of Dataset: Supplementary Material associated with the publication entitled: Foot adaptation to climbing in oven birds and woodcreepers (Furnariida)

Authors list: Killian Leblanc 1,2,*, Romain Pintore 3,4,*, Ana Galvão 5, Ezekiel Heitz 1,2,3, Pauline Provini 1,2,3,#
 
1 Université Paris Cité, Inserm, System Engineering and Evolution Dynamics, F-75004 Paris, France
2 Learning Planet Institute, F-75004 Paris, France
3 Mécanismes adaptatifs et évolution (MECADEV)/UMR 7179, CNRS/Muséum National d’Histoire Naturelle, Paris, France
4 Structure and Motion Laboratory, Department of Comparative Biomedical Sciences, Royal Veterinary College, AL9 7TA, Hatfield, United-Kingdom
5 Laboratório de Ornitologia, Departamento de Zoologia, Universidade Federal do Rio de Janeiro
*: co first authors
 #: corresponding author

Date of data collection (approximate date): 2015-2022
Geographic location of data collection: Universidade Federal do Rio de Janeiro, Rio de Janeiro, Brazil

--------------------------
SHARING/ACCESS INFORMATION
-------------------------- 

Licenses/restrictions placed on the data, or limitations of reuse: public domain, without restriction on use

Recommended citation for the data: 

--------------------
DATA & FILE OVERVIEW
--------------------
File list: 

SM1: Ark identifiers for the specimens available on morphosource
 
SM2: Anatomical description of the 3D landmarks of the 14 phalanges and the tarsometatarsus used in the geometric morphometrics analysis 
 
SM3: Repetition test made with 10 replicates on 3 species morphologically similar
- SM3a:Tarsometatarsus
- SM3b:Phalanx I digit I
- SM3c:Phalanx II digit I
- SM3d:Phalanx I digit II
- SM3e:Phalanx II digit II
- SM3f:Phalanx III digit II
- SM3g:Phalanx I digit III
- SM3h:Phalanx II digit III
- SM3i:Phalanx III digit III
- SM3j:Phalanx IV digit III
- SM3k:Phalanx I digit IV
- SM3l:Phalanx II digit IV
- SM3m:Phalanx III digit IV
- SM3n:Phalanx IV digit IV
- SM3o:Phalanx V digit IV
 
SM4: Classifications of 55 Furnariida used in this study
- SM4a: Locomotion and Perching abilities of each studied species, with the sources used to define them.
- SM4b: Foraging time in different forest strates extracted from Elton traits 
- SM4c: Summary of the syndactyly observed in the studied specimens
 
SM5: SM5 Phylogenetic tree of the species used in the study, based on Moyle et al., 2009
 
SM6: PCA figure of PC1 and PC2 with phylomorphospace
- SM6a:Phalanx II digit II
- SM6b:Phalanx III digit II 
- SM6c:Phalanx II digit III
- SM6d:Phalanx III digit III
- SM6e:Phalanx IV digit III 
- SM6f:Phalanx II digit IV
- SM6g:Phalanx III digit IV
- SM6h:Phalanx IV digit IV 
- SM6i:Phalanx V digit IV 
 
SM7: Theoretical shape deformation
- SM7a:Theoretical shape deformation of phalanx II digit II, at both extremum of the two first axes of the PCA in Hylopezus macularius (A8409) 
- SM7b: Theoretical shape deformation of phalanx III digit II, at the negative side (top) and positive side (bottom) of PC1 (a), and PC2 (b) and (c) in Hylopezus macularius (A8409) 
- SM7c: Theoretical shape deformation of phalanx II digit III, at the negative side (left)and positive side (right) of PC1 (a), and PC2 (b) in Hylopezus macularius (A8409)
- SM7d: Theoretical shape deformation of phalanx III digit III, at the negative side (left) and positive side (right) of PC1 (a), and PC2 (b) in Hylopezus macularius (A8409)
- SM7e: Theoretical shape deformation of phalanx IV digit III, at the negative side (top) and positive side (bottom) of PC1 (a), and PC2 (b) and (c) in Hylopezus macularius (A8409)
- SM7f: Theoretical shape deformation of phalanx II digit IV, at the negative side (left) and positive side (right) of PC1 (a), and PC2 (b) in Hylopezus macularius (A8409)
- SM7g: Theoretical shape deformation of phalanx III digit IV, at the negative side (left) and positive side (right) of PC1 (a), and PC2 (b) in Hylopezus macularius (A8409)
- SM7h: Theoretical shape deformation of phalanx IV digit IV, at the negative side (left) and positive side (right) of PC1 (a), and PC2 (b) in Hylopezus macularius (A8409)
- SM7i: Theoretical shape deformation of phalanx V digit IV, at the negative side (top) and positive side (bottom) of PC1 (a), and PC2 (b) in Hylopezus macularius (A8409)
 
 
SM8: Summary of the statistical tests
- SM7a: Size effect 
- SM7b: Phylogenetic test 
- SM7c: Variance summary of tmt 
- SM7d: Variance summary of hallux 
- SM7e: Variance summary of digit II 
- SM7f: Variance summary of digit III 
- SM7g: Variance summary f digit IV 
 
SM9: R Code used in the analyses
- SM9a Script to load the .pts file and launch the Generalised Procrustes Alignment (GPA) and plot PCA associate
- SM9b Script to visualise the shape deformation along each PC for each PCA
- SM9c Script to launch statistical tests (phylogenetic and size effect).
- SM9d All file requiert to launch the script, including the result of the sliding, the different trees adjusted as for some bones as some specimens had either damaged or absent bones. Finally we add the .ply necessary to visualise shape deformation.


Data derived from other sources:

In SM1, the landmarks defifition is based on previous definitions from Hedrick et al., 2019 and Bjarnason & Benson 2021
In SM4, locomotion categories and perching performance are derived from data in ebird from Wilman et al, 2014.
In SM5, the cladogram is based on the phylogeny of Moyle et al., 2009
In SM9 we included extra variables: climbing index from Claramunt et al., 2012, the foraging behaviours and bodymass from Wilman et al, 2014 and Tobias & Pigot, if the operator wants to vizualise our datasets.


--------------------------
METHODOLOGICAL INFORMATION
--------------------------

Description of methods used for collection/generation of data: 

See Material and Method section in the manuscript.

SM1: We constructed a landmark protocol based on previous works on ungual bones (Hedrick et al. 2019) and tarsometarsus (Bjarnason & Benson 2021). We reused some definitions and create new one in order to capture correctly the shape of each structure.

SM4: We classify species based on the most used locomotor mode from the observation from the liteterature together with photos and videos from ebird. We create a new classification, the perching performance to evaluate the ability of birds to perform acrobatic and climbing abilities. The perching performance is a attemp to refine the hopping categories as hopping hopping and the ground and in branches are not exposed to the same physical stresss. Indeed the birds hopping in branches put their center of mass below the feet more often. 

SM5: In order to construct the phylogenetic tree, we retained the species only present in our study based on the result from Moyle et al., 2009.

