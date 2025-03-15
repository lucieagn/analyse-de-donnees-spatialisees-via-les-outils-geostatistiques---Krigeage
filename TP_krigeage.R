################################################################################
########################## Atelier de B. Iooss #################################
######################## MIM 1 Cote d'Azur 2024 ################################
################################################################################

# Etude de précipitations en Suisse (Swiss rainfall)
# Ref: Christensen, O.F., Diggle, P.J. and Ribeiro Jr, P.J. (2001) Analysing positive-valued spatial data: the transformed Gaussian model. In Monestiez, P., Allard, D. and Froidevaux (eds), GeoENV III - Geostatistics for environmental applications. Quantitative Geology and Geostatistics, Kluwer Series, 11, 287–298

# Ce TP est inspiré de : Manuel d’analyse spatiale (rapport INSEE), chapitre 5 de Jean-Michel Floch

#partie cours : slide 14 
#but du tp faire une cartographie sur les previions pluie
#(
###pour mac utiliser geoR
install.packages('geoR')
install.packages('tcltk')
library(tcltk)
library(geoR)

graphics.off()
rm(list=ls())

# 3 bases dans le package geoR :
# — sic.100 : échantillon de 100 observations qui pourront servir à effectuer les interpolations ;   -Base d'apprentissage (100 stations qui ont enregistrees de la pluie)
# — sic.367 : observations non incluses dans l’échantillon qui permettront de comparer les estimations et les observations ; -Base d'apprentissage (367 stations ....)
# — sic.all : ensemble.

library(geoR)

help(SIC)
summary(sic.all)


#sic.100$coords
#sic.100*data
################################################################################
# 1 Visualisation des donnees

x11() ; points(sic.100, borders=sic.borders,col="green")

x11()
points(sic.100, borders=sic.borders,col="green")
points(sic.367, borders=sic.borders,col="red",add=TRUE)

x11() ; plot.geodata(sic.100,bor=sic.borders) # quartiles in the first plot are divided following blue, green, yellow and red

#y at-il une structuration spatiale ?
################################################################################
# 2 Analyse variographique
#variogramme fct qui depend de h avec la variance . quand h est nul ca vaut 0 et quand h petit le variogramme est faible
#lorsque h augmente le variogramme augmente jusqua une certaine limite vers l'infini plus de corrélation ( la corrélation s'annule)
#donc en l'inifin le variogramme tend vers la variance
#slide 22
#comment est calcule la variogramme : tous les couples de pts dobservations qui sont separes de h
#sur le premier x11() 
#pour chaque couple on peut établior une nuée variographique
help(variog)

#cours :slide 21

# variogram cloud
vario.c <- variog(sic.100, op="cloud") #nuée variographique quand la distance entre deux observation augmente, plus la valeur du variogramme est forte
#binned variogram and stores the cloud.  #"cloud" pour nuée
vario.bc <- variog(sic.100, bin.cloud=TRUE) #boxplot par distance du variogramme juste au dessus.
#ligne noire correspond a la mediane et en fct de la distance cela tend a augmenter
# smoothed variogram
vario.s <- variog(sic.100, op="sm", band=10)#lissage de la moyenne en fct de la distance et quand la distance de sepration augmente le vario augm.
#ca augmente jusqua un certain seuil et puis cela devient constant les valeurs ne sont alors plus siognificative
# binned variogram. #option smoothed = lissage a noyaux= estimation non parametrique de la densite en utilisant lissage de la gaussienne

#tres important
vario.b <- variog(sic.100) # variogramme experimentale
#discretise mes distances avec un certain pas et pour chaque pas 
#classe de distance pour chque classe ont fait la moyenne 
#bin=classe
#petits points =Z(X+h)-Z(h)

#
x11() ; par(mfrow=c(2,2))
plot(vario.c, main="variogram cloud")
plot(vario.bc, bin.cloud=TRUE, main="clouds for binned variogram")
plot(vario.s, main="smoothed variogram", ylim=c(0,25000), xlim=c(0,300))
plot(vario.b, main="binned variogram", ylim=c(0,25000), xlim=c(0,300)) 

# variogrammes directionnels
help(variog4)
vario4 <- variog4(sic.100)
x11() ; plot(vario4,same=FALSE)
#isotropie : la structuration pstaile est la meme = la taille des boule rouge et des boule bleu est similaire verticalemùent et horizontalement
#anistortrope= variabilite spatiale n'est pas la meme horizontalement que verticalement
#plui plus important en longitude uquen largitude

#direction horizontale : evalue la distance horizontalement pour chaque point (a etudier sur le tout premier graphe)
#vario0 :horizontale
#vario90 : vertical
#il faut surtout regarder le debut des variogrammes et la portée 
#ce qui est important cest la vitessse a laquelle on monte verzs un plateau
#vario0 plateau a partir de 80
#vario45 : portée et plus grande car plateau est plus loin
#su vario90 :  plateau a 50

# le phenomene est plus corrélée dans la direction 45 que 135 
#l'angle de portée la plus forte est celle de vario45
#vario135 : creu et du une bande de valurs forte et ensuite plus de correlation puis de nouveau une bande de valeurs tres fortes

################################################################################
# 3 Ajustement d'un modèle de variogramme
#pourquoi besoin u modele : ajouter une fonction prarametrique = modele 
#fonction extremement simple (ex: gaussienne slide 27, loi exponentielle. avec a la portée, plus a est faible pente ardue plus a elevee plus pente faible)
#a = valeur a laqueulle plus bcp de corelation = longueur de correlation
#sigma 2 = valeur du plateau a laquelle je veux que ma courbe arrive

# distance >= a : plus aucune correlation a cette distance
#distance <a : correlation

#on utilise le modele gaussien et le modele spherique
vario.ex<- variog(sic.100,option="bin")
x11() ; plot(vario.ex) # nous donne une idee des parametres de variance et de portee

#exemple de fit
help(variofit)
#methode de descente de gradient pour optimiser la somme residus aux carres entre la courbe que je fit et les points qui sont la
#avec une certaine ponderation selon les points (on choisit de donner du poids aux données qui permettront de lisser les donneesz)
vario.sphe <- variofit(vario.ex,cov.model= "spher", ini.cov.pars=c(15000,100))
#on fit le variogramme experimental on donne le modele = spherique et c(15000,100)=valeurs initiales (on choisit lintervalle 15000 pour semi variance=15000 et portée =100)
#cette fontcion donne les parametres phi=a sigmasg =sigma2
print(vario.sphe)
#3 eme parametre tausq : erreur de mesure = effet de pepite = hypothese qu'il ny ait pas de fluctuaction en 0 est pas valable mais permet de relaxer
#on ne va pas lutiliser
vario.exp <- variofit(vario.ex,cov.model= "exp", ini.cov.pars=c(15000,100))
print(vario.exp)
x11() ; par(mfrow=c(1,2))
plot(vario.ex,main="Sphérique")
lines.variomodel(cov.model="sphe",cov.pars=c(vario.sphe$cov.pars[1],vario.sphe$cov.pars[2]), nugget=0, max.dist=350)
plot(vario.ex,main="Exponentiel")
lines.variomodel(cov.model="exp",cov.pars=c(vario.exp$cov.pars[1],vario.exp$cov.pars[2]), nugget=0, max.dist=350)
#modele sont satisfesant 
#le modele spherique se rapproche plus des 3 premiers pts
#print(vario.exp) -> derniere ligne =la valeur final du le critere qui a ete optimisé
#print(vario.sphe)
#vario.spe est mieux


# autre methode de fit avec likfit (max de vraisemblance)

# ajustement a l'oeil
#quel modele jaurais si je fixais moi-même les paramètres
x11() ; par(mfrow=c(2,2), mar=c(3,3,1,1), mgp =c (2,1,0))
plot(vario.ex,main="Sphérique")
lines.variomodel(cov.model="sphe",cov.pars=c(15000,100), nugget=0, max.dist=350)
plot(vario.ex,main="Exponentiel")
lines.variomodel(cov.model="exp",cov.pars=c(15000,100), nugget=0, max.dist=350)
plot(vario.ex,main="Exponentiel")
lines.variomodel(cov.model="exp",cov.pars=c(15000,30), nugget=0, max.dist=350)
plot(vario.ex,main="Exponentiel avec pépite")
lines.variomodel(cov.model="exp",cov.pars=c(10000,100), nugget=5000, max.dist=350)

################################################################################
# 4 krigeage

#combinaison lineaire de poids de chque observation :poids de krigeage

#C(ui-uj)=covariance données
#C(ui-u)=matrice de covariance données
#systeme de n equation a n inconnues pour trouver les lambda i

vario.ex <- variog(sic.100, bin.cloud=TRUE)
x11() ; plot(vario.ex,main="")
lines.variomodel(cov.model="spher",cov.pars=c(15000,50), nug=0,max.dist=300) #remplacer c par c(14524,68,222) avec les donnes de vario.sphe
 # sur tout le territoire suisse -> on cree une grille
pred.grid <- expand.grid(seq(0,350, l=51), seq(0,220, l=51)) #1re seq = noeuds en absci #2eme seq=noeud en ord et grid va concatener
rgb.palette <- colorRampPalette(c("blue", "lightblue", "orange", "red"),space = "rgb")#creation palette de couleurs

help(krige.conv)
#je krige sic.100 (base dobservation) loc=location ou je veux calculer mes prediction = ma grille
#krig.option
#remplacer c par c(14524,68)

kc <- krige.conv(sic.100, loc = pred.grid, krige=krige.control(cov.model="spherical",cov.pars=c(15000,50)))
x11() ; image(kc, loc = pred.grid,col =rgb.palette(20) ,xlab="Coord X", ylab="Coord Y",borders=sic.borders,main="Estimation") #prediction sur lens de la suisse
x11() ; image(kc, kc$krige.var, loc = pred.grid,col=rgb.palette(20), xlab="Coord X",ylab="Coord Y",borders=sic.borders, main="Variance de krigeage")
#names(kc)=connaitres les attributs dun objet 
#krige.var=variance de krigeage = erreurde preditcion /  cad dispersion possible quand je fais cette prediction
#estimation montre que l'on uit une direction de 45
#variance erreur est plus forte la il n'y a pas d'observation (pas assez d'observation pour predire)


#retrouve t-on ?
#variogramme et lhomogeneite de l'estimation optimale 
#cartographie = lie des echelles d'heterogeneite
#variogramme 
#n'aurait-on pas pu integrer des infos vario ... ? : possible de faire un mix de plusieurs vario 
#dans les données on sait pas servit de l'laltitude, est ce que l'altitude pourrait fournir une info complementaire sur le jeu de données ?
#acces a l'altitude : sic.100$altitude
#rajouter effet de tendance : krigeage effet le plus simple car le moyenne est connue et constante




# validation du krigeage

#sic.367
#predire les données

kc1<- krige.conv(sic.100, loc = sic.100$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,50)))
#prediction de krigeage avec lens dapprentissage pas sur la grille reguliere 
kc2<- krige.conv(sic.100, loc = sic.367$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,50)))
#predit sur l'ensemble test
x11() ; par(mfrow=c(1,2))
plot(sic.100$data,kc1$predict,xlab="Observe",ylab="Predite", main="Echantillon")#residus sur ma base d'appr
#meme pas de surapprentissage car le modele est trs riche et est basee sur les observations 
#tjr le meme resultat avec ce modele
abline(a=0,b=1,col="red")
plot(sic.367$data,kc2$predict, xlab="Observe",ylab="Predite", main="Tests")#residus sur ma bade de test
#abscis les valeurs, ord:les predictions
#sur cette base le graphe est moins bien 
#pas satisfait : pas content car tres dispersee autour de la diagonale (on veut tres pret)
#mais content :quand les valeurs observees sont fortes je predits une valeurs fortes (resp faible)
#certaine pred sont predites 4x plus que les valeurs reels=valeurs tres mal predites

abline(a=0,b=1,col="red")



# erreur de test (mean square error)
#quantifier l'erreur de test 
#je prends mon vecteur de donnes sic.367 et son vecteur de prediction (chaque coord correspond correctement a sa prediction)
MSE <- mean((sic.367$data-kc2$predict)^2)
print(MSE)#3175.393
#comparer moyenne quadratique avec la variance
#si mse=var cela veut dire que mon modele ne me donne rien, 
#
Q2 <- 1 - MSE/var(sic.367$data) # vaut 0 si c mauvais modele constant cela veut dire que je fit parfaitement mes donnees pas derreur impossible/ le modele ne veut rien dire
print(Q2) 
#bu ameliorer Q2
#outil qui permet de choisir le meilleur modele
################################################################################
# 5 Planification d'experiences



#qui je vire : je peux tirer au hasard NON 

#selectionner les 45 stations qui auront une bonne repartition dans lespace sur tout le territoire
#montecarlo = tirage au hasard
#outil pour juger du bon espacement de point dans un espce = distance mpin entre les points dun echantillon
#fonction mindist de la librairie DiceDesign
#45 stations parmi les 457 ou la distance minimale est maximale 
#je ne veux pas quelle soit proc-he les unes des autres

library(DiceDesign)

###########
# Methode a)

R <- 1e4 # nb de repetitions
critere_opt <- 0
for (i in 1:R){ #1 a 10000
  sel <- sample(1:467, 45) # selection de 45 numeros d'observations parmi 467 (taille de sic.all) (tirage sans remise)
  coords <- sic.all$coords[sel,] #jextrais les coord des emplacements geographiques des stations (il y a 2 coord avec 45 lignes)
  critere <- mindist(coords)# japplique la fonction mindist qui prend en entree une matrice renvoie la plus petites lignes 
  if (critere > critere_opt){
    critere_opt <- critere
    sel_opt <- sel
  }
}

# Visualisation des resultats et calcul de l'erreur de test

# sic.sel contient le resultat de la selection des 45 stations
sic.sel <- sic.all
class(sic.sel) <- "geodata"
sic.sel$coords <- sic.all$coords[sel_opt,]
sic.sel$data <- sic.all$data[sel_opt]
sic.sel$covariate <- sic.all$covariate[[1]][sel_opt]

# sic.reste contient toutes les autres stations
sic.reste <- sic.all
class(sic.reste) <- "geodata"
sic.reste$coords <- sic.all$coords[-sel_opt,]
sic.reste$data <- sic.all$data[-sel_opt]
sic.reste$covariate <- sic.all$covariate[[1]][-sel_opt]

x11() ; points(sic.sel, borders=sic.borders,col="red")
x11() ; points(sic.all, borders=sic.borders,col="black")#je revisualise la moditfication de mes donnees

kc<- krige.conv(sic.sel, loc = pred.grid, krige=krige.control(cov.model="spher",cov.pars=c(15000,50)))
x11() ; image(kc, loc = pred.grid,col=rgb.palette(20) ,xlab="Coord X", ylab="Coord Y",borders=sic.borders,main="Estimation")#krigeage sur ces nouvelles donnees
x11() ; image(kc, kc$krige.var, loc = pred.grid,col=rgb.palette(20), xlab="Coord X",ylab="Coord Y",borders=sic.borders, main="Variance de krigeage")#variance plus diversifiées les pt bleus sont mes stations choisies

#je visualise 
kc1<- krige.conv(sic.sel, loc = sic.sel$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,50)))
kc2<- krige.conv(sic.sel, loc = sic.reste$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,50)))
x11() ; par(mfrow=c(1,2))
plot(sic.sel$data,kc1$predict,xlab="Observe",ylab="Predite", main="Echantillon")
abline(a=0,b=1,col="red")
plot(sic.reste$data,kc2$predict, xlab="Observe",ylab="Predite", main="Tests")#est ce je predis bien maintenant avec ma nouvelle base de testa
abline(a=0,b=1,col="red")

#je recaluclue mon erreur quadratiques
MSE <- mean((sic.reste$data-kc2$predict)^2)
print(MSE)
Q2 <- 1 - MSE/var(sic.reste$data)
print(Q2) 
#valeur bcp moins bonnes avec le budget moins bas donc forcement moins dobservations 

#############

# Methode b)
#2eme methode pour choisir les 47 stations 
 
r <- 28 # rayon du cercle de suppression des points (km)
x <- wspDesign(sic.all$coords, dmin = r, init = "center")$design #je tire un pt au miliu ou au hasard je definis un certain rayon et toutes les sattions dans ce cercle je vire 
#init=pt initial = on prend le point au entre
#avec cette methode il faut trouver le bon rayon 
print(dim(x)[1])

# A faire : trouver le r optimal de maniere automatique (avec une boucle 'while' par exemple)

# on recupere les indices dans sic.all de la selection des 45 stations
sel <- NULL
for (i in 1:dim(x)[1]) sel <- c(sel, which(x[i,1] == sic.all$coords[,1] & x[i,2] == sic.all$coords[,2]))
sel_opt <- as.vector(sel)

# Visualisation des resultats et calcul de l'erreur de test

# sic.sel contient le resultat de la selection des 45 stations
sic.sel <- sic.all
class(sic.sel) <- "geodata"
sic.sel$coords <- sic.all$coords[sel_opt,]
sic.sel$data <- sic.all$data[sel_opt]
sic.sel$covariate <- sic.all$covariate[[1]][sel_opt]

# sic.reste contient toutes les autres stations
sic.reste <- sic.all
class(sic.reste) <- "geodata"
sic.reste$coords <- sic.all$coords[-sel_opt,]
sic.reste$data <- sic.all$data[-sel_opt]
sic.reste$covariate <- sic.all$covariate[[1]][-sel_opt]

x11() ; points(sic.sel, borders=sic.borders,col="red")
x11() ; points(sic.all, borders=sic.borders,col="black")

kc<- krige.conv(sic.sel, loc = pred.grid, krige=krige.control(cov.model="spher",cov.pars=c(15000,50)))
x11() ; image(kc, loc = pred.grid,col=rgb.palette(20) ,xlab="Coord X", ylab="Coord Y",borders=sic.borders,main="Estimation")
x11() ; image(kc, kc$krige.var, loc = pred.grid,col=rgb.palette(20), xlab="Coord X",ylab="Coord Y",borders=sic.borders, main="Variance de krigeage")

kc1<- krige.conv(sic.sel, loc = sic.sel$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,50)))
kc2<- krige.conv(sic.sel, loc = sic.reste$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,50)))
x11() ; par(mfrow=c(1,2))
plot(sic.sel$data,kc1$predict,xlab="Observe",ylab="Predite", main="Echantillon")
abline(a=0,b=1,col="red")
plot(sic.reste$data,kc2$predict, xlab="Observe",ylab="Predite", main="Tests")
abline(a=0,b=1,col="red")

MSE <- mean((sic.reste$data-kc2$predict)^2)
print(MSE)
Q2 <- 1 - MSE/var(sic.reste$data)
print(Q2) #0.62
#algorithme deterministe => on choisit les memes stations

#en comparant les deux graphes on observe que cette methode est meilleure car elle selectionne des stations bien reparties dans la carte
###############
# Methode c)
#methode des k plus proche voisin 
#methode du cafe : client choisit sa table le deuxieme sinstalle le plus loin possible
#maximise la distance avec son voisin le plus proche 
#je tire au hasard la premiere station
#je choisis le la station la plus lointaine et ainsi de suite


library(FNN)

# 1er indice
sel <- sample(1:467, 1)
sic.sel <- sic.all
class(sic.sel) <- "geodata"
sic.sel$coords <- matrix(sic.all$coords[sel,], ncol=2)
sic.reste <- sic.all
class(sic.reste) <- "geodata"
sic.reste$coords <- sic.all$coords[-sel,]

#indices suivants
#
for (i in 1:44){
  j <- which.max(knnx.dist(sic.sel$coords, sic.reste$coords, k=1)) # indice dans sic.reste qui maximise entre ces deux coords statiuons selectionnes avec les statiuons pas encore selectionnees 
  sel <- c(sel, which(sic.reste$coords[j,1] == sic.all$coords[,1] & sic.reste$coords[j,2] == sic.all$coords[,2]))
  sic.sel$coords <- sic.all$coords[sel,]
  sic.reste$coords <- sic.all$coords[-sel,]
}
#algorithme sequentiel un a un 
sel_opt <- as.vector(sel)

# Visualisation des resultats et calcul de l'erreur de test

# sic.sel contient le resultat de la selection des 45 stations
sic.sel <- sic.all
class(sic.sel) <- "geodata"
sic.sel$coords <- sic.all$coords[sel_opt,]
sic.sel$data <- sic.all$data[sel_opt]
sic.sel$covariate <- sic.all$covariate[[1]][sel_opt]

# sic.reste contient toutes les autres stations
sic.reste <- sic.all
class(sic.reste) <- "geodata"
sic.reste$coords <- sic.all$coords[-sel_opt,]
sic.reste$data <- sic.all$data[-sel_opt]
sic.reste$covariate <- sic.all$covariate[[1]][-sel_opt]

x11() ; points(sic.sel, borders=sic.borders,col="red")
x11() ; points(sic.all, borders=sic.borders,col="black")
#les stations seleceionnnees sont presentes dans le sextremites les frontieres
kc<- krige.conv(sic.sel, loc = pred.grid, krige=krige.control(cov.model="spher",cov.pars=c(15000,50)))
x11() ; image(kc, loc = pred.grid,col=rgb.palette(20) ,xlab="Coord X", ylab="Coord Y",borders=sic.borders,main="Estimation")
x11() ; image(kc, kc$krige.var, loc = pred.grid,col=rgb.palette(20), xlab="Coord X",ylab="Coord Y",borders=sic.borders, main="Variance de krigeage")

kc1<- krige.conv(sic.sel, loc = sic.sel$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,50)))
kc2<- krige.conv(sic.sel, loc = sic.reste$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,50)))
x11() ; par(mfrow=c(1,2))
plot(sic.sel$data,kc1$predict,xlab="Observe",ylab="Predite", main="Echantillon")
abline(a=0,b=1,col="red")
plot(sic.reste$data,kc2$predict, xlab="Observe",ylab="Predite", main="Tests")
abline(a=0,b=1,col="red")

MSE <- mean((sic.reste$data-kc2$predict)^2)
print(MSE)
Q2 <- 1 - MSE/var(sic.reste$data)
print(Q2)
#meilleur Q2 obtenu

#ne trouvez vous pas quil y a des aspects prob ? : 1) les distances en traversant des frontieres est ce resonnable en prennant en compte les frontieres et les altitudesd qui ne sont pas pris en compte
#pas daspects geographioues pris en compte par rapport au correlation /freontiere pa s pris en compte
#2) choisis 47 stations etape variogrammes : pas faite
#variogramme avec que 47 stations probleme : variogramme a besoin de stations pas trop eloignees (paire de pts de dist de 10Km de rayon )
#separation des pts plus de petites distance on uara pas les memes pts durant le fit : surtout pour le fit 
#3) pas reflechis a la valeur de la sortie lorsque jai fais ma selection => distribution de quantite de pluie pas prise en compte pendant la selection
#distribution de la sortie pas prise en compte 
#uniquement la distance geog est prise en compte  