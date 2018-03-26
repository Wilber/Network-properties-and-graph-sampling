library(ineq)
library(igraph)
library(poweRlaw)
library(ggplot2)
library(reshape)
library(scales)
library(grid)

#setwd("C:/Users/User/Box Sync/WilburStuff/Network Analysis/MulltivariateRegression/modeling/")

#OBJECTIVE: TEST THE HYPOTHESIS THAT THE EXPONENT ALPHA DESCRIBES THE 'FAIRNESS' IN CONNECTIVITY 

##MOTIVATION: 
# (1) https://networkscience.wordpress.com/2012/04/19/power-law-paradox-power-law-exponent-does-not-mean-what-you-think-it-means/ 

# (2) GINI COEFFICIENT: http://www.josechristian.com/programming/inequality-and-lorenz-curve-r/
########################################################################################
#SIMULATION:
#Generate degree distributions with specified alpha:

degreeSimulatedAlpha2<-rpldis(1000,1,2) #found in poweRlaw package!
degreeSimulatedAlpha2_5<-rpldis(1000,1,2.5)
degreeSimulatedAlpha3<-rpldis(1000,1,3)
degreeSimulatedAlpha3_5<-rpldis(1000,1,3.5)
degreeSimulatedAlpha4<-rpldis(1000,1,4)
degreeSimulatedAlpha4_5<-rpldis(1000,1,4.5)

#lapply(1:10, function(i) rpldis(10, 1, 2))

lcolc2<-Lc(degreeSimulatedAlpha2)
lcolc2_5<-Lc(degreeSimulatedAlpha2_5)
lcolc3<-Lc(degreeSimulatedAlpha3)
lcolc3_5<-Lc(degreeSimulatedAlpha3_5)
lcolc4<-Lc(degreeSimulatedAlpha4)
lcolc4_5<-Lc(degreeSimulatedAlpha4_5)


##Use ggplot:
#lcdf <- data.frame(L = lcolc3$L, p = lcolc3$p, Uprob = c(1:length(lcolc3$L)/length(lcolc3$L)))
#create dat frames for each lorenz of each alpha 
lcdf2 <- data.frame(L = lcolc2$L, p = lcolc2$p, alpha = "2")
lcdf2_5 <- data.frame(L = lcolc2_5$L, p = lcolc2_5$p, alpha = "2.5")
lcdf3 <- data.frame(L = lcolc3$L, p = lcolc3$p, alpha = "3")
lcdf3_5 <- data.frame(L = lcolc3_5$L, p = lcolc3_5$p, alpha = "3.5")
lcdf4 <- data.frame(L = lcolc4$L, p = lcolc4$p, alpha = "4")
lcdf4_5 <- data.frame(L = lcolc4_5$L, p = lcolc4_5$p, alpha = "4.5")


#
lcdfAll<-rbind(lcdf2,lcdf2_5,lcdf3,lcdf3_5,lcdf4,lcdf4_5)


##PLOTS:
#par(mar=c(6, 5, 4, 2), mgp=c(2, 0.4, 0), tck=-.01,
#    cex.axis=0.6, las=1, lwd=1.0)

plot(lcdf2$L ~ lcdf2$p, type = "l",lwd=3,
     pch=18, col="grey", panel.first=grid(col="grey80"),
     ylab="", xlab="", main=NA, xaxt='n', yaxt='n')
axis(side=1, at=seq(0,1,0.2), labels=seq(0,1,0.2), cex.axis=1.0)
axis(side=2, at=seq(0,1,0.2), labels=seq(0,1,0.2), cex.axis=1.0)



lines(lcdf2_5$L ~ lcdf2_5$p, col="green", lwd=3)
lines(lcdf3$L ~ lcdf3$p, col="red", lwd=3)
lines(lcdf3_5$L ~ lcdf3_5$p, col="blue", lwd=3)
lines(lcdf4$L ~ lcdf4$p, col="aquamarine4", lwd=3)
lines(lcdf4_5$L ~ lcdf4_5$p, col="burlywood4", lwd=3)

#Add grid lines
abline(h=seq(0,1,0.2),v=seq(0, 1, 0.2), col="gray", lty=3)

#Add diagonal line:
x<-c(0.0,0.25,0.50,0.75,1.0)
y<-c(0.0,0.25,0.50,0.75,1.0)
diagonal<-as.data.frame(cbind(x,y))
lines(diagonal$y ~ diagonal$x, col="black", lwd=3)

legend(0.001,1,legend=c("2","2.5","3","3.5","4","4.5"), 
       pch=15:18, col=c("grey80","green","red","blue",
                        "aquamarine4","burlywood4" ),box.lty=0,cex=0.9)

title(xlab="Proportion of nodes", line=1.6, cex.lab=1.0)
title(ylab="Proportion of edges", line=2.0, cex.lab=1.0)

#export with 5.61 X 5.61


#LORENZ CURVES OF ACTUAL GRGS:
elegansGrid <- read.table("DrosophilaGridChIPY1H", header=F)
elegansGrid <- read.table("elegansGridWS220ChIPY1H", header=F)
elegansGrid <- read.table("yeastGrid08-05-15-version2.txt", header=F) #yeast grid!
elegansGrid <- read.delim("arabidopsisInteractions_lowerCase.txt", header=F)
####Alon's network
elegansGrid <- read.table("~/networks/Grids/yeastGrid08-10-15-Uri-Alon.txt", header=F) #Alon's network, collection of networks here: http://www.weizmann.ac.il/mcb/UriAlon/download/collection-complex-networks
elegansGrid<-elegansGrid[,1:2]
####


el=as.matrix(elegansGrid)
el[,1]=as.character(el[,1])
el[,2]=as.character(el[,2])


g=graph.edgelist(el,directed=TRUE)
g_Simplified<-simplify(g, remove.multiple=TRUE, remove.loops=FALSE)
degree<-degree(g_Simplified, mode="out")
degree<-degree[degree!=0]

#First get degree for each genome
degreeDrosophila<-degree
degreeElegans<-degree
degreeYeast<-degree
degreeArabidopsis<-degree

degreeDrosophila<-(as.data.frame(degreeDrosophila))
degreeElegans<-(as.data.frame(degreeElegans))
degreeYeast<-(as.data.frame(degreeYeast))
degreeArabidopsis<-(as.data.frame(degreeArabidopsis))

#IMPORTANT:
#Fit power laws to get the xmin, so that use degree >=xmin for Lorenz curves
xminDrosophila<-power.law.fit(degreeDrosophila$degreeDrosophila)$xmin
xminElegans<-power.law.fit(degreeElegans$degreeElegans)$xmin
xminYeast<-power.law.fit(degreeYeast$degreeYeast)$xmin
xminArabidopsis<-power.law.fit(degreeArabidopsis$degreeArabidopsis)$xmin

degreeDrosophila<-degreeDrosophila[degreeDrosophila>xminDrosophila]
degreeElegans<-degreeElegans[degreeElegans>xminElegans]
degreeYeast<-degreeYeast[degreeYeast>xminYeast]
degreeArabidopsis<-degreeArabidopsis[degreeArabidopsis>xminArabidopsis]


lorenzDrosophila<-Lc(degreeDrosophila)
lorenzElegans<-Lc(degreeElegans)
lorenzYeast<-Lc(degreeYeast)
lorenzArabidopsis<-Lc(degreeArabidopsis)


##PLOTS:
#par(mar=c(6, 5, 4, 2), mgp=c(2, 0.4, 0), tck=-.01,
#    cex.axis=0.6, las=1, lwd=1.0)


plot(lorenzElegans$L ~ lorenzElegans$p, type = "l",lwd=3,
     pch=18, col="blue", panel.first=grid(col="grey80"),
     ylab="", xlab="", main=NA, xaxt='n', yaxt='n')
axis(side=1, at=seq(0,1,0.2), labels=seq(0,1,0.2), cex.axis=1.0)
axis(side=2, at=seq(0,1,0.2), labels=seq(0,1,0.2), cex.axis=1.0)


lines(lorenzDrosophila$L ~ lorenzDrosophila$p, col="green", lwd=3)
lines(lorenzYeast$L ~ lorenzYeast$p, col="red", lwd=3)

#lines(lorenzArabidopsis$L ~ lorenzArabidopsis$p, col="grey", lwd=2)

#xticks <- seq(0.0, 1.0, 0.2)
#axis(1, at = xticks, labels = xticks)
#Add grid lines
abline(h=seq(0,1,0.2),v=seq(0, 1, 0.2), col="gray", lty=3)

#Add diagonal line:
x<-c(0.0,0.25,0.50,0.75,1.0)
y<-c(0.0,0.25,0.50,0.75,1.0)
diagonal<-as.data.frame(cbind(x,y))
lines(diagonal$y ~ diagonal$x, col="black", lwd=3)


legend(0.001,1,legend=c("C. elegans, 4.12","D. melanogaster, 3.04","S. cerevisiae, 2.0"), 
       pch=15:18, col=c("blue","green","red"),box.lty=0,cex=0.9, text.font=3)

title(xlab="Proportion of TFs", line=1.6, cex.lab=1.0)
title(ylab="Proportion of TGs", line=2.0, cex.lab=1.0)

