require(igraph)
library(doParallel)# for parallel processing!
library(foreach)
library(poweRlaw)

#Objective: Predict number of interactions (PDIs) for complete GRNs 

#Load interactions:
elegansGrid1 <- read.table("DrosophilaGridChIPY1H_observedGrid", header=F) 
#elegansGrid1 <- read.table("elegansGridWS220ChIPY1H_observedGrid", header=F) 
#elegansGrid1 <- read.table("yeastGridChipEvidence.txt", header=F)   
#elegansGrid1 <- read.table("yeastGrid08-05-15-version2_observedGrid", header=F) 
#elegansGrid1 <- read.table("Arabidopsis_GRG_Aug_29_16_2_observedGrid", header=F) 

el=as.matrix(elegansGrid1)
el[,1]=as.character(el[,1])
el[,2]=as.character(el[,2])

g=graph.edgelist(el,directed=TRUE)
g_Simplified<-simplify(g, remove.multiple=TRUE, remove.loops=FALSE)

##Sample 70% of the network:

alpha<-numeric(1000)
#gof<-numeric(1000)
for (i in 1:1000) {
  
  random<-sample(1:nrow(el), 0.7*nrow(el), replace = F) 
  el_sample<-el[random,]
  el_sample[,1]=as.character(el_sample[,1])
  el_sample[,2]=as.character(el_sample[,2])
  g_sample=graph.edgelist(el_sample,directed=TRUE)
  g_sample_simplified<-simplify(g_sample, remove.multiple=TRUE, remove.loops=FALSE)
  
  #Fit a discrete power-law distribution on data out-degree data using "poweRlaw":
  degree_sample<-degree(g_sample_simplified, mode="out")
  degree_sample<-degree_sample[degree_sample!=0]
  
  
  #Fit a discrete power-law distribution on degree, using poweRlaw:
  degree_sample<-degree(subgrid, mode="out")
  degree_sample<-degree_sample[degree_sample!=0]
  degreePowerlawFit_sample = displ$new(degree_sample)
  estPL_sample = estimate_xmin(degreePowerlawFit_sample)
  degreePowerlawFit_sample$setXmin(estPL_sample)
  degreePowerlawFit_sample$setPars(estimate_pars(degreePowerlawFit_sample))
  
  alpha[i]<-estPL_sample$pars
  #  gof[i]<-estPL_sample$gof
  
}


#Calculate the CIs based on the exponent estimates of 70% of the network:
error <- qt(0.975,df=1000-1)*sd(alpha)/sqrt(1000) #alpha is the stimated exponent
alpha_min<-mean(alpha) - error #Lower bound of 95% CI
alpha_max<-mean(alpha) + error #Upper bound of 95% CI

#Sample from power law with exponents alpha_min and alpha_max
#Min alpha

degree<-degree(g_Simplified, mode="out")
degree<-degree[degree!=0]
degreePowerlawFit = displ$new(degree)
estPL = estimate_xmin(degreePowerlawFit)
degreePowerlawFit$setXmin(estPL)
degreePowerlawFit$setPars(estimate_pars(degreePowerlawFit))

#set.seed(2017)
#Use xmin for the fit to power law of GRNs for different organisms 
#Use their corresponging alphas/exponents
xmin=800 #Drosophila
#xmin=1500 # C. elegans 
#xmin=40 #Yeast
#xmin=1 #Arabidopsis
#xmax=17000 # Drosophila Genome size
#xmax=20470 # C. elegans Genome size
#xmax=5000 # Yeast Genome size
#xmax=27655 #TAIR10
exponent = estPL$pars
#exponent = 2
nodeNumber = 166 #Drosophila
#nodeNumber = 1052 #Drosophila TFs in Genome (FlyTF.org)
#nodeNumber = 219 #C. elegans
#nodeNumber = 934 #C. elegans TFs in Genome 
#nodeNumber =    length(degree) # 151 #Yeast 
#nodeNumber = 301 #Yeast TFs in Genome  
#nodeNumber =    length(degree) # 342 #Arabidopsis 
#nodeNumber = 2451 #TFs in Arabidopsis 

#Draw RNs 10000 times, determining the PDIs for each iteration:
PDIs<-numeric(10000)
for (j in 1:10000){
  
  randomNumbers=vector()
  for (i in 1:nodeNumber){
    repeat {
      randomNumber=(xmin*runif(1)^(-1/(exponent-1)))
      if(randomNumber<=xmax){
        break
      }
      
      
    }
    
    randomNumbers[i]<-randomNumber
  }
  randomNumbersDiscrete<-round(randomNumbers)
  
  #PDIs[j]<-sum(randomNumbers)
  PDIs[j]<-sum(randomNumbersDiscrete)
  
}



#Hypothesis testing: estimate the plausibility of th estimate # of PDIs using ... 
#H0: PDI_obs=PDIs_simulated (Use Z-test????)

z_score<-(mean(PDIs) - ecount(g_Simplified))/sd(PDIs)

pvalue2sided=2*pnorm(-abs(z_score))
