#library(ggplot2)
library(igraph)
library(poweRlaw)
library(gridExtra)
library(gsl)
library(numDeriv)
library(doParallel)# for parallel processing!
library(foreach)

#library(lmtest) #For likelihood ratio test- model comparison using lrtest(); but not working well!!!
#library("poweRlaw")


#setwd("/home/wilber/networks/Grids/")

elegansGrid <- read.table("~/networks/Grids/elegansGridWS220ChIPY1H", header=F)
el=as.matrix(elegansGrid)
el[,1]=as.character(el[,1])
el[,2]=as.character(el[,2])
#el<-as.data.frame(el)

g=graph.edgelist(el,directed=TRUE)
#Simplify graph
g_Simplified<-simplify(g, remove.multiple=TRUE, remove.loops=FALSE)


#FIT POWER-LAW FOR THE 100% GRID JUST ONCE!
degree<-degree(g_Simplified, mode="out")
degree<-degree[degree!=0]
degreePowerlawFit = displ$new(degree)
estPL = estimate_xmin(degreePowerlawFit) #estimate parameters! Recommended than "estimate_pars" which uses the smallest Xmin
degreePowerlawFit$setXmin(estPL)
degreePowerlawFit$setPars(estimate_pars(degreePowerlawFit))

############################################################################################

#(I)Determine the KS distance between 100% and EACH sample graph (at a given %), ..
#and then determine the mean of KS difference
#Statistics= determine whether the mean KS difference is ==0 or not!
#Use one-sample t-test: Ho:μx=0, Ha:μx≠0; where 'μx' is the mean KS difference
#Assumption: By the central limit theorem, if the sampling of the parent population .. 
#is independent then the sample 'means' will be approximately normal ..
#especially for large sample size 'n'
#############################################################################################

###use apply instead of for loops!!!

#########################

###Set-up clusters!
#detectCores() detect how many cores are available
# Create cluster with 6 cores
cl <- makeCluster(15) 
# Register cluster
registerDoParallel(cl)
# getDoParWorkers() # Find out how many cores are being used


percentages<-c(0.95,0.90,0.85,0.80,0.75,0.70,0.65,0.60,0.55,0.50,0.45,0.40,0.35,0.30,0.25,0.20,0.15,0.10, 0.05, 0.010, 0.005, 0.0025, 0.001)
for (j in 1:23) { 
  
KS_distances<-numeric(1000)
KS_pValue<-numeric(1000)

x<-foreach (i = 1:1000)  %dopar% {
  
  library(igraph)
  library(poweRlaw)
  library(gridExtra)
  library(gsl)
  library(numDeriv)

  random<-sample(1:nrow(el), percentages[j]*nrow(el), replace = F)
  el_sample<-el[random,]
  el_sample[,1]=as.character(el_sample[,1])
  el_sample[,2]=as.character(el_sample[,2])
  g_sample=graph.edgelist(el_sample,directed=TRUE)
  g_sample_simplified<-simplify(g_sample, remove.multiple=TRUE, remove.loops=FALSE)
  
  #Fit a discrete power-law distribution on data out-degree data using "poweRlaw":
  degree_sample<-degree(g_sample_simplified, mode="out")
  degree_sample<-degree_sample[degree_sample!=0]
  
  degreePowerlawFit_sample = displ$new(degree_sample)
  estPL_sample = estimate_xmin(degreePowerlawFit_sample)
  
  degreePowerlawFit_sample$setXmin(estPL_sample)
  degreePowerlawFit_sample$setPars(estimate_pars(degreePowerlawFit_sample))
  
  #####Likelihood ratio tests to between 100% grid vs 90% grid
  #Option 1: adopt the "compare_distributions()" LR test 
  #Option 2: Use KS statistic: compare distances between the 100% and the sample, and later get the mean distance
  
  #Option 2:
  #ks.test(dist_pdf(degreePowerlawFit), dist_pdf(degreePowerlawFit_sample))
  KS<-ks.test(dist_cdf(degreePowerlawFit), dist_cdf(degreePowerlawFit_sample))
  #KS<-ks.test(degree, degree_sample)
  #KSstatistic<-KS$statistic
  #KSpvalue<-KS$p.value
  
  KS_distances[i]<-KS$statistic #append KS distances to the vector!
  KS_pValue[i]<-KS$p.value #append p value 
  
  results<-list(KS_distances, KS_pValue)
  return(results)
  
  random<-NULL
  
}


#'Assemble' all stats
#KS_D_stats<-numeric(100)
for (i in 1:1000){KS_distances[i]<-(x[[i]][[1]][i])}
for (i in 1:1000){KS_pValue[i]<-(x[[i]][[2]][i])}
#end time
#print(Sys.time()-strt)
#KS_distances_t_test<-t.test(KS_distances) # t-test: test whether the average of all values (differences) is different from zero
#KS_distances_mean<-mean(KS_distances) #get mean KS distance

###Next, determine the number of KS Dstat p values (number of subgraphs) less than 0.01
FractionPvalueGrt0.01<-length(KS_pValue[KS_pValue <=0.01])/length(KS_pValue) #get a cut-off of 95%???

print(c(percentages[j], mean(KS_distances), mean(KS_pValue), FractionPvalueGrt0.01))


}
