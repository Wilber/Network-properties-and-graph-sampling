library(ggplot2)
library(igraph)
library(poweRlaw)
library(gridExtra)
#library(gsl)
library(numDeriv)
library(doParallel)# for parallel processing!
library(foreach)

##OBJECIVE: Sample subnetworks from observed networks, and determine their scale-free properties

#Load interactions (Grid), two columns of TFs and TGs
Grid <- read.table("yeastGridYeasTract.txt", header=F) #Yeast GRG from YeastTract bulk data
el=as.matrix(Grid)
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

###Set-up compute clusters!
detectCores() #detect how many cores are available
# Create cluster with 28 cores
cl <- makeCluster(28) 
# Register cluster
registerDoParallel(cl)
getDoParWorkers() # Find out how many cores are being used

#Define the size (interactions) of the sub-nets 
percentages<-seq(0.02,0.98,0.02)


data<-data.frame()

#For each sub-net size, sample 1000 sub-nets
for (j in 1:length(percentages)) { 
  
  #strt<-Sys.time() #start time
  KS_D_stats<-numeric(1000)      #empty vector to populate KS statistics for 1000 sampling events!
  x<-foreach (i = 1:1000) %dopar% {            #'n' == 1000; 1000 subgraphs for jth %
    
    library(igraph)
    library(poweRlaw)
    library(gridExtra)
    #library(gsl)
    library(numDeriv)  
    
    random<-sample(1:nrow(el), percentages[j]*nrow(el), replace = F) #0.0001 gives sign. diff. 
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
    
    #Get the empirical KS D statistic
    #KS_D_stats[i]<-estPL_sample$KS #append the impirical KS D stat of sample graph to the growing vector os sample D stats!
    
    #Get alpha
    
    KS_D_stats[i]<-estPL_sample$pars
    return (KS_D_stats)
    
    random<-NULL
    
  }
  
  #'Assemble' all stats
  #KS_D_stats<-numeric(100)
  #(x[[i]][[1]][i])
  for (i in 1:1000){KS_D_stats[i]<-(x[[i]][i])}
  dat<-cbind(percentages[j],KS_D_stats)
  data<-rbind(data,dat)
  
}
names(data)<-c("percentage","alpha")

stopCluster(cl)

write.table(data,file="gridSamplingPoweRLawYeasTract", sep = "\t", row.names = FALSE)
