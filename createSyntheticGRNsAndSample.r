require(igraph)
library(doParallel)# for parallel processing!
library(foreach)
require(poweRlaw)

#Objective: create a 'synthetic' yeast GRN, and sample sub-nets
#Compare two fitted distributions

#Yeast:
no_of_interactions = 70830 #Yeast
no_of_nodes = 6275
TFs_genome = 301
out_degree =  2.00
in_degree = 6.41


#set.seed(2016)
#g<-sample_fitness_pl(no_of_nodes,no_of_interactions,out_degree,in_degree)



gridSimulation<-function(grid,iterations){ #processors and interations are numerical
  
  #set.seed(2016)
  cl <- makeCluster(12)
  # Register cluster
  registerDoParallel(cl)
  
  EdgeSelectionSamplingData<-data.frame()
  
  #percentages<-c(0.1,0.50,0.9)
  percentages<-seq(0.02,0.98,0.02) #fraction of the full grid to be sampled, REMEMBER TO REVERSE
  
  for (j in 1:length(percentages)) {
    
    #KS_distances<-numeric(iterations)
    #KS_pValue<-numeric(iterations)
    alpha<-numeric(iterations)
    CD_pValue<- numeric(iterations)
    CD_stat <- numeric(iterations)
    #KS.power.law.fit.pValue<-numeric(iterations)
    #cc_sample<-numeric(iterations)
    
    x<-foreach (i = 1:iterations) %dopar%  {
      
      library(igraph) #require????
      require(ggplot2)
      require(poweRlaw)
      #library(gridExtra)
      #library(gsl)
      library(numDeriv)
      
      #Function for subgrid construction:
      SubgridConstruction<-function (edgeListSubgrid) {
        edgeListSubgrid[,1]=as.character(edgeListSubgrid[,1])
        edgeListSubgrid[,2]=as.character(edgeListSubgrid[,2])
        subgrid=graph.edgelist(edgeListSubgrid,directed=TRUE)
        #g_sample_simplified<-simplify(g_sample, remove.multiple=TRUE, remove.loops=FALSE) #deals with loops
        return(subgrid)
      }
      #Fitting power law with 'igraph':
      gridPowerLawIgraph<-function(grid,mode) { #mode takes options: "in", "out", or "all"
        degree<-degree(grid, mode=mode)
        degree<-degree[degree!=0]
        powerLawFit<-power.law.fit(degree)
        return(powerLawFit)
      }
      
      # degree<-degree[degree>0]
      
      #A function to fit power law dist using the poweRlaw package:
      powerlawfit<-function(degree) {
        poweRlawFit = displ$new(degree)
        estPL = estimate_xmin(poweRlawFit) #estimate parameters! Recommended than "estimate_pars" which uses the smallest Xmin
        poweRlawFit$setXmin(estPL) #set xmin
        poweRlawFit$setPars(estimate_pars(poweRlawFit)) #Set parameters
        return(poweRlawFit)
      }
      
      random<-sample(1:nrow(get.edgelist(grid)), percentages[j]*nrow(get.edgelist(grid)), replace = F) #Generate random #s provide/store the seed!!!
      edgeListSubgrid<-get.edgelist(grid)[random,] #Sample from graph edge list by indexing
      subgrid<-SubgridConstruction(edgeListSubgrid = edgeListSubgrid) #calls the "edgeListSubgrid" function
      
      #Fit a discrete power-law distribution on degree
      #subgridPowerlawStats<-gridPowerLawIgraph(subgrid,"out") #using igraph
      
      
      degree<-degree(subgrid, mode = 'out') #get degree#
      degree<-degree[degree>0]
      #Execute function for fitting a power law using poweRlaw:
      poweRlawFit<-powerlawfit(degree)
      
      ##Fit poisson distribution:
      
      #PoissonFit = dispois$new(degree)
      #PoissonFit$setXmin(poweRlawFit$getXmin())  #To compare distributions, both must have same lower threshold:
      #PoissonFit$setPars(estimate_pars(PoissonFit))
      
      
      # Fit exponential distribution
      ExpFit = disexp$new(degree)
      ExpFit$setXmin(poweRlawFit$getXmin())  #To compare distributions, both must have same lower threshold:
      ExpFit$setPars(estimate_pars(ExpFit))
      
      #Compare the two models:
      compDist = compare_distributions(poweRlawFit, ExpFit)
      
      #get p value
      CD_pValue[i] <- compDist$p_two_sided #Low (less than 0.05?) p value denotes a rejection of null hypothesis, Power law favoured. But large p values denote no model favored!
      
      CD_stat[i] <- compDist$test_statistic   # Vuong's test statistic
      
      #      KS_pValue[i]<-subgridPowerlawStats$KS.p #append KS statistic to the KS_distances vector
      alpha[i]<-poweRlawFit$pars #alpha parameter!
      #cc_sample[i]<-transitivity(subgrid) #Transitivity
      results<-list(alpha,CD_pValue, CD_stat)
      return(results)
      
      
    }
    #Extract results from the list 'x'
    for (i in 1:iterations){alpha[i]<-(x[[i]][[1]][i])}
    for (i in 1:iterations){CD_pValue[i]<-(x[[i]][[2]][i])}
    for (i in 1:iterations){CD_stat[i]<-(x[[i]][[3]][i])}
    
    iterationResults<-data.frame(cbind(percentages[j],alpha,CD_pValue, CD_stat))
    EdgeSelectionSamplingData<-rbind(EdgeSelectionSamplingData,iterationResults)
    
    
  }
  colnames(EdgeSelectionSamplingData)<-c("percentage","alpha","CD_pValue","CD_stat")
  return(EdgeSelectionSamplingData)
  
  stopCluster(cl)
}

#data<-data.frame()
#grids<-list()
#gridsNum=10
#for (i in 1:gridsNum){
#  grid<-sample_fitness_pl(no_of_nodes,no_of_interactions,out_degree,in_degree)
#  dat<-gridSimulation(grid,1000)
#  dat$gridNumber<-i
#  data<-rbind(data,dat)
#}

grid<-sample_fitness_pl(no_of_nodes,no_of_interactions,out_degree,in_degree)
data<-gridSimulation(grid,1000)
write.table(data,file="grid_sampling_synthetic_Comp_Distr_yeast", sep = "\t", row.names = FALSE)


