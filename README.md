# Network-properties-and-graph-sampling
These R codes accompany the Ouma et al (2018) publication in PLoS Computational Biology entitled "Topological and statistical analyses of gene regulatory networks reveal unifying yet quantitatively different emergent properties". 

1) randomEdgeSampling.r: randomly samples subnetworks - and fits a power law function on their out-degree - given an observed gene regulatoru network (GRN)

2) predictPDIs.r: predicts the number of protein-DNA interactions (PDIs) given the power law exponent and the number of transcription factors (TFs) of an organism

3) createSyntheticGRNsAndSample.r: given the predicted number of PDIs for a complete yeast GRN, this code creates a complete synthetic (in silico) GRN and samples sub-nets from the GRN, fitting power law function in each sampling.

4) LorenzCurves.r: Use of Lorenz curves to describe the 'inequality' of TF binding events in GRNs
