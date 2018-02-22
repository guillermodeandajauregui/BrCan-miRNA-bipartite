######################################
#Bipartite network Analysis
#miRNA - mRNA
#in Breast Cancer
######################################

######################################
#libraries
######################################
library(data.table)
library(tidyverse)
library(igraph)
#library(bipartite)
source("NAC.R")
######################################

######################################
#Load networks
######################################
#as edge list
sanos = fread(path_sanos)
casos = fread(path_casos)

sanos = sanos[, c(3,1)]
casos = casos[, c(3,1)]

#make igraph
sanos_nw = igraph::graph_from_edgelist(el = as.matrix(sanos))
casos_nw = igraph::graph_from_edgelist(el = as.matrix(casos))


######################################
#General analysis 
######################################
an_sanos = NetworkAnalyzer(sanos_nw, directed = TRUE)
an_casos = NetworkAnalyzer(casos_nw, directed = TRUE)

######################################
#color nodes
######################################
V(an_sanos$g)[grep(pattern = "hsa-mir", 
                   x = V(an_sanos$g)$name)]$type= "miRNA" 
V(an_sanos$g)[grep(pattern = "hsa-mir", 
                   x = V(an_sanos$g)$name, 
                   invert = TRUE)]$type= "gene" 


V(an_casos$g)[grep(pattern = "hsa-mir", 
                   x = V(an_casos$g)$name)]$type= "miRNA" 
V(an_casos$g)[grep(pattern = "hsa-mir", 
                   x = V(an_casos$g)$name, 
                   invert = TRUE)]$type= "gene" 

########################################
#export for Networkx bipartite analysis
########################################

igraph::write.graph(graph = an_sanos$g, 
                    file = "results/sanos.gml", 
                    format = "gml"
                    )

igraph::write.graph(graph = an_casos$g, 
                    file = "results/casos.gml", 
                    format = "gml"
)

########################################
#outside scripts
########################################

#call the bash script for modification

system(command = "./rgml2nx.sh results/sanos.gml results/sanos_nx.gml")
system(command = "./rgml2nx.sh results/casos.gml results/casos_nx.gml")

#Networkx analysis

system(command = "python BipartiteAnalysis.py results/sanos_nx.gml results/anx_sanos.gml")
system(command = "python BipartiteAnalysis.py results/casos_nx.gml results/anx_casos.gml")
########################################
#Reload Networkx analyzed 
########################################

########################################
#Plot
#degree
#CCdot
#CCribbon
#Redundancy
########################################

