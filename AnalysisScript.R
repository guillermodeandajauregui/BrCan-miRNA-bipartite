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
library(bipartite)
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

########################################
#export for Networkx bipartite analysis
########################################

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

