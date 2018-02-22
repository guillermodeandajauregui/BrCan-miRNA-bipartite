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

#change names of miRNAs 

V(an_sanos$g)$name = gsub(pattern = "-", replacement = "GUIONMEDIO", x = V(an_sanos$g)$name)
V(an_sanos$g)$name = gsub(pattern = ".MIMAT", replacement = "puntoMIMAT", x = V(an_sanos$g)$name)

V(an_casos$g)$name = gsub(pattern = "-", replacement = "GUIONMEDIO", x = V(an_casos$g)$name)
V(an_casos$g)$name = gsub(pattern = ".MIMAT", replacement = "puntoMIMAT", x = V(an_casos$g)$name)

V(an_casos$g)$name
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

anx_sanos = read.graph(file = "results/anx_sanos.gml", "gml")
anx_casos = read.graph(file = "results/anx_casos.gml", "gml")

#change labels back to original names
V(anx_sanos)$label = gsub(pattern = "puntoMIMAT", 
                          replacement = ".MIMAT", 
                          x = V(anx_sanos)$label)

V(anx_sanos)$label = gsub(pattern = "GUIONMEDIO", 
                          replacement = "-", 
                          x = V(anx_sanos)$label)

V(anx_casos)$label = gsub(pattern = "puntoMIMAT", 
                          replacement = ".MIMAT", 
                          x = V(anx_casos)$label)

V(anx_casos)$label = gsub(pattern = "GUIONMEDIO", 
                          replacement = "-", 
                          x = V(anx_casos)$label)

V(anx_sanos)$name = V(anx_sanos)$label
V(anx_casos)$name = V(anx_casos)$label

########################################
#Plot
#degree
#CCdot
#CCribbon
#Redundancy
########################################

