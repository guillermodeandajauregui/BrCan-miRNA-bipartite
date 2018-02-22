#network analyzer clone
#and assorted network functions 

library(igraph)

NetworkAnalyzer = function(g, directed = FALSE, skip.betweenness = FALSE, workaround.betweenness = FALSE){
  if(directed == FALSE){
    V(g)$degree <- igraph::degree(g)
    
    if(workaround.betweenness == TRUE && skip.betweenness == FALSE){
      copy = g
      E(copy)$weight = rep(1, length(E(copy)))
      V(g)$betweenness <- igraph:: betweenness(copy, directed = FALSE)
    }
    
    if(workaround.betweenness == FALSE && skip.betweenness == FALSE){
      V(g)$betweenness <- igraph:: betweenness(g, directed = FALSE)
    }
    
    V(g)$closeness <- igraph:: closeness(g)
    
    average_path_length = average.path.length(g, directed = FALSE)
    ClusteringCoefficient = transitivity(g, type = "global")
    GraphDiameter = diameter(g, directed = FALSE)
    
    component = components(g)
    V(g)$component = component$membership
    
    infomap = infomap.community(g)
    V(g)$infomap = infomap$membership
    
    results = list(g = g, 
                   average_path_length = average_path_length, 
                   ClusteringCoefficient = ClusteringCoefficient, 
                   GraphDiameter = GraphDiameter,
                   component = component,
                   infomap = infomap)
  }
  if(directed == TRUE){
    V(g)$degree <- igraph::degree(g, mode = "all")
    V(g)$indegree <- igraph::degree(g, mode = "in")
    V(g)$outdegree <- igraph::degree(g, mode = "out")
    
    if(workaround.betweenness == TRUE && skip.betweenness == FALSE){
      copy = g
      E(copy)$weight = rep(1, length(E(copy)))
      V(g)$betweenness <- igraph:: betweenness(copy, directed = TRUE)
    }
    
    if(skip.betweenness == FALSE){
      V(g)$betweenness <- igraph:: betweenness(g, directed = TRUE)
    }
    
    V(g)$closeness <- igraph:: closeness(g, mode = "all")
    V(g)$closeness_in <- igraph:: closeness(g, mode = "in")
    V(g)$closeness_out <- igraph:: closeness(g, mode = "out")
    
    average_path_length = average.path.length(g, directed = TRUE)
    ClusteringCoefficient = transitivity(g, type = "global")
    GraphDiameter = diameter(g, directed = TRUE)
    
    component = components(g)
    V(g)$component = component$membership
    
    infomap = infomap.community(g)
    V(g)$infomap = infomap$membership
    
    results = list(g = g, 
                   average_path_length = average_path_length, 
                   ClusteringCoefficient = ClusteringCoefficient, 
                   GraphDiameter = GraphDiameter,
                   component = component,
                   infomap = infomap)
  }
  return(results)
}