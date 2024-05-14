calc_nodal_metric = function(origgraph, community_attr){
  allmetrics <- list()
  allmetrics$local.clustering <- transitivity(origgraph, type = "local") #defaults to read edge weights when available
  allmetrics$strength <- strength(origgraph, v = V(origgraph))
 
  #only compute community-related statistics if the expected community attribute is available
  if (!is.null(vertex_attr(origgraph, community_attr))) {
    allmetrics$community.membership <- vertex_attr(origgraph, community_attr)
    
    ##nodal within vs between modules
    wibw <- wibw_module_degree(origgraph, community_attr=community_attr)
    allmetrics[["within.module.deg"]] <- wibw$Ki
    allmetrics[["between.module.deg"]] <- wibw$Bi
    allmetrics[["part.coeff"]] <- part_coeff(origgraph, allmetrics$community.membership)#, use.parallel = FALSE)
     allmetrics[["efficiency"]] <- local_efficiency(origgraph,directed=F)
    allmetrics[["betweenness.centrality"]] <- betweenness(origgraph,directed=F)
    
  }
  return(allmetrics)
}


#adapted function from igraph
wibw_module_degree <- function (g, community_attr="community") {
  stopifnot(is_igraph(g))
  if (is.null(get.vertex.attribute(g, community_attr))) { 
    warning("Cannot find community vertex attribute: ", community_attr) 
    return(NULL)
  }
  memb <- get.vertex.attribute(g, community_attr)
  modules <- unique(memb)
  N <- length(modules)
  A <- as_adj(g, sparse = FALSE, names = FALSE)
  z_within <- Ki <- rep(0, nrow(A))
  z_between <- Bi <- rep(0, nrow(A))
  Ksi <- sigKsi <- rep(0, nrow(A))
  Bsi <- sigBsi <- rep(0, nrow(A))
  for (S in modules) {
    x_wi <- A[memb == S, ] %*% (memb == S) #from S to S
    x_bw <- A[memb == S, ] %*% (memb != S) #from S to other nodes
    Ki[memb == S] <- x_wi; Ksi[memb==S] <- mean(x_wi); sigKsi[memb==S] <- sd(x_wi)
    Bi[memb == S] <- x_bw; Bsi[memb==S] <- mean(x_bw); sigBsi[memb==S] <- sd(x_bw)
  }
  z_within <- (Ki - Ksi)/sigKsi
  z_within <- ifelse(!is.finite(z_within), 0, z_within)
  z_between <- (Bi - Bsi)/sigBsi
  z_between <- ifelse(!is.finite(z_between), 0, z_between)
  return(data.frame(module=memb, Ki=Ki, Ksi=Ksi, sigKsi=sigKsi, z_within=z_within, Bi=Bi, Bsi=Bsi, sigBsi=sigBsi, z_between=z_between))
}



calc_global_metric <- function(origgraph, community_attr="community", modular_weighted = NULL) {
  #edgelist <- data.frame(get.edgelist(origgraph)) #not used currently
  globmetrics <- list()
  
  globmetrics[["edge_density"]] <- edge_density(origgraph) #does not use edge weights
  globmetrics[["characteristic_path_length"]] <- mean_distance(origgraph) #does not use edge weights
  globmetrics[["clustering_coefficient"]] <- transitivity(origgraph, type = "global") #uses edge weights if the graph has an edge weight attribute
  #globmetrics[["small_worldness"]] <- (mean_distance(origgraph)/mean_distance(rewire(origgraph, with = keeping_degseq(loops = FALSE, niter = 1000))))/(transitivity(origgraph, type = "global")/transitivity(rewire(origgraph, with = keeping_degseq(loops = FALSE, niter = 1000)), type = "global")) 
  globmetrics[["efficiency"]] <- global_efficiency(origgraph,directed=F)
  if (!is.null(vertex_attr(origgraph, community_attr))) {
    if (is.null(modular_weighted)) {
      globmetrics[["modularity"]] <- modularity(origgraph, vertex_attr(origgraph, community_attr))
    } else {
      globmetrics[["modularity"]] <- modularity(origgraph, vertex_attr(origgraph, community_attr), weights = E(origgraph)$weight)
    }
  }
  
  return(globmetrics)
}


assign_communities <- function(allg, comm, attribute="community") {
  stopifnot(is.list(allg))
  stopifnot("communities" %in% class(comm)) #ensure that it's a community object from igraph
  if (is_igraph(allg[[1]])) {
    #assume that we have a list of subjects where each element is a graph (currently used for weighted graphs)
    allg <- lapply(allg, function(subj) {
      subj <- set_vertex_attr(subj, attribute, value=comm$membership)
      return(subj)
    })
  } else if (is.list(allg[[1]]) && is_igraph(allg[[1]][[1]])) {
    #assuming we have a nested list structure: [subjects] x [density thresholds]
    
    allg <- lapply(allg, function(subj) {
      subj_graphs <- lapply(1:length(subj), function(d) {
        subj[[d]] <- set_vertex_attr(subj[[d]], attribute, value=comm$membership)
        return(subj[[d]])
      })
    })
  }
} 

part_coeff <- function(g, memb, A=NULL) {
  stopifnot(is_igraph(g))
  Ki <- check_degree(g)
  N <- max(memb)
  if (is.null(A)) A <- as_adj(g, sparse=FALSE, names=FALSE)
  Kis <- vapply(seq_len(N), function(x) colSums(A * (memb == x)), numeric(length(Ki)))
  PC <- 1 - ((1 / Ki^2) * rowSums(Kis^2))
  
  return(PC)
}

check_degree <- function(g) {
  if ('degree' %in% vertex_attr_names(g)) V(g)$degree else degree(g)
}


