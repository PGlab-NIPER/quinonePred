library("rcdk")
source("bfs.R")
smartcyp_smarts <- read.csv("smartcypenergyvalues.csv")
molecules <- load.molecules("quinone510_obmin.sdf",typing = T,aromaticity = T)
molecules <- lapply(molecules, remove.hydrogens)
adj_mats<- lapply(molecules, get.adjacency.matrix)
#function for calculating accessibility of atoms
span_value <- function(adjmat, atomindex = 0){
                  n <- nrow(adjmat)
                  topo_dist <- matrix(nrow = n,ncol = n)
                  span <- c()
                  if(atomindex >=0 & atomindex <= nrow(adjmat))
                    {
                      for (i in 1:n) 
                        {
                          topo_dist[i,] <- bfs(adjmat,i)
                        }
                      max_span <- max(topo_dist)
                      if(atomindex==0)
                        {
                          for (i in 1:n) 
                            {
                              span[i] <- max(topo_dist[i,])/max_span
                            }
                          return(span)
                        }
                      else
                        {
                          span <- max(topo_dist[atomindex,])/max_span
                          return(span)
                        }
                    }
                  else
                    {
                      print("enter valid atomindex value")
                    }
}


smart_des <- function(mols,adj_mat){
  span <- c()
  cyp_des <- matrix(0,ncol = 71,nrow = length(mols))
  for (i in 1:length(mols)) 
    {
      bits = get.fingerprint(molecule = mols[[i]],type = "substructure",substructure.pattern = smartcyp_smarts$SMARTS)@bits
      cyp_des[i,bits] <- 1
    }
  for (i in 1:length(mols)) {
    bits = which(cyp_des[i,]==1)
    
    for (j in bits) {
      matched_atoms <- c()
      matched_atoms_list <- matches(smartcyp_smarts$SMARTS[j],target = mols[[i]],return.matches = T)[[1]][[2]]
      for (k in 1:length(matched_atoms_list)) {
        matched_atoms <- c(matched_atoms,matched_atoms_list[[k]])
      }
      #start
      if(length(matched_atoms) == 1){
        cyp_des[i,j] <- smartcyp_smarts$E_value[j]*span_value(adjmat =  adj_mat[[i]],atomindex = (matched_atoms+1))
      }
      if(length(matched_atoms) > 1){
        for (k in 1:length(matched_atoms)) {
          span <- c(span,span_value(adjmat = adj_mat[[i]],atomindex = (matched_atoms[k]+1)))
        }
        cyp_des[i,j] <- smartcyp_smarts$E_value[j]*max(span)
      }
      if(length(matched_atoms) < 1 ){
        cyp_des[i,j] <- 0
      }
      
    }
  }
  return(cyp_des)
}


as.data.frame(smart_des(molecules,adj_mat = adj_mats)) -> smartcyp_des




