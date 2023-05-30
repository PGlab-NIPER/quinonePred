Sys.setenv(JAVA_HOME='')
require("randomForest")
require("kernlab")
require("xgboost")
require("rcdk")
require("caret")
source('data/imports/bfs.R')
load("data/bin.RData")
xgb2d <- xgb.load("data/xgb2d.model")
xgbfp <- xgb.load("data/xgbfp.model")
xgbcbd <- xgb.load("data/xgbcbd.model")
smartcyp_smarts <- read.csv("data/imports/smartcypenergyvalues.csv")
smartcyp_smarts <- read.csv("data/imports/smartcypenergyvalues.csv")

infile <- file.choose()
molecules <- rcdk::load.molecules(infile,typing = T,aromaticity = T)
molecules <- lapply(molecules, rcdk::get.largest.component)
rcdk::write.molecules(molecules, 'data/temp/processed.sdf')

molecules <- lapply(molecules, rcdk::remove.hydrogens)
adj_mats<- lapply(molecules, rcdk::get.adjacency.matrix)

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
  cyp_des <- matrix(0,ncol = 9,nrow = length(mols))
  for (i in 1:length(mols)) 
  {
    bits =rcdk::get.fingerprint(molecule = mols[[i]],type = "substructure",substructure.pattern = smartcyp_smarts$SMARTS)@bits
    
    cyp_des[i,bits] <- 1
  }
  for (i in 1:length(mols)) {
    bits = which(cyp_des[i,]==1)
    
    for (j in bits) {
      matched_atoms <- c()
      matched_atoms_list <- rcdk::matches(smartcyp_smarts$SMARTS[j],target = mols[[i]],return.matches = T)[[1]][[2]]
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

smartcyp_des <- as.data.frame(smart_des(mols = molecules,adj_mat = adj_mats))
colnames(smartcyp_des) <- c("V36","V37","V39","V40","V42","V45","V51","V52","V53")

#padel descriptors
system(paste("java -jar data/padel/bin.jar -2d -fingerprints -removesalt -retainorder -detectaromaticity -standardizenitro -descriptortypes data/des_type.xml -dir data/temp/processed.sdf -file data/temp/out.csv"))
des1 <- read.csv('data/temp/out.csv')
des1 <- des1[,-1]
padel_des <- predict(prepro_models[[1]], des1)
padel_des <- padel_des[,feature_names[[1]]]
smartcyp_des <- predict(prepro_models[[2]], smartcyp_des)
padel_fp <- des1[,feature_names$fp]
pubchem_fp <- des1[,1691:2570]

cbd_des <- cbind.data.frame(padel_des,padel_fp,smartcyp_des)[feature_names[[3]]]

#predictions
pred1 <- predict( models$fit.rf,      padel_des , type = 'prob')[,2]
pred2 <- predict( models$fit.rf.fp,   padel_fp , type = 'prob')[,2]
pred3 <- predict( models$fit.rf.cbd,  cbd_des , type = 'prob')[,2]
pred4 <- predict( models$fit.svm,     padel_des , type = 'prob')[,2]
pred5 <- predict( models$fit.svm.fp,  padel_fp , type = 'prob')[,2]
pred6 <- predict( models$fit.svm.cbd, cbd_des , type = 'prob')[,2]
pred7 <- predict( models$fit.knn.cbd, cbd_des , type = 'prob')[,2]
pred8 <- 1-predict( xgb2d,  as.matrix(padel_des) )
pred9 <- 1-predict( xgbfp,  as.matrix(padel_fp) )
pred10 <- 1-predict(xgbcbd, as.matrix(cbd_des) )


probs <- cbind.data.frame(pred1, pred2, pred3, pred4, pred5, pred6, pred7, pred8, pred9, pred10)
ensemble_prob <- apply(probs, 1, mean)
ensemble_pred <- c()
ensemble_pred[ensemble_prob >= 0.5] <- "quinone"
ensemble_pred[ensemble_prob < 0.5] <- "nonquinone"
distance_to_model <- apply(probs, 1, sd)
domain <- ifelse(distance_to_model < 0.3, TRUE,FALSE)
smiles <-  unlist(lapply(molecules, rcdk::get.smiles))
molnames <-  unlist(lapply(molecules, rcdk::get.title))
results <- cbind.data.frame(Name = molnamesSMILES=smiles,Prediction=ensemble_pred, Probability = ensemble_prob, distance_to_model, Applicability_Domain = domain)
write.csv(results, file = "quinone_predictions.csv", row.names = F)


