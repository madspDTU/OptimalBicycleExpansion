#### allPathsInParallel ####
allPathsInParallel <- function(theEdges, theOriPaths){
  allPaths <- foreach(i=1:length(theOriPaths)) %dopar% {
    return(as.vector(theEdges[unlist(theOriPaths[[i]])]));
  }
  return(allPaths);
}


#### allPathLengthsInParallel ####
allPathLengthsInParallel <- function(theOriPaths){
  allPathLengths <- foreach(i=1:length(theOriPaths)) %dopar% {
    return(as.vector(unlist(lapply(theOriPaths[[i]],length))));
  }
  return(allPathLengths);
}


#### allShortestPathsInParallel ####
allShortestPathsInParallel <- function(graph, theCentroidNodes) {
  allOriPaths <- foreach(i=1:dim(theCentroidNodes)[1]) %dopar% {
    fromNode <- igraph::V(graph)[theCentroidNodes$NodeId[i]];
    return(igraph::shortest_paths(graph, from=fromNode, to=igraph::V(graph)[theCentroidNodes$NodeId], 
                                  weight = igraph::E(graph)$weight, mode = "out", output = "epath")$epath);
  }
  return(allOriPaths)
}


#### assignFlows ####
assignFlows <- function(graph, CalculateCorrelations){
  
  start_time <- Sys.time()
  
  invisible(gc());
  
  E(graph)$Benefit <- 0
  E(graph)$Flow <- 0
  
  GeneralisedCostMatrix <- matrix(0,ncol = N_centroids, nrow = N_centroids);
  GeneralisedCostMatrix_Length <- matrix(0,ncol = N_centroids, nrow = N_centroids);
  
  pb <- txtProgressBar(min = 0, max = 100, style = 3, width = 33)
  if(CalculateCorrelations){
    uEffect <- array(0, c(N_centroids,N_centroids,N_segments)); #Three-dimensional array (O, D, Segment)
  }
  
  progr <- 5;
  if(CalculateCorrelations){
    progr <- progr/3; 
  }
  setTxtProgressBar(pb, progr)  
  
  
  #print("Finding all shortest paths")
  allOriPaths <- allShortestPathsInParallel(graph,CentroidNodes);
  progr <- 25;
  if(CalculateCorrelations){
    progr <- progr/3
  }  
  setTxtProgressBar(pb, progr)  
  
  #print("Finding all pathLengths")
  allPathLengths <- allPathLengthsInParallel(allOriPaths);
  progr <- 50;
  if(CalculateCorrelations){
    progr <- progr/3;
  }
  setTxtProgressBar(pb, progr)  
  
  
  #print("Finding all paths")
  allPaths <- allPathsInParallel(E(graph),allOriPaths);
  progr <- 75;
  if(CalculateCorrelations){
    progr <- progr/3;
  }
  setTxtProgressBar(pb, progr)  
  
  
  edgeLengths <- E(graph)$Length;
  edgeLinkIds <- E(graph)$LinkID;
  edgeOnlyNewIDs <- E(graph)$OnlyNewID;
  edgeSegments <- E(graph)$SegmentId;
  edgeWeights <- E(graph)$weight;
  
  for(i in 1:N_centroids){
    if(CalculateCorrelations){
      setTxtProgressBar(pb, 75/3 + 3*25*((i-1)/N_centroids))   
    } else {
      setTxtProgressBar(pb, 75 + 25*((i-1)/N_centroids))  
    }
    
    oriPaths <- allOriPaths[[i]]
    # Extracting the lengths of each path
    pathLengths <- allPathLengths[[i]];
    
    # Paths, one after another as a vector
    paths <- allPaths[[i]];
    
    totalAgg <- 0;
    totalDT <- 0;
    if(CalculateCorrelations){
      pathLengthCumSum <- c(0,cumsum(pathLengths));
      
      # NewIds of the edges of paths
      allTheLengths <- edgeLengths[paths];
      linkIDs <- edgeLinkIds[paths];
      segmentIDs <- edgeSegments[paths];
      for(j in 1:N_centroids){
        if(i != j){ #Ignoring zone-internal stuff
          k <- pathLengthCumSum[j]+1;
          endIndex <- pathLengthCumSum[j+1]; 
          allIds <- linkIDs[k:endIndex];
          
          usedSegmentIDs <- segmentIDs[k:endIndex];
          usedSegmentIndices <- usedSegmentIDs > 0
          usedSegmentIDs <- usedSegmentIDs[usedSegmentIndices];
          if(length(usedSegmentIDs) > 0){
            uniqueUsedSegmentIDs <- unique(usedSegmentIDs);
            N_usedSubProjectsThisCentroid <- length(uniqueUsedSegmentIDs);  
            if(N_usedSubProjectsThisCentroid == 1){
              proportions <- 1;
              affectingSegments <- uniqueUsedSegmentIDs;
            } else {
              segmentLengthPerSegment = allTheLengths[k:endIndex];
              segmentLengthPerSegment <- segmentLengthPerSegment[usedSegmentIndices];
              H3 <- aggregate(list(length=segmentLengthPerSegment), by = list(segment = usedSegmentIDs), FUN = sum);
              affectingSegments <- H3$segment
              proportions <- H3$length / sum(H3$length);
            }
            uEffect[i,j,affectingSegments] <- proportions;
          }
        }
      }
    }
    
    # Costs of the edges of paths
    pathGCs <- edgeWeights[paths];
    pathSPLengths <- edgeLengths[paths]
    # Corresponding destinations
    destinations <- rep(1:N_centroids, times = pathLengths)
    
    # Outputting the generalised costs and lenghts of the shortest path for each destination  
    H2 <- data.table(weight = pathGCs, Destination = destinations);
    H2 <- H2[, list(weight=sum(weight)), by=Destination]; #Equivalent to aggregate
    # And updating table of OD shortest path costs.
    GeneralisedCostMatrix[i,H2$Destination] <- H2$weight;
  }
  
  out <- list("GCs" = GeneralisedCostMatrix);
  if(CalculateCorrelations){
    out[["uEffect"]] <- uEffect;
  }
  
  setTxtProgressBar(pb, 100)  
  end_time <- Sys.time()
  print(end_time-start_time)
  return(out)
}


#### calculateTravelTimeSavings ####
calculateTravelTimeSavings <- function(selectedLinks){
  travelTimeSavings <- 0;
  for(travelerType in 1:N_TravelerTypes){
    OD <- allODs[[travelerType]];
    gcs_Basis <- all_gcs_Basis[[travelerType]];
    
    E(g)$weight                 <- allWeights[[paste0("Normal",travelerType)]];
    E(g)$weight[selectedLinks]  <- allWeights[[paste0("Upgraded",travelerType)]][selectedLinks];
    
    out <- assignFlows(g, FALSE);
    gcs_All <- out$GCs / 60;
    gcs_Dif <- gcs_Basis - gcs_All;
    
    travelTimeSavings <- travelTimeSavings + sum(OD*gcs_Dif);
  }
  return(travelTimeSavings);
}


#### createODMatrix ####
createODMatrix <- function(dat){
  out <- matrix(0,nrow=N_centroids, ncol=N_centroids);
  for( i in 1:dim(dat)[1]){
    val <- dat$value[i];
    if(!is.na(val)){
      fromIndex <- which(CentroidNodes$ZoneID == dat$FromZoneID[i])
      toIndex <- which(CentroidNodes$ZoneID == dat$ToZoneID[i])
      out[fromIndex, toIndex] <- val; 
    }
  }
  return(out);
}


#### determineMaximumSpeed ####
determineMaximumSpeed <- function(vMaxType, bikeType, speedType){
  if(bikeType == 1){
    if(speedType == 1){
      if(vMaxType == 1){
        return(13.6);
      } else if(vMaxType == 2){
        return(15.1);
      } else if(vMaxType == 3){
        return(16.6);
      }
    } else if(speedType == 2){
      if(vMaxType == 1){
        return(16.3)
      } else if(vMaxType == 2){
        return( 17.8) 
      } else if(vMaxType == 3){
        return( 19.3)
      }
    } else if(speedType == 3){
      if(vMaxType == 1){
        return( 19.1) 
      } else if(vMaxType == 2){
        return( 20.8) 
      } else if(vMaxType == 3){
        return( 22.5)
      }
    }
  } else if(bikeType == 2){
    if(speedType == 1){
      if(vMaxType == 1){
        return(15.6)
      } else if(vMaxType == 2){
        return( 17.1)
      } else if(vMaxType == 3){
        return( 18.6) 
      }
    } else if(speedType == 2){
      if(vMaxType == 1){
        return( 18.3) 
      } else if(vMaxType == 2){
        return( 19.8) 
      } else if(vMaxType == 3){
        return( 21.3) 
      }
    } else if(speedType == 3){
      if(vMaxType == 1){
        return( 21.1)
      } else if(vMaxType == 2){
        return( 22.8)
      } else if(vMaxType == 3){
        return( 24.5) 
      }
    }
  } else if (bikeType == 3){
    if(speedType == 1){
      if(vMaxType == 1){
        return( 22.6)
      } else if(vMaxType == 2){
        return( 24.1)
      } else if(vMaxType == 3){
        return( 25.6)
      }
    } else if(speedType == 2){
      if(vMaxType == 1){
        return( 25.3)
      } else if(vMaxType == 2){
        return( 26.8) 
      } else if(vMaxType == 3){
        return( 28.3)
      }
    } else if(speedType == 3){
      if(vMaxType == 1){
        return( 27.3)
      } else if(vMaxType == 2){
        return( 29.8)
      } else if(vMaxType == 3){
        return( 31.5)
      }
    }
  }
}

#### updateNetworkCosts ####
updateNetworkCosts <- function(g, infrastructureTypes, vMax_1, vMax_2, vMax_3){
  weights <- numeric(dim(LinksData)[1])+1000000; # A very high value
  weights <- ifelse(infrastructureTypes == 1, E(g)$Length / (vMax_1 / 3.6) + E(g)$IntersectionDelay, weights);
  weights <- ifelse(infrastructureTypes == 2, E(g)$Length / (vMax_2 / 3.6) + E(g)$IntersectionDelay, weights);
  weights <- ifelse(infrastructureTypes == 3, E(g)$Length / (vMax_3 / 3.6) + E(g)$IntersectionDelay, weights);
  return(weights);
}