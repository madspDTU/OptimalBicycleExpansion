#### Initialization ######

NUMBER_OF_ALLOCATED_CORES <- 4;

#install.packages("data.table")
#install.packages("igraph")
#install.packages("bit64")
#install.packages("stringr")
#install.packages("lpSolve")
#install.packages("profvis")
#install.packages("foreach")
#install.packages("doParallel")
library("data.table")
library("igraph")
library("bit64")
library("stringr")
library("lpSolve")
#library("profvis") #For profiling performance
library("foreach")
library("doParallel")
registerDoParallel(cores=NUMBER_OF_ALLOCATED_CORES);

#Adjust output location here:
outputDir <- "M:/GeoCBA/ModelOutputs"

#Adjust network data location here (must include the files: Connectors.csv, Centroids.csv, Nodes.csv and Links.csv):
NetworkDataDirectory <- "M:/GeoCBA/NetworkTables"

`%notin%` <- Negate(`%in%`)

USE_PREDETERMINED_SEGMENTS <- TRUE;
DIVIDE_BENEFITS_BY_NUMBER_OF_SEGMENTS <- FALSE;
DIVIDE_BENEFITS_PROPORTIONAL_TO_THEIR_INCLUDED_LENGTH <- TRUE;

### Functions #####

allShortestPathsInParallel <- function(graph, theCentroidNodes) {
  allOriPaths <- foreach(i=1:dim(theCentroidNodes)[1]) %dopar% {
    fromNode <- igraph::V(graph)[theCentroidNodes$NodeId[i]];
    return(igraph::shortest_paths(graph, from=fromNode, to=igraph::V(graph)[theCentroidNodes$NodeId], ##TODO Figure out how to skip 0-demand relations
                                  weight = igraph::E(graph)$weight, mode = "out", output = "epath")$epath);
  }
  return(allOriPaths)
}

allPathLengthsInParallel <- function(theOriPaths){
  allPathLengths <- foreach(i=1:length(theOriPaths)) %dopar% {
    return(as.vector(unlist(lapply(theOriPaths[[i]],length))));
  }
  return(allPathLengths);
}

allPathsInParallel <- function(theEdges, theOriPaths){
  allPaths <- foreach(i=1:length(theOriPaths)) %dopar% {
    return(as.vector(theEdges[unlist(theOriPaths[[i]])]));
  }
  return(allPaths);
}

assignFlows <- function(graph, ODDemand, Benefits, CalculateCorrelations, createUniqueSets){
  start_time <- Sys.time()
  
  invisible(gc());
  ODBenefits = ODDemand*Benefits;
  
  E(graph)$Benefit <- 0
  E(graph)$Flow <- 0
  
  GeneralisedCostMatrix <- matrix(0,ncol = N_centroids, nrow = N_centroids);
  GeneralisedCostMatrix_Length <- matrix(0,ncol = N_centroids, nrow = N_centroids);
  
  pb <- txtProgressBar(min = 0, max = 100, style = 3, width = 33)
  if(CalculateCorrelations){
    if(USE_PREDETERMINED_SEGMENTS){
      internalSegmentCorrs <- matrix(0,ncol = N_segments, nrow = N_segments);
      naiveBenefitVec <- numeric(N_segments);
    } else {
      internalCorrStructures <- matrix(0,ncol = N_newLinks, nrow = N_newLinks);
    }
  }
  progr <- 5;
  if(CalculateCorrelations){
    progr <- progr/3; 
  }
  setTxtProgressBar(pb, progr)  
  
  
  if(createUniqueSets){
    print("Initilising stuff")
    internalSubProjects <- data.table(matrix(0,ncol=(N_newLinks+1), nrow=N_centroids*N_centroids));
    subProjectId <- 1;
  }
  
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
  
  edgeBenefits <- E(graph)$Benefit;
  edgeFlows <- E(graph)$Flow;
  edgeLengths <- E(graph)$Length;
  edgeLinkIds <- E(graph)$LinkID;
  edgeOnlyNewIDs <- E(graph)$OnlyNewID;
  edgeSegments <- E(graph)$SegmentId;
  edgeWeights <- E(graph)$weight;
  
  
  #print("Doing all other stuff")
  for(i in 1:N_centroids){
    if(CalculateCorrelations){
      setTxtProgressBar(pb, 75/3 + 3*25*((i-1)/N_centroids))   
    } else {
      setTxtProgressBar(pb, 75 + 25*((i-1)/N_centroids))  
    }
    oriPaths <- allOriPaths[[i]]
    #    fromNode = V(graph)[CentroidNodes$NodeId[i]];
    #  oriPaths <- shortest_paths(graph, from=fromNode, to=V(graph)[CentroidNodes$NodeId], ##TODO Figure out how to skip 0-demand relations
    #                             weight = E(graph)$weight, mode = "out", output = "epath")$epath;
    # Extracting the lengths of each path
    #pathLengths <- as.vector(unlist(lapply(oriPaths,length)));
    pathLengths <- allPathLengths[[i]];
    
    # Paths, one after another as a vector
    #paths <- as.vector(E(graph)[unlist(oriPaths)])
    paths <- allPaths[[i]];
    
    totalAgg <- 0;
    totalDT <- 0;
    if(CalculateCorrelations){
      pathLengthCumSum <- c(0,cumsum(pathLengths));
      
      # Demand for the OD-pairs
      theseBenefits <- ODBenefits[i,];
      # NewIds of the edges of paths
      onlyNewIDs <- edgeOnlyNewIDs[paths];
      allTheLengths <- edgeLengths[paths];
      linkIDs <- edgeLinkIds[paths];
      segmentIDs <- edgeSegments[paths];
      for(j in 1:N_centroids){
        if(i != j){
          k <- pathLengthCumSum[j]+1;
          endIndex <- pathLengthCumSum[j+1]; 
          newIds <- onlyNewIDs[k:endIndex];
          newIds <- newIds[newIds > 0];
          allIds <- linkIDs[k:endIndex];
          
          if(USE_PREDETERMINED_SEGMENTS){
            usedSegmentIDs <- segmentIDs[k:endIndex];
            usedSegmentIndices <- usedSegmentIDs > 0
            usedSegmentIDs <- usedSegmentIDs[usedSegmentIndices];
            if(length(usedSegmentIDs) > 0){
              uniqueUsedSegmentIDs <- unique(usedSegmentIDs);
              N_usedSubProjectsThisCentroid <- length(uniqueUsedSegmentIDs);  
              if(N_usedSubProjectsThisCentroid == 1 | DIVIDE_BENEFITS_BY_NUMBER_OF_SEGMENTS){
                internalSegmentCorrs[uniqueUsedSegmentIDs,uniqueUsedSegmentIDs] <-  
                  internalSegmentCorrs[uniqueUsedSegmentIDs,uniqueUsedSegmentIDs] + 
                  theseBenefits[j] / N_usedSubProjectsThisCentroid;
              } else if(DIVIDE_BENEFITS_PROPORTIONAL_TO_THEIR_INCLUDED_LENGTH){
                segmentLengthPerSegment = allTheLengths[k:endIndex];
                segmentLengthPerSegment <- segmentLengthPerSegment[usedSegmentIndices];
                st5 <- Sys.time();
                H3 <- aggregate(list(length=segmentLengthPerSegment), by = list(segment = usedSegmentIDs), FUN = sum);
                et5 <- Sys.time();
                proportions <- H3$length / sum(H3$length);
                perSegmentBenefit <- theseBenefits[j] * proportions;
                addOnMatrix <- matrix(0,nrow =  N_usedSubProjectsThisCentroid, ncol = N_usedSubProjectsThisCentroid); 
                for(m in 1:N_usedSubProjectsThisCentroid){
                  addOn <- perSegmentBenefit[m] * proportions;
                  #addOn[addOn > addOn[m]] <- addOn[m];   #Never larger than diagonal element 
                  addOnMatrix[m,] <- addOn; 
                }
                matrixSegments <- H3$segment
                internalSegmentCorrs[matrixSegments,matrixSegments] <- 
                  internalSegmentCorrs[matrixSegments,matrixSegments] + addOnMatrix;
              }
              # vector thing
              naiveBenefitVec [uniqueUsedSegmentIDs] <- naiveBenefitVec [uniqueUsedSegmentIDs] + theseBenefits[j]
            }
          } else {
            if(length(newIds) > 0){
              internalCorrStructures[newIds,newIds] <- internalCorrStructures[newIds, newIds] + theseBenefits[j];
              if(createUniqueSets){
                internalSubProjects[subProjectId, 1:= theseBenefits[j]];
                internalSubProjects[subProjectId, (newIds+1) := 1];
                subProjectId <- subProjectId + 1;
              }
            }
          }
        }
        #k <- endIndex + 1;
        if(createUniqueSets){
          print(paste(i,j,subProjectId));
        } 
      }
    }
    
    # Costs of the edges of paths
    pathGCs <- edgeWeights[paths];
    pathSPLengths <- edgeLengths[paths]
    # Corresponding destinations
    destinations <- rep(1:N_centroids, times = pathLengths)
    
    # Benefits for the OD-pairs (not normalised)
    theseBenefits <- ODBenefits[i,destinations];
    
    # Demand for the OD-pairs
    theseDemands <- ODDemand[i,destinations];
    # Binding paths and benefits
    H <- data.table(Benefit = theseBenefits, Flow = theseDemands, Edge = paths);
    H[,list(Benefit=sum(Benefit), Flow =sum(Flow)),by=Edge]
    # Assigning benefits to each of the edges
    edgeBenefits[H$Edge] <- edgeBenefits[H$Edge] + H$Benefit;
    # Assigning flows to each of the edges
    edgeFlows[H$Edge] <- edgeFlows[H$Edge] + H$Flow;
    # Outputting the generalised costs of shortest path for each destination  
    H2 <- data.table(weight = pathGCs, length=pathSPLengths,  Destination = destinations);
    H2 <- H2[, list(weight=sum(weight), length=sum(length)), by=Destination]; #Equivalent to aggregate
    # And updating table of OD shortest path costs.
    GeneralisedCostMatrix[i,H2$Destination] <- H2$weight;
    GeneralisedCostMatrix_Length[i,H2$Destination] <- H2$length;
  }
  
  E(graph)$Benefit <- edgeBenefits;
  E(graph)$Flow <- edgeFlows;
  
  #print("Creating output list")
  
  out <- list("Graph" = graph, "GCs" = GeneralisedCostMatrix, ODLengths = GeneralisedCostMatrix_Length);
  if(CalculateCorrelations){
    if(createUniqueSets){
      out[["subProjects"]] <- internalSubProjects;
    } 
    if(USE_PREDETERMINED_SEGMENTS){
      out[["segmentCorrs"]] <- internalSegmentCorrs;
      out[["naiveBenefitVector"]] <- naiveBenefitVec ;
      print(sum(internalSegmentCorrs))
      print(sum(naiveBenefitVec ))
    } else {
      out[["corrStructures"]] <- internalCorrStructures;
    }
  }
  
  setTxtProgressBar(pb, 100)  
  
  end_time <- Sys.time()
  print(end_time-start_time)
  return(out)
}

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

determineVMax <- function(vMaxType, bikeType, speedType){
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



#131


CONSTRUCTION_COST_PER_KM <- 2.6   ## Fra Hallberg (2021) Bør måske opjusteres
MAINTENANCE_COST_PER_KM_PER_YEAR <- 0.2; # Fra Hallberg (2021) Bør måske opjusteres. 
EVALUATION_PERIOD <- 50;
CYCLIST_UNIT_PRICE_PER_MINUTE  <- 91 / 60;
AADT_FACTOR <- 250;
NAF <- 1.325;
r <- 1.04;
r_2 <- 1.03;
END_OF_FIRST_RATE <- 35;
DISCOUNT_FACTORS <- c(r^(1:END_OF_FIRST_RATE),r^END_OF_FIRST_RATE*r_2^(1:(EVALUATION_PERIOD-END_OF_FIRST_RATE + 1))); #+1 for scrap year
cumDiscFactors <- cumsum(1/DISCOUNT_FACTORS);

CS_PHASEIN <- c((1:10)/10,rep(1,40));

CONSTRUCTION_TIME <- 1; ## [years]

aspSegments = c(3,13,30,50,53,56,63,84,97,108,111,129,132,138,147,157,179,189,195);




discount <- function(totalAmount, discountFactor, years, delay, phaseIn = rep(1,length(years))){
  if(totalAmount == 0){
    return(0);
  } else if(delay == 0 && years == 1) {
    return(totalAmount); 
  } else {
    return(sum(rep(totalAmount / years, years) / discountFactor[delay + 1:years] * phaseIn));
  }
}


## Vores længder er off. Lav dem bedre med længste korteste vej i hvert segment. (alternativt se på priser for hver rute)
## Arbejd på diskontering. (Dette skal gøres.)
## En anden ting: Udvid algoritmen, så vi genberegner benefit-matricen efter hver epoke. 

## Reel længde er 745.8km.
### 575 (nye) + 186.5 (eksisterende) - rute 26 (Køge banesti) 7.14 = 754. Tæt nok på...
### Vores omkorstning 4'072.5 millioner.
### Incentive's: 6'800 millioner, men heraf er en del eksisterende.  Svarende til ca. 5.1 mia. 





### Import network ######

LinksData <- as.data.frame(data.table::fread(paste0(NetworkDataDirectory, "/Links_WithSegments.csv"), dec=","))
NodesData <- as.data.frame(data.table::fread(paste0(NetworkDataDirectory, "/Nodes.csv"), dec=","))
CentroidNodes <- subset(NodesData, ZoneID > -1)
ConnectorsData <- as.data.frame(data.table::fread(paste0(NetworkDataDirectory, "/Connectors.csv")))
CentroidsData <- as.data.frame(data.table::fread(paste0(NetworkDataDirectory, "/Centroids.csv")))
N_segments <- max(LinksData$SegId);


LinksData$SeeminglyUnused <- ifelse(LinksData$Active_Basis == 0 & LinksData$Active_Kommende == 0 & 
                                      LinksData$Active_Plan_18_30 == 0 & LinksData$Active_Plan_30_45 == 0,1, 0)
unimplementedLinks <- LinksData$ID[LinksData$SeeminglyUnused == 1];
unimplementedLinks <- as.data.frame(unimplementedLinks);
colnames(unimplementedLinks) <- "LinkID";
data.table::fwrite(unimplementedLinks, file = paste0(outputDir, "/unimplementedLinks.csv"))

LinksData <- subset(LinksData, SeeminglyUnused == 0);
#LinksData$Active_Plan_30_45[LinksData$SeeminglyUnused == 1] <- 1;

LinksData$InfraType_Basis[LinksData$Active_Basis == 0] <- -1;
LinksData$InfraType_Basis[LinksData$Active_Basis == 1 & LinksData$Super_Basis == 1] <- 3
LinksData$InfraType_Basis[LinksData$Active_Basis == 1 & LinksData$Super_Basis == 0 & LinksData$FreeSpeed >= 22] <- 2
LinksData$InfraType_Basis[LinksData$Active_Basis == 1 & LinksData$Super_Basis == 0 & LinksData$FreeSpeed < 22] <- 1

LinksData$InfraType_Kommende[LinksData$Active_Kommende == 0] <- -1;
LinksData$InfraType_Kommende[LinksData$Active_Kommende == 1 & LinksData$Super_Kommende == 1] <- 3
LinksData$InfraType_Kommende[LinksData$Active_Kommende == 1 & LinksData$Super_Kommende == 0 & LinksData$FreeSpeed >= 22] <- 2
LinksData$InfraType_Kommende[LinksData$Active_Kommende == 1 & LinksData$Super_Kommende == 0 & LinksData$FreeSpeed < 22] <- 1

LinksData$InfraType_Plan_18_30[LinksData$Active_Plan_18_30 == 0] <- -1;
LinksData$InfraType_Plan_18_30[LinksData$Active_Plan_18_30 == 1 & LinksData$Super_Plan_18_30 == 1] <- 3
LinksData$InfraType_Plan_18_30[LinksData$Active_Plan_18_30 == 1 & LinksData$Super_Plan_18_30 == 0 & LinksData$FreeSpeed >= 22] <- 2
LinksData$InfraType_Plan_18_30[LinksData$Active_Plan_18_30 == 1 & LinksData$Super_Plan_18_30 == 0 & LinksData$FreeSpeed < 22] <- 1

LinksData$InfraType_Plan_30_45[LinksData$Active_Plan_30_45 == 0] <- -1;
LinksData$InfraType_Plan_30_45[LinksData$Active_Plan_30_45 == 1 & LinksData$Super_Plan_30_45 == 1] <- 3
LinksData$InfraType_Plan_30_45[LinksData$Active_Plan_30_45 == 1 & LinksData$Super_Plan_30_45 == 0 & LinksData$FreeSpeed >= 22] <- 2
LinksData$InfraType_Plan_30_45[LinksData$Active_Plan_30_45 == 1 & LinksData$Super_Plan_30_45 == 0 & LinksData$FreeSpeed < 22] <- 1


LinksData <- subset(LinksData, select=c("FromNodeId","ToNodeId","FreeSpeed","Shape_Length","IntersectionDelayFor",
                                        "InfraType_Basis","InfraType_Kommende","InfraType_Plan_18_30","InfraType_Plan_30_45","ID","SegId","RouteId"))
colnames(LinksData) <- c("from","to","FreeSpeed","Length","IntersectionDelay",
                         "InfraType_Basis","InfraType_Kommende","InfraType_Plan_18_30","InfraType_Plan_30_45","LinkID","SegmentId","RouteId")
NodesData <- subset(NodesData, select=c("NodeId","X","Y"), by.x = "FromNodeId", by.y = "NodeId");
colnames(NodesData) <- c("name","X","Y")

disappearingLinks <- LinksData$LinkID[LinksData$InfraType_Plan_30_45==-1 & LinksData$InfraType_Basis>0];



N_centroids <- dim(CentroidNodes)[1];
N_edges <- dim(LinksData)[1];

existingSuper <- LinksData$LinkID[LinksData$InfraType_Basis == 3];
existingSuper <- as.data.frame(existingSuper);
colnames(existingSuper) <- "LinkID";
data.table::fwrite(existingSuper, file = paste0(outputDir, "/ExistingSuper.csv"))




### Create correlation structures #####

LinksData$OnlyNewID <- 0;
news <- which(LinksData$InfraType_Basis == -1);
newId <- 0;
for( i in news){
  newId <- newId +1;
  LinksData$OnlyNewID[i] <- newId;
}
N_newLinks <- newId;





# Creates the graph "g"
g <- graph_from_data_frame(LinksData, directed = TRUE, vertices = NodesData)
components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
g <- igraph::induced_subgraph(g, vert_ids);
#g <- delete_vertex_attr(g,"name")   # Deleting the vertex atribute "name". 

N_edges = length(E(g)$OnlyNewID)
linkBenefits = numeric(N_edges);
onlyNewLinksBool <- E(g)$OnlyNewID > 0;

CentroidNodes$Name <- CentroidNodes$NodeId
for(i in 1:dim(CentroidNodes)[1]){
  CentroidNodes$NodeId[i]  <- which(V(g)$name == CentroidNodes$Name[i]) 
}
g <- delete_vertex_attr(g,"name")   # Deleting the vertex atribute "name". 
g <- delete_vertex_attr(g,"X")   # Deleting the vertex atribute "name". 
g <- delete_vertex_attr(g,"Y")   # Deleting the vertex atribute "name". 
g <- delete_edge_attr(g,"OnlyNewID")   # Deleting the vertex atribute "name".
g <- delete_edge_attr(g,"InfraType_Kommende")   # Deleting the vertex atribute "name".
g <- delete_edge_attr(g,"InfraType_Plan_18_30")   # Deleting the vertex atribute "name".




if(FALSE){
  for(whichBase in c("Base","RSP")){
    #for(whichBase in c("RSP")){
    if(whichBase == "RSP"){
      rspSuffix = "_RSP";
      randomStartingPoint <- c(53,108,84,3,63,30,138,179,97,147,132,204,129,198,50,111,56,13,
                               157,167,42,133,26,58,193,54,12,115,165,200,103,98,174)
      freeLinks <- E(g)$LinkID[E(g)$SegmentId %in% randomStartingPoint];
    } else {
      rspSuffix = "";
      freeLinks <- c();
    }
    #profvis({
    bikeTypes <- c(1,2,3);
    speedTypes <- c(1,2,3);
    corrStructures <- matrix(0,ncol=N_newLinks, nrow=N_newLinks);
    segmentCorrs <- matrix(0,ncol=N_segments, nrow=N_segments);
    globalNaiveBenefitVector <- numeric(N_segments);
    totalOD <- matrix(0,nrow=length(CentroidNodes$NodeId), ncol = length(CentroidNodes$NodeId));
    totalOD_B <- matrix(0,nrow=length(CentroidNodes$NodeId), ncol = length(CentroidNodes$NodeId));
    totalOD_B_OD <- matrix(0,nrow=length(CentroidNodes$NodeId), ncol = length(CentroidNodes$NodeId));
    ODData <- as.data.frame(data.table::fread("O:/Public/4233-82676-BIKELONGER-persondata/CBA/OD_BikeType_SpeedType/SC01_MixBase1_MixSc1.csv"));
    ODData <- subset(ODData, !is.na(ODData$value));
    for(bikeType in bikeTypes){
      for(speedType in speedTypes){
        
        OD <- createODMatrix(subset(ODData, BikeType==bikeType & SpeedType == speedType));
        totalOD <- totalOD + OD;
        
        ### Creating benefits ####
        print(c("Biketype, speedtype"))
        print(paste(bikeType,speedType))   
        
        # Initialization
        vMax_1 <- determineVMax(1,bikeType,speedType); # Maximum speed of cyclists used for the assignment (Shared)
        vMax_2 <- determineVMax(2,bikeType,speedType); # Maximum speed of cyclists used for the assignment (Separated)
        vMax_3 <- determineVMax(3,bikeType,speedType); # Maximum speed of cyclists used for the assignment (Super)
        
        ZeroBenefits = matrix(0,nrow = N_centroids, ncol = N_centroids);
        
        print("Running basis scenario")
        theType <- E(g)$InfraType_Basis
        E(g)$weight <- 1000000;
        E(g)$weight <- ifelse(theType == 1, E(g)$Length / (vMax_1 / 3.6) + E(g)$IntersectionDelay, E(g)$weight);
        E(g)$weight <- ifelse(theType == 2, E(g)$Length / (vMax_2 / 3.6) + E(g)$IntersectionDelay, E(g)$weight);
        E(g)$weight <- ifelse(theType == 3, E(g)$Length / (vMax_3 / 3.6) + E(g)$IntersectionDelay, E(g)$weight);
        if(whichBase == "RSP"){
          E(g)$weight <- ifelse(E(g)$LinkID %in% freeLinks, E(g)$Length / (vMax_3 / 3.6) + E(g)$IntersectionDelay, E(g)$weight);
        }
        out <- assignFlows(g, OD, ZeroBenefits, FALSE, FALSE);
        gcs_Basis <- out$GCs / 60;
        ODLengths_Basis <- out$ODLengths;
        g_Basis <- out$Graph
        
        print(mean(gcs_Basis)) 
        saveRDS(gcs_Basis, file=paste0(outputDir,"/gcs_Basis_",bikeType,"_",speedType,rspSuffix,".Rda"))    
        
        print("Running _all_ scenario")
        theType <- E(g)$InfraType_Plan_30_45
        E(g)$weight <- 1000000;
        E(g)$weight <- ifelse(theType == 1, E(g)$Length / (vMax_1 / 3.6) + E(g)$IntersectionDelay, E(g)$weight);
        E(g)$weight <- ifelse(theType == 2, E(g)$Length / (vMax_2 / 3.6) + E(g)$IntersectionDelay, E(g)$weight);
        E(g)$weight <- ifelse(theType == 3, E(g)$Length / (vMax_3 / 3.6) + E(g)$IntersectionDelay, E(g)$weight);
        
        out <- assignFlows(g, OD, ZeroBenefits, FALSE, FALSE);
        
        gcs_All <- out$GCs / 60;
        ODLengths_All <- out$ODLengths;
        g_All <- out$Graph;
        
        
        ## Calculating differences
        OD_B <- gcs_Basis - gcs_All;
        totalOD_B <- totalOD_B + OD_B;
        
        totalOD_B_OD <- totalOD_B_OD + OD_B * OD;
        
        OD_LDif <- ODLengths_Basis - ODLengths_All;
        
        if(min(OD_B) < -1e-12){
          print(paste("Unfortunately, we have",sum(OD_B<0),"OD relations with negative benefits :("));
          print(paste ("For now, setting these to 0"))
          wobb <- which(OD_B<0, arr.ind = TRUE)
        }
        
        print("Assigning benefits (all-base)")
        out <- assignFlows(g, OD, OD_B, TRUE, FALSE);
        
        g_dif <- out$Graph;
        subProjects <- out$subProjects;
        
        linkBenefits <- linkBenefits + E(g_dif)$Benefit
        if(USE_PREDETERMINED_SEGMENTS){
          segmentCorrs <- segmentCorrs + out$segmentCorrs;
          globalNaiveBenefitVector <- globalNaiveBenefitVector + out$naiveBenefitVector;
          
        } else {
          corrStructures <- corrStructures + out$corrStructures;
        }
      }
    }
    print(discount(sum(totalOD_B_OD) * CYCLIST_UNIT_PRICE_PER_MINUTE * AADT_FACTOR * EVALUATION_PERIOD / 1e6
                   ,DISCOUNT_FACTORS,50,0,CS_PHASEIN));
    print(sum(globalNaiveBenefitVector)/sum(totalOD_B_OD))
    
    
    if(USE_PREDETERMINED_SEGMENTS){
      saveRDS(globalNaiveBenefitVector, file=paste0(outputDir,"/globalNaiveBenefitVector_Parallel",rspSuffix,".Rda"))
      if(DIVIDE_BENEFITS_PROPORTIONAL_TO_THEIR_INCLUDED_LENGTH){
        saveRDS(segmentCorrs, file=paste0(outputDir,"/segmentCorrs_Parallel_Per_Length",rspSuffix,".Rda"))
      } else if(DIVIDE_BENEFITS_BY_NUMBER_OF_SEGMENTS){
        saveRDS(segmentCorrs, file=paste0(outputDir,"/segmentCorrs_Parallel",rspSuffix,".Rda"))      
      }
    } else {
      saveRDS(corrStructures, file=paste0(outputDir,"/corrStructures_Parallel",rspSuffix,".Rda"))
    }
    #})
  }
}

if(USE_PREDETERMINED_SEGMENTS){
  globalNaiveBenefitVector <- readRDS(file=paste0(outputDir,"/globalNaiveBenefitVector_Parallel.Rda"))
  globalNaiveBenefitVector_RSP <- readRDS(file=paste0(outputDir,"/globalNaiveBenefitVector_Parallel_RSP.Rda"))
  
  if(DIVIDE_BENEFITS_PROPORTIONAL_TO_THEIR_INCLUDED_LENGTH){
    segmentCorrs <- readRDS(file=paste0(outputDir,"/segmentCorrs_Parallel_Per_Length.Rda"))
    segmentCorrs_RSP <- readRDS(file=paste0(outputDir,"/segmentCorrs_Parallel_Per_Length_RSP.Rda"))
  } else if(DIVIDE_BENEFITS_BY_NUMBER_OF_SEGMENTS){
    segmentCorrsOld <- readRDS(file=paste0(outputDir,"/segmentCorrs_Parallel.Rda"))   
    segmentCorrsOld_RSP <- readRDS(file=paste0(outputDir,"/segmentCorrs_Parallel_RSP.Rda"))   
  }
} else {
  corrStructures <- readRDS(paste0(outputDir,"/corrStructures.Rda"))
  corrStructures_RSP <- readRDS(paste0(outputDir,"/corrStructures_RSP.Rda"))
}






if(FALSE){
  lengthDist <- numeric(60);
  for(i in 1:length(lengthDist)){
    lengthDist[i] <- sum(OD[ODLengths_Basis>=(i-1)*1000 & ODLengths_Basis<i*1000] / sum(OD))
  }
  plot(1:length(lengthDist),cumsum(lengthDist))
  #barplot(lengthDist, names.arg = 1:length(lengthDist), las = 2)
  
  lengthDist_All <- numeric(60);
  for(i in 1:length(lengthDist_All)){
    lengthDist_All[i] <- sum(OD[ODLengths_All>=(i-1)*1000 & ODLengths_All<i*1000] / sum(OD))
  }
  points(1:length(lengthDist_All),cumsum(lengthDist_All), col = 2)
  #barplot(lengthDist_All, names.arg = 1:length(lengthDist), las = 2)
  
  
  crowDist_All <- numeric(300);
  for(i in 1:dim(OD)[1]){
    for( j in 1:dim(OD)[2]){
      if(i != j & OD[i,j]>0){
        crowDist <- sqrt((CentroidNodes$X[i]-CentroidNodes$X[j])^2 + (CentroidNodes$Y[i]-CentroidNodes$Y[i])^2) / 1000;
        ceiledCrowDist <- ceiling(crowDist);
        if(ceiledCrowDist < length(crowDist_All)){
          crowDist_All[ceiledCrowDist] <- crowDist_All[ceiledCrowDist] +  OD[i,j];
        } else {
          print("What??");
        }
      }
    }
  }
  crowDist_All <- crowDist_All / (sum(OD)-sum(diag(OD)));
  points(1:length(crowDist_All),cumsum(crowDist_All), col = 4)
  # plot(1:length(crowDist_All),cumsum(crowDist_All), col = 4, xlim=c(0,20))
}



lengths <- E(g)$Length / 1000;
costs <- lengths * (CONSTRUCTION_COST_PER_KM + EVALUATION_PERIOD * MAINTENANCE_COST_PER_KM_PER_YEAR);
routeIds <- numeric(N_segments) - 1



if(!USE_PREDETERMINED_SEGMENTS){
  # Calculating correlations
  unusedNewIds <- which(diag(corrStructures) == 0);
  print(paste(length(unusedNewIds),"unused new links"));
  usedNewIds <- 1:N_newLinks;
  if(length(unusedNewIds) > 0){
    usedNewIds <- usedNewIds[usedNewIds %notin% unusedNewIds];
  }
  
  usedCorrStructures <- corrStructures[usedNewIds, usedNewIds];
  usedCorr1d <- diag(corrStructures)[usedNewIds];
  newLinkBenefits <- linkBenefits[usedNewIds];
  newCosts <- costs[usedNewIds]
  newLengths <- lengths[usedNewIds]
} else {
  coreCosts <- numeric(N_segments);
  newLengths <- coreCosts;
  for(i in 1:N_segments){
    segmentSelector <- E(g)$SegmentId == i;
    eids <- E(g)[segmentSelector];
    gsub <- subgraph.edges(g,eids=eids)
    if(length(E(gsub))>0){
      newLengths[i] <- farthest_vertices(gsub, directed=FALSE, weights=E(gsub)$Length)$distance/1000;
      coreCosts[i] <- discount(newLengths[i] * CONSTRUCTION_COST_PER_KM, r, CONSTRUCTION_TIME,0) +
        discount(newLengths[i] * MAINTENANCE_COST_PER_KM_PER_YEAR *  (EVALUATION_PERIOD - CONSTRUCTION_TIME), 
                 r, EVALUATION_PERIOD - CONSTRUCTION_TIME,CONSTRUCTION_TIME);
      routeIds[i] <-  unique(E(g)$RouteId[segmentSelector])[1];
    }
    print(paste("segment",i,newLengths[i]))
    print(components(gsub)$csize)
  }
}


sum(newLengths);
sum(coreCosts);


IncentiveDataDirectory <- "M:/GeoCBA/Incentive"
IncentiveData <- as.data.frame(data.table::fread(paste0(IncentiveDataDirectory, "/KeyFigures_Incentive.csv"), dec="."))
IncentiveData$ConstructionCostPerKm <- IncentiveData$ConstructionCost/IncentiveData$Length;
IncentiveData$MaintenanceCostPerKm <- IncentiveData$MaintenanceCost/IncentiveData$Length;


statusDF <- as.data.frame(cbind(1:N_segments, routeIds, newLengths))
colnames(statusDF) <- c("SegmentId", "RouteId","SegmentLength")
routeStatusDF <- aggregate(list(AppliedLength = statusDF$SegmentLength), by = list(RouteId=statusDF$RouteId), FUN=sum)
statusDF <- merge(statusDF,IncentiveData, by = "RouteId", all.x = TRUE);
statusDF$CorrectedConstructionCosts <- statusDF$SegmentLength * statusDF$ConstructionCostPerKm;
statusDF$CorrectedMaintenanceCosts <- statusDF$SegmentLength * statusDF$MaintenanceCostPerKm;
routeStatusDF <- merge(routeStatusDF,IncentiveData, by = "RouteId", all.x = TRUE);
#View(statusDF)
#View(routeStatusDF)
#plot(routeStatusDF$Length,routeStatusDF$AppliedLength)
#outlierDetec <- abs(routeStatusDF$AppliedLength-routeStatusDF$Length) > 1
#points(routeStatusDF$Length[outlierDetec],routeStatusDF$AppliedLength[outlierDetec], col = 2)
#routeStatusDF$Name[outlierDetec]
#routeStatusDF[outlierDetec,]
routeCostsForGIS = statusDF[statusDF$RouteId >=0,c("RouteId","ConstructionCostPerKm","MaintenanceCostPerKm")];
routeCostsForGIS = unique(routeCostsForGIS);
routeCostsForGIS <- routeCostsForGIS[order(routeCostsForGIS$RouteId),];
data.table::fwrite(routeCostsForGIS, file = "M:/GeoCBA/RouteCostsForGIS.csv", 
                   col.names = TRUE, row.names = FALSE, sep = ";", dec = ",");

correctedConstructionCosts <- numeric(N_segments)
correctedConstructionCosts[statusDF$SegmentId] <- statusDF$CorrectedConstructionCosts
correctedMaintenanceCosts <- numeric(N_segments)
correctedMaintenanceCosts[statusDF$SegmentId] <- statusDF$CorrectedMaintenanceCosts;

correctedConstructionCosts[is.na(correctedConstructionCosts)] <- 
  statusDF$SegmentLength[is.na(correctedConstructionCosts)] * CONSTRUCTION_COST_PER_KM;
correctedMaintenanceCosts[is.na(correctedMaintenanceCosts)] <- 
  statusDF$SegmentLength[is.na(correctedMaintenanceCosts)] * MAINTENANCE_COST_PER_KM_PER_YEAR;

#coreCostSum <- 0;
#for(i in 1:length(correctedConstructionCosts)){
#  coreCostSum <- coreCostSum + discount(correctedConstructionCosts[i], r, CONSTRUCTION_TIME,0) +
#      discount(correctedMaintenanceCosts[i] *  (EVALUATION_PERIOD - CONSTRUCTION_TIME), 
#               r, EVALUATION_PERIOD - CONSTRUCTION_TIME,CONSTRUCTION_TIME);
#}
#print(coreCostSum)



#sigmas <- apply(usedCorrStructures, 2, function(x) x/usedCorr1d);

if(!USE_PREDETERMINED_SEGMENTS){
  N_usedNewIds <- length(usedNewIds);
  corrThress <- c(10,25,50,60,70,80,90,100);
  for(corrThres in corrThress){
    th <- usedCorrStructures / diag(usedCorrStructures)
    subProjects <- numeric(N_usedNewIds);
    nextSubProject <- 1;
    tol <- corrThres / 100
    th <- th>=tol;
    
    while(any(subProjects==0)){
      th2 <- rowsums(th)
      nv <- max(th2);
      print(paste(nv,nextSubProject-1,sum(th2>0)))
      hist(th2[th2>0], breaks=0:max(th2))
      if(nv == 1){
        nextOnes <- th2 == nv;
        subProjects[nextOnes] <- nextSubProject:(nextSubProject+sum(nextOnes)-1);
        nextSubProject <- nextSubProject + sum(nextOnes);
      } else {
        nextOne <- which(th2 == nv)[1];
        nextOnes <- th[nextOne,]>=tol;
        subProjects[nextOnes] <- nextSubProject;
        nextSubProject <- nextSubProject + 1;
      }
      th[nextOnes,] <- FALSE;
      th[,nextOnes] <- FALSE;
    }
    plot(ecdf(subProjects))
    if(corrThres == 100){
      subProjects100 <- subProjects
    } else if(corrThres == 90){
      subProjects90 <- subProjects;
    } else if(corrThres == 80){
      subProjects80 <- subProjects
    } else if(corrThres == 70){
      subProjects70 <- subProjects
    } else if(corrThres == 60){
      subProjects60 <- subProjects
    } else if(corrThres == 50){
      subProjects50 <- subProjects
    } else if(corrThres == 25){
      subProjects25 <- subProjects
    } else if(corrThres == 10){
      subProjects10 <- subProjects
    }
  }
  
  plot(ecdf(subProjects100))
  lines(ecdf(subProjects90))
  lines(ecdf(subProjects80))
  lines(ecdf(subProjects70))
  lines(ecdf(subProjects60))
  lines(ecdf(subProjects50))
  lines(ecdf(subProjects25))
  lines(ecdf(subProjects10))
  print(paste(max(subProjects100),max(subProjects90),max(subProjects80),max(subProjects70),
              max(subProjects60),max(subProjects50),max(subProjects25),max(subProjects10)))
  subProjects <- as.data.frame(cbind(subProjects100,subProjects90,subProjects80,subProjects70,
                                     subProjects60,subProjects50,subProjects25,subProjects10))
  colnames(subProjects) <- c("SubProject100","SubProject90","SubProject80","SubProject70",
                             "SubProject60","SubProject50","SubProject25","SubProject10")
  saveRDS(subProjects, file=paste0(outputDir,"/SubProjects.Rda"))
  subProjects <- readRDS(paste0(outputDir, "/SubProjects.Rda"));
} 





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## BLP setup #####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

solutionTypes <- c("LP","LPGreedyStop","LPGreedy","ByActualOrder","ByRandomSegment","ByShortestSegment","ByShortestRoute","ByLongestSegment","ByLongestRoute")

#What about "NettoAfgiftsFaktor", "Skatteforvidningseffekt", osv....

## Generating predetermined implementation orders
{
  set.seed(426992)
  coreSegments <- 1:N_segments;
  #"ByRandomSegment";
  randomSegmentOrder <- sample(coreSegments[newLengths > 0], replace = FALSE);
  #"ByLongestSegment";
  longestSegmentOrder <- coreSegments[newLengths > 0];
  longestSegmentOrder <- longestSegmentOrder[order(-newLengths[newLengths>0])]
  #"ByShortestSegment";
  shortestSegmentOrder <- rev(longestSegmentOrder);
  #"ByLongestRoute";
  sortedRouteIds <- routeStatusDF[routeStatusDF$RouteId > 0, ];
  sortedRouteIds <- sortedRouteIds[order(-sortedRouteIds$AppliedLength), "RouteId"];
  longestRouteOrder <- numeric(0);
  while(length(sortedRouteIds) > 0){
    thisRouteId <- sortedRouteIds[1];
    sortedSegmentIds <- statusDF[statusDF$RouteId == thisRouteId, ];
    sortedSegmentIds <- sortedSegmentIds[order(-sortedSegmentIds$SegmentLength), "SegmentId"];
    longestRouteOrder <- c(longestRouteOrder, sortedSegmentIds);
    sortedRouteIds <- sortedRouteIds[-1];
  }
  ## Including those that do not belong to any route
  longestRouteOrder <- c(longestRouteOrder, as.numeric(statusDF$SegmentId[statusDF$RouteId == -1 & statusDF$SegmentLength > 0])); 
  #"ByShortestRoute";
  shortestRouteOrder <- rev(longestRouteOrder);
  #"ByActual";
  actualOrder <- c(164, 165, 166, 154, 105, 54, 93, 103, 129, 131, 114, 112, 66, 65, 64, 37, 38, 39, 40,
                   28, 34, 35, 50, 51, 52, 53, 55, 57, 91, 88, 84, 82, 80, 175, 190, 49, 47, 43, 42, 144,
                   76, 122, 127, 128, 102, 101, 176, 171, 170, 70, 69, 68, 67, 111, 151)
  actualRealisedOrder <- c(164, 165, 166, 154, 105, 54, 93, 103, 129, 131, 114, 112, 66, 65, 64, 37, 38, 39, 40,
                           28, 34, 35, 50, 51, 52, 53, 55, 57, 91, 88, 84, 82, 80, 175, 190)
}
set.seed(426992)
randomSegmentOrder <- sample(coreSegments[newLengths > 0], replace = FALSE);
randomStartingPoint <- c(53,108,84,3,63,30,138,179,97,147,132,204,129,198,50,111,56,13,157,167,42,133,26,58,193,54,12,115,165,200,103,98,174)


kappaFun <- function(currentT){
  if(currentT == 1){
    kappa <- 1;
  } else {
    kappa <- 1/sum(DISCOUNT_FACTORS[currentT-1])
  } 
  return(kappa);
}
kappaSumFun <- function(firstT,lastT){
  return(kappa <- sum(1/DISCOUNT_FACTORS[(firstT-1):(lastT-1)]))
}

slFun <- function(currentT){
  return(kappaFun(currentT) * sl);
}
SlFun <- function(firstT, lastT){
  if(firstT > lastT){
    return(0*sl);
  }
  return(kappaSumFun(firstT,lastT) * sl);
}
mlFun <- function(currentT){
  return(kappaFun(currentT) * ml);
}
MlFun <- function(firstT, lastT){
  if(firstT > lastT){
    return(0*sl);
  }
  return(kappaSumFun(firstT,lastT) * ml);
}
CCFun <- function(currentT){
  return(kappaFun(currentT) * cl)
}
SVFun <- function(){
  return(kappaFun(finalT+1) * cl);
}

finalT <- 50;
randomFirstT <- 10;

coreB <- segmentCorrs * CYCLIST_UNIT_PRICE_PER_MINUTE * AADT_FACTOR / 1e6; 
slOri <-  apply(coreB,1,sum);
coreB_RSP <- segmentCorrs_RSP * CYCLIST_UNIT_PRICE_PER_MINUTE * AADT_FACTOR / 1e6; 
slOri_RSP <-  apply(coreB_RSP,1,sum);

mlOri <- correctedMaintenanceCosts * NAF;
clOri <- correctedConstructionCosts * NAF;
NVar <- length(clOri);


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## BLP part #####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#for(scenarioType in c("Short","Long","FromRandom","FromRandomUpdate")){
for(scenarioType in c("RandomFirst","Short","Long","FromRandom","FromRandomUpdate")){
  #for(scenarioType in c("RandomFirst")){
  
  
  theMainOutputDir <- paste0(outputDir,"/AdjustedScenario_", scenarioType)
  if(!file.exists(theMainOutputDir)){
    dir.create(theMainOutputDir);
  }
  
  if(scenarioType == "Short"){
    lastYearOfInvestment <- 5;
    startingBudget <- 46.2;    #58.7 sufficient to buy the realized segments within three years (46.2 for also financed within 5)
  } else {
    lastYearOfInvestment <- finalT;
    startingBudget <- 50;
  }
  budgetIncrement <- startingBudget;
  
  budgetArr = numeric(finalT);
  budgetArr[1] <- startingBudget;
  for(i in 2:length(budgetArr)){
    budgetArr[i] <- budgetArr[i-1] + budgetIncrement;
  }  
  
  
  par(mfrow=c(4,2))
  for(solutionType in solutionTypes){
    #for(solutionType in c("ByRandomSegment")){
    
    if((scenarioType != "Short" & solutionType %in% c("ByActualOrder","ByActualRealisedOrder")) | 
       (scenarioType == "Short" & solutionType == "LPGreedyStop")){
      next;
    } 
    
    theOutputDir <- paste0(theMainOutputDir,"/",solutionType);
    if(!file.exists(theOutputDir)){
      dir.create(file.path(theMainOutputDir, solutionType))
    }
    
    if(solutionType == "ByLongestSegment"){
      theOrder <-  longestSegmentOrder;
    } else if(solutionType == "ByLongestRoute"){
      theOrder <- longestRouteOrder
    } else if(solutionType == "ByShortestSegment"){
      theOrder <- shortestSegmentOrder
    } else if(solutionType == "ByShortestRoute"){
      theOrder <- shortestRouteOrder
    } else if(solutionType == "ByRandomSegment"){
      theOrder <- randomSegmentOrder;
    } else if(solutionType == "ByActualOrder"){
      theOrder <- actualOrder;
    } else if(solutionType == "ByActualRealisedOrder"){
      theOrder <- actualRealisedOrder;
    }
    
    
    budget <- startingBudget;
    maxBudgetReached <- FALSE;
    prevSolution <- rep(0,NVar);
    alreadySpent <- 0;
    expectedTotalNPV <- 0;
    expectedCurrentNPV <- 0;
    
    totalCC <- 0;
    currentCC <- 0;
    ecNPVs <- rep(0,finalT);
    ecCSs <- rep(0,finalT);
    remainingBudgets <- rep(0,finalT);
    alreadySpents <- rep(0,finalT);
    oldMaintenances <- rep(0,finalT)
    consistency <- diag(NVar)
    
    print("########################################################")
    print(paste0("########### ", solutionType, " #################"))
    print("########################################################")
    
    if(scenarioType %in% c("FromRandom","FromRandomUpdate")){
      freeSegments = randomStartingPoint;
      theOrder <- theOrder[theOrder %notin% freeSegments];
    } else {
      freeSegments = c();
      if(scenarioType == "RandomFirst"){
        theOtherSolutionType <- solutionType;
        theOtherOrder <- theOrder[theOrder %notin% randomStartingPoint];
        theOrder <- randomStartingPoint;
        solutionType <- "ByRandomSegment"
      }
    }
    
    freeLinks <- as.data.frame(E(g)$LinkID[E(g)$SegmentId %in% freeSegments]);
    colnames(freeLinks) <- "LinkID"
    data.table::fwrite(freeLinks, file = paste0(theOutputDir, "/SelectedLinks_Predetermined_", solutionType, "_T0.csv"))
    
    ml <- mlOri;
    cl <- clOri;
    sl <- slOri;
    if(scenarioType == "FromRandomUpdate"){
      sl <- slOri_RSP;
    }
    
    ml[freeSegments] <- 0;
    cl[freeSegments] <- 0;
    sl[freeSegments] <- 0;
    prevSolution[freeSegments] <- 1;
    
    
    
    for(currentT in 1:finalT){
      #for(currentT in 1:6){
      
      budget <- budgetArr[currentT];
      
      thisSl <- slFun(currentT);
      Sl <- SlFun(currentT+1,finalT)
      thisMl <- mlFun(currentT);
      Ml <- MlFun(currentT+1,finalT);
      CC <- CCFun(currentT);
      SV <- SVFun();
      
      oldMaintenance <- sum(thisMl * prevSolution)
      
      Zt <- Sl - CC - Ml + SV;
      Zt <- Zt*!prevSolution;
      
      CostsT <- matrix(CC, nrow = 1) #In future t's, also the M's until now.
      #CostsT <- (CostsT + Ml) * !prevSolution
      #CostsT <- (CostsT + thisMl) * !prevSolution
      CostsT <- CostsT * !prevSolution
      
      if(currentT <= lastYearOfInvestment){
        if(scenarioType == "RandomFirst" & (currentT == randomFirstT + 1)){
          theOrder <- theOtherOrder;
          solutionType <- theOtherSolutionType;
        }
        
        conLHS <- rbind(CostsT,consistency)
        
        alreadySpent <- alreadySpent + sum(thisMl * prevSolution) ; 
        remainingBudget <- budget - alreadySpent;
        
        consistencyRHS <- matrix(prevSolution,nrow=NVar)
        conRHS <- rbind(remainingBudget,consistencyRHS)
        
        dirs = c("<=",rep(">=",NVar));
        
        if(solutionType == "LP"){
          lpSol <- lp("max", Zt, conLHS, dirs, conRHS, int.vec = 1:NVar, all.bin = TRUE)
          sol <- lpSol$solution;
        } else {
          sol <- prevSolution;
          if(solutionType == "LPGreedy" | solutionType == "LPGreedyStop"){
            remainingCostsOfL = CostsT;
            theOrder <-  ((1:NVar)[sol==0 & newLengths > 0])[order(-Zt[sol==0 & newLengths > 0]/
                                                                     (remainingCostsOfL[sol==0 & newLengths > 0]))]
          }
          while(length(theOrder)>0){
            nextCandidate <- theOrder[1];
            if(solutionType == "LPGreedyStop"){
              if(Zt[nextCandidate] <= 0){
                break;
              }
            }
            if(remainingBudget >= CostsT[nextCandidate]){
              sol[nextCandidate] <- 1;
              remainingBudget <- remainingBudget - CostsT[nextCandidate];
              if(length(theOrder) == 1){
                theOrder <- c();
              } else {
                theOrder <- theOrder[2:length(theOrder)];
              }
            } else {
              break; #Cannot afford it....
            }
          }
        }
      }
      
      
      deltaSol <- sol - prevSolution;
      objval <- sum(Zt * sol);
      
      alreadySpent <- alreadySpent  +  sum(deltaSol * CostsT);
      expectedTotalNPV <- expectedTotalNPV + objval
      expectedCurrentNPV <- expectedCurrentNPV - sum(CC * deltaSol) + sum(SV * deltaSol) +
        sum(thisSl * prevSolution) - oldMaintenance;
      ecNPVs[currentT] <- expectedCurrentNPV;
      ecCSs[currentT] <- sum(thisSl * prevSolution);
      remainingBudgets[currentT] <- remainingBudget;
      alreadySpents[currentT] <- alreadySpent;
      oldMaintenances[currentT] <- oldMaintenance;
      
      print(c(currentT,sum(deltaSol), sum(sol), alreadySpent, budget, budget-alreadySpent, expectedCurrentNPV, expectedTotalNPV))
      if(alreadySpent > budget){
        print(budget - alreadySpent)
      }
      
      totalCC <- totalCC + sum(deltaSol * CC);
      
      prevSolution <- sol;
      
      
      selectedSegments <- which(as.logical(sol));
      selectedLinks <- as.data.frame(E(g)$LinkID[E(g)$SegmentId %in% selectedSegments]);
      colnames(selectedLinks) <- "LinkID";
      data.table::fwrite(selectedLinks, file = paste0(theOutputDir, "/SelectedLinks_Predetermined_", solutionType,
                                                      "_T", currentT,".csv"))
      selectedSegments <- as.data.frame(selectedSegments)
      colnames(selectedSegments) <- c("SegmentId");
      data.table::fwrite(selectedSegments, file = paste0(theOutputDir, "/SelectedSegments_Predetermined_", solutionType,
                                                         "_T", currentT,".csv"))
    }
    
    plot(1:finalT,ecNPVs,main=paste(scenarioType, solutionType,round(expectedTotalNPV,1)), ylim=c(-1600,2600))
    lines(1:finalT, remainingBudgets, col = 2)
    lines(1:finalT, alreadySpents, col = 4)
    abline(h = max(budgetArr), col = 4, lty = 3)
    lines(1:finalT, oldMaintenances, col = 3)
    
    thisRes <- data.frame(cbind(1:finalT, ecCSs, ecNPVs))
    colnames(thisRes) <- c("t","E_CS_t","E_CurrentNPV")
    data.table::fwrite(thisRes, file = paste0(theOutputDir, "/ExpectedEvaluation_Predetermined_", solutionType,
                                              "_T", finalT,".csv"))
  }
  
  totalSV <- sum(SV[as.logical(prevSolution)]);
  totalCS <- expectedTotalNPV + alreadySpent - totalSV
  totalMC <- alreadySpent - totalCC;
  
  print("Predicted CC, SV, MC, CS, NPV, B/C")
  print(c(totalCC,totalSV, totalCS, expectedTotalNPV, totalCS /(totalCC+totalMC-totalSV) ))
  
}

#########################################
###### Real evaluation of solution setup
#########################################

ODData <- as.data.frame(data.table::fread("O:/Public/4233-82676-BIKELONGER-persondata/CBA/OD_BikeType_SpeedType/SC01_MixBase1_MixSc1.csv"));
ODData <- subset(ODData, !is.na(ODData$value));
ZeroBenefits = matrix(0,nrow = N_centroids, ncol = N_centroids);


#########################################
###### Real evaluation of solutions
#########################################


#for( scenarioType in c("FromRandom","Short","Long","FromRandomUpdate","RandomFirst")){
for( scenarioType in c("Long","RandomFirst","Short")){
  mainOutputDir <- paste0(outputDir,"/AdjustedScenario_",scenarioType);
  for(solutionType in solutionTypes){
    theOutputDir <- paste0(mainOutputDir, "/", solutionType)
    
    if((scenarioType != "Short" & solutionType == "ByActualOrder") |
       (scenarioType == "Short" & solutionType == "LPGreedyStop")){
      next;
    } else if(scenarioType == "RandomFirst" & solutionType == "ByRandomSegment"){
      #Simply copy the one from the long run....
      file.copy("M:/GeoCBA/ModelOutputs/AdjustedScenario_Long/ByRandomSegment/Evaluation_Predetermined_ByRandomSegment_T50.csv",
                "M:/GeoCBA/ModelOutputs/AdjustedScenario_RandomFirst/ByRandomSegment/", overwrite = FALSE)
      next;
    }
    
    if(scenarioType != "RandomFirst"){
      freeLinks <- data.table::fread(paste0(theOutputDir, "/SelectedLinks_Predetermined_",solutionType,"_T0.csv"))$LinkID;
    } else {
      freeLinks <- data.table::fread(paste0(theOutputDir, "/SelectedLinks_Predetermined_ByRandomSegment_T0.csv"))$LinkID;
    }
    freeSegments <- logical(N_segments);
    freeSegments[unique(E(g)$SegmentId[E(g)$LinkID %in% freeLinks])] <- TRUE
    
    prevCurrentSegments <- latestEvaluatedSolution <- freeSegments;
    cumulCC <- cumulMC <- cumulSV <- cumulCS <- currentNPV <- currentBCR <-  totalTravelTimeSavings <- 0;
    results <- matrix(nrow=0, ncol = 13)
    colnames(results) <- c("t","CC_t","SV_t","MC_t","CS_t","NPV_t","CumulCC_t","CumulSV_t","CumulMC_t","CumulCS_t","CurrentNPV","CurrentBCR","N_t")
    
    
    ml <- mlOri;
    cl <- clOri;
    sl <- slOri;
    if(scenarioType == "FromRandomUpdate"){
      sl <- slOri_RSP;
    }
    ml[freeSegments] <- 0;
    cl[freeSegments] <- 0;
    sl[freeSegments] <- 0;
    prevCurrentSegments[freeSegments] <- TRUE;
    
    #Trying to load previous half-created file. 
    tFromEarlier <- 1; 
    for(currentT in finalT:2){  
      fname <- paste0(theOutputDir, "/Evaluation_Predetermined_",solutionType,"_T", currentT,".csv");
      if(file.exists(fname) | (scenarioType == "RandomFirst" & currentT == 11)){
        tFromEarlier <- currentT;
        print(paste0("Found evaluation from ", tFromEarlier))
        if(currentT == finalT){
          break;
        }
        if(scenarioType == "RandomFirst" & currentT == 11){
          results <- data.frame(data.table::fread(
            "M:/GeoCBA/ModelOutputs/AdjustedScenario_Long/ByRandomSegment/Evaluation_Predetermined_ByRandomSegment_T50.csv"));
          results <- results[1:currentT,];
        } else {
          results <- data.frame(data.table::fread(paste0(theOutputDir, "/Evaluation_Predetermined_",solutionType,"_T", currentT,".csv")));
        }
        latestCS <- results$CS_t[dim(results)[1]];
        totalTravelTimeSavings <- latestCS * DISCOUNT_FACTORS[currentT] * 1e6 / CYCLIST_UNIT_PRICE_PER_MINUTE / AADT_FACTOR;
        
        results <- results[1:(dim(results)[1]-1),];
        cumulCS <- results$CumulCS_t[dim(results)[1]];
        cumulSV <- results$CumulSV_t[dim(results)[1]];
        cumulCC <- results$CumulCC_t[dim(results)[1]];
        cumulMC <- results$CumulMC_t[dim(results)[1]];
        currentNPV <- results$CurrentNPV[dim(results)[1]]
        
        
        
        if(scenarioType == "RandomFirst" & currentT == 11){
          chosenLinks <- data.table::fread(file =  
                                             "M:/GeoCBA/ModelOutputs/AdjustedScenario_Long/ByRandomSegment/SelectedLinks_Predetermined_ByRandomSegment_T10.csv")$LinkID
        } else {
          chosenLinks <- data.table::fread(file =  paste0(theOutputDir,"/SelectedLinks_Predetermined_",solutionType,"_T", currentT-1,".csv"))$LinkID
        }
        chosenLinks <- chosenLinks[chosenLinks %notin% freeLinks];
        chosenLinksBool <- E(g)$LinkID %in% chosenLinks
        prevCurrentSegments <- logical(N_segments);
        prevCurrentSegments[unique(E(g)$SegmentId[chosenLinksBool])] <- TRUE;
        break;
      }
    }
    
    
    if(tFromEarlier < finalT){
      for(currentT in tFromEarlier:finalT){
        
        fname <- paste0(theOutputDir, "/SelectedLinks_Predetermined_",solutionType,"_T", currentT,".csv");
        if(file.exists(fname)){
          print(fname)
          chosenLinks <- data.table::fread(file = fname)$LinkID
          chosenLinks <- chosenLinks[chosenLinks %notin% freeLinks];
        }
        chosenLinksBool <- E(g)$LinkID %in% chosenLinks
        
        currentSegments <- logical(N_segments);
        currentSegments[unique(E(g)$SegmentId[chosenLinksBool])] <- TRUE;
        
        newlyAddedSegments <- prevCurrentSegments != currentSegments;
        
        MCPart <- sum(mlFun(currentT)[prevCurrentSegments]); # NAF inclued in cl ##Is this currentT or currentT+1-1?
        CCPart <- sum(CCFun(currentT)[newlyAddedSegments]); # NAF included in cl ##Is this currentT or currentT+1-1?
        SVPart <- sum(SVFun()[newlyAddedSegments]); # NAF included in cl
        
        rawConsumerSurplus <- totalTravelTimeSavings / 1e6 * CYCLIST_UNIT_PRICE_PER_MINUTE * AADT_FACTOR
        if(currentT ==1){
          CSPart <- rawConsumerSurplus;
        } else {
          CSPart <- rawConsumerSurplus / DISCOUNT_FACTORS[currentT-1]
        }
        
        
        NPV <- CSPart + SVPart - CCPart - MCPart;
        
        cumulCS <- cumulCS + CSPart;
        cumulSV <- cumulSV + SVPart;
        cumulCC <- cumulCC + CCPart;
        cumulMC <- cumulMC + MCPart;
        
        currentNPV <- currentNPV + NPV;
        currentBCR <- cumulCS / (cumulCC + cumulMC - cumulSV)
        currentN <- sum(currentSegments);
        
        results <- rbind(results, c(currentT, CCPart, SVPart, MCPart, CSPart, NPV, cumulCC, cumulSV, cumulMC, cumulCS, 
                                    currentNPV, currentBCR, currentN))
        data.table::fwrite(data.frame(results), paste0(theOutputDir, "/Evaluation_Predetermined_",solutionType,"_T", currentT,".csv"),
                           quote = FALSE, row.names = FALSE);
        if(file.exists( paste0(theOutputDir, "/Evaluation_Predetermined_",solutionType,"_T", currentT-1,".csv"))){
          file.remove(paste0(theOutputDir, "/Evaluation_Predetermined_",solutionType,"_T", currentT-1,".csv"));
        }
        
        print(data.frame(results))
        
        #Calculatng travel time savings for next period (if needed, i.e. if any changes since last evaluation)
        if(sum(latestEvaluatedSolution != currentSegments) >0){
          chosenType <- ifelse((E(g)$LinkID %in% freeLinks) | chosenLinksBool, E(g)$InfraType_Plan_30_45, E(g)$InfraType_Basis);
          
          
          bikeTypes <- c(1,2,3);
          speedTypes <- c(1,2,3);
          travelTimeSavings <- numeric(length(bikeTypes)*length(speedTypes));
          
          totalOD <- matrix(0,nrow=length(CentroidNodes$NodeId), ncol = length(CentroidNodes$NodeId));
          totalOD_B <- matrix(0,nrow=length(CentroidNodes$NodeId), ncol = length(CentroidNodes$NodeId));
          totalOD_B_OD <- matrix(0,nrow=length(CentroidNodes$NodeId), ncol = length(CentroidNodes$NodeId));
          
          typeIndex <- 0;
          for(bikeType in bikeTypes){
            for(speedType in speedTypes){
              typeIndex <- typeIndex + 1;
              
              
              OD <- createODMatrix(subset(ODData, BikeType==bikeType & SpeedType == speedType));
              
              basisFileName <- paste0(outputDir,"/gcs_Basis_",bikeType,"_",speedType,".Rda");
              if(scenarioType == "FromRandomUpdate"){
                basisFileName <- paste0(outputDir,"/gcs_Basis_",bikeType,"_",speedType,"_RSP.Rda");
              }
              while(!file.exists(basisFileName)){
                print("Temporary connection problems. Waiting 20 seconds to hopefully overcome the issue....")
                Sys.sleep(20); 
              }
              gcs_Basis <- readRDS(file=basisFileName) 
              
              vMax_1 <- determineVMax(1,bikeType, speedType);
              vMax_2 <- determineVMax(2,bikeType, speedType);
              vMax_3 <- determineVMax(3,bikeType, speedType);
              
              print(paste("Running scenario for",bikeType,speedType));
              # All scenario
              E(g)$weight <- 1000000;
              E(g)$weight <- ifelse(chosenType == 1, E(g)$Length / (vMax_1 / 3.6) + E(g)$IntersectionDelay, E(g)$weight);
              E(g)$weight <- ifelse(chosenType == 2, E(g)$Length / (vMax_2 / 3.6) + E(g)$IntersectionDelay, E(g)$weight);
              E(g)$weight <- ifelse(chosenType == 3, E(g)$Length / (vMax_3 / 3.6) + E(g)$IntersectionDelay, E(g)$weight);    
              
              out <- assignFlows(g, OD, ZeroBenefits, FALSE, FALSE);
              gcs_All <- out$GCs / 60;
              
              ODLengths_All <- out$ODLengths;
              ## Calculating differences
              OD_B <- gcs_Basis - gcs_All;
              
              if(min(OD_B) < -1e-12){
                print(paste0("Negative CS :( ", sum(OD_B<0)))
                print(kill + to + kill)
              }
              
              travelTimeSavings[typeIndex] <- sum(OD*OD_B);
            }
          }
          totalTravelTimeSavings <- sum(travelTimeSavings);
          latestEvaluatedSolution <- currentSegments;
        }
        prevCurrentSegments <- currentSegments
        gc();
      }
    }
  }
}


###### # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
## Compare expected and actual CSs #####
###### # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

fs <- 1;
palette(c("black","red2","forestgreen","blue","turquoise2","orchid1", "seagreen1","goldenrod1","gray","chocolate1"));
legendCatsToUse <- c("BLP", "Greedy w/ stop", "Greedy w/o stop","Actual order","Random order","Shorter segments first", "Shorter routes first",
                     "Longer segments first", "Longer routes first");
thisPCHs <- c(20,20,20,20,20,20,20,20,20)
thisCols <- c(1,2,8,8,6,5,4,7,3)
thisCEXs <- c(2,sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2))



for(budgetType in c("Short","Long","RandomFirst")){
  theMainOutputDir <- paste0(outputDir,"/AdjustedScenario_", budgetType);
  #png(paste0(theMainOutputDir,"/EvaluationComparisonPlot_",budgetType,".png"), width = 700, height = 1000, pointsize = 20);
  
  fs <- fs * 2;
  allEvalSums <- NULL;
  for(solutionType in solutionTypes){
    theOutputDir <- paste0(theMainOutputDir,"/",solutionType);
    if((budgetType != "Short" & solutionType == "ByActualOrder") |
       (budgetType == "Short" & solutionType == "LPGreedyStop")){
      skipI <- which(solutionType == solutionTypes)
      next;
    }
    
    if(file.exists(paste0(theOutputDir, "/Evaluation_Predetermined_", solutionType,"_T", finalT,".csv"))){
      thisExpRes <- data.frame(data.table::fread(paste0(theOutputDir, "/ExpectedEvaluation_Predetermined_", solutionType,
                                                        "_T", finalT,".csv")))
      thisRes <- data.frame(data.table::fread(paste0(theOutputDir, "/Evaluation_Predetermined_", solutionType,
                                                     "_T", finalT,".csv")))
      thisRes <- merge(thisRes, thisExpRes);
      pdf(paste0(theMainOutputDir,"/CSComparison_",solutionType,"_",budgetType,".pdf"), width=12.5, height = 7.0, pointsize = 14);
      par(oma = c(0, 0, 0, 0), mar = c(fs*2, fs*2.3, 0.1, 0.1))
      plot(NA,xlim=c(0.8,49.2),ylim=c(-5,84),
           xlab = "Year", ylab = "Consumer Surplus [million DKK]", 
           mgp=c(sqrt(2)*fs,1,0), cex.lab = fs, cex.axis = fs);
      abline(v=seq(-5,55,1), col = rgb(0.975,0.975,0.975), lwd = 0.05)
      abline(v=seq(-5,55,5), col = rgb(0.9,0.9,0.9), lwd = 0.1)
      abline(h=seq(-10,90,2), col = rgb(0.975,0.975,0.975), lwd = 0.05)
      abline(h=seq(-20,90,10), col = rgb(0.9,0.9,0.9), lwd = 0.1)
      
      points(1:finalT, thisRes$CS_t, col = 1, pch = 20, cex = 2*fs)
      points(1:finalT, thisRes$E_CS_t, col = 2, pch = 20, cex = sqrt(2)*fs)
      #points(1:finalT, thisRes$E_CS_t - thisRes$CS_t, pch = 20, col = 4)
      legend("topright",c("Proper evaluation", "Linear approximation"), pch=c(20,20),col=c(1,2), pt.cex = c(2*fs,sqrt(2)*fs), cex = fs, bg = rgb(1,1,1));
      dev.off();
      
      thisRes$Config <- solutionType
      if(is.null(allEvalSums)){
        allEvalSums <- thisRes;
      } else {
        allEvalSums <- rbind(allEvalSums, thisRes)
      }
    }
  }
  
  
  fs <- fs / 2;
  allEvalSums$Config <- as.character(allEvalSums$Config);
  configs <- unique(allEvalSums$Config)
  
  if(length(configs)>0){
    
    pdf(paste0(theMainOutputDir,"/NPV_",budgetType,".pdf"), width=12.5, height = 7.0, pointsize = 14);
    par(oma = c(0, 0, 0, 0), mar = c(fs*3.5, fs*3.5, 0.1, 0.1))
    plot(NA,xlim=c(1,50),ylim=c(-1500,1500),
         xlab = "Year", ylab = "Net Present Value [million DKK]", 
         mgp=c(2.4*fs,1,0), cex.lab = fs, cex.axis = fs);
    abline(h=seq(-2500,2500,100), col = rgb(0.975,0.95,0.975), lwd = 0.05)
    abline(h=seq(-2500,2500,500), col = rgb(0.9,0.9,0.9), lwd = 0.1)
    abline(v=seq(-10,100,1), col = rgb(0.975,0.975,0.975), lwd = 0.05)
    abline(v=seq(-10,100,5), col = rgb(0.9,0.9,0.9), lwd = 0.1)
    for(i in 1:length(solutionTypes)){
      if(i == skipI){
        next;
      }
      conf = solutionTypes[i];
      col = thisCols[i];
      pch = thisPCHs[i];
      cex = thisCEXs[i];
      thisEvalSum = subset(allEvalSums, Config == conf);
      points(thisEvalSum$t, thisEvalSum$CurrentNPV, col = col, pch = pch, lwd = 2, cex=cex)
      lines(thisEvalSum$t, thisEvalSum$CurrentNPV, col = col, pch = pch, lty = 3)
    }
    legend("bottomleft",legendCatsToUse[1:length(solutionTypes) != skipI], 
           pch = thisPCHs[1:length(solutionTypes) != skipI], col = thisCols[1:length(solutionTypes) != skipI], 
           lty = 3, y.intersp=1, lwd = 2, text.width = 9.5, pt.cex = thisCEXs[1:length(solutionTypes) != skipI], 
           bg = rgb(1,1,1))
    dev.off();
    
    pdf(paste0(theMainOutputDir,"/BCR_",budgetType,".pdf"), width=12.5, height = 7.0, pointsize = 14);
    par(oma = c(0, 0, 0, 0), mar = c(fs*3.5, fs*3.5, 0.1, 0.1))
    plot(NA,xlim=c(1,50),ylim=c(0.03,2.97),
         xlab = "Year", ylab = "Benefit/cost ratio",
         mgp=c(2.4*fs,1,0), cex.lab = fs, cex.axis = fs);
    abline(h=seq(-2.500,6.500,0.100), col = rgb(0.975,0.95,0.975), lwd = 0.05)
    abline(h=seq(-2.500,6.500,0.500), col = rgb(0.9,0.9,0.9), lwd = 0.1)
    abline(v=seq(-10,100,1), col = rgb(0.975,0.975,0.975), lwd = 0.05)
    abline(v=seq(-10,100,5), col = rgb(0.9,0.9,0.9), lwd = 0.1)
    for(i in 1:length(solutionTypes)){
      if(i == skipI){
        next;
      }
      conf = solutionTypes[i];
      col = thisCols[i];
      pch = thisPCHs[i];
      cex = thisCEXs[i];
      thisEvalSum = subset(allEvalSums, Config == conf);
      points(thisEvalSum$t, thisEvalSum$CurrentBCR, col = col, pch = pch, lwd = 2, cex = cex)
      lines(thisEvalSum$t, thisEvalSum$CurrentBCR, col = col, pch = pch, lty = 3)
    }
    legend("topleft",legendCatsToUse[1:length(solutionTypes) != skipI], 
           pch = thisPCHs[1:length(solutionTypes) != skipI], col = thisCols[1:length(solutionTypes) != skipI], 
           lty = 3, y.intersp=1, lwd = 2, text.width = 9.5, pt.cex = thisCEXs[1:length(solutionTypes) != skipI],
           bg = rgb(1,1,1))
    dev.off();
    
  }
  
}

##############################################################################################################
############## Temporary NAF corrector (and for correcting other mishabs #####################################
##############################################################################################################
if(FALSE){
  for( scenarioType in c("Long","Short")){
    mainOutputDir <- paste0(outputDir,"/AdjustedScenario_",scenarioType);
    for(solutionType in solutionTypes){
      theOutputDir <- paste0(mainOutputDir, "/", solutionType)
      if(file.exists(paste0(theOutputDir, "/Evaluation_Predetermined_",solutionType,"_T", finalT,".csv"))){
        results <- data.frame(data.table::fread(paste0(theOutputDir, "/Evaluation_Predetermined_",solutionType,"_T", finalT,".csv")));
        results[c("CC_t","SV_t","MC_t","CumulCC_t","CumulSV_t","CumulMC_t")] <- 
          results[c("CC_t","SV_t","MC_t","CumulCC_t","CumulSV_t","CumulMC_t")] * NAF
        results$NPV_t <- results$CS_t + results$SV_t - results$CC_t - results$MC_t;
        results$CurrentNPV <- results$CumulCS_t + results$CumulSV_t - results$CumulCC_t - results$CumulMC_t;
        results$CurrentBCR <- results$CumulCS_t /(-results$CumulSV_t + results$CumulCC_t + results$CumulMC_t);
        data.table::fwrite(results, paste0(theOutputDir, "/Evaluation_Predetermined_",solutionType,"_T", finalT,"_WithNAF.csv"),
                           quote = FALSE, row.names = FALSE);
      }
    }
  }
}

if(FALSE){
  for( scenarioType in c("Long","Short")){
    mainOutputDir <- paste0(outputDir,"/AdjustedScenario_",scenarioType);
    for(solutionType in solutionTypes){
      theOutputDir <- paste0(mainOutputDir, "/", solutionType)
      if(file.exists(paste0(theOutputDir, "/Evaluation_Predetermined_",solutionType,"_T", finalT,".csv"))){
        results <- data.frame(data.table::fread(paste0(theOutputDir, "/Evaluation_Predetermined_",solutionType,"_T", finalT,".csv")));
        if("N_t" %notin% colnames(results)){
          print(solutionType)
          Nts <- numeric(finalT);
          for(currentT in 1:finalT){
            selectedSegments <- data.frame(data.table::fread(paste0(theOutputDir,  "/SelectedSegments_Predetermined_",solutionType,"_T", currentT,".csv")));
            Nts[currentT] <- dim(selectedSegments)[1]
          }
          results$N_t <- Nts;
          data.table::fwrite(results, paste0(theOutputDir,"/Evaluation_Predetermined_",solutionType,"_T", finalT,".csv"),
                             quote = FALSE, row.names = FALSE);
        }
      }
    }
  }
}




######################################################
######## Entering the old part of the code ###########
######################################################




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
### Plotting different scenarios #####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

coreDir = "M:/GeoCBA/ModelOutputs";
plotDir <- "M:/GeoCBA/ModelOutputs/Plots/Analysis"
palette(c("black","red2","forestgreen","blue","turquoise2","orchid1", "seagreen1","goldenrod1","gray","chocolate1"));
fs <- 1;

plotDir <- "M:/GeoCBA/Figures"


printKPIs <- TRUE;

alternativeStartingPoint = FALSE;

selectedConfigsOuter <- c("LP", "LPGreedy", "ByLongestRoute","ByShortestRoute",
                          "ByLongestSegment","ByShortestSegment","ByRandomSegment","ByActual");

forecastVoTByGDP <- FALSE;

efficiencyLowerBound <- 1;


maintenanceFactors <- c(1,2,4);
efficiencyLowerBounds <- c("0","05","08","09","1")  
#includeNAFs <- c(TRUE,FALSE);  ##Evaluation summaries are done without NAF. To include NAF, change this boolean. 
includeNAF <- TRUE; #Already included in the real-life evaluations!





firstParts1 <- c("1_1_1", "1_1_1", "ByLongestRoute","ByShortestRoute",
                 "ByLongestSegment","ByShortestSegment","ByRandomSegment");
lastParts1 <- c("_0_SemiNaive", "_1_SemiNaive", "","","","","");
legendCats1 <- c("Heuristic 2", "Heuristic 2 w/ int. stop",
                 "Longer routes first", "Shorter routes first",
                 "Longer segments first", "Shorter segments first", "Random order")
cols1 <- c(2,2,3,10,7,8,9);
PCHs1 <- c(1,2,1,1,1,1,1);
firstParts2 <- c("1_1_1","1_1_1","1_1_1","1_1_1","1_1_1","1_1_1","1_1_1","1_1_1")
lastParts2 <- c("_0_Naive", "_0_SemiNaive", "_0_Pessimistic1", "_0_Pessimistic3", 
                "_1_Naive","_1_SemiNaive", "_1_Pessimistic1", "_1_Pessimistic3");
legendCats2 <- c("Heuristic 1", "Heuristic 2", "Heuristic 3", "Heuristic 4",
                 "Heuristic 1 w/ int. stop", "Heuristic 2 w/ int. stop", "Heuristic 3 w/ int. stop", "Heuristic 4 w/ int. stop")
cols2 <- c(6,2,5,4,6,2,5,4);
PCHs2 <- c(1,1,1,1,2,2,2,2);

cols3 <- c(cols2, 3,10,7,8,9)
PCHs3 <- c(PCHs2, c(4,4,4,4,4))
firstParts3 <- c(firstParts2,"ByLongestRoute","ByShortestRoute",
                 "ByLongestSegment","ByShortestSegment","ByRandomSegment")
lastParts3 <- c(lastParts2, "","","","","","")
legendCats3 <- c(legendCats2,  "Longer routes first", "Shorter routes first",
                 "Longer segments first", "Shorter segments first", "Random order");

cols4 <- c(cols2[1:4], 3,10,7,8,9,1)
PCHs4 <- c(PCHs2[1:4], c(4,4,4,4,4),20)
firstParts4 <- c(firstParts2[1:4],"ByLongestRoute","ByShortestRoute",
                 "ByLongestSegment","ByShortestSegment","ByRandomSegment","ByActual")
lastParts4 <- c(lastParts2[1:4], "","","","","","")
legendCats4 <- c(legendCats2[1:4],  "Longer routes first", "Shorter routes first",
                 "Longer segments first", "Shorter segments first", "Random order","Actual order");


#plotType = 1: Best compared to common, plotType = 2: Comparing optis. 3: All,  4: Vs Actual
#for(PLOT_TYPE in c(1,2,3,4)){
for(PLOT_TYPE in c(4)){
  
  if(PLOT_TYPE == 1){
    thisFirstParts <- firstParts1
    thisLastParts <- lastParts1
    legendCats <- legendCats1
    cols <- cols1;
    PCHs <- PCHs1;
  } else if(PLOT_TYPE == 2){
    thisFirstParts <- firstParts2
    thisLastParts <- lastParts2
    legendCats <- legendCats2
    cols <- cols2;
    PCHs <- PCHs2;
  } else if(PLOT_TYPE == 3){
    thisFirstParts <- firstParts3
    thisLastParts <- lastParts3
    legendCats <- legendCats3
    cols <- cols3;
    PCHs <- PCHs3;
  } else if(PLOT_TYPE == 4){
    thisFirstParts <- firstParts4
    thisLastParts <- lastParts4
    legendCats <- legendCats4
    cols <- cols4;
    PCHs <- PCHs4;
  }
  
  
  if(PLOT_TYPE == 4){
    maxBudgets <- c(500);
    asps <- c(FALSE)
  } else {
    maxBudgets <- c(900000)
    asps <- c(FALSE,TRUE)
  }
  
  for(maxBudget in maxBudgets){
    for(alternativeStartingPoint in asps){
      for(maintenanceFactor in c(1)){#c(4,2,1)){
        identifier <- paste0(maintenanceFactor, "_", efficiencyLowerBound, "_", includeNAF);
        
        
        dirs = list.files("M:/GeoCBA/ModelOutputs", pattern = "*Configuration*")
        allEvalSums = NULL;
        
        selectedConfigs <- paste(thisFirstParts,maintenanceFactor, sep ="_")
        selectedConfigs <- paste0(selectedConfigs,thisLastParts)
        
        #selectedConfigs <- selectedConfigsOuter
        if(alternativeStartingPoint){
          selectedConfigs <- paste(selectedConfigs, "ASP", sep = "_");
          aspSuffix = "_ASP"
        } else {
          aspSuffix = "";
        }
        if(maxBudget == 500){
          selectedConfigs <- paste(selectedConfigs, maxBudget, sep = "_");
          mbSuffix = "_500"
        } else {
          mbSuffix = "";
        }
        
        for(d in dirs){
          thisDir = paste0(coreDir,"/",d);
          f = list.files(thisDir, pattern = "*EvaluationSummary*");
          if(length(f) > 0){
            f = paste0(thisDir,"/",f)
            
            found <- FALSE; 
            for(sc in selectedConfigs){
              #    if(!grepl("By",sc,fixed=TRUE)){ 
              #      sc <- paste(sc,efficiencyLowerBound, sep = "_")
              #    }
              if(grepl(paste0("EvaluationSummary_",sc,".csv"),f,fixed=TRUE)){ found <- TRUE; break; } }
            
            if(found){
              thisEvalSum <- read.csv(f);
              
              if(forecastVoTByGDP){
                thisEvalSum$ConsumerSurplus =  thisEvalSum$ConsumerSurplus * 1.009^thisEvalSum$Year
              }
              
              if(includeNAF){
                thisEvalSum[,c("ConstructionCosts","MaintenanceCosts","FutureScrapValue","ScrapValue")] = 
                  NAF * thisEvalSum[,c("ConstructionCosts","MaintenanceCosts","FutureScrapValue","ScrapValue")];
                thisEvalSum$DiscountedNetCashFlow <- thisEvalSum$ConsumerSurplus + thisEvalSum$FutureScrapValue - 
                  thisEvalSum$ConstructionCosts - thisEvalSum$MaintenanceCosts;
              }
              
              thisEvalSum$NetPresentValue <- cumsum(thisEvalSum$DiscountedNetCashFlow)
              thisEvalSum$BCR <- cumsum(thisEvalSum$ConsumerSurplus) / 
                cumsum(thisEvalSum$ConstructionCosts - thisEvalSum$FutureScrapValue + thisEvalSum$MaintenanceCosts);
              thisEvalSum$IRR <- (thisEvalSum$BCR^(1/thisEvalSum$Year) - 1 ) * 100;
              
              if(printKPIs){
                print(sc)
                print(c("Construction Costs", "Scrap Value", "Maintenance Costs","Time Benefits", "Net Present Value"))
                print(c(sum(thisEvalSum$ConstructionCosts), sum(thisEvalSum$ScrapValue), sum(thisEvalSum$MaintenanceCosts),
                        sum(thisEvalSum$ConsumerSurplus), thisEvalSum$NetPresentValue[dim(thisEvalSum)[1]]));
                print(paste0("B/C: ", thisEvalSum$BCR[dim(thisEvalSum)[1]]))
                
              }
              
              if(is.null(allEvalSums)){
                allEvalSums = thisEvalSum;
              } else {
                allEvalSums = rbind(allEvalSums,thisEvalSum);
              }
              
            }
          }
        }
        
        allEvalSums$Config <- as.character(allEvalSums$Config);
        foundConfigs <- unique(allEvalSums$Config)
        configs <- selectedConfigs
        configsBl <- logical(length(configs));
        for(i in 1:length(configs)){
          if( configs[i] %in% foundConfigs){
            configsBl[i] <- TRUE;
          }
        }
        configs <- configs[configsBl];
        legendCatsToUse <- legendCats[configsBl];
        thisPCHs <- PCHs[configsBl];
        thisCols <- cols[configsBl];
        print(configs)
        
        if(length(configs)< length(selectedConfigs)){
          identifier <- paste0("INCOMPLETE_", identifier);
        }
        
        #png(paste0(plotDir,"/NetPresentValueComparisons_", identifier, ".png"), width=1200, height = 900, res = 120);
        png(paste0(plotDir,"/NetPresentValueComparisons_Final", maintenanceFactor,
                   "_Type", PLOT_TYPE, aspSuffix, mbSuffix, ".png"), width=1200, height = 900, res = 120);
        par(oma = c(0, 0, 0, 0), mar = c(fs*3.5, fs*3.5, 0.1, 0.1))
        plot(NA,xlim=c(1,50),ylim=c(-1433,1200),
             xlab = "Year", ylab = "Net Present Value [million DKK]", 
             mgp=c(2.4*fs,1,0), cex.lab = fs, cex.axis = fs);
        abline(h=seq(-2500,2500,100), col = rgb(0.975,0.95,0.975), lwd = 0.05)
        abline(h=seq(-2500,2500,500), col = rgb(0.9,0.9,0.9), lwd = 0.1)
        abline(v=seq(-10,100,1), col = rgb(0.975,0.975,0.975), lwd = 0.05)
        abline(v=seq(-10,100,5), col = rgb(0.9,0.9,0.9), lwd = 0.1)
        for(i in 1:length(configs)){ ## TODO corrects formats
          conf = configs[i];
          col = thisCols[i];
          pch = thisPCHs[i];
          thisEvalSum = subset(allEvalSums, Config == conf);
          points(thisEvalSum$Year, thisEvalSum$NetPresentValue, col = col, pch = pch, lwd = 2, cex=sqrt(2))
          lines(thisEvalSum$Year, thisEvalSum$NetPresentValue, col = col, pch = pch, lty = 3)
        }
        legend("bottomleft",legendCatsToUse, pch = thisPCHs, col = thisCols, 
               lty = 3, y.intersp=1, lwd = 2, text.width = 9.5, pt.cex = sqrt(2))
        dev.off();
        
        
        png(paste0(plotDir,"/BCR_Final", maintenanceFactor,
                   "_Type", PLOT_TYPE, aspSuffix, mbSuffix, ".png"), width=1200, height = 900, res = 120);
        par(oma = c(0, 0, 0, 0), mar = c(fs*3.5, fs*3.5, 0.1, 0.1))
        plot(NA,xlim=c(1,50),ylim=c(-0.02,5.02),
             xlab = "Year", ylab = "Benefit/cost ratio",
             mgp=c(2.4*fs,1,0), cex.lab = fs, cex.axis = fs);
        abline(h=seq(-2.500,6.500,0.100), col = rgb(0.975,0.95,0.975), lwd = 0.05)
        abline(h=seq(-2.500,6.500,0.500), col = rgb(0.9,0.9,0.9), lwd = 0.1)
        abline(v=seq(-10,100,1), col = rgb(0.975,0.975,0.975), lwd = 0.05)
        abline(v=seq(-10,100,5), col = rgb(0.9,0.9,0.9), lwd = 0.1)
        for(i in 1:length(configs)){ ## TODO corrects formats
          conf = configs[i];
          col = thisCols[i];
          pch = thisPCHs[i];
          thisEvalSum = subset(allEvalSums, Config == conf);
          points(thisEvalSum$Year, thisEvalSum$BCR, col = col, pch = pch, lwd = 2, cex=sqrt(2))
          lines(thisEvalSum$Year, thisEvalSum$BCR, col = col, pch = pch, lty = 3)
        }
        legend("topleft",legendCatsToUse, pch = thisPCHs, col = thisCols,
               lty = 3, y.intersp=1, lwd = 2, text.width = 9.5)
        dev.off();
        
        
        png(paste0(plotDir,"/IRR_Final", maintenanceFactor,
                   "_Type", PLOT_TYPE, aspSuffix, mbSuffix, ".png"), width=1200, height = 900, res = 120);
        par(oma = c(0, 0, 0, 0), mar = c(fs*3.5, fs*3.5, 0.1, 0.1))
        plot(NA,xlim=c(1,50), ylim=c(-16,11),
             xlab = "Year", ylab = "Internal rate of return [%]",
             mgp=c(2.4*fs,1,0), cex.lab = fs, cex.axis = fs);
        abline(h=seq(-30,30,1), col = rgb(0.975,0.95,0.975), lwd = 0.05)
        abline(h=seq(-30,30,5), col = rgb(0.9,0.9,0.9), lwd = 0.1)
        abline(v=seq(-10,100,1), col = rgb(0.975,0.975,0.975), lwd = 0.05)
        abline(v=seq(-10,100,5), col = rgb(0.9,0.9,0.9), lwd = 0.1)
        for(i in 1:length(configs)){ ## TODO corrects formats
          conf = configs[i];
          col = thisCols[i];
          pch = thisPCHs[i];
          thisEvalSum = subset(allEvalSums, Config == conf);
          points(thisEvalSum$Year, thisEvalSum$IRR, col = col, pch = pch, lwd = 2, cex=sqrt(2))
          lines(thisEvalSum$Year, thisEvalSum$IRR, col = col, pch = pch, lty = 3)
        }
        legend("bottomright",legendCatsToUse, pch = thisPCHs, col = thisCols, 
               lty = 3, y.intersp=1, lwd = 2, text.width = 9.5)
        dev.off();
        
      }
    }
  }
}





#### ######
#Functions for running the whole damn thing

##Sparsitity may be utilised?

updateAuxMat <- function(funAuxMat, newSelections){
  funAuxMat[,newSelections] <- 0;
  return(funAuxMat);
}
calculateZ <- function(deltas){
  V <- sum(deltas*(diag(usedCorrStructures)  - apply(usedCorrStructures * (deltas-1),2,max))) * CYCLIST_UNIT_PRICE_PER_MINUTE * AADT_FACTOR * EVALUATION_PERIOD / 1e6;
  return(sum(V));
}
calculateC <- function(deltas){
  return(sum(newCosts*deltas))
}
calculateInvestedLength <- function(deltas){
  return(sum(newLengths*deltas))
}
calculatePotentialOld <- function(auxMat){
  P <- newLengths*apply(auxMat,2, 
                        function(x){
                          m <- max(x);
                          w <- x == m;
                          if(any(w==0)){
                            cost <- sum(newCosts[w]);
                            m2 <- max(x[!w]);
                            return(m - m2 + cost)
                          } else {
                            return(0);
                          }
                        });
  return(P);
}
takeStepOld <- function(funDeltas, n){
  ##Chooses the link(s) that upgrade(s) a single link the most
  auxMat <- calculateAuxMat(funDeltas);
  P <- calculatePotentialOld(auxMat);
  w <- logical(length(P));
  nSelected <- 0;
  while(nSelected < n){
    i <- which.max(P);
    P[i] <- -Inf
    m <- max(auxMat[,i]);
    wTemp <- auxMat[,i] == m;
    w[wTemp] <- TRUE;
    nSelected <- nSelected + sum(wTemp);
  }
  funDeltas[w] <- 1;
  return(funDeltas);
}

updatePotentialMatrixSemiOld <- function(P,funAuxMat,indices){
  st <- Sys.time();
  if(sum(indices)>1){
    P[indices, indices] <- apply(funAuxMat[indices,indices], 1, calculatePotentialAux);
  } 
  et <- Sys.time();
  #  print(paste("Time spent updating potential matrix:", et-st));
  return(P);
}
calculatePotentialAux <- function(x, subProjs){
  out <- numeric(length(x));
  m <- max(x)
  w <- x == m;
  if(any(!w)){
    w <- which(w)[1];
    subP <- subProjs[w];
    w <- subProjs == subP;
    m2 <- max(x[!w]);
    out[w] <- m - m2;
  }
  return(out);
}
calculatePotentialAuxSegments <- function(x, internalDeltas){
  out <- numeric(length(x));
  if(length(out) == 0){
    return(out);
  }
  m <- max(x, na.rm=TRUE)
  # Sometimes m is not meaningful, because x has length 0
  w <- x == m;
  if(sum(w)==1){ ## w is the only maximum element. 
    m2 <- max(x[!w], na.rm=TRUE);
    out[w] <- m - m2;
  }
  return(out);
}
calculateWinOpti <- function(x, w){
  maxInSet <- max(x[w]);
  maxOutsideSet <- max(x[!w]);
  if(maxInSet > maxOutsideSet){
    return(maxInSet - maxOutsideSet)
  } else {
    return(0)
  }
}

updatePotentialMatrix <- function(P,funAuxMat,indices,subProjs){
  st <- Sys.time();
  if(sum(indices)>1){
    P[indices, indices] <- apply(funAuxMat[indices,indices], 1, calculatePotentialAux, subProjs[indices]);
  } else {
    P[indices, indices] <- calculatePotentialAux(funAuxMat[indices,indices],subProjs[indices]);
  }
  et <- Sys.time();
  #  print(paste("Time spent updating potential matrix:", et-st));
  return(P);
}

updatePotentialMatrixSegments <- function(funAuxMat, P, funDeltas){
  if(NAIVENESS_LEVEL == "Pessimistic"){
    rowMaxs = apply(funAuxMat, 1, max);
    P = matrix(0, nrow = dim(funAuxMat)[1], ncol = dim(funAuxMat)[2]);
    for(m in 1:dim(funAuxMat)[1]){
      maxVal <- rowMaxs[m];
      if(maxVal > 0){
        w = funAuxMat[m,] == maxVal;
        wIndex <- which(w);
        if(SETTINGNO == 1){
          if(length(wIndex) == 1){
            if(wIndex == m | funDeltas[m] == 1){ #if diagonal element is maximum or row already selected
              maxVal2 <- max(funAuxMat[m,!w], na.rm=TRUE);
              P[i,wIndex] <- maxVal - maxVal2
            }
          }
        } else if(SETTINGNO == 2){
          if(funDeltas[m] == 0){
            P[m,m] <- funAuxMat[m,m];
          } else {
            if(length(wIndex) == 1){
              maxVal2 <- max(funAuxMat[m,!w], na.rm=TRUE);
              P[i,wIndex] <- maxVal - maxVal2
            }
          }
        } else if(SETTINGNO == 3){
          if(funDeltas[m] == 0){
            P[m,m] <- rowMaxs[m];
          } else {
            if(length(wIndex) == 1){
              maxVal2 <- max(funAuxMat[m,!w], na.rm=TRUE);
              P[i,wIndex] <- maxVal - maxVal2
            }
          }
        }
      }
    }
  }
  return(P);
}

updatePotentialVector <- function(P, pVec){
  if(NAIVENESS_LEVEL == "Pessimistic"){
    pVec <- colSums(P); 
  }
  return(pVec);
}

takeStep <- function(funDeltas, funAuxMat, P, pVec, relevantIndices, subProjs){
  
  ##Chooses the link(s) that provide the largest total utility across all links
  P <- updatePotentialMatrix(P, funAuxMat, relevantIndices, subProjs);
  pVec <- updatePotentialVector(P, pVec, relevantIndices);
  
  hist.data = hist(pVec, plot=F, breaks = 100)
  hist.data$counts = log10(hist.data$counts+1)
  plot(hist.data, xlab = "Predicted NPV [mio. DKK] (probably overestimated)", ylab = "Frequency (log10)");
  
  w <- logical(length(pVec));
  if("An old apprach" == "True") {
    w[pVec == max(pVec)] <- TRUE;
  } else {
    subPDeterminant <- which(pVec == max(pVec))[1];
    subP <- subProjs[subPDeterminant];
    ws <- subProjs == subP;
    w[ws] <- TRUE;
  }
  nSelected <- sum(as.numeric(w));
  maxPVec <- max(pVec, na.rm = TRUE);
  #m <- maxP;
  #while(maxP > 0  && (maxP - m) / maxP < threshold){
  #  wTemp <- P == m;
  #  P[wTemp] <- -Inf
  #  w[wTemp] <- TRUE;
  #  nSelected <- nSelected + sum(wTemp);
  #  m <- max(P);
  #}
  funConverged <- FALSE;
  if(any(funDeltas[w]==1)){
    print("Link already chosen :/");
    funConverged <- TRUE;
  } else if(maxPVec < 1){
    print("Potential no longer positive")
    funConverged <- TRUE;
  } else {
    pVec[w] <- -Inf;
    P[,w] <- -Inf;
  }
  #print(paste0("Chosen links: ", which(w)));
  funDeltas[w] <- 1;
  
  out <- list(deltas = funDeltas, converged = funConverged, P = P, pVec = pVec);
  return(out);
}

takeStepSegments <- function(funDeltas, funAuxMat, P, pVec){
  ##Chooses the segment that provides the largest total utility across all links
  P <- updatePotentialMatrixSegments(funAuxMat, P, funDeltas);
  pVec <- updatePotentialVector(P, pVec);
  out <- list(P = P, pVec = pVec);
  return(out)
}

createPotentialVectorBasedOnLinkBenefits <- function(funDeltas, linkBenefits_b){
  pVec <- linkBenefits_b;
  pVec[funDeltas==1] = 0;
  return(pVec);
}

takeNaiveStepSegments <- function(funDeltas, linkBenefits_b){
  ##Chooses the segment that provides the largest total utility across all links
  return(list(pVec = createPotentialVectorBasedOnLinkBenefits(funDeltas, linkBenefits_b)));
}

selectFromPotentialVector <- function(pVec, funDeltas, nToSelect, efficiencyLowerBound){
  funConverged <- FALSE;
  
  effVec <- pVec / newCosts;
  effVec[!is.finite(effVec)] <- 0;
  
  #hist.data = hist(effVec, plot=F, breaks = 100)
  #hist.data$counts = log10(hist.data$counts+1)
  #plot(hist.data, xlab = "B/C ratio", ylab = "Frequency (log10)");
  
  w <- logical(length(effVec));
  nSelected <- 0;
  while(nSelected < nToSelect){
    shorterEffVec = effVec[!w];
    maxPVec <- max(effVec[!w], na.rm = TRUE);
    if(maxPVec>efficiencyLowerBound){
      w[effVec == maxPVec] <- TRUE;
    } else {
      print(paste0("      Efficiency does not exceed lower bound: "  , maxPVec , " <= ", efficiencyLowerBound))
      if(efficiencyLowerBound > 0){
        break;
      }
      selectedOne <- which(!funDeltas & !w & statusDF$SegmentLength > 0)[1];
      if(!is.na(selectedOne)){
        w[selectedOne] <- TRUE; ## Adding the first none-selected...
      } else {
        break;
      }
    } 
    nSelected <- sum(w);
  }
  
  
  funConverged <- FALSE;
  if(nSelected == 0){
    #    print("No more benefit to achieve : (")
    funConverged <- TRUE;
  } else if(any(funDeltas[w]==1)){ ##For debugging
    print("Segment already chosen :/");
    print(which(funDeltas[w]==1))
    print(maxPVec);
    print(which(w));
    funConverged <- TRUE;
  } else {
    if(maxPVec < 1){
      #print(paste0("WARNING!   B/c ratio < 1: ", round(maxPVec,2)));
      #funConverged <- TRUE;
    }
  }
  #print(paste0("Chosen links: ", which(w)));
  out <- list(Candidates = w, converged = funConverged);
  return(out);
}

findOptimalCombination <- function(combinations, originalAuxMat, newCosts, evaluationFunction) {
  triv <- logical(N_segments);
  allBCs <- foreach(i=1:dim(combinations)[2]) %dopar% {
    w <- triv;
    w[combinations[,i]] <- TRUE;
    val <- sum(apply(originalAuxMat, 1, FUN=evaluationFunction, w));
    return(val/sum(newCosts[w]));
  }
  allBCs <- unlist(allBCs);
  maxBC <- max(allBCs);
  whichOne <- which(allBCs == maxBC)[1]; #For ties, take the first
  w <- triv;
  w[combinations[,whichOne]] <- TRUE;
  maxB <- sum(apply(originalAuxMat, 1, FUN=calculateWinOpti, w));
  out = list(bestIndex = whichOne, bestValue = maxBC, benefitOfBestValue = maxB);
  return(out)
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Run the whole damn thing #####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

ACTUAL_MAX_BUDGET <- sum(sapply(correctedConstructionCosts, function(x){ discount(x, DISCOUNT_FACTORS, 1, 0)})) +
  sum(sapply(correctedMaintenanceCosts, function(x){discount(x * EVALUATION_PERIOD, DISCOUNT_FACTORS, EVALUATION_PERIOD,0)}))


## Generating predetermined implementation orders
{
  set.seed(426992)
  coreSegments <- 1:N_segments;
  #"ByRandomSegment";
  randomSegmentOrder <- sample(coreSegments[newLengths > 0], replace = FALSE);
  #"ByLongestSegment";
  longestSegmentOrder <- coreSegments[newLengths > 0];
  longestSegmentOrder <- longestSegmentOrder[order(-newLengths[newLengths>0])]
  #"ByShortestSegment";
  shortestSegmentOrder <- rev(longestSegmentOrder);
  #"ByLongestRoute";
  sortedRouteIds <- routeStatusDF[routeStatusDF$RouteId > 0, ];
  sortedRouteIds <- sortedRouteIds[order(-sortedRouteIds$AppliedLength), "RouteId"];
  longestRouteOrder <- numeric(0);
  while(length(sortedRouteIds) > 0){
    thisRouteId <- sortedRouteIds[1];
    sortedSegmentIds <- statusDF[statusDF$RouteId == thisRouteId, ];
    sortedSegmentIds <- sortedSegmentIds[order(-sortedSegmentIds$SegmentLength), "SegmentId"];
    longestRouteOrder <- c(longestRouteOrder, sortedSegmentIds);
    sortedRouteIds <- sortedRouteIds[-1];
  }
  ## Including those that do not belong to any route
  longestRouteOrder <- c(longestRouteOrder, as.numeric(statusDF$SegmentId[statusDF$RouteId == -1 & statusDF$SegmentLength > 0])); 
  #"ByShortestRoute";
  shortestRouteOrder <- rev(longestRouteOrder);
  #"ByActual";
  actualOrder <- c(164, 165, 166, 154, 105, 54, 93, 103, 129, 131, 114, 112, 66, 65, 64, 37, 38, 39, 40,
                   38, 34, 35, 50, 51, 52, 53, 55, 57, 91, 88, 84, 82, 80, 175, 190, 49, 47, 43, 42, 144,
                   76, 122, 127, 127, 102, 101, 176, 171, 170, 70, 69, 68, 67, 111, 151)
  }


byOptimalStrategies <- c(TRUE,FALSE); ## Never FALSE before TRUE!
#byOptimalStrategies <- c(FALSE); ## Never FALSE before TRUE!
#byOptimalStrategies <- c(TRUE); ## Never FALSE before TRUE!
#byOptimalStrategies <- c(FALSE); ## Never FALSE before TRUE!


efficiencyLowerBounds <- c(0,0.5,0.8,0.9,1);
efficiencyLowerBounds <- c(0,1);
#efficiencyLowerBounds <- c(0);

segmentPoolCapacities <- c(10,20);
segmentPoolCapacities <- c(1);
segmentsAddedPerIterations<- c(1,2,5,10,20);
segmentsAddedPerIterations<- c(1,2,5,10);
segmentsAddedPerIterations<- c(1);
segmentChunkSizes <- c(1,2,3,4);
segmentChunkSizes <- c(1);
maintenanceFactors <- c(1,2,4);
maintenanceFactors <- c(1);


byWhatsOuter <- c("ByActual","ByRandomSegment","ByLongestSegment","ByShortestSegment","ByLongestRoute","ByShortestRoute")
#byWhatsOuter <- c("ByRandomSegment")#,"ByLongestSegment","ByShortestSegment","ByLongestRoute","ByShortestRoute")

#byWhatsOuter <- c("ByActual")

alternativeStartingPoints = c(FALSE,TRUE);
#alternativeStartingPoints = c(FALSE);
NAIVENESS_LEVELs = c("Naive","SemiNaive", "Pessimistic");
#NAIVENESS_LEVELs = c("SemiNaive");


SETTINGNOs <- c(1,3);
#SETTINGNOs <- c(1);

MAX_BUDGETs <- c(ACTUAL_MAX_BUDGET,500);



for(MAX_BUDGET in MAX_BUDGETs){
  for(byOptimalStrategy in byOptimalStrategies){
    for(NAIVENESS_LEVEL in NAIVENESS_LEVELs){
      for(alternativeStartingPoint in alternativeStartingPoints){
        for(SETTINGNO in SETTINGNOs){
          if(byOptimalStrategy){
            byWhats <- c("Doesn't matter");
            if(NAIVENESS_LEVEL != "Pessimistic" & SETTINGNO != SETTINGNOs[1]){
              next;
            }
          } else {
            byWhats <- byWhatsOuter;
            segmentPoolCapacities <- c(1);
            segmentsAddedPerIterations <- c(1);
            segmentChunkSizes <- c(1);
            efficiencyLowerBounds <- c(0); 
            if(SETTINGNO != SETTINGNOs[1]){
              next;
            }
            if(NAIVENESS_LEVEL != NAIVENESS_LEVELs[1]){
              next;
            }
          } 
          print(NAIVENESS_LEVEL)
          for(efficiencyLowerBound in efficiencyLowerBounds){
            ELBAsString <- gsub("[[:punct:]]", "", round(efficiencyLowerBound,2))
            for(maintenanceFactor in maintenanceFactors){
              for(byWhat in byWhats){
                if (byWhat == "ByActual"){
                  thisSegmentOrder <- actualOrder;
                } else if (byWhat == "ByRandomSegment"){
                  thisSegmentOrder <- randomSegmentOrder;
                } else if (byWhat == "ByLongestSegment"){
                  thisSegmentOrder <- longestSegmentOrder;
                } else if (byWhat == "ByShortestSegment"){
                  thisSegmentOrder <- shortestSegmentOrder;    
                } else if (byWhat == "ByLongestRoute"){
                  thisSegmentOrder <- longestRouteOrder;    
                } else if (byWhat == "ByShortestRoute"){
                  thisSegmentOrder <- shortestRouteOrder;    
                } else {
                  thisSegmentOrder <- numeric(0); ## Doesn't matter
                }
                if(alternativeStartingPoint){
                  thisSegmentOrder <- thisSegmentOrder[thisSegmentOrder %notin% aspSegments];
                }
                for(segmentPoolCapacity in segmentPoolCapacities){
                  for(segmentsAddedPerIteration in segmentsAddedPerIterations){
                    for(segmentChunkSize in segmentChunkSizes){
                      
                      SEGMENT_POOL_CAPACITY <- segmentPoolCapacity;
                      SEGMENTS_ADDED_PER_ITERATION <- segmentsAddedPerIteration;
                      SEGMENT_CHUNK_SIZE <- segmentChunkSize;
                      MAINTENANCE_FACTOR <- maintenanceFactor;
                      
                      if(SEGMENTS_ADDED_PER_ITERATION > SEGMENT_POOL_CAPACITY){
                        next; #Obviously skipping this one
                      }
                      if(SEGMENT_CHUNK_SIZE == 1){
                        if(SEGMENT_POOL_CAPACITY  == segmentPoolCapacities[1] & 
                           SEGMENTS_ADDED_PER_ITERATION == segmentsAddedPerIterations[1]){
                          SEGMENT_POOL_CAPACITY <- 1;
                          SEGMENTS_ADDED_PER_ITERATION <- 1;
                        } else {
                          next; # If not first time, skip it, because all the other times are redundant
                        }
                      }
                      
                      configu <- paste0( SEGMENT_POOL_CAPACITY,"_",SEGMENTS_ADDED_PER_ITERATION,"_",SEGMENT_CHUNK_SIZE, "_",
                                         MAINTENANCE_FACTOR);
                      
                      
                      if(!byOptimalStrategy){
                        configu <- paste0(byWhat,"_",MAINTENANCE_FACTOR)
                      } else {
                        configu <- paste0(configu, "_", ELBAsString);
                        configu <- paste0(configu,"_",NAIVENESS_LEVEL)
                        if(NAIVENESS_LEVEL == "Pessimistic"){
                          configu <- paste0(configu, SETTINGNO);
                        } 
                      }
                      
                      if(alternativeStartingPoint){
                        configu <- paste0(configu,"_ASP");
                      }
                      
                      if(MAX_BUDGET == 500){
                        configu <- paste0(configu,"_",MAX_BUDGET);
                      }
                      
                      
                      
                      fullStop <- FALSE;
                      
                      coreB <- segmentCorrs * CYCLIST_UNIT_PRICE_PER_MINUTE * AADT_FACTOR / 1e6; 
                      if(NAIVENESS_LEVEL == "Naive"){
                        coreBenefitVector <- globalNaiveBenefitVector * CYCLIST_UNIT_PRICE_PER_MINUTE * AADT_FACTOR / 1e6;
                      } else if(NAIVENESS_LEVEL == "SemiNaive"){
                        coreBenefitVector <- apply(coreB,1,sum); ## TODO
                      }
                      usedCost <- 0;
                      totalExpectedBenefit <- 0;
                      nextCost <- 0;
                      nowYear <- 0;
                      EPOCH_TIME_SPAN <- 2;
                      deltas <- numeric(N_segments);
                      
                      itCounter <- 0;
                      newIndices <- logical(N_segments);
                      nextAuxMat <- segmentCorrs
                      potentialVector <- numeric(N_segments)
                      if(alternativeStartingPoint){
                        deltas[aspSegments] <- 1;
                      }
                      freeSegments <- deltas;
                      
                      budgets <- seq(100,MAX_BUDGET,100);
                      #budgets <- c(5,10,15,20,25,30,40,50,75,100);
                      loopSt <- Sys.time();
                      if(SEGMENT_POOL_CAPACITY < SEGMENT_CHUNK_SIZE | SEGMENT_POOL_CAPACITY < SEGMENTS_ADDED_PER_ITERATION){
                        print("Illegal combination")
                        next;
                      } else {
                        thisOutputDir <- paste0(outputDir,"/Configuration_",configu)
                        dir.create(file.path(thisOutputDir), showWarnings = FALSE)
                        selectFilesInDir <- list.files(thisOutputDir, pattern = "Selected*");
                        for(f in selectFilesInDir){
                          file.remove(paste0(thisOutputDir,"/",f));
                        }
                        print("-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=")
                        print(paste0("=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~ ", configu, " =-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~"));
                        print("~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-~~=-")
                      }
                      for(budget in budgets){
                        if(fullStop){
                          print("Terminated")
                        } else {
                          if(budget %% 500 == 100){
                            print(cat(paste("It.", "time", "C", "Budget", "B","B/C", "Length", "#Segs", "Year","", sep="\t")))
                          }
                          
                          if(budget != budgets[1]){
                            nowYear <- nowYear + EPOCH_TIME_SPAN;
                          }
                          implementedYears <- EVALUATION_PERIOD - CONSTRUCTION_TIME - nowYear;
                          if(NAIVENESS_LEVEL %in% c("Naive", "SemiNaive")){
                            benefitVector <- coreBenefitVector
                            benefitVector <- sapply(benefitVector, function(x){(discount(x * implementedYears, DISCOUNT_FACTORS, implementedYears, CONSTRUCTION_TIME + nowYear))});
                          }
                          newCosts <- sapply(correctedConstructionCosts, function(x){ discount(x, DISCOUNT_FACTORS, CONSTRUCTION_TIME, nowYear)}) + #construction
                            sapply(correctedMaintenanceCosts, function(x){discount(x * implementedYears, DISCOUNT_FACTORS, implementedYears,CONSTRUCTION_TIME + nowYear)}) /
                            MAINTENANCE_FACTOR -  # maintenance
                            sapply(correctedConstructionCosts, function(x){ discount(x, DISCOUNT_FACTORS, 1, EVALUATION_PERIOD)}); #scrap value
                          newCosts <- newCosts * NAF;
                          nextAuxMat <- coreB;
                          nextAuxMat <- updateAuxMat(nextAuxMat, which(deltas==1))
                          nextAuxMat <- matrix(sapply(nextAuxMat, function(x){(discount(x * implementedYears, DISCOUNT_FACTORS, implementedYears, CONSTRUCTION_TIME + nowYear))}), 
                                               nrow=N_segments, ncol=N_segments);
                          originalAuxMat <- nextAuxMat;
                          potentialMatrix <- matrix(0, nrow=N_segments, ncol=N_segments); 
                          potentialMatrix <- updatePotentialMatrixSegments(nextAuxMat,potentialMatrix,numeric(N_segments))
                          originalPotentialMatrix <- potentialMatrix;
                          potentialVector <- updatePotentialVector(potentialMatrix,potentialVector);
                          originalPotentialVector <- potentialVector;
                          originalDeltas <- deltas;
                          
                          activeChunkSize <- SEGMENT_CHUNK_SIZE;
                          actuallyChosenSomethingThisIteration <- FALSE;
                          budgetExceeded <- FALSE;
                          segmentPool <- logical(N_segments); 
                        }
                        while(TRUE){
                          
                          if(fullStop){
                            break;
                          }
                          auxMat <- nextAuxMat;
                          
                          et <- Sys.time();
                          if(actuallyChosenSomethingThisIteration){
                            print(cat(paste(itCounter, round(et-loopSt,2), round(usedCost,1), budget, round(totalExpectedBenefit,1), round(totalExpectedBenefit/usedCost,2), round(calculateInvestedLength(deltas),1),  sum(deltas == 1),  nowYear ,"",  sep="\t")))
                            #print(which(deltas == 1));
                          }
                          
                          actuallyChosenSomethingThisIteration <- FALSE;
                          cantFindMoreCandidates <- FALSE;
                          if(byOptimalStrategy){
                            if(NAIVENESS_LEVEL == "Pessimistic"){
                              stepOut <- takeStepSegments(deltas, auxMat, potentialMatrix, potentialVector);
                              potentialMatrixTemp <- stepOut$P;
                            } else { #SemiNaive or Naive
                              stepOut <- takeNaiveStepSegments(deltas, benefitVector);
                            }
                            potentialVectorTemp <- stepOut$pVec;
                            selectOut <- selectFromPotentialVector(potentialVectorTemp, deltas, 
                                                                   SEGMENTS_ADDED_PER_ITERATION,
                                                                   efficiencyLowerBound);
                            segmentCandidates <- selectOut$Candidates;
                            cantFindMoreCandidates <- selectOut$converged
                            
                            segmentPool[segmentCandidates] <- TRUE;  
                            
                          } else { ## Using predetermined order
                            if(length(thisSegmentOrder) > 0){
                              segmentPool[thisSegmentOrder[1]] <- TRUE;
                            } else {
                              cantFindMoreCandidates <- TRUE;
                            }
                          }
                          
                          segmentsInPool <- sum(segmentPool);
                          if(segmentsInPool == 0){
                            print("Converged");
                            fullStop <- TRUE;
                            break;
                          } else if(segmentsInPool >= SEGMENT_POOL_CAPACITY | cantFindMoreCandidates){
                            deltas <- originalDeltas;
                            potentialVector <- originalPotentialVector;
                            
                            while(TRUE){
                              nextDeltas <- originalDeltas;
                              nToSelect <- activeChunkSize;
                              
                              if(segmentsInPool > nToSelect){
                                ## make triplets
                                combinations <- combn(which(segmentPool), nToSelect);
                                
                                outOptiFind <- findOptimalCombination(combinations, originalAuxMat, newCosts, calculateWinOpti);
                                bestIndex <- outOptiFind$bestIndex;
                                bestValue <- outOptiFind$bestValue;
                                benefitOfBestValue <- outOptiFind$benefitOfBestValue;
                                
                                newDeltas <- logical(N_segments);
                                newDeltas[combinations[,bestIndex]] <- TRUE;
                                
                              } else {
                                newDeltas <- segmentPool;
                                benefitOfBestValue <- sum(apply(originalAuxMat, 1, FUN=calculateWinOpti, newDeltas));
                              }
                              
                              nextCost <- usedCost + sum(newCosts[newDeltas]);
                              
                              if(nextCost <= budget){
                                usedCost <- nextCost;
                                break;
                              } else if( activeChunkSize == 1){
                                budgetExceeded <- TRUE;
                                break;
                              } # else try again with one fewer to select
                              activeChunkSize <- activeChunkSize - 1;
                            }
                            
                            
                            if(budgetExceeded){
                              #print("Budget exceeded")
                              break;
                            }
                            
                            nextDeltas <- originalDeltas + newDeltas;
                            
                            potentialMatrixTemp <- originalPotentialMatrix;
                            potentialMatrixTemp[,newDeltas] <- -Inf;
                            if(NAIVENESS_LEVEL %in% c("Naive","SemiNaive")){
                              potentialVectorTemp <- createPotentialVectorBasedOnLinkBenefits(nextDeltas,benefitVector); 
                            } else {
                              potentialVectorTemp <- updatePotentialVector(potentialMatrixTemp,originalPotentialVector);
                            }
                            
                            nextAuxMat <- updateAuxMat(originalAuxMat, newDeltas);
                            if(sum(newDeltas)>1){
                              rSums <- rowSums(originalAuxMat[,newDeltas]);
                            } else {
                              rSums <- originalAuxMat[,newDeltas];
                            }
                            
                            originalAuxMat <- nextAuxMat; 
                            originalPotentialVector <- potentialVectorTemp;
                            originalPotentialMatrix <- potentialMatrixTemp;
                            originalDeltas <- nextDeltas;
                            
                            totalExpectedBenefit <- totalExpectedBenefit + benefitOfBestValue;
                            
                            segmentPool <- logical(N_segments); 
                            actuallyChosenSomethingThisIteration <- TRUE;
                            
                            if(!byOptimalStrategy){
                              thisSegmentOrder <- thisSegmentOrder[-1];
                            }
                            
                          } else {
                            selectOut <- selectFromPotentialVector(potentialVectorTemp, deltas, 1, efficiencyLowerBound);
                            nextDeltas <- deltas + selectOut$Candidates
                            ### nextDeltas... and tempDeltas. How to do this?
                          }
                          
                          
                          if(budgetExceeded){
                            break;
                          }
                          
                          
                          
                          
                          newIndices <- nextDeltas != deltas;
                          potentialVectorTemp[newIndices] <- -Inf;
                          potentialMatrixTemp[,newIndices] <- -Inf;
                          
                          
                          potentialMatrix <- potentialMatrixTemp;
                          potentialVector <- potentialVectorTemp;
                          
                          
                          deltas <- nextDeltas;
                          
                          if(!any(deltas==0)){
                            print("")
                            print("All links selected")
                            print("")
                            break;
                          }
                          itCounter <- itCounter + 1;
                          
                          if(!actuallyChosenSomethingThisIteration){ ## otherwise done is separate block
                            nextAuxMat <- updateAuxMat(auxMat, newIndices);
                            if(sum(newIndices)>1){
                              rSums <- rowSums(auxMat[,newIndices]);
                            } else {
                              rSums <- auxMat[,newIndices];
                            }
                          }
                          #Feed this into next iteration, so that only the relevant rows/columns of auxMat are updated
                          
                        }
                        selectedSegments <- which(as.logical(deltas));
                        selectedLinks <- as.data.frame(E(g)$LinkID[E(g)$SegmentId %in% selectedSegments]);
                        colnames(selectedLinks) <- "LinkID";
                        data.table::fwrite(selectedLinks, file = paste0(thisOutputDir, "/SelectedLinks_Predetermined_", configu,
                                                                        "_Budget", budget,".csv"))
                        selectedSegments <- as.data.frame(selectedSegments)
                        colnames(selectedSegments) <- c("SegmentId");
                        data.table::fwrite(selectedSegments, file = paste0(thisOutputDir, "/SelectedSegments_Predetermined_", configu,
                                                                           "_Budget", budget,".csv"))
                        if(budget == budgets[1]){
                          if(alternativeStartingPoint){
                            selectedSegments <- aspSegments;
                            selectedLinks <- as.data.frame(E(g)$LinkID[E(g)$SegmentId %in% selectedSegments]);
                            colnames(selectedLinks) <- "LinkID";
                            selectedSegments <- as.data.frame(selectedSegments)
                            colnames(selectedSegments) <- c("SegmentId");
                            data.table::fwrite(selectedLinks, file = paste0(thisOutputDir, "/SelectedLinks_Predetermined_Budget0.csv"))
                            data.table::fwrite(selectedSegments, file = paste0(thisOutputDir, "/SelectedSegments_Predetermined_Budget0.csv"))
                          } else {
                            data.table::fwrite(subset(selectedLinks,FALSE), file = paste0(thisOutputDir, "/SelectedLinks_Predetermined_Budget0.csv"))
                            data.table::fwrite(subset(selectedSegments, FALSE), file = paste0(thisOutputDir, "/SelectedSegments_Predetermined_Budget0.csv"))
                          }
                        }
                        if(fullStop){
                          break;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## True evaluation of plans ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

ACTUAL_MAX_BUDGET <- sum(sapply(correctedConstructionCosts, function(x){ discount(x, DISCOUNT_FACTORS, 1, 0)})) +
  sum(sapply(correctedMaintenanceCosts, function(x){discount(x * EVALUATION_PERIOD, DISCOUNT_FACTORS, EVALUATION_PERIOD,0)}));


alternativeStartingPoints <- c(FALSE,TRUE);
#alternativeStartingPoints <- c(TRUE);



overwriteExistingEvaluation <- TRUE;
byOptimalStrategies <- c(TRUE,FALSE);
byOptimalStrategies <- c(TRUE);
#byOptimalStrategies <- c(FALSE);


NAIVENESS_LEVELs <- c("Pessimistic", "SemiNaive", "Naive") 
NAIVENESS_LEVELs <- c("Naive","SemiNaive") 
NAIVENESS_LEVELs <- c("Pessimistic") 

SETTINGNOs <- c(3,1);
#SETTINGNOs <- c(3);

efficiencyLowerBounds <- c(1,0);

MAX_BUDGETs <- c(ACTUAL_MAX_BUDGET,500);

for(MAX_BUDGET in MAX_BUDGETs){
  for(alternativeStartingPoint in alternativeStartingPoints){
    for(byOptimalStrategy in byOptimalStrategies){
      for(NAIVENESS_LEVEL in NAIVENESS_LEVELs){
        for(SETTINGNO in SETTINGNOs){
          if(byOptimalStrategy){
            if(NAIVENESS_LEVEL != "Pessimistic" & SETTINGNO != SETTINGNOs[1]){
              next;
            }
          } else if(!(NAIVENESS_LEVEL == NAIVENESS_LEVELs[1] & SETTINGNO == SETTINGNOs[1])){
            next;
          }
          efficiencyLowerBoundsOuter <- c(0,0.5,0.8,0.9,1);
          efficiencyLowerBoundsOuter <- c(0,1);
          segmentPoolCapacities <- c(10,20);
          segmentPoolCapacities <- c(1);
          segmentsAddedPerIterations <- c(1,2,5,10,20);
          segmentsAddedPerIterations<- c(1,5);
          segmentsAddedPerIterations<- c(5);
          segmentsAddedPerIterations<- c(1);
          segmentChunkSizes <- c(1,2,3,4);
          segmentChunkSizes <- c(1,4);
          segmentChunkSizes <- c(1,4);
          segmentChunkSizes <- c(1);
          maintenanceFactors <- c(1,2,4);
          maintenanceFactors <- c(1);
          byWhatsOuter <- c("ByActual","ByRandomSegment","ByLongestSegment","ByShortestSegment","ByLongestRoute","ByShortestRoute")
          #byWhatsOuter <- c("ByActual");
          
          EPOCH_TIME_SPAN <- 2;
          budgets <- seq(100,sum(sapply(correctedConstructionCosts, function(x){ discount(x, DISCOUNT_FACTORS, 1, 0)})) +
                           sum(sapply(correctedMaintenanceCosts, function(x){discount(x * EVALUATION_PERIOD, DISCOUNT_FACTORS, EVALUATION_PERIOD,0)}))
                         ,100);
          ZeroBenefits = matrix(0,nrow = N_centroids, ncol = N_centroids);
          
          
          bikeTypes <- c(1,2,3);
          speedTypes <- c(1,2,3);
          totalOD <- matrix(0,nrow=length(CentroidNodes$NodeId), ncol = length(CentroidNodes$NodeId));
          totalOD_B <- matrix(0,nrow=length(CentroidNodes$NodeId), ncol = length(CentroidNodes$NodeId));
          totalOD_B_OD <- matrix(0,nrow=length(CentroidNodes$NodeId), ncol = length(CentroidNodes$NodeId));
          ODData <- as.data.frame(data.table::fread("O:/Public/4233-82676-BIKELONGER-persondata/CBA/OD_BikeType_SpeedType/SC01_MixBase1_MixSc1.csv"));
          ODData <- subset(ODData, !is.na(ODData$value));
          
          
          if(byOptimalStrategy){
            byWhats <- c("Doesn't matter");
            efficiencyLowerBounds <- efficiencyLowerBoundsOuter;
          } else {
            byWhats <- byWhatsOuter;
            segmentPoolCapacities <- c(1);
            segmentsAddedPerIterations <- c(1);
            segmentChunkSizes <- c(1);
            efficiencyLowerBounds <- c(0);
          } 
          
          for(efficiencyLowerBound in efficiencyLowerBounds){
            ELBAsString <- gsub("[[:punct:]]", "", round(efficiencyLowerBound,2))
            for(maintenanceFactor in maintenanceFactors){
              for(byWhat in byWhats){
                for(segmentPoolCapacity in segmentPoolCapacities){
                  for(segmentsAddedPerIteration in segmentsAddedPerIterations){
                    for(segmentChunkSize in segmentChunkSizes){
                      
                      SEGMENT_POOL_CAPACITY <- segmentPoolCapacity;
                      SEGMENTS_ADDED_PER_ITERATION <- segmentsAddedPerIteration;
                      SEGMENT_CHUNK_SIZE <- segmentChunkSize;
                      MAINTENANCE_FACTOR <- maintenanceFactor;
                      
                      if(SEGMENTS_ADDED_PER_ITERATION > SEGMENT_POOL_CAPACITY){
                        next; #Obviously skipping this one
                      }
                      if(SEGMENT_CHUNK_SIZE == 1){
                        if(SEGMENT_POOL_CAPACITY  == segmentPoolCapacities[1] & 
                           SEGMENTS_ADDED_PER_ITERATION == segmentsAddedPerIterations[1]){
                          SEGMENT_POOL_CAPACITY <- 1;
                          SEGMENTS_ADDED_PER_ITERATION <- 1;
                        } else {
                          next; # If not first time, skip it, because all the other times are redundant
                        }
                      }
                      
                      
                      
                      configu <- paste0( SEGMENT_POOL_CAPACITY,"_",SEGMENTS_ADDED_PER_ITERATION,"_",SEGMENT_CHUNK_SIZE, "_",
                                         MAINTENANCE_FACTOR);
                      if(!byOptimalStrategy){
                        configu <- paste0(byWhat,"_",MAINTENANCE_FACTOR)
                      } else {
                        configu <- paste0(configu, "_", ELBAsString);
                        configu <- paste0(configu,"_",NAIVENESS_LEVEL)
                        if(NAIVENESS_LEVEL == "Pessimistic"){
                          configu <- paste0(configu, SETTINGNO);
                        } 
                      }
                      
                      if(alternativeStartingPoint){
                        configu <- paste0(configu, "_ASP");
                      }
                      if(MAX_BUDGET == 500){
                        configu <- paste0(configu,"_",MAX_BUDGET);
                      }
                      
                      
                      thisOutputDir <- paste0(outputDir,"/Configuration_",configu)
                      freeLinks <- data.table::fread(file = paste0(thisOutputDir, "/SelectedLinks_Predetermined_Budget0.csv"))$LinkID
                      
                      
                      if(!overwriteExistingEvaluation & file.exists(paste0(thisOutputDir, "/EvaluationSummary_",configu,".csv"))){
                        next; 
                      }
                      
                      prevCurrentSegments <- logical(N_segments);
                      container = list("Approach" = character(0), "Budget" = numeric(0), "Cost" = numeric(0), "ConsumerSurplus" = numeric(0), "Gain" = numeric(0));
                      
                      constructionCostsPerYear <- numeric(EVALUATION_PERIOD+1)
                      maintenanceCostsPerYear <- numeric(EVALUATION_PERIOD+1);
                      futureScrapValuePerYear <- numeric(EVALUATION_PERIOD+1);
                      consumerSurplusPerYear <- numeric(EVALUATION_PERIOD+1);
                      totalTravelTimeSavings <- 0;
                      for(thisYear in 0:EVALUATION_PERIOD){
                        budget <- budgets[ceiling((thisYear+1)/EPOCH_TIME_SPAN)];
                        yearIndex <- thisYear + 1;
                        
                        print(paste0(configu, ": Year ", thisYear));
                        
                        if(thisYear %% EPOCH_TIME_SPAN == 0){ ##In the first year, and every 2nd year after this, we select (construct) new links
                          fname <- paste0(thisOutputDir, "/SelectedLinks_Predetermined_",configu,"_Budget", budget,".csv");
                          if(file.exists(fname)){
                            print(fname)
                            chosenLinks <- data.table::fread(file = fname)$LinkID
                            chosenLinks <- chosenLinks[chosenLinks %notin% freeLinks];
                          }
                          chosenLinksBool <- E(g)$LinkID %in% chosenLinks
                          
                          currentSegments <- logical(N_segments);
                          currentSegments[unique(E(g)$SegmentId[chosenLinksBool])] <- TRUE;
                          
                          newlyAddedSegments <- prevCurrentSegments != currentSegments;
                          
                          benefitsOfThisEpochNotYetCalculated <- sum(newlyAddedSegments) > 0;
                        } else {
                          newlyAddedSegments <- logical(N_segments);
                        }
                        
                        
                        
                        ## In the year after finishing construction, we calculate benefits 
                        if(thisYear %% EPOCH_TIME_SPAN == (CONSTRUCTION_TIME %% EPOCH_TIME_SPAN) ){ 
                          
                          if( benefitsOfThisEpochNotYetCalculated){ #Otherwise investments were the same as the previous time around
                            benefitsOfThisEpochNotYetCalculated <- FALSE;
                            
                            chosenType <- ifelse(chosenLinksBool, E(g)$InfraType_Plan_30_45, E(g)$InfraType_Basis);
                            
                            travelTimeSavings <- numeric(length(bikeTypes)*length(speedTypes));
                            typeIndex <- 0;
                            for(bikeType in bikeTypes){
                              for(speedType in speedTypes){
                                typeIndex <- typeIndex + 1;
                                
                                OD <- createODMatrix(subset(ODData, BikeType==bikeType & SpeedType == speedType));
                                
                                basisFileName <- paste0(outputDir,"/gcs_Basis_",bikeType,"_",speedType,".Rda");
                                while(!file.exists(basisFileName)){
                                  print("Temporary connection problems. Waiting 20 seconds to hopefully overcome the issue....")
                                  Sys.sleep(20); 
                                }
                                gcs_Basis <- readRDS(file=basisFileName) 
                                
                                vMax_1 <- determineVMax(1,bikeType, speedType);
                                vMax_2 <- determineVMax(2,bikeType, speedType);
                                vMax_3 <- determineVMax(3,bikeType, speedType);
                                
                                VMax_1 <- determineVMax(1,bikeType, speedType);
                                VMax_2 <- determineVMax(2,bikeType, speedType);
                                VMax_3 <- determineVMax(3,bikeType, speedType);
                                
                                print(paste("Running scenario for",bikeType,speedType));
                                # All scenario
                                E(g)$weight <- 1000000;
                                E(g)$weight <- ifelse(chosenType == 1, E(g)$Length / (vMax_1 / 3.6) + E(g)$IntersectionDelay, E(g)$weight);
                                E(g)$weight <- ifelse(chosenType == 2, E(g)$Length / (vMax_2 / 3.6) + E(g)$IntersectionDelay, E(g)$weight);
                                E(g)$weight <- ifelse(chosenType == 3, E(g)$Length / (vMax_3 / 3.6) + E(g)$IntersectionDelay, E(g)$weight);    
                                
                                out <- assignFlows(g, OD, ZeroBenefits, FALSE, FALSE);
                                gcs_All <- out$GCs / 60;
                                
                                ODLengths_All <- out$ODLengths;
                                ## Calculating differences
                                OD_B <- gcs_Basis - gcs_All;
                                
                                if(min(OD_B) < -1e-12){
                                  print(paste0("Negative CS :( ", sum(OD_B<0)))
                                  print(kill + to + kill)
                                }
                                
                                travelTimeSavings[typeIndex] <- sum(OD*OD_B);
                                
                                ODBDist <- numeric(61)
                                for(i in 1:length(ODBDist)){
                                  ODBDist[i] <- sum(OD[OD_B >=(i-1)/6 & OD_B<i/6]/sum(OD) )
                                }
                                plot((0:60)/6, cumsum(ODBDist), ylim=c(ODBDist[1]*0.99,1), typ = "l", xlim = c(0,5), lwd=2,
                                     main = paste0(bikeType, "_", speedType, ":   Avg. min. saved: ", round(sum(OD_B*OD)/sum(OD),3)), 
                                     xlab = "Minutes saved", ylab = "Proportion")
                                abline(v=sum(OD_B*OD)/sum(OD), col = 8, lty = 2)
                                abline(h=1, col =8)
                                
                                # ## Check if distance distributions are correct.....
                                # {
                                #   lengthDist <- numeric(60);
                                #   for(i in 1:length(lengthDist)){
                                #     lengthDist[i] <- sum(OD[ODLengths_Basis>=(i-1)*1000 & ODLengths_Basis<i*1000] / sum(OD))
                                #   }
                                #   plot(1:length(lengthDist),cumsum(lengthDist), xlim = c(0,20),
                                #        main = paste0("CS for ", bikeType, "_", speedType, ": " ,
                                #                      round(discount(travelTimeSavings[typeIndex] / 1e6 * CYCLIST_UNIT_PRICE_PER_MINUTE * AADT_FACTOR,
                                #                        DISCOUNT_FACTORS,1,thisYear),3)))
                                #   #barplot(lengthDist, names.arg = 1:length(lengthDist), las = 2)
                                #   
                                #   lengthDist_All <- numeric(60);
                                #   for(i in 1:length(lengthDist_All)){
                                #     lengthDist_All[i] <- sum(OD[ODLengths_All>=(i-1)*1000 & ODLengths_All<i*1000] / sum(OD))
                                #   }
                                #   points(1:length(lengthDist_All),cumsum(lengthDist_All), col = 2)
                                # }
                              }
                            }    
                            totalTravelTimeSavings <- sum(travelTimeSavings);
                          }
                        }
                        
                        # Taken from the detailed calculations
                        thisConsumerSurplus <- discount(totalTravelTimeSavings / 1e6 * CYCLIST_UNIT_PRICE_PER_MINUTE * AADT_FACTOR,
                                                        DISCOUNT_FACTORS,1,thisYear);
                        #Newly added only
                        baseConstructionCosts <- sum(correctedConstructionCosts[newlyAddedSegments]);
                        thisConstructionCost <- discount(baseConstructionCosts, DISCOUNT_FACTORS, CONSTRUCTION_TIME, thisYear)
                        thisScrap <- discount(baseConstructionCosts, DISCOUNT_FACTORS, 1, EVALUATION_PERIOD);
                        
                        # Applies to all links currently implemented
                        thisMaintenanceCost <- discount(sum(correctedMaintenanceCosts[currentSegments]), DISCOUNT_FACTORS,1, thisYear) /
                          MAINTENANCE_FACTOR;
                        
                        
                        # Adding to arrays
                        consumerSurplusPerYear[yearIndex] <- thisConsumerSurplus;
                        maintenanceCostsPerYear[yearIndex] <- thisMaintenanceCost;
                        constructionCostsPerYear[yearIndex] <- thisConstructionCost;
                        futureScrapValuePerYear[yearIndex] <- thisScrap
                        
                        #Printing status
                        print(paste("B:",round(thisConsumerSurplus,2), "C:", round(thisConstructionCost+thisMaintenanceCost,2)));
                        print(paste0("NPV " , thisYear, ": ", round(sum(consumerSurplusPerYear[1:yearIndex]) - 
                                                                      sum(constructionCostsPerYear[1:yearIndex]) - sum(maintenanceCostsPerYear[1:yearIndex]),2)));
                        
                        #Preparing next iteration
                        prevCurrentSegments <- currentSegments; 
                      }
                      
                      realScrap <- c(rep(0,EVALUATION_PERIOD+1),sum(futureScrapValuePerYear));
                      consumerSurplusPerYear <- c(consumerSurplusPerYear,0);
                      maintenanceCostsPerYear <- c(maintenanceCostsPerYear,0);
                      constructionCostsPerYear <- c(constructionCostsPerYear,0);
                      futureScrapValuePerYear <- c(futureScrapValuePerYear,0);
                      discountedNetCashFlows <- consumerSurplusPerYear + realScrap - constructionCostsPerYear - maintenanceCostsPerYear;
                      NPVs <- cumsum(discountedNetCashFlows)
                      
                      print(paste("Scrap:",round(sum(realScrap),2)));
                      
                      evaluationSummary <- 
                        as.data.frame(cbind(rep(configu,EVALUATION_PERIOD+2),  0:(EVALUATION_PERIOD+1), consumerSurplusPerYear,
                                            constructionCostsPerYear, maintenanceCostsPerYear, futureScrapValuePerYear, realScrap,
                                            discountedNetCashFlows, NPVs));
                      colnames(evaluationSummary) <- 
                        c("Config","Year","ConsumerSurplus","ConstructionCosts","MaintenanceCosts","FutureScrapValue","ScrapValue",
                          "DiscountedNetCashFlow","NetPresentValue")
                      data.table::fwrite(evaluationSummary, file = paste0(thisOutputDir, "/EvaluationSummary_",configu,".csv"))
                      
                      print(paste0("Final result (Net Present Value (NPV)) of ", configu, " is: ", NPVs[EVALUATION_PERIOD+2]))
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}




####### Feasible solutions #####################

hejdu <- hejdu[order(hejdu)]
feasibleSolutions <- 1;
for(M1 in length(hejdu):0){
  feasibleSolutions <- feasibleSolutions + 1; #that of M1
  print(paste0(M1,": ", feasibleSolutions))
  remBudget1 <- 100 - hejdu[M1];
  if(M1 == 1){
    print(paste(M1,log10(feasibleSolutions), sep = ", "))
    next;
  }
  belowBudget <- which(hejdu[1:(M1-1)] <= remBudget1);
  if(length(belowBudget) == 0){
    print(paste(M1,log10(feasibleSolutions), sep = ", "))
    next;
  }
  largestM2 <- max(belowBudget);
  for(M2 in largestM2:1){
    feasibleSolutions <- feasibleSolutions + 1; # (M1,M2)
    remBudget2 <- remBudget1 - hejdu[M2];
    if(M2 == 1){
      print(paste(M1,M2,log10(feasibleSolutions), sep = ", "))
      next;
    }
    belowBudget <- which(hejdu[1:(M2-1)] <= remBudget2);
    if(length(belowBudget) == 0){
      print(paste(M1,M2,log10(feasibleSolutions), sep = ", "))
      next;
    }
    largestM3 <- max(belowBudget);
    for(M3 in largestM3:1){
      feasibleSolutions <- feasibleSolutions + 1; # (M1,M2,M3)
      remBudget3 <- remBudget2 - hejdu[M3];
      if(M3 == 1){
        print(paste(M1,M2,M3,log10(feasibleSolutions), sep = ", "))
        next;
      }
      belowBudget <- which(hejdu[1:(M3-1)] <= remBudget3);
      if(length(belowBudget) == 0){
        print(paste(M1,M2,log10(feasibleSolutions), sep = ", "))
        next;
      }
      largestM4 <- max(belowBudget);
      for(M4 in largestM4:1){
        feasibleSolutions <- feasibleSolutions + 1; # (M1,M2,M3,M4)
        remBudget4 <- remBudget3 - hejdu[M4];
        if(M4 == 1){
          print(paste(M1,M2,M3,M4,log10(feasibleSolutions), sep = ", "))
          next;
        }
        belowBudget <- which(hejdu[1:(M4-1)] <= remBudget4);
        if(length(belowBudget) == 0){
          print(paste(M1,M2,M3,M4,log10(feasibleSolutions), sep = ", "))
          next;
        }
        largestM5 <- max(belowBudget);
        for(M5 in largestM5:1){
          feasibleSolutions <- feasibleSolutions + 1; # (M1,M2,M3,M4)
          remBudget5 <- remBudget4 - hejdu[M5];
          if(M5 == 1){
            next;
          }
          belowBudget <- which(hejdu[1:(M5-1)] <= remBudget5);
          if(length(belowBudget) == 0){
            next;
          }
          largestM6 <- max(belowBudget);
          for(k in largestM6:1){
            sols <- combn(largestM6,k);
            for(col in 1:dim(sols)[2]){
              theSum = sum(hejdu[sols[,col]]);
              if(theSum <= 100){
                feasibleSolutions <- feasibleSolutions + 1;
              }
            }
          }
        }
        print(paste(M1,M2,M3,M4,log10(feasibleSolutions), sep = ", "));
      }
    }
  }
}
print(feasibleSolutions)



feasibleSolutions <- 0;
N <- 1000000;
highestPossible <- max(which(cumsum(hejdu)<100))
for(j in 1:highestPossible){
  successes <- 0;
  for(n in 1:N){
    theSample = sample(1:length(hejdu),j, replace = FALSE);
    if(sum(hejdu[theSample]) < 100){
      successes <- successes + 1;
    }
  }
  p <- successes / N;
  numberOf = choose(length(hejdu),j)
  feasibleSolutions <- feasibleSolutions + numberOf * p;
  print(paste(j,feasibleSolutions))
}
print(feasibleSolutions)



