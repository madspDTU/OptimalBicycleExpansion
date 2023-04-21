#### Initialization ######
NUMBER_OF_ALLOCATED_CORES <- 4;

library("data.table")
library("igraph")
library("bit64")
library("stringr")
library("lpSolve")
library("foreach")
library("doParallel")
registerDoParallel(cores=NUMBER_OF_ALLOCATED_CORES);

#Adjust output location here:
outputDir <- "C:/workAtHome/PYNT/Output"
if(!file.exists(outputDir)){
  dir.create(outputDir);
}

#Adjust network data location here (must include the files: Connectors.csv, Centroids.csv, Nodes.csv and Links.csv):
NetworkDataDirectory <- "C:/workAtHome/PYNT/NetworkTables"

`%notin%` <- Negate(`%in%`)

source("c:/workAtHome/PYNT/RScripts/BP_HelperFunctions.R") #Load helper functions. Adjust this to point at the helper functions file

#### CONSTANTS ################

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
CONSTRUCTION_TIME <- 1; ## [years]



### Import network ######

# NodesData <- as.data.frame(data.table::fread(paste0(NetworkDataDirectory, "/Nodes.csv"), dec=","))
# CentroidNodes <- subset(NodesData, ZoneID > -1, select = c("NodeId","ZoneID"));
# NodesData <- subset(NodesData, select=c("NodeId"));
# colnames(NodesData) <- c("name")
# N_centroids <- dim(CentroidNodes)[1];
# 
# 
# LinksData <- as.data.frame(data.table::fread(paste0(NetworkDataDirectory, "/Links_WithSegments.csv"), dec=","))
# N_segments <- max(LinksData$SegId);
# 
# LinksData$SeeminglyUnused <- ifelse(LinksData$Active_Basis == 0 & LinksData$Active_Kommende == 0 & 
#                                       LinksData$Active_Plan_18_30 == 0 & LinksData$Active_Plan_30_45 == 0,1, 0)
# unimplementedLinks <- LinksData$ID[LinksData$SeeminglyUnused == 1];
# unimplementedLinks <- as.data.frame(unimplementedLinks);
# colnames(unimplementedLinks) <- "LinkID";
# data.table::fwrite(unimplementedLinks, file = paste0(outputDir, "/unimplementedLinks.csv"))
# 
# LinksData <- subset(LinksData, SeeminglyUnused == 0);
# #LinksData$Active_Plan_30_45[LinksData$SeeminglyUnused == 1] <- 1;
# 
# LinksData$InfraType_Basis[LinksData$Active_Basis == 0] <- -1;
# LinksData$InfraType_Basis[LinksData$Active_Basis == 1 & LinksData$Super_Basis == 1] <- 3
# LinksData$InfraType_Basis[LinksData$Active_Basis == 1 & LinksData$Super_Basis == 0 & LinksData$FreeSpeed >= 22] <- 2
# LinksData$InfraType_Basis[LinksData$Active_Basis == 1 & LinksData$Super_Basis == 0 & LinksData$FreeSpeed < 22] <- 1
# 
# LinksData$InfraType_Kommende[LinksData$Active_Kommende == 0] <- -1;
# LinksData$InfraType_Kommende[LinksData$Active_Kommende == 1 & LinksData$Super_Kommende == 1] <- 3
# LinksData$InfraType_Kommende[LinksData$Active_Kommende == 1 & LinksData$Super_Kommende == 0 & LinksData$FreeSpeed >= 22] <- 2
# LinksData$InfraType_Kommende[LinksData$Active_Kommende == 1 & LinksData$Super_Kommende == 0 & LinksData$FreeSpeed < 22] <- 1
# 
# LinksData$InfraType_Plan_18_30[LinksData$Active_Plan_18_30 == 0] <- -1;
# LinksData$InfraType_Plan_18_30[LinksData$Active_Plan_18_30 == 1 & LinksData$Super_Plan_18_30 == 1] <- 3
# LinksData$InfraType_Plan_18_30[LinksData$Active_Plan_18_30 == 1 & LinksData$Super_Plan_18_30 == 0 & LinksData$FreeSpeed >= 22] <- 2
# LinksData$InfraType_Plan_18_30[LinksData$Active_Plan_18_30 == 1 & LinksData$Super_Plan_18_30 == 0 & LinksData$FreeSpeed < 22] <- 1
# 
# LinksData$InfraType_Plan_30_45[LinksData$Active_Plan_30_45 == 0] <- -1;
# LinksData$InfraType_Plan_30_45[LinksData$Active_Plan_30_45 == 1 & LinksData$Super_Plan_30_45 == 1] <- 3
# LinksData$InfraType_Plan_30_45[LinksData$Active_Plan_30_45 == 1 & LinksData$Super_Plan_30_45 == 0 & LinksData$FreeSpeed >= 22] <- 2
# LinksData$InfraType_Plan_30_45[LinksData$Active_Plan_30_45 == 1 & LinksData$Super_Plan_30_45 == 0 & LinksData$FreeSpeed < 22] <- 1
# 
# disappearingLinks <- LinksData$LinkID[LinksData$InfraType_Plan_30_45==-1 & LinksData$InfraType_Basis>0];
# 
# 
# LinksData <- subset(LinksData, select=c("FromNodeId","ToNodeId","Shape_Length","IntersectionDelayFor",
#                                         "InfraType_Basis","InfraType_Plan_30_45","ID","SegId","RouteId"))
# colnames(LinksData) <- c("from","to","Length","IntersectionDelay",
#                          "NormalInfraType","UpgradedInfraType","LinkID","SegmentId","RouteId")
# N_edges <- dim(LinksData)[1];
# 
# data.table::fwrite(LinksData, paste0(outputDir, "/LinksData.csv"));
# data.table::fwrite(NodesData, paste0(outputDir, "/NodesData.csv"));
# data.table::fwrite(CentroidsData, paste0(outputDir, "/CentroidsData.csv"));
# 




#Load network data
LinksData <- data.table::fread(paste0(outputDir,"/LinksData.csv"));
NodesData <- data.table::fread(paste0(outputDir,"/NodesData.csv"));
CentroidsData <- data.table::fread(paste0(outputDir,"/CentroidsData.csv"));

# Creates the graph "g"
g <- graph_from_data_frame(LinksData, directed = TRUE, vertices = NodesData)
components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
g <- igraph::induced_subgraph(g, vert_ids);


CentroidNodes$Name <- CentroidNodes$NodeId
for(i in 1:dim(CentroidNodes)[1]){
  CentroidNodes$NodeId[i]  <- which(V(g)$name == CentroidNodes$Name[i]) 
}
g <- delete_vertex_attr(g,"name")   # Deleting the vertex atribute "name". 


N_TravelerTypes <- 9;
allWeights = list();
for(travelerType in 1:N_TravelerTypes){
  bikeType = (travelerType-1) %/% 3 + 1;
  speedType = (travelerType-1) %% 3 + 1;
  print(paste(bikeType,speedType))
  
  vMax_1 <- determineMaximumSpeed(1,bikeType,speedType); # Maximum speed of cyclists used for the assignment (Shared)
  vMax_2 <- determineMaximumSpeed(2,bikeType,speedType); # Maximum speed of cyclists used for the assignment (Separated)
  vMax_3 <- determineMaximumSpeed(3,bikeType,speedType); # Maximum speed of cyclists used for the assignment (Super)
  
  allWeights[[paste0("Normal",travelerType)]] <- updateNetworkCosts(g, E(g)$NormalInfraType, vMax_1, vMax_2, vMax_3);
  allWeights[[paste0("Upgraded",travelerType)]] <- updateNetworkCosts(g, E(g)$UpgradedInfraType, vMax_1, vMax_2, vMax_3);
}
g <- delete_edge_attr(g,"NormalInfraType")   # Deleting the vertex atribute "name". 
g <- delete_edge_attr(g,"UpgradedInfraType")   # Deleting the vertex atribute "name". 




# ODData <- as.data.frame(data.table::fread("O:/Public/4233-82676-BIKELONGER-persondata/CBA/OD_BikeType_SpeedType/SC01_MixBase1_MixSc1.csv"));
# ODData <- subset(ODData, !is.na(ODData$value));
# ODData$TravelerType <- 0;
# allODs = list();
# for(travelerType in 1:N_TravelerTypes){
#   bikeType = (travelerType-1) %/% 3 + 1;
#   speedType = (travelerType-1) %% 3 + 1;
#   ODData$TravelerType[ODData$SpeedType == speedType & ODData$BikeType == bikeType] <- travelerType;
# }
# data.table::fwrite(ODData[,c("FromZoneID","ToZoneID","TravelerType","value")], paste0(outputDir,"/ODData.csv"));


ODData <- data.table::fread(paste0(NetworkDataDirectory,"/ODData.csv"));
for(travelerType in 1:N_TravelerTypes){
  print(travelerType);
  allODs[[travelerType]] <- createODMatrix(subset(ODData, TravelerType == travelerType)); #OD for this traveller type
}



freeLinks <- E(g)$LinkID[E(g)$SegmentId %in% startingPoint];
globalNaiveBenefitVector <- numeric(N_segments);
all_gcs_Basis <- list();
for(travelerType in 1:N_TravelerTypes){
  print(paste("Traveler type",travelerType))
  
  print("Running basis scenario")
  E(g)$weight <- allWeights[[paste0("Normal",travelerType)]]
  out <- assignFlows(g, FALSE);
  gcs_Basis <- out$GCs / 60 * CYCLIST_UNIT_PRICE_PER_MINUTE; #Converting from seconds to minutes, and from minutes to DKK
  all_gcs_Basis[[travelerType]] <- gcs_Basis;
  
  print("Running full scenario")
  E(g)$weight <- allWeights[[paste0("Upgraded",travelerType)]]
  out <-  assignFlows(g, TRUE);
  gcs_All <- out$GCs / 60 * CYCLIST_UNIT_PRICE_PER_MINUTE; #Converting from seconds to minutes, and from minutes to DKK
  uEffect <- out$uEffect;
  
  ## Calculating cost differences
  Delta_gcs <- gcs_Basis - gcs_All;
  
  ## Multiplying by demand consumer surpluses
  thisCS <- Delta_gcs * allODs[[travelerType]];
  
  ## Calculating the consumer surplus attributed to each segment
  uEffect_CS <- numeric(N_segments);
  for(i in 1:N_segments){
    uEffect_CS[i] = sum(thisCS * uEffect[,,i]);
  }
  
  #Summing across traveler types
  globalNaiveBenefitVector <- globalNaiveBenefitVector + uEffect_CS;
}
saveRDS(globalNaiveBenefitVector, paste0(outputDir,"/BenefitVector.Rda"));
globalNaiveBenefitVector <- readRDS(file=paste0(outputDir,"/BenefitVector.Rda"))





lengths <- E(g)$Length / 1000;
costs <- lengths * (CONSTRUCTION_COST_PER_KM + EVALUATION_PERIOD * MAINTENANCE_COST_PER_KM_PER_YEAR);
routeIds <- numeric(N_segments) - 1

#Calculating lengths
newLengths <- numeric(N_segments);
for(i in 1:N_segments){
  segmentSelector <- E(g)$SegmentId == i;
  eids <- E(g)[segmentSelector];
  gsub <- subgraph.edges(g,eids=eids)
  if(length(E(gsub))>0){
    newLengths[i] <- farthest_vertices(gsub, directed=FALSE, weights=E(gsub)$Length)$distance/1000;
    routeIds[i] <-  unique(E(g)$RouteId[segmentSelector])[1];
  }
  print(paste("segment",i,newLengths[i]))
  print(components(gsub)$csize)
}



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
routeCostsForGIS = statusDF[statusDF$RouteId >=0,c("RouteId","ConstructionCostPerKm","MaintenanceCostPerKm")];
routeCostsForGIS = unique(routeCostsForGIS);
routeCostsForGIS <- routeCostsForGIS[order(routeCostsForGIS$RouteId),];
data.table::fwrite(routeCostsForGIS, file = "M:/GeoCBA/RouteCostsForGIS.csv", 
                   col.names = TRUE, row.names = FALSE, sep = ";", dec = ",");


statusDF$CorrectedConstructionCosts[statusDF$RouteId < 0] <- statusDF$SegmentLength[statusDF$RouteId < 0] * CONSTRUCTION_COST_PER_KM;
statusDF$CorrectedMaintenanceCosts[statusDF$RouteId < 0] <- statusDF$SegmentLength[statusDF$RouteId < 0] * MAINTENANCE_COST_PER_KM_PER_YEAR;


statusDF <- statusDF[,c("SegmentId","SegmentLength","CorrectedConstructionCosts","CorrectedMaintenanceCosts")];
colnames(statusDF) <- c("SegmentId","Length","ConstructionCosts","MaintenanceCosts");
statusDF <- statusDF[order(statusDF$SegmentId),]
data.table::fwrite(statusDF, paste0(outputDir, "/StatusDF.csv"));



statusDF <- data.table::fread(paste0(outputDir, "/StatusDF.csv"))
correctedConstructionCosts <- statusDF$ConstructionCosts
correctedMaintenanceCosts <- statusDF$MaintenanceCosts



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

slOri <-  globalNaiveBenefitVector * AADT_FACTOR / 1e6; #Annual benefits in million DKK
#slOri_RSP <-  apply(coreB_RSP,1,sum);

mlOri <- correctedMaintenanceCosts * NAF;
clOri <- correctedConstructionCosts * NAF;
NVar <- length(clOri);

constraintDirections = c("<=",rep(">=",NVar));


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## BLP part #####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


fetchOrder <- function(solutionType){
  theOrder <- c();
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
  return(theOrder);
}


for(scenarioType in c("RandomFirst","Short","Long")){
  
  theMainOutputDir <- paste0(outputDir,"/", scenarioType)
  if(!file.exists(theMainOutputDir)){
    dir.create(theMainOutputDir);
  }
  
  if(scenarioType == "Short"){
    lastYearOfInvestment <- 5;
    startingBudget <- 46.2; 
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
    
    if((scenarioType != "Short" & solutionType %in% c("ByActualOrder","ByActualRealisedOrder")) | 
       (scenarioType == "Short" & solutionType == "LPGreedyStop")){
      next;
    } 
    
    theOutputDir <- paste0(theMainOutputDir,"/",solutionType);
    if(!file.exists(theOutputDir)){
      dir.create(file.path(theMainOutputDir, solutionType))
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
    
    freeSegments = c();
    freeLinks <- as.data.frame(E(g)$LinkID[E(g)$SegmentId %in% freeSegments]);
    colnames(freeLinks) <- "LinkID"
    data.table::fwrite(freeLinks, file = paste0(theOutputDir, "/SelectedLinks", solutionType, "_T0.csv"))
    
    theOrder <- fetchOrder(solutionType); # In case of predetermined orders;
    if(scenarioType == "RandomFirst"){
      theOtherSolutionType <- solutionType;
      theOtherOrder <- theOrder[theOrder %notin% randomStartingPoint];
      theOrder <- randomStartingPoint;
      solutionType <- "ByRandomSegment"
    }
    
    
    ml <- mlOri;
    cl <- clOri;
    sl <- slOri;
    
    ml[freeSegments] <- 0;
    cl[freeSegments] <- 0;
    sl[freeSegments] <- 0;
    prevSolution[freeSegments] <- 1;  #We optimize wrt u and not DeltaU in the implemnetation!
    
    
    for(currentT in 1:finalT){
      
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
      CostsT <- CostsT * !prevSolution
      
      if(currentT <= lastYearOfInvestment){
        if(scenarioType == "RandomFirst" & (currentT == randomFirstT + 1)){
          theOrder <- theOtherOrder;
          solutionType <- theOtherSolutionType;
        }
        
        alreadySpent <- alreadySpent + sum(thisMl * prevSolution); #Add what is needed for maintenance in this time period 
        remainingBudget <- budget - alreadySpent;
        
        #If solving the optimization problem
        if(solutionType == "LP"){
          #Left-handside constraints. 
          conLHS <- rbind(CostsT,consistency)  #1st part: Construction cost, 2nd part: unit weights
          
          #ight-handside constraints
          consistencyRHS <- matrix(prevSolution,nrow=NVar) 
          conRHS <- rbind(remainingBudget,consistencyRHS); #1st part: Remaining budget, 2nd part: Already implemented (cannot downgrade)  
          
          lpSol <- lp("max", Zt, conLHS, constraintDirections, conRHS, int.vec = 1:NVar, all.bin = TRUE)
          sol <- lpSol$solution;
        } else { #For other strategies, keep investing until budget is reached.
          sol <- prevSolution;
          if(solutionType == "LPGreedy" | solutionType == "LPGreedyStop"){
            remainingCostsOfL = CostsT;
            #Ordering by highest benefit per construction cost. 
            theOrder <-  ((1:NVar)[sol==0 & newLengths > 0])[order(-Zt[sol==0 & newLengths > 0]/
                                                                     (remainingCostsOfL[sol==0 & newLengths > 0]))]
          }
          while(length(theOrder)>0){
            #Investing in segments as long as there is budget available or no more segments to select
            nextCandidate <- theOrder[1]; 
            if(solutionType == "LPGreedyStop"){
              if(Zt[nextCandidate] <= 0){ #If greedy w. stop, and segment not profitable, then break
                break;
              }
            }
            if(remainingBudget >= CostsT[nextCandidate]){ #If this segment can be afforded
              sol[nextCandidate] <- 1; #Invest in the segment
              remainingBudget <- remainingBudget - CostsT[nextCandidate]; #Adjust remaining budget
              #Remove segment from the order
              if(length(theOrder) == 1){
                theOrder <- c();
              } else {
                theOrder <- theOrder[2:length(theOrder)];
              }
            } else {
              break; #Cannot afford it.... (budget exceeded);
            }
          }
        }
      }
      
      deltaSol <- sol - prevSolution;
      objval <- sum(Zt * sol);
      
      alreadySpent <- alreadySpent  +  sum(deltaSol * CostsT); #Update already spent with construction costs of actions
      expectedTotalNPV <- expectedTotalNPV + objval
      expectedCurrentNPV <- expectedCurrentNPV - sum(CC * deltaSol) + sum(SV * deltaSol) +
        sum(thisSl * prevSolution) - oldMaintenance;
      ecNPVs[currentT] <- expectedCurrentNPV;
      ecCSs[currentT] <- sum(thisSl * prevSolution);
      
      remainingBudgets[currentT] <- remainingBudget;
      alreadySpents[currentT] <- alreadySpent;
      oldMaintenances[currentT] <- oldMaintenance;
      
      totalCC <- totalCC + sum(deltaSol * CC);
      prevSolution <- sol;
      
      #Print status
      print(c(currentT,sum(deltaSol), sum(sol), alreadySpent, budget, budget-alreadySpent, expectedCurrentNPV, expectedTotalNPV))
      
      
      #Outputting selected links
      selectedSegments <- which(as.logical(sol));
      selectedLinks <- as.data.frame(E(g)$LinkID[E(g)$SegmentId %in% selectedSegments]);
      colnames(selectedLinks) <- "LinkID";
      data.table::fwrite(selectedLinks, file = paste0(theOutputDir, "/SelectedLinks", solutionType,
                                                      "_T", currentT,".csv"))
      #Outputting selected segments
      selectedSegments <- as.data.frame(selectedSegments)
      colnames(selectedSegments) <- c("SegmentId");
      data.table::fwrite(selectedSegments, file = paste0(theOutputDir, "/SelectedSegments", solutionType,
                                                         "_T", currentT,".csv"))
    }
    
    thisRes <- data.frame(cbind(1:finalT, ecCSs, ecNPVs))
    colnames(thisRes) <- c("t","E_CS_t","E_CurrentNPV")
    data.table::fwrite(thisRes, file = paste0(theOutputDir, "/ExpectedEvaluation", solutionType,"_T", finalT,".csv"))
  }
  
  totalSV <- sum(SV[as.logical(prevSolution)]);
  totalCS <- expectedTotalNPV + alreadySpent - totalSV
  totalMC <- alreadySpent - totalCC;
  
  print("Predicted CC, SV, MC, CS, NPV, B/C")
  print(c(totalCC,totalSV, totalCS, expectedTotalNPV, totalCS /(totalCC+totalMC-totalSV) ))
  
}


#########################################
###### Real evaluation of solutions
#########################################

for( scenarioType in c("Long","RandomFirst","Short")){
  mainOutputDir <- paste0(outputDir,"/",scenarioType);
  for(solutionType in solutionTypes){
    theOutputDir <- paste0(mainOutputDir, "/", solutionType)
    
    if((scenarioType != "Short" & solutionType == "ByActualOrder") |
       (scenarioType == "Short" & solutionType == "LPGreedyStop")){
      next;
    } else if(scenarioType == "RandomFirst" & solutionType == "ByRandomSegment"){
      #Simply copy the one from the long run....
      file.copy("M:/GeoCBA/ModelOutputs/Long/ByRandomSegment/EvaluationByRandomSegment_T50.csv",
                "M:/GeoCBA/ModelOutputs/RandomFirst/ByRandomSegment/", overwrite = FALSE)
      next;
    }
    
    if(scenarioType != "RandomFirst"){
      freeLinks <- data.table::fread(paste0(theOutputDir, "/SelectedLinks",solutionType,"_T0.csv"))$LinkID;
    } else {
      freeLinks <- data.table::fread(paste0(theOutputDir, "/SelectedLinksByRandomSegment_T0.csv"))$LinkID;
    }
    freeSegments <- logical(N_segments);
    freeSegments[unique(E(g)$SegmentId[E(g)$LinkID %in% freeLinks])] <- TRUE
    
    prevCurrentSegments <- latestEvaluatedSolution <- freeSegments;
    cumulCC <- cumulMC <- cumulSV <- cumulCS <- currentNPV <- currentBCR <-  travelTimeSavings <- 0;
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
    
    for(currentT in 1:finalT){
      
      fname <- paste0(theOutputDir, "/SelectedLinks",solutionType,"_T", currentT,".csv");
      if(file.exists(fname)){
        print(paste(solutionType,currentT));
        chosenLinks <- data.table::fread(file = fname)$LinkID
        chosenLinks <- chosenLinks[chosenLinks %notin% freeLinks];
      }
      chosenLinksBool <- E(g)$LinkID %in% chosenLinks
      
      currentSegments <- logical(N_segments);
      currentSegments[unique(E(g)$SegmentId[chosenLinksBool])] <- TRUE;
      
      newlyAddedSegments <- prevCurrentSegments != currentSegments;
      
      MCPart <- sum(mlFun(currentT)[prevCurrentSegments]);  
      CCPart <- sum(CCFun(currentT)[newlyAddedSegments]);  
      SVPart <- sum(SVFun()[newlyAddedSegments]);
      
      #From travel time savings to actual, scaled consumer surplus
      rawConsumerSurplus <- travelTimeSavings / 1e6 * CYCLIST_UNIT_PRICE_PER_MINUTE * AADT_FACTOR 
      #Discounting the consumer surplus
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
      currentN <- sum(currentSegments); #Number of selected segments (cumulative)
      
      results <- rbind(results, c(currentT, CCPart, SVPart, MCPart, CSPart, NPV, cumulCC, cumulSV, cumulMC, cumulCS, 
                                  currentNPV, currentBCR, currentN))
      #Write the latest results
      data.table::fwrite(data.frame(results), paste0(theOutputDir, "/Evaluation",solutionType,"_T", currentT,".csv"),
                         quote = FALSE, row.names = FALSE);
      #Remove the previous results
      if(file.exists( paste0(theOutputDir, "/Evaluation",solutionType,"_T", currentT-1,".csv"))){
        file.remove(paste0(theOutputDir, "/Evaluation",solutionType,"_T", currentT-1,".csv"));
      }
      #Print results
      print(data.frame(results))
      
      #Performing traffic assignment for next time period (if needed, i.e. if any changes since last evaluation),
      #The output is sum of travel time savings. 
      if(any(latestEvaluatedSolution != currentSegments)){
        travelTimeSavings <- calculateTravelTimeSavings(chosenLinksBool);
        latestEvaluatedSolution <- currentSegments;
      }
      
      prevCurrentSegments <- currentSegments
      gc();
    }
  }
}










### The rest can be deleted.... ###



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
  theMainOutputDir <- paste0(outputDir,"/", budgetType);
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
    
    if(file.exists(paste0(theOutputDir, "/Evaluation", solutionType,"_T", finalT,".csv"))){
      thisExpRes <- data.frame(data.table::fread(paste0(theOutputDir, "/ExpectedEvaluation", solutionType,
                                                        "_T", finalT,".csv")))
      thisRes <- data.frame(data.table::fread(paste0(theOutputDir, "/Evaluation", solutionType,
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
    mainOutputDir <- paste0(outputDir,"/",scenarioType);
    for(solutionType in solutionTypes){
      theOutputDir <- paste0(mainOutputDir, "/", solutionType)
      if(file.exists(paste0(theOutputDir, "/Evaluation",solutionType,"_T", finalT,".csv"))){
        results <- data.frame(data.table::fread(paste0(theOutputDir, "/Evaluation",solutionType,"_T", finalT,".csv")));
        results[c("CC_t","SV_t","MC_t","CumulCC_t","CumulSV_t","CumulMC_t")] <- 
          results[c("CC_t","SV_t","MC_t","CumulCC_t","CumulSV_t","CumulMC_t")] * NAF
        results$NPV_t <- results$CS_t + results$SV_t - results$CC_t - results$MC_t;
        results$CurrentNPV <- results$CumulCS_t + results$CumulSV_t - results$CumulCC_t - results$CumulMC_t;
        results$CurrentBCR <- results$CumulCS_t /(-results$CumulSV_t + results$CumulCC_t + results$CumulMC_t);
        data.table::fwrite(results, paste0(theOutputDir, "/Evaluation",solutionType,"_T", finalT,"_WithNAF.csv"),
                           quote = FALSE, row.names = FALSE);
      }
    }
  }
}

if(FALSE){
  for( scenarioType in c("Long","Short")){
    mainOutputDir <- paste0(outputDir,"/",scenarioType);
    for(solutionType in solutionTypes){
      theOutputDir <- paste0(mainOutputDir, "/", solutionType)
      if(file.exists(paste0(theOutputDir, "/Evaluation",solutionType,"_T", finalT,".csv"))){
        results <- data.frame(data.table::fread(paste0(theOutputDir, "/Evaluation",solutionType,"_T", finalT,".csv")));
        if("N_t" %notin% colnames(results)){
          print(solutionType)
          Nts <- numeric(finalT);
          for(currentT in 1:finalT){
            selectedSegments <- data.frame(data.table::fread(paste0(theOutputDir,  "/SelectedSegments",solutionType,"_T", currentT,".csv")));
            Nts[currentT] <- dim(selectedSegments)[1]
          }
          results$N_t <- Nts;
          data.table::fwrite(results, paste0(theOutputDir,"/Evaluation",solutionType,"_T", finalT,".csv"),
                             quote = FALSE, row.names = FALSE);
        }
      }
    }
  }
}
