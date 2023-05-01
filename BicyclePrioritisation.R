#### Initialization ######
NUMBER_OF_ALLOCATED_CORES <- 4; #Number of cares in parallel computing

library("data.table")
library("igraph")
library("bit64")
library("stringr")
library("lpSolve")
library("foreach")
library("doParallel")
registerDoParallel(cores=NUMBER_OF_ALLOCATED_CORES); # For parallel processing


baseFolder <- "/path/to/somewhere"; #Adjust primary folder here <------
#Adjust output location here:
inputDir <- paste0(baseFolder,"/Input");

outputDir <- paste0(baseFolder,"/Output")
if(!file.exists(outputDir)){
  dir.create(outputDir);
}


`%notin%` <- Negate(`%in%`) #Function that is the inverse of the %in% function.

source(paste0(baseFolder,"/BP_HelperFunctions.R")) #Load helper functions. Adjust this to point at the helper functions file

#### CONSTANTS ################

CONSTRUCTION_COST_PER_KM <- 2.6   ## From Hallberg (2021). Only used as a last resort.
MAINTENANCE_COST_PER_KM_PER_YEAR <- 0.2; # From Hallberg (2021). Only used as a last resort. 
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



#Load network data ####
LinksData <- data.table::fread(paste0(inputDir,"/LinksData.csv"));
N_segments <- max(LinksData$SegmentId);
NodesData <- data.table::fread(paste0(inputDir,"/NodesData.csv"));
CentroidsData <- data.table::fread(paste0(inputDir,"/CentroidsData.csv"));
N_centroids <- dim(CentroidsData)[1];

# Creates the graph "g" ####
g <- graph_from_data_frame(LinksData, directed = TRUE, vertices = NodesData)
components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
g <- igraph::induced_subgraph(g, vert_ids);

CentroidsData$Name <- CentroidsData$NodeId
for(i in 1:dim(CentroidsData)[1]){
  CentroidsData$NodeId[i]  <- which(V(g)$name == CentroidsData$Name[i]) 
}
g <- delete_vertex_attr(g,"name")   # Deleting the vertex atribute "name". 

#Calculaiting link costs (weights) for all traveller types (normal and upgraded)
N_TravelerTypes <- 9;
allWeights = list();
for(travelerType in 1:N_TravelerTypes){
  bikeType = (travelerType-1) %/% 3 + 1;
  speedType = (travelerType-1) %% 3 + 1;
  
  vMax_1 <- determineMaximumSpeed(1,bikeType,speedType); # Maximum speed of cyclists used for the assignment (Shared)
  vMax_2 <- determineMaximumSpeed(2,bikeType,speedType); # Maximum speed of cyclists used for the assignment (Separated)
  vMax_3 <- determineMaximumSpeed(3,bikeType,speedType); # Maximum speed of cyclists used for the assignment (Super)
  
  allWeights[[paste0("Normal",travelerType)]] <- updateNetworkCosts(g, E(g)$NormalInfraType, vMax_1, vMax_2, vMax_3);
  allWeights[[paste0("Upgraded",travelerType)]] <- updateNetworkCosts(g, E(g)$UpgradedInfraType, vMax_1, vMax_2, vMax_3);
}
g <- delete_edge_attr(g,"NormalInfraType")   # Deleting the vertex atribute "name". 
g <- delete_edge_attr(g,"UpgradedInfraType")   # Deleting the vertex atribute "name". 


## Loading segments costs ####
statusDF <- as.data.frame(data.table::fread(paste0(inputDir, "/SegmentCosts.csv")));
correctedConstructionCosts <- statusDF$ConstructionCosts
correctedMaintenanceCosts <- statusDF$MaintenanceCosts
newLengths <- statusDF$Length;

routeStatusDF <- as.data.frame(data.table::fread(paste0(inputDir, "/RouteData.csv")))

## Reading demand (OD) data ####
ODData <- data.table::fread(paste0(inputDir,"/ODData.csv"));
allODs <- list();
for(travelerType in 1:N_TravelerTypes){
  allODs[[travelerType]] <- createODMatrix(subset(ODData, TravelerType == travelerType)); #OD for this traveller type
}


# Calculating the extreme scenarios (basis and all upgraded) + storing the costs + determining the linear effects. 
# You also read it directly if already previously calculated. (see 35 lines below)
csEffectVector <- numeric(N_segments);
all_gcs_Basis <- list();
for(travelerType in 1:N_TravelerTypes){
  print(paste("Traveler type:",travelerType))
  
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
  csEffectVector <- csEffectVector + uEffect_CS;
}
# Storing data, so it can be loaded quite for later uses. 
saveRDS(csEffectVector, paste0(outputDir,"/BenefitVector.Rda"));
saveRDS(all_gcs_Basis, paste0(outputDir,"/All_GCs_Basis.Rda"));
#Reading the stored copies
csEffectVector <- readRDS(file=paste0(outputDir,"/BenefitVector.Rda"))
all_gcs_Basis <- readRDS(paste0(outputDir,"/All_GCs_Basis.Rda"));



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# BLP setup # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

solutionTypes <- c("LP","LPGreedyStop","LPGreedy","ByActualOrder","ByRandomSegment","ByShortestSegment","ByShortestRoute","ByLongestSegment","ByLongestRoute")

## Generating predetermined implementation orders (for ad-hoc reference strategies)
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
    sortedSegmentIds <- sortedSegmentIds[order(-sortedSegmentIds$Length), "SegmentId"];
    longestRouteOrder <- c(longestRouteOrder, sortedSegmentIds);
    sortedRouteIds <- sortedRouteIds[-1];
  }
  ## Including those that do not belong to any route
  sortedNonRouteSegments <- statusDF[statusDF$RouteId == -1 & statusDF$Length > 0,];
  sortedNonRouteSegments <- sortedNonRouteSegments$SegmentId[order(-sortedNonRouteSegments$Length)]
  longestRouteOrder <- c(longestRouteOrder, sortedNonRouteSegments); 
  #"ByShortestRoute";
  shortestRouteOrder <- rev(longestRouteOrder);
  #"ByActual";
  actualOrder <- c(164, 165, 166, 154, 105, 54, 93, 103, 129, 131, 114, 112, 66, 65, 64, 37, 38, 39, 40,
                   28, 34, 35, 50, 51, 52, 53, 55, 57, 91, 88, 84, 82, 80, 175, 190, 49, 47, 43, 42, 144,
                   76, 122, 127, 128, 102, 101, 176, 171, 170, 70, 69, 68, 67, 111, 151)
}
set.seed(426992)
randomSegmentOrder <- sample(coreSegments[newLengths > 0], replace = FALSE);
randomStartingPoint <- c(53,108,84,3,63,30,138,179,97,147,132,204,129,198,50,111,56,13,157,167,42,133,26,58,193,54,12,115,165,200,103,98,174)

#Discount factor function
kappaFun <- function(currentT){
  if(currentT == 1){
    kappa <- 1;
  } else {
    kappa <- 1/sum(DISCOUNT_FACTORS[currentT-1])
  } 
  return(kappa);
}
# Function for sum of discount factors
kappaSumFun <- function(firstT,lastT){
  return(kappa <- sum(1/DISCOUNT_FACTORS[(firstT-1):(lastT-1)]))
}

#Benefit vector approximation function
sbFun <- function(currentT){
  return(kappaFun(currentT) * sb);
}
#Cumulative benefit vector approximation function
SbFun <- function(firstT, lastT){
  if(firstT > lastT){
    return(0*sb);
  }
  return(kappaSumFun(firstT,lastT) * sb);
}
#Maintenance cost vector function
mbFun <- function(currentT){
  return(kappaFun(currentT) * mb);
}
#Cumulative maintenance vector function
MbFun <- function(firstT, lastT){
  if(firstT > lastT){
    return(0*sb);
  }
  return(kappaSumFun(firstT,lastT) * mb);
}
#Construction cost vector function
CCFun <- function(currentT){
  return(kappaFun(currentT) * cb)
}
#Scrap value vector function
SVFun <- function(){
  return(kappaFun(finalT+1) * cb);
}

finalT <- 50;
randomFirstT <- 10; #In scenario with alternative staring strategy, how long does random order strategy last.

sbOri <-  csEffectVector * AADT_FACTOR / 1e6; #Annual benefits in million DKK
mbOri <- correctedMaintenanceCosts * NAF;
cbOri <- correctedConstructionCosts * NAF;
NVar <- length(cbOri);
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
  } 
  return(theOrder);
}


for(scenarioType in c("Short","Long","RandomFirst")){

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
    
    if((scenarioType != "Short" & solutionType =="ByActualOrder") | 
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
    print(c("Year","#Segments","Cumul#Segments","Spent","Budget","Unspent","NPV_t","NPV"));
    
    theOrder <- fetchOrder(solutionType); # In case of predetermined orders;
    if(scenarioType == "RandomFirst"){
      theOtherSolutionType <- solutionType;
      theOtherOrder <- theOrder[theOrder %notin% randomStartingPoint];
      theOrder <- randomStartingPoint;
      solutionType <- "ByRandomSegment"
    }
    
    
    mb <- mbOri;
    cb <- cbOri;
    sb <- sbOri;
     
    for(currentT in 1:finalT){
      
      budget <- budgetArr[currentT];
      
      thisSb <- sbFun(currentT);
      Sb <- SbFun(currentT+1,finalT)
      thisMb <- mbFun(currentT);
      Mb <- MbFun(currentT+1,finalT);
      CC <- CCFun(currentT);
      SV <- SVFun();
      
      oldMaintenance <- sum(thisMb * prevSolution)
      
      Zt <- Sb - CC - Mb + SV;
      Zt <- Zt*!prevSolution;
      
      CostsT <- matrix(CC, nrow = 1) #In future t's, also the M's until now.
      CostsT <- CostsT * !prevSolution
      
      if(currentT <= lastYearOfInvestment){
        if(scenarioType == "RandomFirst" & (currentT == randomFirstT + 1)){
          theOrder <- theOtherOrder;
          solutionType <- theOtherSolutionType;
        }
        alreadySpent <- alreadySpent + sum(thisMb * prevSolution); #Add what is needed for maintenance in this time period 
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
        sum(thisSb * prevSolution) - oldMaintenance;
      ecNPVs[currentT] <- expectedCurrentNPV;
      ecCSs[currentT] <- sum(thisSb * prevSolution);
      
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
  print(c(totalCC,totalSV, totalMC, totalCS, expectedTotalNPV, totalCS /(totalCC+totalMC-totalSV) ))
  
}


#########################################
###### Real evaluation of solutions
#########################################

for( scenarioType in c("Long","Short","RandomFirst")){
  mainOutputDir <- paste0(outputDir,"/",scenarioType);
  for(solutionType in solutionTypes){
    theOutputDir <- paste0(mainOutputDir, "/", solutionType)
    
    if((scenarioType != "Short" & solutionType == "ByActualOrder") |
       (scenarioType == "Short" & solutionType == "LPGreedyStop")){
      next;
    } else if(scenarioType == "RandomFirst" & solutionType == "ByRandomSegment"){
      #Simply copy the one from the long run, as they will be identical....
      file.copy(paste0(outputDir,"/Long/ByRandomSegment/EvaluationByRandomSegment_T50.csv"),
                paste0(outputDir,"/RandomFirst/ByRandomSegment/EvaluationByRandomSegment_T50.csv"), overwrite = FALSE)
      next;
    }

    prevCurrentSegments <- latestEvaluatedSolution <- logical(N_segments);
    cumulCC <- cumulMC <- cumulSV <- cumulCS <- currentNPV <- currentBCR <-  travelTimeSavings <- 0;
    results <- matrix(nrow=0, ncol = 13)
    colnames(results) <- c("t","CC_t","SV_t","MC_t","CS_t","NPV_t","CumulCC_t","CumulSV_t","CumulMC_t","CumulCS_t","CurrentNPV","CurrentBCR","N_t")
    
    mb <- mbOri;
    cb <- cbOri;
    sb <- sbOri;
    
    if(scenarioType == "RandomFirst"){
      startT <- randomFirstT + 1;
      
      fname <- paste0(theOutputDir, "/SelectedLinksByRandomSegment_T", randomFirstT,".csv");
      chosenLinks <- data.table::fread(file = fname)$LinkID
      chosenLinksBool <- E(g)$LinkID %in% chosenLinks
      currentSegments <- logical(N_segments);
      currentSegments[unique(E(g)$SegmentId[chosenLinksBool])] <- TRUE;
      prevCurrentSegments <- currentSegments;
      
      results <- as.data.frame(data.table::fread(paste0(outputDir, "/Long/ByRandomSegment/EvaluationByRandomSegment_T", finalT,".csv")));
      cumulCC <- results$CumulCC_t[randomFirstT];
      cumulSV <- results$CumulSV_t[randomFirstT];
      cumulMC <- results$CumulMC_t[randomFirstT];
      cumulCS <- results$CumulCS_t[randomFirstT];
      currentNPV <- results$CurrentNPV[randomFirstT];
      travelTimeSavings <- results$CS_t[randomFirstT+1] * 1e6 / AADT_FACTOR; #Reverse engineering so that the CSPart will be correct
      results <- results[1:randomFirstT,];
    } else {
      startT <- 1;
    }
    
    for(currentT in startT:finalT){
      print(paste(solutionType,currentT));
      
      fname <- paste0(theOutputDir, "/SelectedLinks",solutionType,"_T", currentT,".csv");
     
      chosenLinks <- data.table::fread(file = fname)$LinkID
      chosenLinksBool <- E(g)$LinkID %in% chosenLinks
      
      currentSegments <- logical(N_segments);
      currentSegments[unique(E(g)$SegmentId[chosenLinksBool])] <- TRUE;
      
      newlyAddedSegments <- prevCurrentSegments != currentSegments;
      
      MCPart <- sum(mbFun(currentT)[prevCurrentSegments]);  
      CCPart <- sum(CCFun(currentT)[newlyAddedSegments]);  
      SVPart <- sum(SVFun()[newlyAddedSegments]);
      
      #From travel time savings to actual, scaled consumer surplus
      rawConsumerSurplus <- travelTimeSavings / 1e6 * AADT_FACTOR 
      #Discounting the consumer surplus
      if(currentT == startT){
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
