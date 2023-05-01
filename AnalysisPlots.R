baseFolder <- "/path/to/somewhere"; #Adjust primary folder here <------
inputDir <- paste0(baseFolder,"/Input")
outputDir <- paste0(baseFolder,"/Output")
solutionTypes <- c("LP","LPGreedyStop","LPGreedy","ByActualOrder","ByRandomSegment","ByShortestSegment","ByShortestRoute","ByLongestSegment","ByLongestRoute")

###### # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
## Plot NPV and BCR development #####
###### # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

fs <- 1;
palette(c("black","red2","forestgreen","blue","turquoise2","orchid1", "seagreen1","goldenrod1","gray","chocolate1"));
legendCatsToUse <- c("BLP", "Greedy w/ stop", "Greedy w/o stop","Actual order","Random order","Shorter segments first", "Shorter routes first",
                     "Longer segments first", "Longer routes first");
thisPCHs <- c(20,20,20,20,20,20,20,20,20)
thisCols <- c(1,2,8,8,6,5,4,7,3)
thisCEXs <- c(2,sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2))

for(budgetType in c("Long","Short","RandomFirst")){
  theMainOutputDir <- paste0(outputDir,"/", budgetType);
  #png(paste0(theMainOutputDir,"/EvaluationComparisonPlot_",budgetType,".png"), width = 700, height = 1000, pointsize = 20);
  
  allEvalSums <- NULL;
  for(solutionType in solutionTypes){
    theOutputDir <- paste0(theMainOutputDir,"/",solutionType);
    if((budgetType != "Short" & solutionType == "ByActualOrder") |
       (budgetType == "Short" & solutionType == "LPGreedyStop")){
      skipI <- which(solutionType == solutionTypes)
      next;
    }
    
    if(file.exists(paste0(theOutputDir, "/Evaluation", solutionType,"_T", finalT,".csv"))){
      thisRes <- data.frame(data.table::fread(paste0(theOutputDir, "/Evaluation", solutionType,
                                                     "_T", finalT,".csv")))
      thisRes$Config <- solutionType
      if(is.null(allEvalSums)){
        allEvalSums <- thisRes;
      } else {
        allEvalSums <- rbind(allEvalSums, thisRes)
      }
    }
  }
  
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