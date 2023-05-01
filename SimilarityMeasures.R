### Initialization #####
library(data.table)

baseFolder <- "/path/to/somewhere"; #Adjust primary folder here <------
finalT <- 50; #Length of the evaluation period

pairedTTest <- function(r,a){
  d <- r - a;
  mu_d <- mean(d);
  s_d <- sd(d);
  SE_d <- s_d / sqrt(length(d));
  z_d <- mu_d / SE_d;
  p_d <- pnorm(-abs(z_d))*2
  return(c(mu_d,s_d,SE_d,z_d,p_d));
}

reFun <- function(r,a){
  return(mean(abs(r-a)/abs(r)));
}

nmaeFun <- function(r,a){
  return(sum(abs(r-a))/length(r) / mean(r));
}

nrmseFun <- function(r,a){
  return(sqrt(mseFun(r,a)) / mean(r));
}

## Calculating measures ####
v1s <- v2s <- v3s <- v4s<- v5s <- c();
for(budgetType in c("Long")){
  for(solutionType in c("LP","LPGreedyStop","LPGreedy","ByActualOrder","ByRandomSegment","ByShortestSegment","ByShortestRoute","ByLongestSegment","ByLongestRoute")){
    
    if(!file.exists(paste0(baseFolder, "/Output/", budgetType, 
                           "/", solutionType,"/Evaluation", solutionType,"_T50.csv"))){
      next;
    }    
    trueCS <- data.frame(data.table::fread(paste0(baseFolder, "/Output/", budgetType, 
                                                  "/", solutionType,"/Evaluation", solutionType,"_T50.csv")))
    Ns <- trueCS$N_t;
    trueCS <- trueCS$CS_t;
    approxCS <- data.frame(data.table::fread(paste0(baseFolder, "/Output/", budgetType, 
                                                    "/", solutionType,"/ExpectedEvaluation", solutionType,"_T50.csv")))$E_CS_t
    
    #Plotting the two consumer surpluses
    fs <- 2;
    pdf(paste0(baseFolder,"/Output/CSComparison_",solutionType,"_",budgetType,".pdf"), width=12.5, height = 7.0, pointsize = 14);
    par(oma = c(0, 0, 0, 0), mar = c(fs*2, fs*2.3, 0.1, 0.1))
    plot(NA,xlim=c(0.8,49.2),ylim=c(-5,84),
         xlab = "Year", ylab = "Consumer Surplus [million DKK]", 
         mgp=c(sqrt(2)*fs,1,0), cex.lab = fs, cex.axis = fs);
    abline(v=seq(-5,55,1), col = rgb(0.975,0.975,0.975), lwd = 0.05)
    abline(v=seq(-5,55,5), col = rgb(0.9,0.9,0.9), lwd = 0.1)
    abline(h=seq(-10,90,2), col = rgb(0.975,0.975,0.975), lwd = 0.05)
    abline(h=seq(-20,90,10), col = rgb(0.9,0.9,0.9), lwd = 0.1)
    points(1:finalT, trueCS, col = 1, pch = 20, cex = 2*fs)
    points(1:finalT, approxCS, col = 2, pch = 20, cex = sqrt(2)*fs)
    legend("topright",c("Proper evaluation", "Linear approximation"), pch=c(20,20),col=c(1,2), pt.cex = c(2*fs,sqrt(2)*fs), cex = fs, bg = rgb(1,1,1));
    dev.off();
    
    considerBool <- c(FALSE, Ns[2:length(Ns)] > Ns[1:(length(Ns)-1)]);
    
    trueCS <- trueCS[considerBool];
    approxCS <- approxCS[considerBool];
    
    z <- abs(pairedTTest(trueCS,approxCS)[4])
    nmae <- nmaeFun(trueCS, approxCS)
    nrmse <- nrmseFun(trueCS, approxCS)
    re <- reFun(trueCS, approxCS)
    zSigned <- pairedTTest(trueCS,approxCS)[4];
    
    print(solutionType);
    print(c(z,nmae,nrmse,re, zSigned))
    v1s <- c(v1s, z);
    v2s <- c(v2s, nmae)
    v3s <- c(v3s, nrmse)
    v4s <- c(v4s, re)
    v5s <- c(v5s, zSigned)
  }
}
print(c(min(v1s),min(v2s),min(v3s),min(v4s),min(v5s)))
print(c(mean(v1s),mean(v2s),mean(v3s),mean(v4s),mean(v5s)))
print(c(max(v1s),max(v2s),max(v3s),max(v4s),max(v5s)))