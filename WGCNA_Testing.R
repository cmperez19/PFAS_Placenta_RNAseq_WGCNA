library(WGCNA)
options(stringsAsFactors = FALSE)

load("raw_norm_combat.RData")
#setwd("~/PFAS")
#load("PFAS_Placenta_RNAseq.RData")

combat.adj = t(combat.adj)


# Choose a set of soft-thresholding powers
#powers = c(c(1:10), seq(from = 12, to=20, by=2))

# powers = 1:20
# #Call the network topology analysis function 
# #Analysis of scale free topology for multiple hard thresholds 
# sft = pickSoftThreshold(data = combat.adj, dataIsExpr = T, powerVector = powers)
# # Plot the results:
# sizeGrWindow(9, 5)
# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# 
# 
# k=softConnectivity(datE=combat.adj,power=18)
# # Plot a histogram of k and a scale free topology plot
# sizeGrWindow(10,5)
# par(mfrow=c(1,2))
# hist(k)
# scaleFreePlot(k, main="Check scale free topology\n")

powers1=seq(from=1,to=20,by=1) 
pst=pickSoftThreshold(combat.adj,powerVector=powers1, moreNetworkConcepts=T)[[2]] 
cex1=1.2
cex.axis=1.5
cex.lab=1.5
cex.main=1.5 
par(mfrow=c(2,2));
plot(powers1,-sign(pst[,3])*pst[,2],type="n",
     xlab="Soft Threshold",
     ylab="SFT,signed Rˆ2",cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab, main="Scale Free Fit Index Rˆ2");
text(powers1,-sign(pst[,3])*pst[,2], labels=powers1,cex=cex1,col="red")
# this line corresponds to using an Rˆ2 cut-off of h 
abline(h=0.95,col="red")
plot(powers1,pst$Density,
     type="n", xlab="Soft Threshold",ylab="Density", cex.axis=cex.axis, cex.main=cex.main,cex.lab=cex.lab, main="Density") 
text(powers1,pst$Density,labels=powers1,cex=cex1,col="red") 
plot(powers1,pst$Heterogeneity,type="n",xlab="Soft Threshold", ylab="Heterogeneity",cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, main="Heterogeneity") 
text(powers1,pst$Heterogeneity,labels=powers1,cex=cex1,col="red") 
plot(powers1,pst$Centralization,type="n",xlab="Soft Threshold", ylab="Centralization",cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab, main="Centralization") 
text(powers1, pst$Centralization, labels=powers1,cex=cex1, col="red")


