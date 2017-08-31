#!/usr/bin/env Rscript
### 2017-08-28 by Jiyoung Lee
### USAGE: Rscript ./scripts/CalPCC_CoexpNetwork.R
### To get Co-expression Networks from gene expression profiles. 
### Input: [GMX|ATH]_FPKM_Subset.csv
### Output: CoexpNetwork_[GMX|ATH]_FPKM.[csv|Rdata]

### Setting a project path.
#HOME_PATH=Sys.getenv("HOME")
#PROJ_PATH=paste0(HOME_PATH,"/ATH_GMA")
#setwd(PROJ_PATH)
PROJ_PATH=getwd()

# https://github.com/jasdumas/data-exploreR/blob/master/www/flattenCorrMatrix.R
flattenCorrMatrix <- function(cormat, pmat) {
	ut <- upper.tri(cormat)
	data.frame(row=rownames(cormat)[row(cormat)[ut]], 
	column=rownames(cormat)[col(cormat)[ut]], cor=(cormat)[ut], p=pmat[ut])
}

### Install matrixStats package
install.packages("matrixStats")


### For loop to calculate Co-expression Network for each gene expression profile.
setwd(paste(PROJ_PATH,"/result/"))

Inputs = c("GMX_FPKM.csv", "ATH_FPKM.csv")
for (Input in Inputs){
		
	### Read input file 
	GE=read.table(Input, sep = ",", header = T,stringsAsFactors= F, row.names=1)


	### Filtering some genes out
	### SD != 0
	library(matrixStats)
	GE.tmp=as.matrix(GE)
	GE.rowSds = transform(GE.tmp, SD=rowSds(GE.tmp, na.rm=TRUE))
	GE.rowSds = GE.rowSds[which(GE.rowSds [,8]!=0),]
	GE.rowSds = GE.rowSds[order(GE.rowSds [,8]),]
	GE.tmp = GE.rowSds[, (1:7)] 
		
	### any row >= 0.5 
	row_sub = apply(GE.tmp, 1, function(row) any(row >=0.5 ))
	GE.any.5= GE.tmp[row_sub,]
	GE.tmp=GE.any.5
	
	### rowVars >0.5
	GE.tmp=as.matrix(GE.tmp)
	GE.rowVars = transform(GE.tmp, Vars=rowVars(GE.tmp, na.rm=TRUE))	
	GE.rowVars.5 = GE.rowVars[which(GE.rowVars[,8] > 0.5),]
	GE.rowVars.5 = GE.rowVars.5[order(row.names(GE.rowVars.5)), ]
	GE.tmp=GE.rowVars.5[, -8]
	
	
	### Calculating correlation between genes with corr & P-val
	install.packages("Hmisc")
	library("Hmisc")	
	GE_gene_rcor <- rcorr(as.matrix(t(GE.tmp)))
	
	### Converting cor matrix to cor table by flattenCorrMatrix function
	GE_gene_rcor_list =flattenCorrMatrix(GE_gene_rcor$r, GE_gene_rcor$P)


	### Filtering rcor result
	pvalCut=0.001
	corCut=0.99

	### p-val cutoff
	GE_gene_rcor_list_tmp = GE_gene_rcor_list
	GE_gene_rcor_list_filterP= GE_gene_rcor_list_tmp[which(GE_gene_rcor_list_tmp$p<=pvalCut),]
	GE_gene_rcor_list_tmp = GE_gene_rcor_list_filterP
	
	### cor cutoff
	GE_gene_rcor_list_filterCor = GE_gene_rcor_list_tmp[which(GE_gene_rcor_list_tmp$cor>=corCut),]
	GE_gene_rcor_list_tmp= GE_gene_rcor_list_filterCor
	GE_gene_rcor_list_Out = cbind(GE_gene_rcor_list_tmp[,1:2], cor=round(GE_gene_rcor_list_tmp[,3],3), pval=signif(GE_gene_rcor_list_tmp[,4],3))
	
	
	### Writing filtered rcor result into csv file
	Output=paste0("CoexpNetwork_", Input)
	write.csv(GE_gene_rcor_list_Out, file=Output)

	
	### Saving co-expression network into RData
	RDataFile=paste0("CoexpNetwork_", gsub(".csv", "", Input),".RData")
	print("RDataFile name:"); print(RDataFile)
	save(GE, GE.tmp, GE_gene_rcor_list_Out, file = RDataFile)
}

