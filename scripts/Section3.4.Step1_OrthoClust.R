#!/usr/bin/env Rscript
### 2017-08-28 by Jiyoung Lee
### USAGE: Rscript ./scripts/Section3.4.Step1_OrthoClust.R
### To get Co-expression Networks from gene expression profiles. 
### Input: CoexpNetwork_ATH_FPKM.csv: this is a list of all co-expression edges in Arabidopsis data.
###        CoexpNetwork_GMX_FPKM.csv: this is a list fo all co-expression edges in soybean data.
###        ARATH2GLYMA.RBH.example.txt: this is a list of RBH genes between two species.
###
### Output: Orthoclust_Results.csv, 
###         Orthoclust_Results.RData, 
###         Orthoclust_Results_Summary.csv, 


### Setting a project path.
PROJ_PATH=getwd()
setwd(paste(PROJ_PATH,"/results/"))


### Loading GMX network data
GMX_GE_gene_rcor_list = read.table("CoexpNetwork_GMX_FPKM.csv",header=TRUE,row.names=1,as.is=T)

### Loading ATH network data
ATH_GE_gene_rcor_list = read.table("CoexpNetwork_ATH_FPKM.csv",header=TRUE,row.names=1,as.is=T)

### Co-expression input files for OrthoClust R
GMX_edgelist= as.matrix(GMX_GE_gene_rcor_list[, 1:2])
ATH_edgelist= as.matrix(ATH_GE_gene_rcor_list[, 1:2])

### Loading RHB gene pairs for OrthoClust R
GmxAthOrthoFile="ARATH2GLYMA.RBH_extracted.txt"
GAortho=read.table(GmxAthOrthoFile, sep = "\t", header = F,stringsAsFactors= F)
GA_orthologs=GAortho[, 1:2]


### Run OrthoClust:  cluster two simulated networks using OrthoClust
library(OrthoClust)
KAPPA=3
orthoclust_results=OrthoClust2(Eg1=GMX_edgelist, Eg2=ATH_edgelist, list_orthologs=GA_orthologs, kappa=KAPPA)


### OrthoClust result into csv file
GMX_modules=as.matrix(orthoclust_results[[1]])
ATH_modules=as.matrix(orthoclust_results[[2]])
GMX_modules=as.data.frame(cbind(GMX_modules, "GMX", rownames(GMX_modules)),stringsAsFactors=FALSE)
ATH_modules=as.data.frame(cbind(ATH_modules, "ATH", rownames(ATH_modules)),stringsAsFactors=FALSE)
Merged_modules=rbind(GMX_modules, ATH_modules)
colnames(Merged_modules) = c("ModuleID", "Species", "GeneID")
Merged_modules[, 1] <- as.numeric(Merged_modules[, 1])
Merged_modules = Merged_modules[order(Merged_modules$ModuleID),]
rownames(Merged_modules) = NULL


### Organizing module info into a module number table 
NumModule =as.numeric(rownames(as.matrix(tail(table(orthoclust_results$species1_clust), n=1))));
ortho_module_list = matrix(0, nrow = NumModule, ncol = 4)
ortho_module_list = cbind(as.matrix(1:NumModule), ortho_module_list) 
ortho_tb1=table(orthoclust_results$species1_clust)
ortho_tb2=table(orthoclust_results$species2_clust)
ortho_tb1_list = cbind(Module=as.numeric(rownames(as.matrix(ortho_tb1))), NoGenes=as.numeric(as.matrix(ortho_tb1)))
ortho_tb2_list = cbind(Module=as.numeric(rownames(as.matrix(ortho_tb2))), NoGenes=as.numeric(as.matrix(ortho_tb2)))

### Filling in the module number table 
for (i in 1:NumModule)
{
	idmx1 = match(ortho_module_list[i,1], ortho_tb1_list[,1])
	modNum1= ortho_tb1_list[idmx1,1]
	numGen1= ortho_tb1_list[idmx1,2]
	ortho_module_list[i,2]=modNum1
	ortho_module_list[i,3]=numGen1
	
	idmx2 = match(ortho_module_list[i,1], ortho_tb2_list[,1])
	modNum2= ortho_tb2_list[idmx2,1]
	numGen2= ortho_tb2_list[idmx2,2]
	ortho_module_list[i,4]=modNum2
	ortho_module_list[i,5]=numGen2
}
ortho_module_list=ortho_module_list[, c(1,3,5)]
ortho_module_list[is.na(ortho_module_list)] <- 0

### Calculating percentages for the module number table
#install.packages("matrixStats")
library(matrixStats)
tmp = ortho_module_list
ortho_module_list=cbind(ortho_module_list,sum=ortho_module_list[,2]+ortho_module_list[,3])
ortho_module_list=cbind(ortho_module_list,ratio1=round(ortho_module_list[,2]/ortho_module_list[,4],3), ratio2=round(ortho_module_list[,3]/ortho_module_list[,4],3))
colnames(ortho_module_list)=c("ModuleID","NumOfGenesFromSpe1","NumOfGenesFromSpe2","TotalNumOfGenesFromModule","RatioOfGenesFromSpe1","RatioOfGenesFromSpe2")
ortho_module_list_sumOrdered <- ortho_module_list[order(ortho_module_list[,4],decreasing = TRUE),] 


### Saving OrthoClust result into csv
OrthoOutFile=paste0("Orthoclust_Results.csv")
write.csv(Merged_modules , OrthoOutFile)

### Saving a table of OrthoClust result into csv file
TableFile=gsub(".csv", "_Summary.csv", OrthoOutFile)
write.csv(ortho_module_list_sumOrdered , TableFile)

### Saving OrthoClust results into RData
OrthoRDataFile=gsub(".csv", ".RData", OrthoOutFile)
save(ortho_module_list_sumOrdered, orthoclust_results, GMX_edgelist, ATH_edgelist, GA_orthologs, KAPPA, file = OrthoRDataFile)



