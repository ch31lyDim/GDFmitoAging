library(GENIE3)
library(data.table)
exprMatr = read.table("/home/data3/dingyanheng/COP/data/TPM_name_mat",sep = '\t',
                    header = T)
rownames(exprMatr) = exprMatr$Name
exprMatr = as.matrix(exprMatr[,-1])
TF = as.data.frame(fread("/home/data3/dingyanheng/COP/TF_list/Canis_lupus_familiaris_TF.txt",na.strings = ''))$Symbol
TF = unique(na.omit(TF))
TF = TF[TF%in%rownames(exprMatr)]
weightMat <- GENIE3(exprMatr,treeMethod = "RF",nCores = 63,verbose = T,nTrees = 2500)
linkList <- getLinkList(weightMat)
