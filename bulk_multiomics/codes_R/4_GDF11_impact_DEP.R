library(data.table)
library(tidyverse)
library(DEP)
library(uwot)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(ggplot2)
library(ggrepel)
library(circlize)
library(DEP)
prot = read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/data/all_sample.xls",sep = '\t',header = T)
prot = prot[,-2]
data_unique <- make_unique(prot,"Gene","Protein",delim = ";")
# data_unique = data_unique[data_unique$Gene!='',]
LFQ_cols = 3:20
Design_mat = data.frame(Sample = colnames(data_unique)[LFQ_cols],
                        Genotype = c(rep("Control",6),rep("sh",6),rep("OverExpr",6)),
                        Condition = rep(c(rep("Normal",3),rep("H2O2",3)),3)
                        )

LFQ_norm = LFQ_cols[Design_mat$Condition=='Normal']
prot_design_norm = data.frame(label = colnames(data_unique)[LFQ_norm],
                         condition = c(rep("Control",3),rep("sh",3),rep("OverExpr",3)),
                         replicate = rep(1:3,3))

data_se = make_se(data_unique,LFQ_norm,prot_design_norm)
data_filter <- filter_missval(data_se,thr = 3)
data_norm = normalize_vsn(data_se)
data_imp <- impute(data_filter,fun = "min",q=0.01)
data_diff <- test_diff(data_imp, type = "manual",test = c("sh_vs_Control","OverExpr_vs_Control","OverExpr_vs_sh"))
dep <- add_rejections(data_diff,alpha = 0.05,lfc = log2(2))
data_results <- get_results(dep)

assay(data_imp)["GDF11",]
write.table(data_results,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_normal",
            col.names = T,row.names = F,sep = '\t',quote = F)


LFQ_H2O2 = LFQ_cols[Design_mat$Condition=='H2O2']
prot_design_H2O2 = data.frame(label = colnames(data_unique)[LFQ_H2O2],
                         condition = c(rep("Control",3),rep("sh",3),rep("OverExpr",3)),
                         replicate = rep(1:3,3))

data_se = make_se(data_unique,LFQ_H2O2,prot_design_H2O2)
data_filter <- filter_missval(data_se,thr = 3)
data_norm = normalize_vsn(data_se)
data_imp <- impute(data_filter,fun = "min",q = 0.01)
data_diff <- test_diff(data_imp, type = "manual",test = c("sh_vs_Control","OverExpr_vs_Control","OverExpr_vs_sh"))
dep <- add_rejections(data_diff,alpha = 0.05,lfc = log2(1.2))
data_results <- get_results(dep)
write.table(data_results,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_H2O2",
            col.names = T,row.names = F,sep = '\t',quote = F)







