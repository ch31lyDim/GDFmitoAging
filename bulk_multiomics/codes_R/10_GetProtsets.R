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
library(ggrepel)
library(ggforce)
prot = as.data.frame(fread("/home/data5/dingyanheng/COP/GDF11/Proteome/data/all_sample.xls"))
prot = prot[,-2]
data_unique <- make_unique(prot,"Gene","Protein",delim = ";")
data_unique = data_unique[data_unique$Gene!='',]
LFQ_cols = 3:20
prot_design = data.frame(label = colnames(data_unique)[LFQ_cols],
                         condition = substr(colnames(data_unique)[LFQ_cols],1,nchar(colnames(data_unique)[LFQ_cols])-3),
                         replicate = rep(1:3,6))
data_se = make_se(data_unique,LFQ_cols,prot_design)
data_filter <- filter_missval(data_se,thr = 3)
data_norm = normalize_vsn(data_se)
data_imp <- impute(data_norm,fun = "MinProb",q=0.01)
data_diff <- test_diff(data_imp, type = "manual",
                       test = c("C_H_vs_C_NC","SH_H_vs_SH_NC","OE_H_vs_OE_NC"))
dep <- add_rejections(data_diff,alpha = 0.05,lfc = log2(2))
data_results <- get_results(dep)

DEP_CH_vs_CN_UP = data_results$name[data_results$C_H_vs_C_NC_p.val<0.05 & data_results$C_H_vs_C_NC_ratio >1]
DEP_CH_vs_CN_DN = data_results$name[data_results$C_H_vs_C_NC_p.val<0.05 & data_results$C_H_vs_C_NC_ratio < -1]

DEP_OEH_vs_OEN_UP = data_results$name[data_results$SH_H_vs_SH_NC_p.val<0.05 & data_results$SH_H_vs_SH_NC_ratio >1]
DEP_OEH_vs_OEN_DN = data_results$name[data_results$SH_H_vs_SH_NC_p.val<0.05 & data_results$SH_H_vs_SH_NC_ratio < -1]

DEP_shH_vs_shN_UP = data_results$name[data_results$OE_H_vs_OE_NC_p.val<0.05 & data_results$OE_H_vs_OE_NC_ratio >1]
DEP_shH_vs_shN_DN = data_results$name[data_results$OE_H_vs_OE_NC_p.val<0.05 & data_results$OE_H_vs_OE_NC_ratio < -1]


write.table(DEP_CH_vs_CN_UP,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/DEP_CH_vs_CN_UP",
            col.names = F,row.names = F,quote = F)
write.table(DEP_CH_vs_CN_DN,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/DEP_CH_vs_CN_DN",
            col.names = F,row.names = F,quote = F)

write.table(DEP_OEH_vs_OEN_UP,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/DEP_OEH_vs_OEN_UP",
            col.names = F,row.names = F,quote = F)
write.table(DEP_OEH_vs_OEN_DN,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/DEP_OEH_vs_OEN_DN",
            col.names = F,row.names = F,quote = F)

write.table(DEP_shH_vs_shN_UP,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/DEP_shH_vs_shN_UP",
            col.names = F,row.names = F,quote = F)
write.table(DEP_shH_vs_shN_DN,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/DEP_shH_vs_shN_DN",
            col.names = F,row.names = F,quote = F)



prot = read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/data/all_sample.xls",sep = '\t',header = T)
prot = prot[,-2]
data_unique <- make_unique(prot,"Gene","Protein",delim = ";")
data_unique = data_unique[data_unique$Gene!='',]
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



DEP_OEN_vs_CN_UP = data_results$name[data_results$OverExpr_vs_Control_p.val < 0.05 & data_results$OverExpr_vs_Control_ratio > 1]
DEP_OEN_vs_CN_DN = data_results$name[data_results$OverExpr_vs_Control_p.val < 0.05 & data_results$OverExpr_vs_Control_ratio < -1]

DEP_shN_vs_CN_UP = data_results$name[data_results$OverExpr_vs_Control_p.val < 0.05 & data_results$OverExpr_vs_Control_ratio > 1]
DEP_shN_vs_CN_DN = data_results$name[data_results$OverExpr_vs_Control_p.val < 0.05 & data_results$OverExpr_vs_Control_ratio < -1]

DEP_OEN_vs_shN_UP = data_results$name[data_results$OverExpr_vs_Control_p.val < 0.05 & data_results$OverExpr_vs_Control_ratio > 1]
DEP_OEN_vs_shN_DN = data_results$name[data_results$OverExpr_vs_Control_p.val < 0.05 & data_results$OverExpr_vs_Control_ratio < -1]

write.table(DEP_OEN_vs_CN_UP,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/DEP_OEN_vs_CN_UP",
            col.names = F,row.names = F,quote = F)
write.table(DEP_OEN_vs_CN_DN,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/DEP_OEN_vs_CN_DN",
            col.names = F,row.names = F,quote = F)

write.table(DEP_shN_vs_CN_UP,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/DEP_shN_vs_CN_UP",
            col.names = F,row.names = F,quote = F)
write.table(DEP_shN_vs_CN_DN,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/DEP_shN_vs_CN_DN",
            col.names = F,row.names = F,quote = F)

write.table(DEP_OEN_vs_shN_UP,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/DEP_OEN_vs_shN_UP",
            col.names = F,row.names = F,quote = F)
write.table(DEP_OEN_vs_shN_DN,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/DEP_OEN_vs_shN_DN",
            col.names = F,row.names = F,quote = F)


