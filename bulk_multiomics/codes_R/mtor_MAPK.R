library(clusterProfiler)
library(ComplexHeatmap)
library(data.table)
library(rtracklayer)
library(dplyr)
library(ggplotify)
library(ggplot2)
library(cowplot)
library(MetBrewer)
library(ggpubr)
library(org.Cf.eg.db)
DEP = fread("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP.txt")
DEP_norm <- fread("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_normal")
DEP_H2O2 <- fread("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_H2O2")

DEP_df <- data.frame(ID = DEP$name,
                     CH_vs_CN = DEP$C_H_vs_C_NC_ratio,
                     OEH_vs_OEN = DEP$OE_H_vs_OE_NC_ratio,
                     ShH_vs_ShN = DEP$SH_H_vs_SH_NC_ratio,
                     OEN_vs_CN = DEP_norm$OverExpr_vs_Control_ratio,
                     ShN_vs_CN = DEP_norm$sh_vs_Control_ratio,
                     OEN_vs_ShN = DEP_norm$OverExpr_vs_sh_ratio,
                     OEH_vs_CH = DEP_H2O2$OverExpr_vs_Control_ratio,
                     ShH_vs_CH = DEP_H2O2$sh_vs_Control_ratio,
                     OEH_vs_ShH = DEP_H2O2$OverExpr_vs_sh_ratio)

rownames(DEP_df) <- DEP_df$ID
DEP_df <- DEP_df[,2:ncol(DEP_df)]



GeneSets <- lapply(1:9,FUN = function(i){
    FC <- DEP_df[,i]
    names(FC) <- rownames(DEP_df)
    FC <- sort(FC,decreasing = T)
    return(FC)
})
names(GeneSets) <- colnames(DEP_df)



gsea_go_res<- compareCluster(GeneSets, 
                    "gseGO",
                    ont = "ALL",
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    verbose = T,
                    OrgDb = org.Cf.eg.db,
                    keyType = "SYMBOL")
dotplot(gsea_go_res)

GeneSets_ENTREZID <- lapply(1:9,function(i){
    tmp <- GeneSets[[i]]
    names(tmp) <- bitr(names(tmp),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Cf.eg.db,drop = T)[,2]
    tmp <- tmp[!is.na(names(tmp))]
    return(tmp)
})

names(GeneSets_ENTREZID) <- colnames(DEP_df)
gsea_kegg_res<- compareCluster(GeneSets_ENTREZID, 
                    "gseKEGG",
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    verbose = T,
                    organism = "cfa",
                    pAdjustMethod = "fdr")
gsea_kegg_res_sig<- compareCluster(GeneSets_ENTREZID, 
                    "gseKEGG",
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    verbose = T,
                    organism = "cfa",
                    pAdjustMethod = "fdr")
dotplot(gsea_kegg_res_sig,split = ".sign",size = "NES")+
    facet_grid(~.sign,scales = "free")
gsea_kegg_res_focused_df <- gsea_kegg_res@compareClusterResult[gsea_kegg_res@compareClusterResult$ID%in%c("cfa04150","cfa04010"),]



library(org.Cf.eg.db)
keggrest <- KEGGREST::keggGet("cfa04150")
genes_mtor<-unlist(lapply(keggrest[[1]]$GENE,function(x) strsplit(x,';')))
genes_mtor <- genes_mtor[genes_mtor%in%as.data.frame(org.Cf.egSYMBOL2EG)[,"symbol"]]

keggrest <- KEGGREST::keggGet("cfa04010")
genes_MAPK <- unlist(lapply(keggrest[[1]]$GENE,function(x) strsplit(x,';')))
genes_MAPK <- genes_MAPK[genes_MAPK%in%as.data.frame(org.Cf.egSYMBOL2EG)[,"symbol"]]

gene_lists <- list(mTOR_Signaling = genes_mtor,MAPK_Signaling = genes_MAPK)

library(GSVA)
prot_expr <- read.table("/data1/dyh/GDF11/prot")

gsva_param <- ssgseaParam(as.matrix(prot_expr),gene_lists)
GSVA_res <- gsva(gsva_param)
GSVA_res <- as.data.frame(GSVA_res)
GSVA_res_long <- reshape2::melt(as.matrix(GSVA_res))
colnames(GSVA_res_long) <- c("Pathway","Sample","ssGSEA_Score")
meta_df <- data.frame(Sample = colnames(prot_expr),
                      Group = c(rep("C",6),rep("sh",6),rep("OE",6)),
                      Treatment = rep(c(rep("Norm",3),rep("H2O2",3)),3))
plt_data <- merge(GSVA_res_long,meta_df,by = "Sample")

mtor_expr <- prot_expr[genes_mtor,]
mtor_expr <- as.data.frame(t(mtor_expr))
mtor_expr$Group <- c(rep("C",6),rep("sh",6),rep("OE",6))
mtor_expr$Treatment <- rep(c(rep("Norm",3),rep("H2O2",3)),3)
mtor_expr_long <- reshape2::melt(as.matrix(mtor_expr))
colnames(mtor_expr_long) <- c("Sample","Prot","Expr")
mtor_expr_long <- mtor_expr_long[mtor_expr_long$Prot%in%genes_mtor,]
mtor_expr_long <- merge(mtor_expr_long,meta_df,by = "Sample")

mtor_expr_long$Expr <- as.numeric(mtor_expr_long$Expr)
mtor_expr_long$Group <- factor(mtor_expr_long$Group,levels = c("sh","C","OE"))

#heatmap
#long df 2 wide df
mtor_expr_wide <- reshape2::dcast(mtor_expr_long[mtor_expr_long$Treatment=='Norm',], Sample ~ Prot, value.var = "Expr")
rownames(mtor_expr_wide) <- mtor_expr_wide$Sample
mtor_expr_wide <- mtor_expr_wide[,-1]
ComplexHeatmap::Heatmap(mtor_expr_wide)

ggviolin(data = mtor_expr_long,x = "Group",y = "Expr",color = "Treatment",size = 1.2,facet.by = 'Treatment')+
       xlab("")+
       theme(legend.title = element_blank())+
       theme(aspect.ratio = 1/1)


library(ggplot2)
library(ggpubr)

plt_data$Group <- factor(plt_data$Group,levels = c("sh","C","OE"))
plt_data <- plt_data[plt_data$Group%in%c("C","OE"),]
plt_data <- plt_data[plt_data$Pathway=='mTOR_Signaling',]
p<-ggline(data = plt_data,x = "Group",y = "ssGSEA_Score",color = "Treatment",add =c("mean_se","jitter"),scales = "free")+
       xlab("")+
       theme(legend.title = element_blank())+
       theme(aspect.ratio = 1/1)
p + stat_compare_means(aes(group = Treatment,color = Treatment), 
                      label = "p.format",
                      method = "t.test",  # 或者使用"wilcox.test"
                      label.y = max(plt_data$ssGSEA_Score) + 0.1) 
ggsave("/data1/dyh/GDF11/plots/mtor.pdf",width = 4,height = 4)
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/mtor_MAPK.png",width = 8,height = 4)



selected_Genes_mTOR <- c("EIF4B","FZD2","GRB2","MAP2K1","RHOA","RPS6KA1","RRAGB")
selected_Genes_MAPK <- c("FGFR1","GRB2","IRAK4","LOC474850","MAP2K1","MAP3K20","MAX","NFATC1","RELB",
                         "RPS6KA1","STMN1","TGFBR2")

expr_plt <- reshape2::melt(as.matrix(prot_expr))
colnames(expr_plt) <- c("Prot","Sample","Expr")
expr_plt <- expr_plt[expr_plt$Prot%in%c(selected_Genes_mTOR,selected_Genes_MAPK),]
expr_plt <- merge(expr_plt,meta_df)

ggbarplot(data = expr_plt,x = "Prot",y = "Expr",color = "Group",add =c("mean_se","jitter"),facet.grid = TRUE,scales = "free",size = 1.2,
         position = position_dodge())+
       xlab("")+
       theme(legend.title = element_blank())

expr_plt <- expr_plt[expr_plt$Group%in%c("C","OE"),]
library(dplyr)
y_max <- max(expr_plt$Expr, na.rm = TRUE)
y_range <- max(expr_plt$Expr, na.rm = TRUE) - min(expr_plt$Expr, na.rm = TRUE)
ggline(data = expr_plt, x = "Group", y = "Expr", color = "Treatment", 
       add = c("mean_se", "jitter"),
       facet.by = c("Prot"), facet.grid = TRUE, scales = "free", size = 1.2) +
       xlab("") +
       theme(aspect.ratio = 1/1) +
       theme(legend.title = element_blank()) +
       # Norm组的p值
       stat_compare_means(
         aes(group = Group),
         method = "t.test",
         label = "p.format",
         color = treatment_colors[1],
         data = ~ .x[.x$Treatment == "Norm",],
         label.y = y_max + 0.05*y_range  # 较低位置
       ) +
       # H2O2组的p值
       stat_compare_means(
         aes(group = Group),
         method = "t.test",
         label = "p.format",
         color = treatment_colors[2],
         data = ~ .x[.x$Treatment == "H2O2",],
         label.y = y_max + 0.15*y_range  # 较高位置
       )
ggsave("/data1/dyh/GDF11/plots/mtor_MAPK_expr.pdf",width = 10,height = 10)


dat <- read.csv("/home/data5/dingyanheng/csv-jianyang-set.csv")
jur <- dat$Journal.Book


