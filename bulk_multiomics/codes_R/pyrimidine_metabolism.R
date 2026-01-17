library(org.Cf.eg.db)
keggrest <- KEGGREST::keggGet("cfa00240")
genes_pm<-unlist(lapply(keggrest[[1]]$GENE,function(x) strsplit(x,';')))
genes_pm <- genes_pm[genes_pm%in%as.data.frame(org.Cf.egSYMBOL2EG)[,"symbol"]]

library(data.table)
TPM <- fread("/data1/dyh/GDF11/TPM_name_mat",sep = '\t',header = T,data.table = FALSE)
TPM <- TPM[TPM$Name %in% genes_pm,]
rownames(TPM) <- TPM$Name
TPM <- TPM[,-1]

library(ComplexHeatmap)
library(circlize)
# 设置颜色映射函数
tpm_matrix <- as.matrix(TPM)
tpm_scaled <- t(scale(t(tpm_matrix))) 
col_fun <- colorRamp2(c(min(tpm_scaled), 0, max(tpm_scaled)), 
                      c("navy", "white", "firebrick3"))
# 为NT5E行创建特殊的行标签
row_labels <- rownames(tpm_scaled)
row_label_colors <- ifelse(row_labels == "NT5E", "red", "black")
# 创建热图
ht <- Heatmap(tpm_scaled,
              name = "Z-score",
              col = col_fun,
              cluster_rows = TRUE,
              cluster_columns = F,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_names_gp = gpar(col = row_label_colors, fontface = ifelse(row_labels == "NT5E", "bold", "plain")),
              column_names_gp = gpar(fontsize = 10),
              row_names_side = "left",
              row_dend_side = "right",
              heatmap_legend_param = list(title = "Z-score"))


drawModHeatMap = function(Expr,show=F){
    Group = as.data.frame(substr(colnames(Expr),1,nchar(colnames(Expr))-1))
    rownames(Group) = colnames(Expr)
    colnames(Group) = "Group"
    top_an = HeatmapAnnotation(df = Group,name="Group",
        col = list("Group" = c("CH"=Anno_col[1],"CN"=Anno_col[2],"OEN"=Anno_col[3],
                               "OEH"=Anno_col[4],"shH"=Anno_col[5],"shN"=Anno_col[6])))
    bot_an = HeatmapAnnotation(Expr = anno_boxplot(Expr,gp = gpar(fill = sort(rep(1:6,3)))),
                               meanExpr = anno_lines(colMeans(Expr),gp = gpar(fill = sort(rep(1:6,3)))))
    colFunc = circlize::colorRamp2(c(min(Expr),mean(unlist(Expr)),max(Expr)),
                                   c(Heat_col[10],"white",Heat_col[1]))
    p = Heatmap(Expr,name = "Z-score",
                top_annotation = top_an,col = colFunc,show_row_names = T,cluster_rows = T,
                bottom_annotation = bot_an,cluster_columns = F,show_column_names = T,
                row_names_gp = gpar(col = row_label_colors, fontface = ifelse(row_labels == "NT5E", "bold", "plain")))
    return(p)
}
Anno_col = met.brewer('Derain',6)
Heat_col = met.brewer("Benedictus",10)
drawModHeatMap((as.matrix(tpm_scaled)))

pdf("/data1/dyh/GDF11/plots/pyrimidine_NT5E_highlighted.pdf", width = 5, height = 10)
draw(drawModHeatMap((as.matrix(tpm_scaled))), heatmap_legend_side = "right")
dev.off()

CD73_expr <- as.data.frame(unlist(TPM["NT5E",]))
CD73_expr$Sample <- rownames(CD73_expr)
colnames(CD73_expr) <- c("Expr","Sample")
CD73_expr$Group <- substr(CD73_expr$Sample,1,nchar(CD73_expr$Sample)-1)
CD73_expr$Group <- factor(CD73_expr$Group,levels = c("shN","shH","CN","CH","OEN","OEH"))
CD73_expr$Group1 <- ifelse(CD73_expr$Group %in% c("shN","shH"),"sh",
                            ifelse(CD73_expr$Group %in% c("CN","CH"),"CN","OE"))
CD73_expr$Group1 <- factor(CD73_expr$Group1,levels = c("sh","CN","OE"))
CD73_expr$Group2 <- ifelse(CD73_expr$Group %in% c("shN","CN","OEN"),"N","H")
CD73_expr$Group2 <- factor(CD73_expr$Group2,levels = c("N","H"))

library(ggplot2)
library(ggpubr)
ggboxplot(data = CD73_expr,x = "Group",y = "Expr",color = "Group1",add =c("mean_se","jitter"))+
       xlab("")+
       theme(legend.title = element_blank())+
       theme(aspect.ratio = 1/1)+
       stat_compare_means(comparison = list(c("CH","CN"),c("shN","shH"),c("OEN","OEH")),label="p.signif",size=6,method = "wilcox.test")
ggsave("/data1/dyh/GDF11/plots/NT5E_expr.pdf",height = 4,width = 6)


mean(CD73_expr$Expr[1:3])/mean(CD73_expr$Expr[4:6])
mean(CD73_expr$Expr[7:9])/mean(CD73_expr$Expr[10:12])
mean(CD73_expr$Expr[13:15])/mean(CD73_expr$Expr[16:18])
