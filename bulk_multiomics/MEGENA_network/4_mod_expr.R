library(igraph)
library(ggplot2)
library(MetBrewer)
library(ComplexHeatmap)
datExpr = read.table("/home/data3/dingyanheng/COP/data/TPM_name_mat",sep = '\t',header = T)
rownames(datExpr) = datExpr$Name
datExpr = datExpr[,-1]
datExpr_sd = apply(datExpr,MARGIN = 1,sd)
datExpr = datExpr[which(datExpr_sd!=0),]
datExpr = log2(datExpr+1)
mod = readRDS("/home/data3/dingyanheng/COP/MEGENA/module.rds")
GDF11_mode = lapply(mod[[1]],function(x){return("GDF11"%in%x)})
GDF11_mode = names(GDF11_mode)[which(GDF11_mode==TRUE)]

mod_gene = sapply(GDF11_mode,function(x){return(mod$modules[[x]])})
ExprMod = lapply(mod_gene,function(x){return(datExpr[x,])})

Anno_col = met.brewer('Derain',6)
Heat_col = met.brewer("Benedictus",10)
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
    p = Heatmap(Expr,name = "log2(TPM+1)",
                top_annotation = top_an,col = colFunc,show_row_names = show,
                bottom_annotation = bot_an,cluster_columns = F)
    return(p)
}
pdf("/home/data3/dingyanheng/COP/MEGENA/modeexpr/mod1_expr.pdf",height = 5,width = 5)
p_mod1 = drawModHeatMap(ExprMod[[1]])
print(p_mod1)
dev.off()

pdf("/home/data3/dingyanheng/COP/MEGENA/modeexpr/mod2_expr.pdf",height = 5,width = 5)
p_mod2 = drawModHeatMap(ExprMod[[2]])
print(p_mod2)
dev.off()

pdf("/home/data3/dingyanheng/COP/MEGENA/modeexpr/mod3_expr.pdf",height = 6,width = 6)
p_mod3 = drawModHeatMap(ExprMod[[3]],show=T)
print(p_mod3)
dev.off()


GetModMeanExpr = function(Expr){
    
    return(colMeans(Expr))
}

child_mod = setdiff(names(mod$modules),mod$module.table$module.parent)
ExprMod = lapply(mod[[1]][child_mod],function(x){return(datExpr[x,])})

ExprMean = lapply(ExprMod,GetModMeanExpr)
ExprMeanDf = do.call(rbind,ExprMean)

Group = as.data.frame(substr(colnames(datExpr),1,nchar(colnames(datExpr))-1))
    rownames(Group) = colnames(datExpr)
    colnames(Group) = "Group"
top_an = HeatmapAnnotation(df = Group,name="Group",
    col = list("Group" = c("CH"=Anno_col[1],"CN"=Anno_col[2],"OEN"=Anno_col[3],
                        "OEH"=Anno_col[4],"shH"=Anno_col[5],"shN"=Anno_col[6])))
colFunc = circlize::colorRamp2(c(min(ExprMeanDf),mean(unlist(ExprMeanDf)),max(ExprMeanDf)),
                                   c(Heat_col[10],"white",Heat_col[1]))
bot_an = HeatmapAnnotation(ExprMean = anno_lines(t(ExprMeanDf),gp = gpar(col=1:100),height = unit(4,"cm")))
p = Heatmap(ExprMeanDf,name = "ModuleMeanExpr",col = colFunc,
                top_annotation = top_an,show_row_names = F,
                cluster_columns = F,bottom_annotation = bot_an)

pdf("/home/data3/dingyanheng/COP/MEGENA/modeexpr/mode_expr.pdf",width = 5,height = 5)
print(p)
dev.off()


same_level = mod$modules[which(mod$module.table$module.parent=="c1_165")]
same_level$c1_165 = mod$modules$c1_165

ExprMod = lapply(same_level,function(x){return(datExpr[x,])})
ExprModLong = lapply(ExprMod,function(x){return(as.data.frame(as.table(as.matrix(x))))})

ExprModLong = do.call(rbind,ExprModLong)
ExprModLong$module = substr(rownames(ExprModLong),1,6)
colnames(ExprModLong) = c("Gene","Group","MeanExpr","module")
ExprModLong$Design = substr(ExprModLong$Group,1,nchar(as.character(ExprModLong$Group))-1)

ExprMean = lapply(ExprMod,GetModMeanExpr)
ExprMeanDf = do.call(rbind,ExprMean)

pdf("/home/data3/dingyanheng/COP/MEGENA/modeexpr/mod_corr_same_parent.pdf",height = 5,width = 5)
corrplot::corrplot(cor(t(ExprMeanDf)),method = "color",
                    order = "hclust",addrect = 2,addCoef.col = "white",number.font = 4,
                    tl.col = "black")
dev.off()

mean_se_ <- ggpubr::mean_se_

ExprModLong$level = "child"
ExprModLong$level[ExprModLong$module=="c1_165"] = "parent"
ggpubr::ggline(ExprModLong,x="Design",y="MeanExpr",plot_type = "l",
                color="module",add=c("mean"),point.size=0,
                position = position_dodge(0.5),alpha = 0.5)+
                ggsci::scale_color_nejm()+
                ggpubr::theme_pubr()
ggsave("/home/data3/dingyanheng/COP/MEGENA/modeexpr/mod_expr_same_parent.pdf",height = 4,width = 4.5)

