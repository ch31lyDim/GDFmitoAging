library(ggvenn)
library(ggpubr)
library(ggsci)
library(ggrepel)
GetDiffProt = function(DEP,Log2FCThres,pvalThres,Group1,Group2)
{
    FC_col = paste(Group1,"vs",Group2,"ratio",sep = '_')
    Pval_col = paste(Group1,"vs",Group2,"p.val",sep = '_')
    DEP = DEP[,c("name","ID",FC_col,Pval_col)]
    colnames(DEP) = c("name","ID","log2FC","pval")
    # DEP$pval = p.adjust(DEP$pval,"fdr")
    UP = DEP$name[DEP$log2FC>Log2FCThres&DEP$pval<pvalThres]
    DN = DEP$name[DEP$log2FC< -Log2FCThres&DEP$pval<pvalThres]
    return(list(UP = UP, DN = DN))
}

DEP_norm = read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_normal",header = T)
DEP_norm[DEP_norm$name=="GDF11",]
OverExpr_vs_Control_norm = GetDiffProt(DEP_norm,log2(2),0.05,"OverExpr","Control")
sh_vs_Control_norm = GetDiffProt(DEP_norm,log2(2),0.05,"sh","Control")
OverExpr_vs_sh_norm = GetDiffProt(DEP_norm,log2(2),0.05,"OverExpr","sh")

norm_UP = list(OEN_vs_CN_UP = OverExpr_vs_Control_norm$UP,
         shN_vs_CN_UP = sh_vs_Control_norm$UP,
         OEN_vs_shN_UP = OverExpr_vs_sh_norm$UP)
ggvenn(norm_UP,set_name_size = 8,text_size = 6,stroke_linetype = 0,fill_alpha = 0.2)+
scale_fill_nejm()
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_norm_UP_venn.pdf")
norm_DN = list(OEN_vs_CN_DN = OverExpr_vs_Control_norm$DN,
         shN_vs_CN_DN = sh_vs_Control_norm$DN,
         OEN_vs_shN_DN = OverExpr_vs_sh_norm$DN)
ggvenn(norm_DN,set_name_size = 8,text_size = 6,stroke_linetype = 0,fill_alpha = 0.2)+
scale_fill_nejm()
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_norm_DN_venn.pdf")

norm_OvsC_VS_SvsC = intersect(unlist(OverExpr_vs_Control_norm),unlist(sh_vs_Control_norm))
norm_OvsC_VS_SvsC_df = DEP_norm[DEP_norm$name%in%norm_OvsC_VS_SvsC,c("name","OverExpr_vs_Control_ratio","sh_vs_Control_ratio")]
colnames(norm_OvsC_VS_SvsC_df) = c("name","log2FC_O_vs_C","log2FC_S_vs_C")

norm_OvsC_VS_SvsC_df$sig = ifelse(norm_OvsC_VS_SvsC_df$log2FC_O_vs_C*norm_OvsC_VS_SvsC_df$log2FC_S_vs_C>0,"Same","Diff")
ggscatter(norm_OvsC_VS_SvsC_df,x="log2FC_O_vs_C",y='log2FC_S_vs_C',color="sig",add = 'reg.line')+
        geom_vline(xintercept = 0,col = "gray", linetype = 'dashed')+
        geom_hline(yintercept = 0,col = "gray", linetype = 'dashed')+
        scale_color_manual(values = c("Same"="#021f68","Diff"="#7d0d0d"))+
        geom_text_repel(mapping = aes(label=name,color=sig))+
        stat_cor(method = 'pearson', aes(x = log2FC_O_vs_C, y = log2FC_O_vs_C,color=sig),size = 4)
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/norm_OvsC_VS_SvsC_scatter.pdf",width = 6,height = 6)


DEP_H2O2 = read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_H2O2",header = T)

OverExpr_vs_Control_H2O2 = GetDiffProt(DEP_H2O2,log2(2),0.05,"OverExpr","Control")
sh_vs_Control_H2O2 = GetDiffProt(DEP_H2O2,log2(2),0.05,"sh","Control")
OverExpr_vs_sh_H2O2 = GetDiffProt(DEP_H2O2,log2(2),0.05,"OverExpr","sh")

H2O2_UP = list(OEH_vs_CH_UP = OverExpr_vs_Control_H2O2$UP,
         shH_vs_CH_UP = sh_vs_Control_H2O2$UP,
         OEH_vs_shH_UP = OverExpr_vs_sh_H2O2$UP)
ggvenn(H2O2_UP,set_name_size = 8,text_size = 6,stroke_linetype = 0,fill_alpha = 0.2)+
scale_fill_nejm()
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_H2O2_UP_venn.pdf")

H2O2_DN = list(OEH_vs_CH_DN = OverExpr_vs_Control_H2O2$DN,
         shH_vs_CH_DN = sh_vs_Control_H2O2$DN,
         OEH_vs_shH_DN = OverExpr_vs_sh_H2O2$DN)
ggvenn(H2O2_DN,set_name_size = 8,text_size = 6,stroke_linetype = 0,fill_alpha = 0.2)+
scale_fill_nejm()
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_H2O2_DN_venn.pdf")



H2O2_OvsC_VS_SvsC = intersect(unlist(OverExpr_vs_Control_H2O2),unlist(sh_vs_Control_H2O2))
H2O2_OvsC_VS_SvsC_df = DEP_H2O2[DEP_H2O2$name%in%H2O2_OvsC_VS_SvsC,c("name","OverExpr_vs_Control_ratio","sh_vs_Control_ratio")]
colnames(H2O2_OvsC_VS_SvsC_df) = c("name","log2FC_O_vs_C","log2FC_S_vs_C")

H2O2_OvsC_VS_SvsC_df$sig = ifelse(H2O2_OvsC_VS_SvsC_df$log2FC_O_vs_C*H2O2_OvsC_VS_SvsC_df$log2FC_S_vs_C>0,"Same","Diff")
ggscatter(H2O2_OvsC_VS_SvsC_df,x="log2FC_O_vs_C",y='log2FC_S_vs_C',color="sig",add = 'reg.line')+
        geom_vline(xintercept = 0,col = "gray", linetype = 'dashed')+
        geom_hline(yintercept = 0,col = "gray", linetype = 'dashed')+
        scale_color_manual(values = c("Same"="#021f68","Diff"="#7d0d0d"))+
        geom_text_repel(mapping = aes(label=name,color=sig))+
        stat_cor(method = 'pearson', aes(x = log2FC_O_vs_C, y = log2FC_O_vs_C,color=sig),size = 4)
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/H2O2_OvsC_VS_SvsC_scatter.pdf",width = 6,height = 6)



x = list(H2O2_DEP = H2O2_OvsC_VS_SvsC,
         norm_DEP = norm_OvsC_VS_SvsC)
ggvenn(x,set_name_size = 8,text_size = 6,stroke_linetype = 0,fill_alpha = 0.2)+
scale_fill_nejm()
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/OvsC_VS_SvsC_norm_vs_H2O2.pdf")

OvsC_VS_SvsC_only_in_H2O2 = setdiff(H2O2_OvsC_VS_SvsC,norm_OvsC_VS_SvsC)
OvsC_VS_SvsC_only_in_norm = setdiff(norm_OvsC_VS_SvsC,H2O2_OvsC_VS_SvsC)
OvsC_VS_SvsC_both_in_H2O2_and_norm = intersect(norm_OvsC_VS_SvsC,H2O2_OvsC_VS_SvsC)

OvsC_VS_SvsC_both_in_H2O2_and_norm_df = rbind(
        DEP_H2O2[DEP_H2O2$name%in%OvsC_VS_SvsC_both_in_H2O2_and_norm,c("name","OverExpr_vs_Control_ratio","sh_vs_Control_ratio")],
        DEP_norm[DEP_norm$name%in%OvsC_VS_SvsC_both_in_H2O2_and_norm,c("name","OverExpr_vs_Control_ratio","sh_vs_Control_ratio")])
OvsC_VS_SvsC_both_in_H2O2_and_norm_df$sig = c(rep("H2O2",length(OvsC_VS_SvsC_both_in_H2O2_and_norm)),
                                              rep("Norm",length(OvsC_VS_SvsC_both_in_H2O2_and_norm)))
colnames(OvsC_VS_SvsC_both_in_H2O2_and_norm_df) = c("name","log2FC_O_vs_C","log2FC_S_vs_C","sig")
ggscatter(OvsC_VS_SvsC_both_in_H2O2_and_norm_df,x="log2FC_O_vs_C",y='log2FC_S_vs_C',color="sig",add = 'reg.line')+
        geom_vline(xintercept = 0,col = "gray", linetype = 'dashed')+
        geom_hline(yintercept = 0,col = "gray", linetype = 'dashed')+
        scale_color_manual(values = c("Norm"="#021f68","H2O2"="#7d0d0d"))+
        geom_text_repel(mapping = aes(label=name,color=sig))+
        stat_cor(method = 'pearson', aes(x = log2FC_O_vs_C, y = log2FC_O_vs_C,color=sig),size = 4)
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/OvsC_VS_SvsC_shared.pdf",width = 6,height = 6)

OvsC_VS_SvsC_both_in_H2O2_and_norm_dist = data.frame(name = unique(OvsC_VS_SvsC_both_in_H2O2_and_norm_df$name))
dis = vector()
for(i in 1:nrow(OvsC_VS_SvsC_both_in_H2O2_and_norm_dist))
{
    name = OvsC_VS_SvsC_both_in_H2O2_and_norm_dist$name[i]
    sub_df = OvsC_VS_SvsC_both_in_H2O2_and_norm_df[OvsC_VS_SvsC_both_in_H2O2_and_norm_df$name==name,]
    x1 = sub_df[1,2]
    x2 = sub_df[2,2]
    y1 = sub_df[1,3]
    y2 = sub_df[2,3]
    dist = sqrt((x2-x1)^2+(y2-y1)^2)
    dis[i] = dist
}
OvsC_VS_SvsC_both_in_H2O2_and_norm_dist$dist = dis
OvsC_VS_SvsC_both_in_H2O2_and_norm_dist = arrange(OvsC_VS_SvsC_both_in_H2O2_and_norm_dist,dis)
ggbarplot(OvsC_VS_SvsC_both_in_H2O2_and_norm_dist,y="dist",x="name",fill='name')+
    rotate()+
    guides(fill=FALSE)
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Dist_OvsC_VS_SvsC.pdf",width = 7,height = 7)





