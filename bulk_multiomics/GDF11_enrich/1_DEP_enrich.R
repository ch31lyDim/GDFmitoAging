library(org.Cf.eg.db)
library(clusterProfiler)
library(ggplot2)
library(MetBrewer)
protExpression <- read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/data/prot")

protExpression$ID <- rownames(protExpression)
protExpression$ID <- sapply(protExpression$ID,FUN = function(x){strsplit(x,"\\.")[[1]][1]})

protExpression <- aggregate(. ~ ID, data = protExpression, FUN = max)
rownames(protExpression) <- protExpression$ID
protExpression <- protExpression[,-1]

DEP <- read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP.txt",header = T)

FC_C <- DEP[,c("name","C_H_vs_C_NC_ratio")]
FC_OE <- DEP[,c("name","OE_H_vs_OE_NC_ratio")]
FC_SH <- DEP[,c("name","SH_H_vs_SH_NC_ratio")]

C_ProtLis <- FC_C$C_H_vs_C_NC_ratio
names(C_ProtLis) <- FC_C$name
C_ProtLis <- sort(C_ProtLis,decreasing = T)

gsego_C <- gseGO(
    C_ProtLis ,
    ont = "ALL",#可以替换为BP CC MF，分别富集
    OrgDb = org.Cf.eg.db,
    keyType = "SYMBOL",
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE,
    eps = 1e-20,
    pAdjustMethod = "BH"
)
dotplot(gsego_C,showCategory = 40,split=".sign")+facet_wrap(~.sign,scales = "free")



#_-------------------------------------------
DEP = na.omit(read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP.txt",header = T))

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
C = GetDiffProt(DEP,log2(2),0.05,"C_H","C_NC")
O = GetDiffProt(DEP,log2(2),0.05,"OE_H","OE_NC")
S = GetDiffProt(DEP,log2(2),0.05,"SH_H","SH_NC")

only_in_o = list(UP = setdiff(O$UP,union(C$UP,S$UP)),
                 DN = setdiff(O$DN,union(C$DN,S$DN)))
only_in_s = list(UP = setdiff(S$UP,union(C$UP,O$UP)),
                 DN = setdiff(S$DN,union(C$DN,O$DN)))
OxS_shared = list(BothUP = intersect(S$UP,O$UP),
                  BothDN = intersect(S$DN,O$DN),
                  O_DN_S_UP = intersect(S$UP,O$DN))
os_overlap_DEP = list(OE_unique_UP=only_in_o$UP,
                      OE_unique_DN=only_in_o$DN,
                      sh_unique_UP=only_in_s$UP,
                      sh_unique_DN=only_in_s$DN,
                      OE_sh_both_UP = OxS_shared$BothUP,
                      OE_sh_both_DN = OxS_shared$BothDN,
                      OE_DN_S_UP = OxS_shared$O_DN_S_UP)

GO_os_overlap_DEP = compareCluster(os_overlap_DEP,"enrichGO",
                                OrgDb = org.Cf.eg.db,
                                keyType="SYMBOL",
                                ont = "ALL",
                                pvalueCutoff = 0.01,
                                qvalueCutoff = 1,
                                pAdjustMethod = "none")

clusterProfiler::dotplot(GO_os_overlap_DEP,
        showCategory = 3,
        includeAll=T,
        label_format = 50,
        split="ONTOLOGY",
        x = "Cluster",
        color = "-log10(p.adjust)")+
        facet_grid(.~ONTOLOGY, scales = "free", space = "free",
                margins = unit(c(0,0,0,0),"cm"))+
        theme_bw()+
        # scale_y_discrete(expand=expansion(mult=c(.1,.1)))+
        # scale_x_discrete(expand=expansion(mult=c(.1,.1)))+
        theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
        scale_fill_met_c("Cassatt1",direction = -1)+
        theme(panel.border = element_rect(colour = "black", size = 0.5))+
        theme(panel.grid = element_blank())+
        theme(panel.grid.major.y = element_line(color = "grey",linetype  = 'dashed'))+
        #set facet background color to white
        theme(strip.background = element_rect(fill = "white"))+
        theme(strip.text  = element_text(size = 10,face = "bold"))+
        theme(legend.title = element_text(face = "bold"))
ggsave("/home/data5/dingyanheng/COP/pub_branch/GDF11/DEP_set_GO.pdf",width = 5,height = 10)



#--------------------------------------------------------------------------
#norm:
DEP_norm = na.omit(read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_normal",header = T))
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
O_vs_C_norm = GetDiffProt(DEP_norm,log2(2),0.05,"OverExpr","Control")
S_vs_C_norm = GetDiffProt(DEP_norm,log2(2),0.05,"sh","Control")
O_vs_S_norm = GetDiffProt(DEP_norm,log2(2),0.05,"OverExpr","sh")

DEP_norm = list(DEP_OEN_vs_CN_UP = O_vs_C_norm$UP,
                DEP_OEN_vs_CN_DN = O_vs_C_norm$DN,
                DEP_shN_vs_CN_UP = S_vs_C_norm$UP,
                DEP_shN_vs_CN_DN = S_vs_C_norm$DN,
                DEP_OEN_vs_shN_UP = O_vs_S_norm$UP,
                DEP_OEN_vs_shN_DN = O_vs_S_norm$DN)
DEP_norm_EN = lapply(DEP_norm,function(x){return(bitr(x,"SYMBOL","ENTREZID",org.Cf.eg.db)[,2])})
KEGG_DEP_norm = compareCluster(DEP_norm_EN,"enrichKEGG",
                                organism = "cfa",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 1,
                                pAdjustMethod = "BH")
KEGG_DEP_norm@compareClusterResult$C1 = substr(KEGG_DEP_norm@compareClusterResult$Cluster,
                                                 nchar(as.character(KEGG_DEP_norm@compareClusterResult$Cluster))-1,
                                                 nchar(as.character(KEGG_DEP_norm@compareClusterResult$Cluster)))
KEGG_DEP_norm@compareClusterResult$C2 = substr(KEGG_DEP_norm@compareClusterResult$Cluster,
                                                 5,
                                                 nchar(as.character(KEGG_DEP_norm@compareClusterResult$Cluster))-3)
dotplot(KEGG_DEP_norm,x="C2",showCategory = 5,includeAll=T,label_format = 50,
        color = "-log10(p.adjust)")+
            facet_grid(.~C1)+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            ggtitle("Young\nKEGG Enrichment of DEPs")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            scale_x_discrete(limits = c("OEN_vs_CN","shN_vs_CN","OEN_vs_shN"))+
            theme_bw()+
            theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
            scale_fill_met_c("Cassatt1",direction = -1)+
            theme(panel.border = element_rect(colour = "black", linewidth = 0.5))+
            theme(panel.grid = element_blank())+
            theme(panel.grid.major.y = element_line(color = "grey",linetype  = 'dashed'))+
            #set facet background color to white
            theme(strip.background = element_rect(fill = "white"))+
            theme(strip.text  = element_text(size = 10,face = "bold"))+
            theme(legend.title = element_text(face = "bold"))+
            theme(plot.title = element_text(face = "bold"))
ggsave("/home/data5/dingyanheng/COP/pub_branch/GDF11/KEGG_DEP_Young.pdf",height = 5,width = 6)


#H2O2:
DEP_H2O2 = na.omit(read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_H2O2",header = T))
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
O_vs_C_H2O2 = GetDiffProt(DEP_H2O2,log2(2),0.05,"OverExpr","Control")
S_vs_C_H2O2 = GetDiffProt(DEP_H2O2,log2(2),0.05,"sh","Control")
O_vs_S_H2O2 = GetDiffProt(DEP_H2O2,log2(2),0.05,"OverExpr","sh")

DEP_H2O2 = list(DEP_OEH_vs_CH_UP = O_vs_C_H2O2$UP,
                DEP_OEN_vs_CH_DN = O_vs_C_H2O2$DN,
                DEP_shH_vs_CH_UP = S_vs_C_H2O2$UP,
                DEP_shH_vs_CH_DN = S_vs_C_H2O2$DN,
                DEP_OEH_vs_shH_UP = O_vs_S_H2O2$UP,
                DEP_OEH_vs_shH_DN = O_vs_S_H2O2$DN)
DEP_H2O2_EN = lapply(DEP_H2O2,function(x){return(bitr(x,"SYMBOL","ENTREZID",org.Cf.eg.db)[,2])})
KEGG_DEP_H2O2 = compareCluster(DEP_H2O2,"enrichGO",
                                OrgDb = org.Cf.eg.db,
                                keyType="SYMBOL",
                                ont = "ALL",
                                pvalueCutoff = .05,
                                qvalueCutoff = 1,
                                pAdjustMethod = "BH")
KEGG_DEP_H2O2@compareClusterResult$C1 = substr(KEGG_DEP_H2O2@compareClusterResult$Cluster,
                                                 nchar(as.character(KEGG_DEP_H2O2@compareClusterResult$Cluster))-1,
                                                 nchar(as.character(KEGG_DEP_H2O2@compareClusterResult$Cluster)))
KEGG_DEP_H2O2@compareClusterResult$C2 = substr(KEGG_DEP_H2O2@compareClusterResult$Cluster,
                                                 5,
                                                 nchar(as.character(KEGG_DEP_H2O2@compareClusterResult$Cluster))-3)
dotplot(KEGG_DEP_H2O2,x="C2",showCategory = 5,includeAll=T,label_format = 50,
        color = "-log10(p.adjust)")+
            facet_grid(.~C1)+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            ggtitle("Senescence\nGO Enrichment of DEPs")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            scale_x_discrete(limits = c("OEH_vs_CH","shH_vs_CH","OEH_vs_shH"))+
            theme_bw()+
            theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
            scale_fill_met_c("Cassatt1",direction = -1)+
            theme(panel.border = element_rect(colour = "black", linewidth = 0.5))+
            theme(panel.grid = element_blank())+
            theme(panel.grid.major.y = element_line(color = "grey",linetype  = 'dashed'))+
            #set facet background color to white
            theme(strip.background = element_rect(fill = "white"))+ 
            theme(strip.text  = element_text(size = 10,face = "bold"))+
            theme(legend.title = element_text(face = "bold"))+
            theme(plot.title = element_text(face = "bold"))
ggsave("/home/data5/dingyanheng/COP/pub_branch/GDF11/GO_DEP_Senescence.pdf",height = 4,width = 6.5)