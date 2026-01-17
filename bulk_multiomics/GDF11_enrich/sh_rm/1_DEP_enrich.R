library(org.Cf.eg.db)
library(clusterProfiler)
library(ggplot2)
library(MetBrewer)
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
# S = GetDiffProt(DEP,log2(2),0.05,"SH_H","SH_NC")

only_in_o = list(UP = setdiff(O$UP,(C$UP)),
                 DN = setdiff(O$DN,(C$DN)))
only_in_c = list(UP = setdiff(C$UP,(O$UP)),
                 DN = setdiff(C$DN,(O$DN)))
OxC_shared = list(BothUP = intersect(C$UP,O$UP),
                  BothDN = intersect(C$DN,O$DN),
                  O_DN_C_UP = intersect(C$UP,O$DN),
                  O_UP_C_DN = intersect(C$DN,O$UP))

os_overlap_DEP = list(OE_unique_UP=only_in_o$UP,
                      OE_unique_DN=only_in_o$DN,
                      C_unique_UP=only_in_c$UP,
                      C_unique_DN=only_in_c$DN,
                      C_OE_both_UP = OxC_shared$BothUP,
                      C_OE_both_DN = OxC_shared$BothDN,
                      OE_DN_C_UP = OxC_shared$O_DN_C_UP,
                      OE_UP_C_DN = OxC_shared$O_UP_C_DN)

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
        label_format = 500,
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
        theme(legend.title = element_text(face = "bold"))+
        coord_flip()
ggsave("/home/data5/dingyanheng/COP/pub_branch/GDF11/sh_rm/DEP_set_GO.pdf",width = 15,height = 6)



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
# S_vs_C_norm = GetDiffProt(DEP_norm,log2(2),0.05,"sh","Control")
# O_vs_S_norm = GetDiffProt(DEP_norm,log2(2),0.05,"OverExpr","sh")

DEP_norm = list(DEP_OEN_vs_CN_UP = O_vs_C_norm$UP,
                DEP_OEN_vs_CN_DN = O_vs_C_norm$DN)
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
dotplot(KEGG_DEP_norm,x="Count",showCategory = 5,includeAll=T,label_format = 50,
        color = "-log10(p.adjust)")+
            facet_grid(.~C1)+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            ggtitle("Young\nKEGG Enrichment of DEPs(OEGDF vs Control)")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("Count")+
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
            theme(plot.title = element_text(face = "bold",size = 10))+
            theme(panel.spacing = unit(1, "lines"))
ggsave("/home/data5/dingyanheng/COP/pub_branch/GDF11/sh_rm/KEGG_DEP_Young.pdf",height = 3,width = 7)


#H2O2:
DEP_H2O2 = na.omit(read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_H2O2",header = T))


DEP_H2O2_fc <- DEP_H2O2$OverExpr_vs_Control_ratio
names(DEP_H2O2_fc) <- DEP_H2O2$name
DEP_H2O2_fc <- sort(DEP_H2O2_fc,decreasing = T)
KEGG_DEP_H2O2 = gseGO(OrgDb = org.Cf.eg.db,
                                geneList = DEP_H2O2_fc,
                                keyType="SYMBOL",
                                ont = "CC",
                                pvalueCutoff = 1,
                                pAdjustMethod = "BH")

dotplot(KEGG_DEP_H2O2)
library(GseaVis)
gseaNb(
    object= KEGG_DEP_H2O2,
    geneSetID= "GO:0099023",
    subPlot= 3,
    addPval= T,
    pvalX= 0.95,
    pvalY= 0.8,
    curveCol = met.brewer("Cassatt1",20)[c(20,8,5,1)],
    htCol = met.brewer("Cassatt1",20)[c(20,8,5,1)],
    rankCol = met.brewer("Cassatt1",20)[c(20,10,1)],
    base_size = 9
)

ggsave("/home/data5/dingyanheng/COP/pub_branch/GDF11/sh_rm/GO_DEP_Senescence.pdf",width = 5,height = 4)
