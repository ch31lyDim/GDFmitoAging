library(clusterProfiler)
library(org.Cf.eg.db)
library(MetBrewer)
library(ggpubr)
library(ggplot2)
library(ggsci)
init = function()
{
    Log2FCThres <<- log2(2)
    pvalThres <<- 0.05
}
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
DrawRRvgo = function(GO_res,Out){
    for(ONT in c("BP","CC","MF")){
        simMatrix <- calculateSimMatrix(GO_res$ID,
                                orgdb="org.Cf.eg.db",
                                ont=ONT,
                                method="Rel")
        if(anyNA(simMatrix)){next()}
        if('numeric'%in%class(simMatrix)){next()}
        scores <- setNames(-log10(GO_res$qvalue), GO_res$ID)
        reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores=scores,
                                    threshold=0.7,
                                    orgdb="org.Cf.eg.db")
        pdf(paste0(Out,"_",ONT,"_Heat.pdf"),height = 12,width = 17)
        pht = heatmapPlot(simMatrix,
                reducedTerms,
                annotateParent=TRUE,
                annotationLabel="parentTerm",
                fontsize=10)
        print(pht)
        dev.off()

        if(nrow(simMatrix) > 3){
            pdf(paste0(Out,"_",ONT,"_Scatter.pdf"))
            pst = scatterPlot(simMatrix, reducedTerms)
            print(pst)
            dev.off()
        }
        
        pdf(paste0(Out,"_",ONT,"_TreeMap.pdf"))
        ptm = treemapPlot(reducedTerms)
        print(ptm)
        dev.off()

        pdf(paste0(Out,"_",ONT,"_WordCloud.pdf"))
        pwc = wordcloudPlot(reducedTerms, min.freq=1, colors="black")
        print(pwc)
        dev.off()
    } 
}
EnrichCore = function(GeneSet,Out,suffix="")
{
     GO_enrich = compareCluster(GeneSet,fun="enrichGO",OrgDb = org.Cf.eg.db,
                         keyType = "SYMBOL",
                         ont = "ALL",pvalueCutoff = 0.05,
                         qvalueCutoff = 1,
                         pAdjustMethod = "none")
    GO_enrich@compareClusterResult$Description[which(GO_enrich@compareClusterResult$Description=="oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen")] = "oxidoreductase activity"
    GO_enrich@compareClusterResult$Description[which(GO_enrich@compareClusterResult$Description=="oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen")] = "oxidoreductase activity"
    print(is.null(GO_enrich))
    if(!is.null(GO_enrich))
	{
        print("enter")
		dotplot(GO_enrich,x="Count",color = "-log10(p.adjust)",split = "ONTOLOGY")+
            facet_grid(ONTOLOGY~Cluster, scale="free",space = 'fixed')+
            theme_bw(base_size = 18)+
            scale_color_gradientn(colors=(rev(met.brewer("Benedictus"))))+
            theme(panel.grid.major.y = element_line(linewidth = 1.5,linetype = 'dashed'))
		OutGO = paste0("/data1/dyh/GDF11/plots/",Out,suffix,"_GO_enrich.pdf")
		ggsave(OutGO,height = 15,width = 12)
		OutGO = paste0("/data1/dyh/GDF11/plots/",Out,suffix,"_GO_enrich")
        write.table(GO_enrich@compareClusterResult,OutGO,
        col.names = T,row.names = F,quote = F,sep = '\t')

        UP = GO_enrich@compareClusterResult[GO_enrich@compareClusterResult$Cluster=="UP",]
        # if(nrow(UP)>0){
        #     UP_DIR = paste0("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/GO/",Out,suffix,"_UP")
        #     DrawRRvgo(UP,UP_DIR)
        # }
        # DN = GO_enrich@compareClusterResult[GO_enrich@compareClusterResult$Cluster=="DN",]
        # if(nrow(DN)>0){
        #     DN_DIR = paste0("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/GO/",Out,suffix,"_DN")
        #     DrawRRvgo(DN,DN_DIR)
        }
    GO_all = compareCluster(GeneSet,fun="enrichGO",OrgDb = org.Cf.eg.db,
                         keyType = "SYMBOL",
                         ont = "ALL",pvalueCutoff = 1,
                         qvalueCutoff = 1,
                         pAdjustMethod = "none")
	if(!is.null(GO_all))
	{
		 OutGO = paste0("/data1/dyh/GDF11/plots/",Out,suffix,"_GO_in")
    	 write.table(GO_all@compareClusterResult,OutGO,
            col.names = T,row.names = F,quote = F,sep = '\t') 
	}               
    GeneSet_ENTREZ = lapply(GeneSet,function(x){
                                          res = bitr(x,
                                          fromType = "SYMBOL",
                                          toType = 'ENTREZID',
                                          OrgDb = org.Cf.eg.db)
                                          return(res[,2])})
    KEGG_enrich = compareCluster(GeneSet_ENTREZ,fun = "enrichKEGG",organism = 'cfa',
                               pvalueCutoff = 0.05,
                                qvalueCutoff = 1,
                                pAdjustMethod = "BH")
	if(!is.null(KEGG_enrich))
	{
		dotplot(KEGG_enrich,x="Count",color = "-log10(p.adjust)")+
			theme_bw(base_size = 18)+
            facet_grid(.~Cluster, scale="free",space = 'fixed')+
			scale_color_gradientn(colors=(rev(met.brewer("Benedictus"))))+
			theme(panel.grid.major.y = element_line(linewidth = 1.5,linetype = 'dashed'))
		OutKEGG = paste0("/data1/dyh/GDF11/plots/",Out,suffix,"_KEGG_enrich.pdf")
		ggsave(OutKEGG)
		DEP_KEGG_res = KEGG_enrich@compareClusterResult
		OutDir = paste0("/data1/dyh/GDF11/plots/",Out,suffix,'_KEGG_enrich')
		write.table(DEP_KEGG_res,OutDir,
				col.names = T,row.names = F,quote = F,sep = '\t')
	}
    KEGG_invlove = compareCluster(GeneSet_ENTREZ,fun = "enrichKEGG",organism = 'cfa',
                                pvalueCutoff = 1,
                                qvalueCutoff = 1,
                                pAdjustMethod = "BH")
	if(!is.null(KEGG_invlove))
	{
		DEP_KEGG_res = KEGG_invlove@compareClusterResult
		OutDir = paste0("/data1/dyh/GDF11/plots/",Out,suffix,'_KEGG_involve')
		write.table(DEP_KEGG_res,OutDir,
				col.names = T,row.names = F,quote = F,sep = '\t')
	}
}
Enrich = function(DEP,Log2FCThres,pvalThres,Group1,Group2,suffix="")
{
    GeneSet = GetDiffProt(DEP,Log2FCThres,pvalThres,Group1,Group2)
    Out = paste(Group1,Group2,sep = '_')
    EnrichCore(GeneSet,Out,suffix)
}

init()
library(readxl)
#aged_vs_normal:
DEP = na.omit(as.data.frame(read_xlsx("/data1/dyh/GDF11/difftable/DEP.xlsx")))
Enrich(DEP,Log2FCThres,pvalThres,"C_H","C_NC","")
Enrich(DEP,Log2FCThres,pvalThres,"OE_H","OE_NC","")
Enrich(DEP,Log2FCThres,pvalThres,"SH_H","SH_NC","")

C = GetDiffProt(DEP,log2(2),0.05,"C_H","C_NC")
O = GetDiffProt(DEP,log2(2),0.05,"OE_H","OE_NC")
S = GetDiffProt(DEP,log2(2),0.05,"SH_H","SH_NC")

aged_vs_norm_dep = list(DEP_CH_vs_CN_UP = C$UP,
                        DEP_CH_vs_CN_DN = C$DN,
                        DEP_OEH_vs_OEN_UP = O$UP,
                        DEP_OEH_vs_OEN_DN = O$DN,
                        DEP_shH_vs_shN_UP = S$UP,
                        DEP_shH_vs_shN_DN = S$DN)
GO_aged_vs_norm_dep = compareCluster(aged_vs_norm_dep,"enrichGO",
                                OrgDb = org.Cf.eg.db,
                                keyType="SYMBOL",
                                ont = "ALL",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 1,
                                pAdjustMethod = "BH")
GO_aged_vs_norm_dep@compareClusterResult$C1 = substr(GO_aged_vs_norm_dep@compareClusterResult$Cluster,
                                                 nchar(as.character(GO_aged_vs_norm_dep@compareClusterResult$Cluster))-1,
                                                 nchar(as.character(GO_aged_vs_norm_dep@compareClusterResult$Cluster)))
GO_aged_vs_norm_dep@compareClusterResult$C2 = substr(GO_aged_vs_norm_dep@compareClusterResult$Cluster,
                                                 5,
                                                 nchar(as.character(GO_aged_vs_norm_dep@compareClusterResult$Cluster))-3)
dotplot(GO_aged_vs_norm_dep,x="C2",showCategory = 5,includeAll=T,label_format = 50,
        color = "-log10(p.adjust)")+
            theme_classic2(base_size = 18)+
            # theme(aspect.ratio = 1/3)+
            scale_color_gradientn(colors=rev(met.brewer("Cassatt1")))+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            facet_grid(ONTOLOGY~C1, scale="free_y",space = 'free_y',drop = F,shrink = F)+
            # theme(plot.margin = margin(l=1,r=0.3,unit = "cm"))+
            ggtitle("GO Enrichment of DEPs")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/GO/GO_aged_vs_norm.pdf",height = 10,width = 12)

aged_vs_norm_dep_EN = lapply(aged_vs_norm_dep,function(x){return(bitr(x,"SYMBOL","ENTREZID",org.Cf.eg.db)[,2])})

KEGG_aged_vs_norm_dep = compareCluster(aged_vs_norm_dep_EN,"enrichKEGG",
                                organism = 'cfa',
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 1,
                                pAdjustMethod = "BH")
KEGG_aged_vs_norm_dep@compareClusterResult$C1 = substr(KEGG_aged_vs_norm_dep@compareClusterResult$Cluster,
                                                 nchar(as.character(KEGG_aged_vs_norm_dep@compareClusterResult$Cluster))-1,
                                                 nchar(as.character(KEGG_aged_vs_norm_dep@compareClusterResult$Cluster)))
KEGG_aged_vs_norm_dep@compareClusterResult$C2 = substr(KEGG_aged_vs_norm_dep@compareClusterResult$Cluster,
                                                 5,
                                                 nchar(as.character(KEGG_aged_vs_norm_dep@compareClusterResult$Cluster))-3)
dotplot(KEGG_aged_vs_norm_dep,x="C2",showCategory = 5,includeAll=T,label_format = 50,
        color = "-log10(p.adjust)")+
            theme_classic2(base_size = 18)+
            # theme(aspect.ratio = 1/3)+
            scale_color_gradientn(colors=rev(met.brewer("Cassatt1")))+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            facet_grid(.~C1, scale="free_y",space = 'free_y',drop = F,shrink = F)+
            # theme(plot.margin = margin(l=1,r=0.3,unit = "cm"))+
            ggtitle("KEGG Enrichment of DEPs")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            scale_x_discrete(limits = c("CH_vs_CN","OEH_vs_OEN","shH_vs_shN"))
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/KEGG/KEGG_aged_vs_norm.pdf",height = 10,width = 12)

only_in_o = list(UP = setdiff(O$UP,union(C$UP,S$UP)),
                 DN = setdiff(O$DN,union(C$DN,S$DN)))
EnrichCore(only_in_o,"only_in_O_","")

only_in_s = list(UP = setdiff(S$UP,union(C$UP,O$UP)),
                 DN = setdiff(S$DN,union(C$DN,O$DN)))
EnrichCore(only_in_s,"only_in_S_","")

OxS_shared = list(BothUP = intersect(S$UP,O$UP),
                  BothDN = intersect(S$DN,O$DN),
                  O_DN_S_UP = intersect(S$UP,O$DN))
EnrichCore(OxS_shared,"OxS_shared_","")

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
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 1,
                                pAdjustMethod = "BH")

dotplot(GO_os_overlap_DEP,showCategory = 20,label_format = 3,color = "-log10(p.adjust)")+
            scale_x_discrete(limits = c("OE_unique_DN","OE_sh_both_DN"))+
            scale_y_discrete(limits = c("mitochondrial ribosome",
                                        "mitochondrial large ribosomal subunit",
                                        "mitochondrial protein-containing complex",
                                        "mitochondrial matrix"))+
            theme_classic2(base_size = 18)+
            # theme(aspect.ratio = 1/3)+
            scale_fill_gradientn(colors=rev(met.brewer("Cassatt1")))+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            # facet_grid(.~C1, scale="free_y",space = 'free_y',drop = F,shrink = F)+
            # theme(plot.margin = margin(l=1,r=0.3,unit = "cm"))+
            ggtitle("GO Enrichment of DEPs")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(aspect.ratio = 1/1)
ggsave("/data1/dyh/GDF11/plots/1G.pdf",height = 8,width = 8)
p<-enrichplot::cnetplot(GO_os_overlap_DEP)+
    scale_fill_manual(values = c("OE_unique_UP" = "#E63946", 
                                        "OE_unique_DN" = "#E63946",
                                        "sh_unique_UP" = "#A8DADC",
                                        "sh_unique_DN" = "#457B9D",
                                        "OE_sh_both_UP" = "#2A9D8F",
                                        "OE_sh_both_DN" = "#1D3557",
                                        "OE_DN_S_UP" = "#F1FAEE",
                                        "OE_DN_S_DN" = "#F1FAEE"))
ggsave("/data1/dyh/GDF11/plots/DEP_os_overlap_GO_cnet.pdf",height = 6.5,width = 6)
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/GO/GO_os_overlap_DEP.pdf")

os_overlap_DEP_EN = lapply(os_overlap_DEP,function(x){return(bitr(x,"SYMBOL","ENTREZID",org.Cf.eg.db)[,2])})

KEGG_os_overlap_DEP = compareCluster(os_overlap_DEP_EN,"enrichKEGG",
                                organism = 'cfa',
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 1,
                                pAdjustMethod = "BH")
dotplot(KEGG_os_overlap_DEP,x="Cluster",showCategory = 20,label_format = 50,color = "-log10(p.adjust)",)+
            # scale_x_discrete(limits = c("OE_unique_UP","OE_unique_DN",
            #                             "sh_unique_UP","sh_unique_DN",
            #                             "OE_sh_both_UP","OE_sh_both_DN",
            #                             "OE_DN_S_UP","OE_DN_S_DN"))+
            theme_classic2(base_size = 18)+
            # theme(aspect.ratio = 1/3)+
            scale_color_gradientn(colors=rev(met.brewer("Cassatt1")))+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            # facet_grid(.~C1, scale="free_y",space = 'free_y',drop = F,shrink = F)+
            # theme(plot.margin = margin(l=1,r=0.3,unit = "cm"))+
            ggtitle("KEGG Enrichment of DEPs")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/KEGG/KEGG_os_overlap_DEP.pdf")






DEP_norm = na.omit(read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_normal",header = T))
Enrich(DEP_norm,Log2FCThres,pvalThres,"sh","Control","_norm")
Enrich(DEP_norm,Log2FCThres,pvalThres,"OverExpr","Control","_norm")
Enrich(DEP_norm,Log2FCThres,pvalThres,"OverExpr","sh","_norm")

DEP_H2O2 = na.omit(read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_H2O2",header = T))
Enrich(DEP_H2O2,Log2FCThres,pvalThres,"sh","Control","_H2O2")
Enrich(DEP_H2O2,Log2FCThres,pvalThres,"OverExpr","Control","_H2O2")
Enrich(DEP_H2O2,Log2FCThres,pvalThres,"OverExpr","sh","_H2O2")

O_vs_C_norm = GetDiffProt(DEP_norm,log2(2),0.05,"OverExpr","Control")
S_vs_C_norm = GetDiffProt(DEP_norm,log2(2),0.05,"sh","Control")
O_vs_S_norm = GetDiffProt(DEP_norm,log2(2),0.05,"OverExpr","sh")


DEP_norm = list(DEP_OEN_vs_CN_UP = O_vs_C_norm$UP,
                DEP_OEN_vs_CN_DN = O_vs_C_norm$DN,
                DEP_shN_vs_CN_UP = S_vs_C_norm$UP,
                DEP_shN_vs_CN_DN = S_vs_C_norm$DN,
                DEP_OEN_vs_shN_UP = O_vs_S_norm$UP,
                DEP_OEN_vs_shN_DN = O_vs_S_norm$DN)

GO_DEP_norm = compareCluster(DEP_norm,"enrichGO",
                                OrgDb = org.Cf.eg.db,
                                keyType="SYMBOL",
                                ont = "ALL",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 1,
                                pAdjustMethod = "BH")

GO_DEP_norm@compareClusterResult$C1 = substr(GO_DEP_norm@compareClusterResult$Cluster,
                                                 nchar(as.character(GO_DEP_norm@compareClusterResult$Cluster))-1,
                                                 nchar(as.character(GO_DEP_norm@compareClusterResult$Cluster)))
GO_DEP_norm@compareClusterResult$C2 = substr(GO_DEP_norm@compareClusterResult$Cluster,
                                                 5,
                                                 nchar(as.character(GO_DEP_norm@compareClusterResult$Cluster))-3)
dotplot(GO_DEP_norm,x="C2",showCategory = 5,includeAll=T,label_format = 50,
        color = "-log10(p.adjust)")+
            theme_classic2(base_size = 18)+
            # theme(aspect.ratio = 1/3)+
            scale_color_gradientn(colors=rev(met.brewer("Cassatt1")))+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            facet_grid(ONTOLOGY~C1, scale="free_y",space = 'free_y',drop = F,shrink = F)+
            # theme(plot.margin = margin(l=1,r=0.3,unit = "cm"))+
            ggtitle("GO Enrichment of DEPs")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            scale_x_discrete(limits = c("OEN_vs_CN","shN_vs_CN","OEN_vs_shN"))
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/GO/GO_DEP_norm.pdf")

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
            theme_classic2(base_size = 18)+
            # theme(aspect.ratio = 1/3)+
            scale_color_gradientn(colors=rev(met.brewer("Cassatt1")))+
            facet_grid(.~C1)+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            # theme(plot.margin = margin(l=1,r=0.3,unit = "cm"))+
            ggtitle("KEGG Enrichment of DEPs")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            scale_x_discrete(limits = c("OEN_vs_CN","shN_vs_CN","OEN_vs_shN"))
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/KEGG/KEGG_DEP_norm.pdf")


O_vs_C_H2O2 = GetDiffProt(DEP_H2O2,log2(2),0.05,"OverExpr","Control")
S_vs_C_H2O2 = GetDiffProt(DEP_H2O2,log2(2),0.05,"sh","Control")
O_vs_S_H2O2 = GetDiffProt(DEP_H2O2,log2(2),0.05,"OverExpr","sh")


DEP_H2O2 = list(DEP_OEH_vs_CH_UP = O_vs_C_H2O2$UP,
                DEP_OEH_vs_CH_DN = O_vs_C_H2O2$DN,
                DEP_shH_vs_CH_UP = S_vs_C_H2O2$UP,
                DEP_shH_vs_CH_DN = S_vs_C_H2O2$DN,
                DEP_OEH_vs_shH_UP = O_vs_S_H2O2$UP,
                DEP_OEH_vs_shH_DN = O_vs_S_H2O2$DN)

GO_DEP_H2O2 = compareCluster(DEP_H2O2,"enrichGO",
                                OrgDb = org.Cf.eg.db,
                                keyType="SYMBOL",
                                ont = "ALL",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 1,
                                pAdjustMethod = "BH")

GO_DEP_H2O2@compareClusterResult$C1 = substr(GO_DEP_H2O2@compareClusterResult$Cluster,
                                                 nchar(as.character(GO_DEP_H2O2@compareClusterResult$Cluster))-1,
                                                 nchar(as.character(GO_DEP_H2O2@compareClusterResult$Cluster)))
GO_DEP_H2O2@compareClusterResult$C2 = substr(GO_DEP_H2O2@compareClusterResult$Cluster,
                                                 5,
                                                 nchar(as.character(GO_DEP_H2O2@compareClusterResult$Cluster))-3)
dotplot(GO_DEP_H2O2,x="C2",showCategory = 5,includeAll=T,label_format = 50,
        color = "-log10(p.adjust)")+
            theme_classic2(base_size = 18)+
            # theme(aspect.ratio = 1/3)+
            scale_color_gradientn(colors=rev(met.brewer("Cassatt1")))+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            facet_grid(ONTOLOGY~C1, scale="free_y",space = 'free_y',drop = F,shrink = F)+
            # theme(plot.margin = margin(l=1,r=0.3,unit = "cm"))+
            ggtitle("GO Enrichment of DEPs")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            scale_x_discrete(limits = c("OEH_vs_CH","shH_vs_CH","OEH_vs_shH"))
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/GO/GO_DEP_H2O2.pdf")

DEP_H2O2_EN = lapply(DEP_H2O2,function(x){return(bitr(x,"SYMBOL","ENTREZID",org.Cf.eg.db)[,2])})
KEGG_DEP_H2O2 = compareCluster(DEP_H2O2_EN,"enrichKEGG",
                                organism = "cfa",
                                pvalueCutoff = 0.05,
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
            theme_classic2(base_size = 18)+
            # theme(aspect.ratio = 1/3)+
            scale_color_gradientn(colors=rev(met.brewer("Cassatt1")))+
            facet_grid(.~C1)+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            # theme(plot.margin = margin(l=1,r=0.3,unit = "cm"))+
            ggtitle("KEGG Enrichment of DEPs")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            scale_x_discrete(limits = c("OEH_vs_CH","shH_vs_CH","OEH_vs_shH"))
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/KEGG/KEGG_DEP_H2O2.pdf")




















only_in_O_vs_C_norm = list(UP = setdiff(O_vs_C_norm$UP,S_vs_C_norm$UP),
                           DN = setdiff(O_vs_C_norm$DN,S_vs_C_norm$DN))
only_in_S_vs_C_norm = list(UP = setdiff(S_vs_C_norm$UP,O_vs_C_norm$UP),
                           DN = setdiff(S_vs_C_norm$DN,O_vs_C_norm$DN))
O_vs_C_and_S_vs_C_shared_norm = list(UP = intersect(S_vs_C_norm$UP,O_vs_C_norm$UP),
                                     DN = intersect(S_vs_C_norm$DN,O_vs_C_norm$DN))

only_in_O_vs_C_H2O2 = list(UP = setdiff(O_vs_C_H2O2$UP,S_vs_C_H2O2$UP),
                           DN = setdiff(O_vs_C_H2O2$DN,S_vs_C_H2O2$DN))
only_in_S_vs_C_H2O2 = list(UP = setdiff(S_vs_C_H2O2$UP,O_vs_C_H2O2$UP),
                           DN = setdiff(S_vs_C_H2O2$DN,O_vs_C_H2O2$DN))
O_vs_C_and_S_vs_C_shared_H2O2 = list(UP = intersect(S_vs_C_H2O2$UP,O_vs_C_H2O2$UP),
                                     DN = intersect(S_vs_C_H2O2$DN,O_vs_C_H2O2$DN))

EnrichCore(only_in_O_vs_C_norm,"only_in_O_vs_C","_norm")
EnrichCore(only_in_S_vs_C_norm,"only_in_S_vs_C","_norm")
EnrichCore(O_vs_C_and_S_vs_C_shared_norm,"O_vs_C_and_S_vs_C_shared","_norm")

EnrichCore(only_in_O_vs_C_H2O2,"only_in_O_vs_C","_H2O2")
EnrichCore(only_in_S_vs_C_H2O2,"only_in_S_vs_C","_H2O2")
EnrichCore(O_vs_C_and_S_vs_C_shared_H2O2,"O_vs_C_and_S_vs_C_shared","_H2O2")


