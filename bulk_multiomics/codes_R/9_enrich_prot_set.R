library(clusterProfiler)
library(org.Cf.eg.db)
library(MetBrewer)
library(ggpubr)
library(ggplot2)
library(ggsci)
library(rrvgo)
library(ggvenn)
EnrichCore = function(GeneSet,Out,suffix="")
{
     GO_enrich = enrichGO(GeneSet,OrgDb = org.Cf.eg.db,
                         keyType = "SYMBOL",
                         ont = "ALL",pvalueCutoff = 0.05,
                         qvalueCutoff = 1,
                         pAdjustMethod = "none")
    GO_enrich@result$Description[which(GO_enrich@result$Description=="oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen")] = "oxidoreductase activity"
    GO_enrich@result$Description[which(GO_enrich@result$Description=="oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen")] = "oxidoreductase activity"
    if(!is.null(GO_enrich))
	{
		dotplot(GO_enrich,split = "ONTOLOGY")+
            facet_grid(ONTOLOGY~., scale="free",space = 'free')+
            theme_bw(base_size = 18)+
            scale_color_gradientn(colors=(met.brewer("Benedictus")))+
            theme(panel.grid.major.y = element_line(linewidth = 1.5,linetype = 'dashed'))
		OutGO = paste0("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/GO/",Out,suffix,"_GO_enrich.pdf")
		ggsave(OutGO,height = 15,width = 12)
		OutGO = paste0("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/GO/",Out,suffix,"_GO_enrich")
        write.table(GO_enrich@result,OutGO,
        col.names = T,row.names = F,quote = F,sep = '\t')
    }
    GeneSet_ENTREZ = res = bitr(GeneSet,
                                fromType = "SYMBOL",
                                toType = 'ENTREZID',
                                OrgDb = org.Cf.eg.db)
                                          
    KEGG_enrich = enrichKEGG(GeneSet_ENTREZ$ENTREZID,organism = 'cfa',
                               pvalueCutoff = 0.05,
                                qvalueCutoff = 1,
                                pAdjustMethod = "none")
	if(!is.null(KEGG_enrich))
	{
		dotplot(KEGG_enrich)+
			theme_bw(base_size = 18)+
			scale_color_gradientn(colors=(rev(met.brewer("Benedictus"))))+
			theme(panel.grid.major.y = element_line(linewidth = 1.5,linetype = 'dashed'))
		OutKEGG = paste0("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/KEGG/",Out,suffix,"_KEGG_enrich.pdf")
		ggsave(OutKEGG)
		DEP_KEGG_res = KEGG_enrich@result
		OutDir = paste0("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/KEGG/",Out,suffix,'_KEGG_enrich')
		write.table(DEP_KEGG_res,OutDir,
				col.names = T,row.names = F,quote = F,sep = '\t')
	}
    KEGG_invlove = enrichKEGG(GeneSet_ENTREZ$ENTREZID,organism = 'cfa',
                               pvalueCutoff = 1,
                                qvalueCutoff = 1,
                                pAdjustMethod = "none")
	if(!is.null(KEGG_invlove))
	{
		DEP_KEGG_res = KEGG_invlove@result
		OutDir = paste0("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/KEGG/",Out,suffix,'_KEGG_involve')
		write.table(DEP_KEGG_res,OutDir,
				col.names = T,row.names = F,quote = F,sep = '\t')
	}
}
O_vs_C_DEP_diff = unlist(read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/O_vs_C_DEP_diff"))
O_vs_C_DEP_same = unlist(read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/O_vs_C_DEP_same"))
S_vs_C_DEP_diff = unlist(read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/S_vs_C_DEP_diff"))
S_vs_C_DEP_same = unlist(read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/S_vs_C_DEP_same"))

EnrichCore(O_vs_C_DEP_diff,"O_vs_C_DEP_diff","")
EnrichCore(O_vs_C_DEP_same,"O_vs_C_DEP_same","")
EnrichCore(S_vs_C_DEP_diff,"S_vs_C_DEP_diff","")
EnrichCore(S_vs_C_DEP_same,"S_vs_C_DEP_same","")

O_UP_C_DN = "SEPTIN5"
O_DN_C_UP = setdiff(O_vs_C_DEP_diff,O_UP_C_DN)

S_UP_C_DN = c("CXXC1","SLC7A11","TMEM254")
S_DN_C_UP = setdiff(S_vs_C_DEP_diff,S_UP_C_DN)


diff_reg_dep = list()
diff_reg_dep = list(O_UP_C_DN=O_UP_C_DN,
                    O_DN_C_UP=O_DN_C_UP,
                    S_UP_C_DN=S_UP_C_DN,
                    S_DN_C_UP=S_DN_C_UP)
ggvenn(diff_reg_dep,text_size = 6,fill_alpha = 0.2,stroke_size = 0)+
scale_fill_lancet()
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/diff_reg_dep.png")


GO_diff_reg_dep = compareCluster(diff_reg_dep,"enrichGO",
                                OrgDb = org.Cf.eg.db,
                                keyType="SYMBOL",
                                ont = "ALL",
                                pvalueCutoff = 1,
                                qvalueCutoff = 1,
                                pAdjustMethod = "none")
dotplot(GO_diff_reg_dep,color = "-log10(pvalue)",size='Count')+
            theme_classic2(base_size = 18)+
            # theme(aspect.ratio = 1/3)+
            scale_color_gradientn(colors=rev(met.brewer("Cassatt1")))+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            # facet_grid(.~Cluster)
            # theme(plot.margin = margin(l=1,r=0.3,unit = "cm"))+
            ggtitle("GO Annotation of Diff regulated DEPs")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/GO/GO_diff_reg_dep.pdf")

diff_reg_dep_EN = lapply(diff_reg_dep,function(x){return(bitr(x,"SYMBOL","ENTREZID",org.Cf.eg.db)[,2])})
KEGG_diff_reg_dep = compareCluster(diff_reg_dep_EN,"enrichKEGG",
                                organism = 'cfa',
                                pvalueCutoff = 1,
                                qvalueCutoff = 1,
                                pAdjustMethod = "none")
dotplot(KEGG_diff_reg_dep,color = "-log10(pvalue)",size='Count')+
            theme_classic2(base_size = 18)+
            # theme(aspect.ratio = 1/3)+
            scale_color_gradientn(colors=rev(met.brewer("Cassatt1")))+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            # facet_grid(.~Cluster)
            # theme(plot.margin = margin(l=1,r=0.3,unit = "cm"))+
            ggtitle("KEGG Annotation of Diff regulated DEPs")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/Enrich/KEGG/KEGG_diff_reg_dep.pdf")
