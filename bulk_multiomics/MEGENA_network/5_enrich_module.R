library(clusterProfiler)
library(org.Cf.eg.db)
library(MetBrewer)
library(ggplot2)
library(ggpubr)
library(enrichplot)
mod = readRDS("/data1/dyh/GDF11/MEGENA/module.rds")
same_level = mod$modules[which(mod$module.table$module.parent=="c1_165")]
# same_level$c1_165 = mod$modules$c1_165
GO_en = compareCluster(same_level,"enrichGO",
                       ont = "ALL",
                        keyType = "SYMBOL",
                       OrgDb = org.Cf.eg.db,
                       qvalueCutoff = 1,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.1)
dotplot(GO_en,showCategory = 5,includeAll=T,color = "-log10(p.adjust)",label_format = 50)+
            theme_classic2(base_size = 18)+
            # theme(aspect.ratio = 1/3)+
            scale_color_gradientn(colors=rev(met.brewer("Cassatt1")))+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            facet_grid(ONTOLOGY~., scale="free_y",space = 'free_y',drop = F,shrink = F)+
            # theme(plot.margin = margin(l=1,r=0.3,unit = "cm"))+
            ggtitle("GO Enrichment of Modules")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            scale_x_discrete(limits = names(same_level))
ggsave("/home/data3/dingyanheng/COP/MEGENA/module_enrich_GO.pdf")


GO_en = pairwise_termsim(GO_en)
enrichplot::emapplot(GO_en,alpha=0.2,lwd=0,cluster.params = list(cluster = T,legend = F,label_format = 5),
            pie.params = list(legend_n = 2),shadowtext = F,node_label = 'group',
            cex.params = list(label_group = 2,line = 0.2,category_node=3),layout.params = list(layout = "kk"),
            hilight.params = list(category = "c1_475"))+
    scale_fill_manual(values = met.brewer("Redon",6))+
    theme(legend.position = "bottom")+
    theme(legend.title = element_blank(),legend.text = element_text(size = 18))
ggsave("/home/data3/dingyanheng/COP/MEGENA/module_enrich_GO_ema.pdf")














same_level_en = lapply(same_level,function(x){return(bitr(x,"SYMBOL","ENTREZID",org.Cf.eg.db)[,2])})
KEGG_en = compareCluster(same_level_en,"enrichKEGG",
                        organism  = "cfa",
                        qvalueCutoff = 1,
                        pAdjustMethod = "none",
                        pvalueCutoff = 0.01)
dotplot(KEGG_en,showCategory = 60,includeAll=T,color = "-log10(p.adjust)",label_format = 50)+
scale_color_gradientn(colors=rev(met.brewer("Cassatt1")))+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            # theme(plot.margin = margin(l=1,r=0.3,unit = "cm"))+
            ggtitle("KEGG Enrichment of Modules")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/data3/dingyanheng/COP/MEGENA/module_enrich_KEGG.pdf",height = 12,width = 10)
