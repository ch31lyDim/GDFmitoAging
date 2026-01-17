
library(pbmcapply)
library(DiscriMiner)
library(ggsci)
library(ggrepel)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(clusterProfiler)
library(MetBrewer)
library(ggpubr)
PreProcess = function(mode){
    intensity = t(read.table(paste0("/home/data5/dingyanheng/COP/GDF11/Metabolome/data/data_norm_shift_",mode)))
    # intensity = intensity[,-(1:3)]
    Design <<- data.frame(sample = colnames(intensity),
                        Group = substr(colnames(intensity),1,nchar(colnames(intensity))-2),
                        Genotype = substr(colnames(intensity),1,nchar(colnames(intensity))-3),
                        Condition = substr(colnames(intensity),nchar(colnames(intensity))-2,nchar(colnames(intensity))-2),
                        Replicate = substr(colnames(intensity),nchar(colnames(intensity)),nchar(colnames(intensity))))
    MetaboName = as.data.frame(fread(paste0("/home/data5/dingyanheng/COP/GDF11/Metabolome/data/metabo_",mode)))
    MetaboName$Name = iconv(MetaboName$Name, from = 'UTF-8', to = 'ASCII//TRANSLIT')
    MetaboName = MetaboName[!is.na(MetaboName$Name),]
    colnames(MetaboName)[1] = 'ID'
    intensity = intensity[rownames(intensity)%in%MetaboName$ID,]
    return(list(MetaboName,intensity))
}
intensity_neg = PreProcess("neg")[[2]]
intensity_pos = PreProcess("pos")[[2]]
MetaboName_neg = PreProcess("neg")[[1]]
MetaboName_pos = PreProcess("pos")[[1]]

TestCore = function(i,df1,df2){
    Exp1 = df1[i,]
    Exp2 = df2[i,]
    test_res = t.test(Exp1,Exp2)
    pval = test_res$p.value
    log2fc = mean(Exp1)-mean(Exp2)
    return(c(log2fc,pval))
}
TestGroup = function(Group1,Group2,test_sample1,test_sample2,mode){
    intensity = get(paste0("intensity_",mode))
    test_df1 = intensity[sort(rownames(intensity)),test_sample1]
    test_df2 = intensity[sort(rownames(intensity)),test_sample2]
    TestRes = t(sapply(1:nrow(test_df1),FUN=TestCore,df1=test_df1,df2=test_df2))
    res = data.frame(ID = rownames(test_df1),log2FC = TestRes[,1],pval = TestRes[,2])
    plsda = plsDA(t(cbind(test_df1,test_df2)),c(rep(Group1,length(test_sample1)),rep(Group2,length(test_sample2))),
                  cv = "LOO")
    res$VIP = plsda$VIP[,ncol(plsda$VIP)]
    MetaboName = get(paste0("MetaboName_",mode))
    res4Print = full_join(res,MetaboName[,c("ID","Name","Kegg_ID")])
    return(list(res,res4Print))
}
TestDiff = function(Group1,Group2){
    test_sample1 = Design$sample[Design$Group==Group1]
    test_sample2 = Design$sample[Design$Group==Group2]
    res_neg = TestGroup(Group1,Group2,test_sample1,test_sample2,'neg')[[1]]
    res_pos = TestGroup(Group1,Group2,test_sample1,test_sample2,'pos')[[1]]
    res = rbind(res_neg,res_pos)
    Out = paste0("/home/data5/dingyanheng/COP/GDF11/Metabolome/res/DiffMets/",Group1,"_",Group2,"_DEM")
    res4Print_neg = TestGroup(Group1,Group2,test_sample1,test_sample2,'neg')[[2]]
    res4Print_pos = TestGroup(Group1,Group2,test_sample1,test_sample2,'pos')[[2]]
    res4Print = rbind(res4Print_neg,res4Print_pos)
    write.table(res4Print,Out,row.names = F,col.names = T,quote = F,sep = '\t')
    return(res)
}

GetDiffMet = function(Group1,Group2,Log2FCThres,pvalThres,VIPThres)
{
    diff_df = inner_join(TestDiff(Group1,Group2),MetaboName)
    UP = na.omit(diff_df[diff_df$log2FC>Log2FCThres & diff_df$pval<pvalThres & diff_df$VIP>VIPThres,'Kegg_ID'])
    UP = gsub("cpd:","",UP)
    DN = na.omit(diff_df[diff_df$log2FC< -Log2FCThres & diff_df$pval<pvalThres & diff_df$VIP>VIPThres,'Kegg_ID'])
    DN = gsub("cpd:","",DN)
    Diff_list = list(UP = UP,DN = DN)
    return(Diff_list)
}
KEGG_background = as.data.frame(fread("/home/data5/dingyanheng/COP/GeneAnno/dog_KEGGtoMet"))
TERM2GENE = KEGG_background[,c(1,3)]
TERM2NAME = unique(KEGG_background[,c(1,2)])
MetaboName = rbind(MetaboName_neg,MetaboName_pos)
MetaboName = MetaboName[,c("ID","Name","Kegg_ID")]

#norm
MetSet_OEN_vs_CN = GetDiffMet("OEN","CN",log2(1.2),0.05,1)
MetSet_shN_vs_CN = GetDiffMet("shN","CN",log2(1.2),0.05,1)
MetSet_OEN_vs_shN = GetDiffMet("OEN","shN",log2(1.2),0.05,1)

DEM_norm = list(DEM_OEN_vs_CN_UP = MetSet_OEN_vs_CN$UP,
                DEM_OEN_vs_CN_DN = MetSet_OEN_vs_CN$DN,
                DEM_shN_vs_CN_UP = MetSet_shN_vs_CN$UP,
                DEM_shN_vs_CN_DN = MetSet_shN_vs_CN$DN,
                DEM_OEN_vs_shN_UP = MetSet_OEN_vs_shN$UP,
                DEM_OEN_vs_shN_DN = MetSet_OEN_vs_shN$DN)

DEM_norm_KEGG = compareCluster(DEM_norm,"enricher",
            TERM2GENE=TERM2GENE,TERM2NAME=TERM2NAME,
            pvalueCutoff = 0.05,qvalueCutoff = 1,pAdjustMethod = 'BH')


DEM_norm_KEGG@compareClusterResult$C1 = substr(DEM_norm_KEGG@compareClusterResult$Cluster,
                                                 nchar(as.character(DEM_norm_KEGG@compareClusterResult$Cluster))-1,
                                                 nchar(as.character(DEM_norm_KEGG@compareClusterResult$Cluster)))
DEM_norm_KEGG@compareClusterResult$C2 = substr(DEM_norm_KEGG@compareClusterResult$Cluster,
                                                 5,
                                                 nchar(as.character(DEM_norm_KEGG@compareClusterResult$Cluster))-3)
p_norm <- dotplot(DEM_norm_KEGG,x="C2",showCategory = 5,includeAll=T,label_format = 50,
        color = "-log10(p.adjust)")+
        
            theme_bw()+
            scale_fill_gradientn(colors=rev(met.brewer("Cassatt1")))+
            facet_grid(.~C1,shrink = F,drop = F)+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            ggtitle("Young\nKEGG Enrichment of DEMs")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(face = 'bold'))+
            theme(legend.title = element_text(face = 'bold'))+
            theme(strip.text = element_text(face = 'bold'))+
            theme(aspect.ratio = 3/1)

#Senescence;
MetSet_OEH_vs_CH = GetDiffMet("OEH","CH",log2(1.2),0.05,1)
MetSet_shH_vs_CH = GetDiffMet("shH","CH",log2(1.2),0.05,1)
MetSet_OEH_vs_shH = GetDiffMet("OEH","shH",log2(1.2),0.05,1)

DEM_H2O2 = list(DEM_OEH_vs_CH_UP = MetSet_OEH_vs_CH$UP,
                DEM_OEH_vs_CH_DN = MetSet_OEH_vs_CH$DN,
                DEM_shH_vs_CH_UP = MetSet_shH_vs_CH$UP,
                DEM_shH_vs_CH_DN = MetSet_shH_vs_CH$DN,
                DEM_OEH_vs_shH_UP = MetSet_OEH_vs_shH$UP,
                DEM_OEH_vs_shH_DN = MetSet_OEH_vs_shH$DN)

DEM_H2O2_KEGG = compareCluster(DEM_H2O2,"enricher",
            TERM2GENE=TERM2GENE,TERM2NAME=TERM2NAME,
            pvalueCutoff = 0.05,qvalueCutoff = 1,pAdjustMethod = 'BH')


DEM_H2O2_KEGG@compareClusterResult$C1 = substr(DEM_H2O2_KEGG@compareClusterResult$Cluster,
                                                 nchar(as.character(DEM_H2O2_KEGG@compareClusterResult$Cluster))-1,
                                                 nchar(as.character(DEM_H2O2_KEGG@compareClusterResult$Cluster)))
DEM_H2O2_KEGG@compareClusterResult$C2 = substr(DEM_H2O2_KEGG@compareClusterResult$Cluster,
                                                 5,
                                                 nchar(as.character(DEM_H2O2_KEGG@compareClusterResult$Cluster))-3)
p_h2o2<-dotplot(DEM_H2O2_KEGG,x="C2",showCategory = 5,includeAll=T,label_format = 50,
        color = "-log10(p.adjust)")+
            theme_bw()+
            scale_fill_gradientn(colors=rev(met.brewer("Cassatt1")))+
            facet_grid(.~C1,shrink = F,drop = F)+
            theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1))+
            ggtitle("Senescence\nKEGG Enrichment of DEMs")+
            theme(strip.background = element_rect(fill = "White"))+
            xlab("")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(face = 'bold'))+
            theme(strip.text = element_text(face = 'bold'))+
            theme(aspect.ratio = 3/1)

DEMPlot <- ggpubr::ggarrange(p_norm,p_h2o2,ncol = 2,nrow = 1,common.legend = T)+
    theme(legend.position = 'bottom')

ggsave("/home/data5/dingyanheng/COP/pub_branch/GDF11/DEM_KEGG.pdf",width = 13,height = 5)
