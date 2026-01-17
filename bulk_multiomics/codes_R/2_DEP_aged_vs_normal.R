library(data.table)
library(tidyverse)
library(DEP)
library(uwot)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(ggplot2)
library(ggrepel)
library(circlize)
library(ggrepel)
library(ggforce)
prot = as.data.frame(fread("/home/data5/dingyanheng/COP/GDF11/Proteome/data/all_sample.xls"))
# prot = prot[,-2]
data_unique <- make_unique(prot,"Gene","Protein",delim = ";")
# data_unique = data_unique[data_unique$Gene!='',]
LFQ_cols = 4:21
prot_design = data.frame(label = colnames(data_unique)[LFQ_cols],
                         condition = substr(colnames(data_unique)[LFQ_cols],1,nchar(colnames(data_unique)[LFQ_cols])-3),
                         replicate = rep(1:3,6))
data_se = make_se(data_unique,LFQ_cols,prot_design)
data_filter <- filter_missval(data_se,thr = 3)
data_norm = normalize_vsn(data_se)
data_imp <- impute(data_norm,fun = "MinProb",q=0.01)
data_diff <- test_diff(data_imp, type = "manual",
                       test = c("C_H_vs_C_NC","SH_H_vs_SH_NC","OE_H_vs_OE_NC"))
dep <- add_rejections(data_diff,alpha = 0.05,lfc = log2(2))
data_results <- get_results(dep)

DEP_CH_vs_CN_UP = data_results$name[data_results$C_H_vs_C_NC_p.val<0.05 & data_results$C_H_vs_C_NC_ratio >1]
DEP_CH_vs_CN_DN = data_results$name[data_results$C_H_vs_C_NC_p.val<0.05 & data_results$C_H_vs_C_NC_ratio < -1]

DEP_CH_vs_CN_UP = data_results$name[data_results$C_H_vs_C_NC_p.val<0.05 & data_results$C_H_vs_C_NC_ratio >1]
DEP_CH_vs_CN_DN = data_results$name[data_results$C_H_vs_C_NC_p.val<0.05 & data_results$C_H_vs_C_NC_ratio < -1]

DEP_CH_vs_CN_UP = data_results$name[data_results$C_H_vs_C_NC_p.val<0.05 & data_results$C_H_vs_C_NC_ratio >1]
DEP_CH_vs_CN_DN = data_results$name[data_results$C_H_vs_C_NC_p.val<0.05 & data_results$C_H_vs_C_NC_ratio < -1]

prot_mat = assay(data_imp)
write.table(data_results,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP.txt",row.names = F,col.names = T,sep = '\t',quote = F)
write.table(prot_mat,"/home/data5/dingyanheng/COP/GDF11/Proteome/data/prot",
            col.names = T,row.names = T,quote = F,sep = '\t')
#UMAP
UmapRes = umap(t(prot_mat),n_threads = 78,n_sgd_threads = 78,n_neighbors = 5,n_components = 2)
rownames(UmapRes) = gsub('_',"",rownames(UmapRes))
rownames(UmapRes) = gsub('NC',"N",rownames(UmapRes))
rownames(UmapRes) = gsub('SH',"sh",rownames(UmapRes))
UmapPlotDf = data.frame(UMAP1 = UmapRes[,1], UMAP2 = UmapRes[,2],
                        Sample = rownames(UmapRes), 
                        Group = substr(rownames(UmapRes),1,nchar(rownames(UmapRes))-1),
                        Replicate = substr(rownames(UmapRes),nchar(rownames(UmapRes)),nchar(rownames(UmapRes))))

ggplot(data=UmapPlotDf,mapping = aes(x = UMAP1, y = UMAP2, color = Group))+
    geom_point(size = 6 )+
    theme_bw()+
    theme(aspect.ratio = 1/1)+
    scale_color_nejm()+
    scale_fill_nejm()+
    theme(axis.title = element_text(size=18))+
    theme(legend.text = element_text(size=18))+
    theme(axis.text = element_text(size=18))+
    theme(legend.title =  element_text(size=18))+
    guides(color = guide_legend(order = 1) , shape = guide_legend(order = 0))+
    # theme(panel.grid = element_blank())
    geom_text_repel(mapping = aes(label = Sample),max.overlaps = Inf,na.rm = T,verbose = T,color='black',size=4)+
    geom_mark_ellipse(mapping = aes(color = Group,fill = Group),show.legend = F,size=0,linetype = 0)
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/UMAP.pdf",height = 8,width = 10)

#PCA
PCA_res = prcomp(prot_mat)
rownames(PCA_res$rotation) = gsub('_',"",rownames(PCA_res$rotation))
rownames(PCA_res$rotation) = gsub('NC',"N",rownames(PCA_res$rotation))
rownames(PCA_res$rotation) = gsub('SH',"sh",rownames(PCA_res$rotation))
PCAPlotDf = data.frame(PC1 = PCA_res$rotation[,1], PC2 = PCA_res$rotation[,2],
                       Sample = rownames(PCA_res$rotation), 
                       Group = substr(rownames(PCA_res$rotation),1,nchar(rownames(PCA_res$rotation))-1),
                       Replicate = substr(rownames(PCA_res$rotation),nchar(rownames(PCA_res$rotation)),nchar(rownames(PCA_res$rotation))))
variance <- PCA_res$sdev^2
variance_ratio <- paste0(round((100* variance / sum(variance)),2),"%")

ggplot(data=PCAPlotDf,mapping = aes(x = PC1, y = PC2, color = Group))+
    geom_point(size = 6 )+
    theme_bw()+
    theme(aspect.ratio = 1/1)+
    scale_color_nejm()+
    scale_fill_nejm()+
    theme(axis.title = element_text(size=18))+
    theme(legend.text = element_text(size=18))+
    theme(axis.text = element_text(size=18))+
    theme(legend.title =  element_text(size=18))+
    guides(color = guide_legend(order = 1) , shape = guide_legend(order = 0))+
    # theme(panel.grid = element_blank())+
    xlab(paste0("PC1 ","(",variance_ratio[1],")"))+
    ylab(paste0("PC2 ","(",variance_ratio[2],")"))+
    geom_text_repel(mapping = aes(label = Sample),max.overlaps = Inf,na.rm = T,verbose = T,color='black',size=4)+
    geom_mark_ellipse(mapping = aes(color = Group,fill = Group),show.legend = F,size=0,linetype = 0)
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/PCA.pdf",height = 8,width = 10)
#volcano:
DEP = read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP.txt",sep = '\t',header = T)
DEP_c = dplyr::select(DEP,name,C_H_vs_C_NC_ratio,C_H_vs_C_NC_p.val)
DEP_o = dplyr::select(DEP,name,OE_H_vs_OE_NC_ratio,OE_H_vs_OE_NC_p.val)
DEP_s = dplyr::select(DEP,name,SH_H_vs_SH_NC_ratio,SH_H_vs_SH_NC_p.val)

DrawVolcano = function(DEP,fcshres,pvalThres,Out)
{
    colnames(DEP) = c("GeneID","log2FC","p.adj")
    # DEP$p.adj = p.adjust(DEP$p.adj,'fdr')
    DEP$diffexpressed <- "NO"
    DEP$diffexpressed[DEP$log2FC > log2(fcshres) & DEP$p.adj < pvalThres] <- "UP"
    DEP$diffexpressed[DEP$log2FC < -log2(fcshres) & DEP$p.adj < pvalThres] <- "DOWN"
    DEP$p.adj[which(DEP$p.adj==0)] = .Machine$double.xmin
    DEP$delabel <- ifelse(DEP$GeneID %in% c(head(DEP[DEP$diffexpressed=="DOWN",][order(DEP$p.adj[DEP$diffexpressed=="DOWN"]), "GeneID"], 10),
                                            head(DEP[DEP$diffexpressed=="UP",][order(DEP$p.adj[DEP$diffexpressed=="UP"]), "GeneID"], 10)),DEP$GeneID, NA)
    DEP$diffexpressed = factor(DEP$diffexpressed,levels = c("UP","DOWN","NO"))

    ggplot(data = DEP, aes(x = log2FC, y = -log10(p.adj), color = diffexpressed,label = delabel)) +
        geom_vline(xintercept = c(-log2(fcshres), log2(fcshres)), col = "gray", linetype = 'dashed') +
        geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
        geom_point(size = 1.5,alpha=0.5)+
        theme_classic(base_size = 20)+
        theme(axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
            axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
            plot.title = element_text(hjust = 0.5))+
        theme(aspect.ratio = 1/1)+
        scale_color_manual(values = c(pal_nejm()(2)[1], pal_nejm()(2)[2],"grey"),
                        labels = c(paste("Upregulated"," (",sum(DEP$diffexpressed=="UP"),")",sep = ""),
                                    paste("Downregulated"," (",sum(DEP$diffexpressed=="DOWN"),")",sep = ""),
                                    paste("Not significant"," (",sum(DEP$diffexpressed=="NO"),")",sep = "")))+
        # coord_cartesian(ylim = c(0, 320), xlim = c(-8, 11))+
        labs(color = '',
        x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"))+
        # scale_x_continuous(breaks = seq(-10, 10, 2))+
        # theme(legend.position = c(1.25, .5),legend.box.margin = margin(r=200))+
        geom_text_repel(max.overlaps = Inf,na.rm = T,verbose = T,max.time = 2.5)
    ggsave(Out,width = 10,height = 10)
}
DrawVolcano(DEP_c,2,0.05,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/Volcano_c.pdf")
DrawVolcano(DEP_o,2,0.05,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/Volcano_o.pdf")
DrawVolcano(DEP_s,2,0.05,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/Volcano_s.pdf")

#DEP_heat
DrawHeat = function(Identifier,DEP,fcshres,Out)
{
    colnames(DEP) = c("GeneID","log2FC","p.adj")
    # DEP$p.adj = p.adjust(DEP$p.adj,'fdr')
    DEP$diffexpressed <- "NO"
    DEP$diffexpressed[DEP$log2FC > log2(fcshres) & DEP$p.adj < 0.05] <- "UP"
    DEP$diffexpressed[DEP$log2FC < -log2(fcshres) & DEP$p.adj < 0.05] <- "DOWN"
    DEP_list = DEP$GeneID[DEP$diffexpressed!="NO"]
    if(Identifier == 'C'){DEP_intensity = assay(data_imp)[DEP_list,1:6]}
    if(Identifier == 'S'){DEP_intensity = assay(data_imp)[DEP_list,7:12]}
    if(Identifier == 'O'){DEP_intensity = assay(data_imp)[DEP_list,13:18]}
    ColFunc = colorRamp2(c(min(DEP_intensity), median(unlist(DEP_intensity)), max(DEP_intensity)), c("#377EB8", "white", "#E41A1C"),space = "sRGB")
    col_ha = HeatmapAnnotation(Group=c(rep("Normal",3),rep("H2O2",3)),Replicate = as.character(c(1,2,3,1,2,3)),
                            col = list(Group = c("H2O2"="#377EB8","Normal"="#E41A1C"),
                                        Replicate = c("1"="#DD5129","2"="#0F7BA2","3"="#43B284"))
                                )
    htp = Heatmap(DEP_intensity,name = "log2 Intensity",show_row_names = FALSE,col = ColFunc,top_annotation = col_ha,
                heatmap_legend_param = list(legend_direction = "vertical",
                                            legend_height = unit(2.7, "cm"),
                                            legend_size = 20),
                show_column_names = F,clustering_method_columns = "ward.D2")
    pdf(Out,height = 6,width = 7)
    print(htp)
    dev.off()
}
DrawHeat("C",DEP_c,2,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/Heat_c.pdf")
DrawHeat("S",DEP_s,2,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/Heat_s.pdf")
DrawHeat("O",DEP_o,2,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/Heat_o.pdf")

