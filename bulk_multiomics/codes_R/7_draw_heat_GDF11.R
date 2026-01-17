library(ComplexHeatmap)
library(data.table)
prot = read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/data/prot")
colnames(prot) = gsub("SH","sh",gsub("_","",gsub("NC","N",colnames(prot))))

DEP_normal = read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_normal",header = T)
DEP_H2O2 = read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_H2O2",header = T)

DEP_O_vs_C_norm = DEP_normal[,c(1,2,4,14)]
DEP_S_vs_C_norm = DEP_normal[,c(1,2,3,13)]
DEP_O_vs_S_norm = DEP_normal[,c(1,2,5,15)]

DEP_O_vs_C_H2O2 = DEP_H2O2[,c(1,2,4,14)]
DEP_S_vs_C_H2O2 = DEP_H2O2[,c(1,2,3,13)]
DEP_O_vs_S_H2O2 = DEP_H2O2[,c(1,2,5,15)]

prot_design = data.frame(label = colnames(prot),
                         condition = substr(colnames(prot),1,nchar(colnames(prot))-1),
                         replicate = rep(1:3,6))

DrawHeat = function(DEP_str,Log2FCThres,PvalThres,OutDir)
{
    ID1 = substr(DEP_str,5,5)
    if(ID1=="S"){ID1='sh'}
    if(ID1=="O"){ID1="OE"}
    ID2 = substr(DEP_str,10,10)
    if(ID2=="S"){ID2='sh'}
    if(ID2=="O"){ID2="OE"}
    ID3 = substr(DEP_str,nchar(DEP_str)-3,nchar(DEP_str))
    if(ID3=="H2O2"){ID3="H"}
    if(ID3=='norm'){ID3="N"}
    Group1 = paste0(paste0(ID1,ID3),1:3)
    Group2 = paste0(paste0(ID2,ID3),1:3)
    G1 = paste0(ID1,ID3)
    G2 = paste0(ID2,ID3)
    DEP = get(DEP_str)
    colnames(DEP) = c("Name","ID","p.adj","log2FC")
    # DEP$p.adj = p.adjust(DEP$p.adj,'fdr')
    DEP$diffexpressed <- "NO"
    DEP$diffexpressed[DEP$log2FC > Log2FCThres & DEP$p.adj < PvalThres] <- "UP"
    DEP$diffexpressed[DEP$log2FC < -Log2FCThres & DEP$p.adj < PvalThres] <- "DOWN"
    DEP$p.adj[which(DEP$p.adj==0)] = .Machine$double.xmin
    DEP_list = DEP$Name[DEP$diffexpressed!="NO"]
    DEP_intensity = prot[DEP_list,c(Group1,Group2)]
    
    collis = list(Group = c(G1="#377EB8",G2="#E41A1C"),
                                        Replicate = c("1"="#DD5129","2"="#0F7BA2","3"="#43B284"))
    names(collis$Group) = c(G1,G2)
    ColFunc = colorRamp2(c(min(DEP_intensity), median(unlist(DEP_intensity)), max(DEP_intensity)), c("#377EB8", "white", "#E41A1C"),space = "sRGB")
    col_ha = HeatmapAnnotation(Group=c(rep(G1,3),rep(G2,3)),Replicate = as.character(c(1,2,3,1,2,3)),
                                col = collis)
    htp = Heatmap(DEP_intensity,name = "log2 Intensity",show_row_names = F,col = ColFunc,top_annotation = col_ha,
                heatmap_legend_param = list(legend_direction = "vertical",
                                            legend_height = unit(2.7, "cm"),
                                            legend_size = 20),
                show_column_names = F)
    pdf(OutDir,height = 6,width = 7)
    print(htp)
    dev.off()
}

DrawHeat("DEP_O_vs_C_norm",Log2FCThres = log2(2),PvalThres = 0.05,
            OutDir = "/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_O_vs_C_norm_Heat.pdf")
DrawHeat("DEP_S_vs_C_norm",Log2FCThres = log2(2),PvalThres = 0.05,
            OutDir = "/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_S_vs_C_norm_Heat.pdf")
DrawHeat("DEP_O_vs_S_norm",Log2FCThres = log2(2),PvalThres = 0.05,
            OutDir = "/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_O_vs_S_norm_Heat.pdf")
DrawHeat("DEP_O_vs_C_H2O2",Log2FCThres = log2(2),PvalThres = 0.05,
            OutDir = "/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_O_vs_C_H2O2_Heat.pdf")
DrawHeat("DEP_S_vs_C_H2O2",Log2FCThres = log2(2),PvalThres = 0.05,
            OutDir = "/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_S_vs_C_H2O2_Heat.pdf")
DrawHeat("DEP_O_vs_S_H2O2",Log2FCThres = log2(2),PvalThres = 0.05,
            OutDir = "/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_O_vs_S_H2O2_Heat.pdf")