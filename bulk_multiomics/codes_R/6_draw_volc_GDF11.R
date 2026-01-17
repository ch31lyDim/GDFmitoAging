library(data.table)
prot = read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/data/prot")
DEP_normal = read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_normal",header = T)
DEP_H2O2 = read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_H2O2",header = T)

DEP_o_vs_c_norm = DEP_normal[,c(1,2,3,13)]
DEP_s_vs_c_norm = DEP_normal[,c(1,2,5,15)]
DEP_o_vs_s_norm = DEP_normal[,c(1,2,4,14)]

DEP_o_vs_c_H2O2 = DEP_H2O2[,c(1,2,3,13)]
DEP_s_vs_c_H2O2 = DEP_H2O2[,c(1,2,5,15)]
DEP_o_vs_s_H2O2 = DEP_H2O2[,c(1,2,4,14)]
DrawVolcano = function(DEP,Log2FCThres,PvalThres,OutDir)
{
    colnames(DEP) = c("Name","ID","pval","log2FC")
    # DEP$pval = p.adjust(DEP$pval,'fdr')
    DEP$diffexpressed <- "NO"
    DEP$diffexpressed[DEP$pval < PvalThres & DEP$log2FC > Log2FCThres] <- "UP"
    DEP$diffexpressed[DEP$pval < PvalThres  & DEP$log2FC < -Log2FCThres] <- "DOWN"
    DEP$delabel <- ifelse(DEP$ID %in% c(head(DEP[DEP$diffexpressed=="DOWN",][order(DEP$pval[DEP$diffexpressed=="DOWN"]), "ID"], 10),
                                            head(DEP[DEP$diffexpressed=="UP",][order(DEP$pval[DEP$diffexpressed=="UP"]), "ID"], 10)),DEP$Name, NA)
    DEP$delabel = ifelse(nchar(DEP$delabel)>20,NA,DEP$delabel)
    DEP$diffexpressed = factor(DEP$diffexpressed,levels = c("UP","DOWN","NO"))

    DEP$FoldChange = 2^DEP$log2FC
    ggplot(data = DEP, aes(x = log2FC, y = -log10(pval), color = diffexpressed,label = delabel))+
        geom_vline(xintercept = c(-Log2FCThres,Log2FCThres), col = "gray", linetype = 'dashed')+
        geom_hline(yintercept = -log10(PvalThres), col = "gray", linetype = 'dashed')+
        geom_point(alpha=0.5)+
        theme_classic(base_size = 20)+
        theme(axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
                axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
                plot.title = element_text(hjust = 0.5))+
        theme(aspect.ratio = 1/1)+
        scale_color_manual(values = c(pal_nejm()(2)[1], pal_nejm()(2)[2],"grey"),
                        labels = c(paste("Upregulated"," (",sum(DEP$diffexpressed=="UP"),")",sep = ""),
                                    paste("Downregulated"," (",sum(DEP$diffexpressed=="DOWN"),")",sep = ""),
                                    paste("Not significant"," (",sum(DEP$diffexpressed=="NO"),")",sep = "")))+
        labs(color = '',
            x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"))+
        geom_text_repel(max.overlaps = Inf,na.rm = T,verbose = T)
ggsave(OutDir,width = 7.5,height =7.5)
}
DrawVolcano(DEP_o_vs_c_norm,Log2FCThres = log2(2),PvalThres = 0.05,
            OutDir = "/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_o_vs_c_norm_volcano.pdf")
DrawVolcano(DEP_s_vs_c_norm,Log2FCThres = log2(2),PvalThres = 0.05,
            OutDir = "/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_s_vs_c_norm_volcano.pdf")
DrawVolcano(DEP_o_vs_s_norm,Log2FCThres = log2(2),PvalThres = 0.05,
            OutDir = "/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_o_vs_s_norm_volcano.pdf")
DrawVolcano(DEP_o_vs_c_H2O2,Log2FCThres = log2(2),PvalThres = 0.05,
            OutDir = "/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_o_vs_c_H2O2_volcano.pdf")
DrawVolcano(DEP_s_vs_c_H2O2,Log2FCThres = log2(2),PvalThres = 0.05,
            OutDir = "/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_s_vs_c_H2O2_volcano.pdf")
DrawVolcano(DEP_o_vs_s_H2O2,Log2FCThres = log2(2),PvalThres = 0.05,
            OutDir = "/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_o_vs_s_H2O2_volcano.pdf")

