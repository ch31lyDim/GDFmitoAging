library(rtracklayer)
library(GSVA)
library(ggpubr)
library(ggsci)
library(dplyr)
expr = read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/data/prot",sep = '\t',header = T)
colnames(expr) <- c("CN1","CN2","CN3","CH1","CH2","CH3",
                    "shN1","shN2","shN3","shH1","shH2","shH3",
                    "OEN1","OEN2","OEN3","OEH1","OEH2","OEH3")
sample_info <- data.frame(sample = colnames(expr),
                          genotype = c(rep("C",6),rep("sh",6),rep("OE",6)),
                          condition = rep(c(rep("norm",3),rep("H2O2",3)),3),
                          treatment = substr(colnames(expr),1,nchar(colnames(expr))-1))
rownames(sample_info) <- sample_info$sample


expr = expr[,which(sample_info$genotype!='sh')]
sample_info <- sample_info[sample_info$genotype!='sh',]


Anno <- as.data.frame(import("/home/data5/dingyanheng/COP/GeneAnno/Canis_lupus_familiaris.ROS_Cfam_1.0.110.chr.gtf"))
Anno <- Anno[Anno$type=='gene',]
Anno <- Anno[,c("gene_id","gene_name")]
# Anno <- na.omit(Anno)
aging_markers <- c("ENSCAFG00845009606","CDKN2B","TP53","LMNB1","MMP3","IL6","TNF","KAT7","YAP1","ENSCAFG00845000688")

aging_markers <- Anno[Anno$gene_id%in%aging_markers | Anno$gene_name%in%aging_markers,]

DEG_P3P6 = read.table("/home/data5/dingyanheng/COP/P3P6/Transcriptome/results/DEG.txt")
P3P6_MARKER = read.table("/home/data5/dingyanheng/COP/P3P6/TxPxM/res/KeyGene1")
colnames(P3P6_MARKER) = "gene_name"
DEG_P3P6$gene_id <- rownames(DEG_P3P6)
DEG_P3P6 <- inner_join(DEG_P3P6,Anno,by="gene_id")
DEG_P3P6 <- DEG_P3P6[,c("gene_id","gene_name","log2FoldChange","padj")]

P3P6_MARKER <- inner_join(P3P6_MARKER,DEG_P3P6,by="gene_name")
P3P6_MARKER_UP <- P3P6_MARKER$gene_name[P3P6_MARKER$log2FoldChange>0]
P3P6_MARKER_DN <- P3P6_MARKER$gene_name[P3P6_MARKER$log2FoldChange<0]

gsea_res <- t(gsva(as.matrix(expr),list(AGING_MARKERS = aging_markers$gene_name,
                                         P3P6_MARKER_UP = P3P6_MARKER_UP,
                                         P3P6_MARKER_DN = P3P6_MARKER_DN),
                method="ssgsea",kcdf="Gaussian",abs.ranking=F,ssgsea.norm=F,tau = 0.5))

sample_info$condition[sample_info$condition=="norm"] <- "Young"
sample_info$condition[sample_info$condition=="H2O2"] <- "Senescence"

gsea_res <- cbind(gsea_res,sample_info)
gsea_res <- tidyr::pivot_longer(gsea_res,cols=AGING_MARKERS:P3P6_MARKER_DN,names_to="GeneSet",values_to = "ssGSEA_score")
gsea_res <- as.data.frame(gsea_res)
# gsea_res$ssGSEA_score <- scale(gsea_res$ssGSEA_score)
gsea_res$condition <- factor(gsea_res$condition,levels = c("Young","Senescence"))
gsea_res$genotype <- factor(gsea_res$genotype,levels = c("C","OE"))


gsea_res_AGING_MARKERS <- gsea_res[gsea_res$GeneSet=='AGING_MARKERS',]
gsea_res_P3P6_MARKER_UP <- gsea_res[gsea_res$GeneSet=='P3P6_MARKER_UP',]
gsea_res_P3P6_MARKER_DN <- gsea_res[gsea_res$GeneSet=='P3P6_MARKER_DN',]

p1<-ggline(gsea_res_AGING_MARKERS,x="genotype",y="ssGSEA_score",color = "condition",Group = "condition",
          add = c("jitter","mean_se"),size = 2,add.params = list(size = 3))+
    facet_grid(.~condition)+
    scale_color_manual(values = met.brewer("Cassatt1",2))+
    theme_classic(base_size = 26)+
    ggtitle("AGING_MARKERS")+
    theme(aspect.ratio = 2/1)+
    theme(legend.position = "none")+
    theme(axis.title.x = element_blank())+
    theme(axis.text.y = element_blank())+
    theme(axis.ticks.y = element_blank())
p2<-ggline(gsea_res_P3P6_MARKER_UP,x="genotype",y="ssGSEA_score",color = "condition",Group = "condition",
          add = c("jitter","mean_se"),size = 2,add.params = list(size = 3))+
    facet_grid(.~condition)+
    scale_color_manual(values = met.brewer("Cassatt1",2))+
    theme_classic(base_size = 26)+
    ggtitle("P3P6_MARKER_UP")+
    theme(aspect.ratio = 2/1)+
    theme(legend.position = "none")+
    theme(axis.title.x = element_blank())+
    theme(axis.text.y = element_blank())+
    theme(axis.ticks.y = element_blank())
    
p3<-ggline(gsea_res_P3P6_MARKER_DN,x="genotype",y="ssGSEA_score",color = "condition",Group = "condition",
          add = c("jitter","mean_se"),size = 2,add.params = list(size = 3))+
    facet_grid(.~condition)+
    scale_color_manual(values = met.brewer("Cassatt1",2))+
    theme_classic(base_size = 26)+
    ggtitle("P3P6_MARKER_DN")+
    theme(aspect.ratio = 2/1)+
    theme(legend.position = "none")+
    theme(axis.title.x = element_blank())+
    theme(axis.text.y = element_blank())+
    theme(axis.ticks.y = element_blank())
p <- ggarrange(p1,p2,p3,ncol = 3,nrow = 1)
p
ggsave("/home/data5/dingyanheng/COP/pub_branch/GDF11/sh_rm/ssPSEA.pdf",width = 14,height = 7)










