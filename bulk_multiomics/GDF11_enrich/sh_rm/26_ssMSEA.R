library(data.table)
aging_atlas = fread("/home/data5/dingyanheng/COP/P3P6/TxPxM/res/aging_atlas.txt",na.strings = '')
aging_atlas = aging_atlas[aging_atlas$V3=='hMSCs',]

c_in_aging_atlas_UP = na.omit(aging_atlas$V13[aging_atlas$V5<0])
c_in_aging_atlas_DN = na.omit(aging_atlas$V13[aging_atlas$V5>0])


meta_net <- fread("/home/data5/dingyanheng/COP/GDF11/Metabolome/data/metabo_neg")
meta_pos <- fread("/home/data5/dingyanheng/COP/GDF11/Metabolome/data/metabo_pos")
meta <- rbind(meta_net,meta_pos)
meta$Kegg_ID <- gsub("cpd:","",meta$Kegg_ID)
meta <- meta[,c(3:11,15:24)]

meta <- aggregate(.~Kegg_ID,sum,data = meta)
rownames(meta) <- meta$Kegg_ID
meta <- meta[,2:ncol(meta)]
colnames(meta) <- c("CN1","CN2","CN3","shN1","shN2","shN2","OEN1","OEN2","OEN3",
                    "CH1","CH2","CH3","shH1","shH2","shH3","OEH1","OEH2","OEH3")
expr <- meta

sample_info <- data.frame(sample = colnames(expr),
                          genotype = rep(c(rep("C",3),rep("sh",3),rep("OE",3)),2),
                          condition = rep(c(rep("norm",9),rep("H2O2",9)),1),
                          treatment = substr(colnames(expr),1,nchar(colnames(expr))-1))

expr = expr[,which(sample_info$genotype!='sh')]
sample_info <- sample_info[sample_info$genotype!='sh',]


gsea_res <- t(gsva(as.matrix(expr),list(AGING_ATLAS_UP = c_in_aging_atlas_UP,
                                        AGING_ATLAS_DN = c_in_aging_atlas_DN),
                method="ssgsea",kcdf="Gaussian",abs.ranking=F,ssgsea.norm=F,tau = 0.5))

sample_info$condition[sample_info$condition=="norm"] <- "Young"
sample_info$condition[sample_info$condition=="H2O2"] <- "Senescence"
gsea_res <- cbind(gsea_res,sample_info)
gsea_res <- tidyr::pivot_longer(gsea_res,cols=AGING_ATLAS_UP:AGING_ATLAS_DN,names_to="GeneSet",values_to = "ssGSEA_score")
gsea_res <- as.data.frame(gsea_res)
# gsea_res$ssGSEA_score <- scale(gsea_res$ssGSEA_score)
gsea_res$condition <- factor(gsea_res$condition,levels = c("Young","Senescence"))
gsea_res$genotype <- factor(gsea_res$genotype,levels = c("C","OE"))

gsea_res_AGING_ATLAS_UP <- gsea_res[gsea_res$GeneSet=='AGING_ATLAS_UP',]
gsea_res_AGING_ATLAS_DN <- gsea_res[gsea_res$GeneSet=='AGING_ATLAS_DN',]

p1<-ggline(gsea_res_AGING_ATLAS_UP,x="genotype",y="ssGSEA_score",color = "condition",Group = "condition",
          add = c("jitter","mean_se"),size = 2,add.params = list(size = 3))+
    facet_grid(.~condition)+
    scale_color_manual(values = met.brewer("Cassatt1",2))+
    theme_classic(base_size = 26)+
    ggtitle("AGING_ATLAS_UP")+
    theme(aspect.ratio = 2/1)+
    theme(legend.position = "none")+
    theme(axis.title.x = element_blank())+
    theme(axis.text.y = element_blank())+
    theme(axis.ticks.y = element_blank())
p2<-ggline(gsea_res_AGING_ATLAS_DN,x="genotype",y="ssGSEA_score",color = "condition",Group = "condition",
          add = c("jitter","mean_se"),size = 2,add.params = list(size = 3))+
    facet_grid(.~condition)+
    scale_color_manual(values = met.brewer("Cassatt1",2))+
    theme_classic(base_size = 26)+
    ggtitle("AGING_ATLAS_DN")+
    theme(aspect.ratio = 2/1)+
    theme(legend.position = "none")+
    theme(axis.title.x = element_blank())+
    theme(axis.text.y = element_blank())+
    theme(axis.ticks.y = element_blank())
    
p <- ggarrange(p1,p2,ncol = 2,nrow = 1)
p
ggsave("/home/data5/dingyanheng/COP/pub_branch/GDF11/sh_rm/ssMSEA.pdf",width = 14,height = 7)
