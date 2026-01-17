library(mixOmics)
TPM <- read.table("/home/data5/dingyanheng/COP/GDF11/Transcriptome/data/TPM")
protMat <- read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/data/prot")
colnames(protMat) <- c("CN1","CN2","CN3","CH1","CH2","CH3",
                    "shN1","shN2","shN3","shH1","shH2","shH3",
                    "OEN1","OEN2","OEN3","OEH1","OEH2","OEH3")


metMat_neg <- t(read.table("/home/data5/dingyanheng/COP/GDF11/Metabolome/data/data_norm_shift_neg"))
metMat_pos <- t(read.table("/home/data5/dingyanheng/COP/GDF11/Metabolome/data/data_norm_shift_pos"))
metMat <- rbind(metMat_neg,metMat_pos)
metMat <- metMat[,!colnames(metMat)%in%c("QC1","QC2","QC3")]
colnames(metMat) <- c("CN1","CN2","CN3","shN1","shN2","shN3","OEN1","OEN2","OEN3",
                    "CH1","CH2","CH3","shH1","shH2","shH3","OEH1","OEH2","OEH3")


protMat <- protMat[,colnames(TPM)]
metMat <- metMat[,colnames(TPM)]

TPM<-TPM[,1:12]
protMat<-protMat[,1:12]
metMat<-metMat[,1:12]

normIndex <- 1:6
senIndex <- 7:12

TPM <- t(TPM)
TPM <- log2(TPM+1)
nzvTPM <- apply(TPM,2,var)
TPM <- TPM[,nzvTPM!=0]
list.keepX <- list(RNA = c(50,50),Prot = c(30, 30),Met = c(30,30))

X_norm <- list(RNA = TPM[normIndex,], Prot = t(protMat[,normIndex]), Met = t(metMat[,normIndex]))
Y_norm <- c(rep("Sen",3),rep("Norm",3))
MyResult.diablo <- block.splsda(X_norm, Y_norm,keepX = list.keepX, ncomp = 2,near.zero.var = T)
selected_res <- selectVar(MyResult.diablo)
normRNA <- selected_res$RNA$name
normProt <- selected_res$Prot$name
normMet <- selected_res$Met$name


X_sen <- list(RNA = TPM[senIndex,], Prot = t(protMat[,senIndex]), Met = t(metMat[,senIndex]))
Y_sen <- c(rep("Sen",3),rep("sen",3))
MyResult.diablo <- block.splsda(X_sen, Y_sen,keepX = list.keepX, ncomp = 2,near.zero.var = T)
selected_res <- selectVar(MyResult.diablo)
senRNA <- selected_res$RNA$name
senProt <- selected_res$Prot$name
senMet <- selected_res$Met$name


negAnno <- fread("/home/data5/dingyanheng/COP/GDF11/Metabolome/data/metabo_neg",sep = '\t')
posAnno <- fread("/home/data5/dingyanheng/COP/GDF11/Metabolome/data/metabo_pos",sep = '\t')

metAnno <- rbind(negAnno,posAnno)[,c("Compound_ID","Kegg_ID")]
metAnno$Kegg_ID <- gsub("cpd:","",metAnno$Kegg_ID)

normMet_kegg <- na.omit(metAnno$Kegg_ID[metAnno$Compound_ID%in%normMet])
senMet_kegg <- na.omit(metAnno$Kegg_ID[metAnno$Compound_ID%in%senMet])

KEGG_background = as.data.frame(data.table::fread("/home/data5/dingyanheng/COP/GeneAnno/dog_KEGGtoMet"))

TERM2GENE = KEGG_background[,c(1,3)]
TERM2NAME = unique(KEGG_background[,c(1,2)])



library(clusterProfiler)
library(org.Cf.eg.db)
RNAnormKEGGRes <- enrichKEGG(bitr(normRNA,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Cf.eg.db)[,2],organism = "cfa",
                            pvalueCutoff = 0.05,qvalueCutoff = 1,pAdjustMethod = 'none')
RNAsenKEGGRes <- enrichKEGG(bitr(senRNA,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Cf.eg.db)[,2],organism = "cfa",
                            pvalueCutoff = 0.05,qvalueCutoff = 1,pAdjustMethod = 'none')

ProtNormKEGGRes <- enrichKEGG(bitr(normProt,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Cf.eg.db)[,2],organism = "cfa",
                            pvalueCutoff = 0.05,qvalueCutoff = 1,pAdjustMethod = 'none')
ProtSenKEGGRes <- enrichKEGG(bitr(senProt,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Cf.eg.db)[,2],organism = "cfa",
                            pvalueCutoff = 0.05,qvalueCutoff = 1,pAdjustMethod = 'none')

MetNormKEGGRes <- enricher(normMet_kegg,TERM2GENE=TERM2GENE,TERM2NAME=TERM2NAME,
                    pvalueCutoff = 0.05,qvalueCutoff = 1,pAdjustMethod = 'none')
MetSenKEGGRes <- enricher(senMet_kegg,TERM2GENE=TERM2GENE,TERM2NAME=TERM2NAME,
                    pvalueCutoff = 0.05,qvalueCutoff = 1,pAdjustMethod = 'none')

RNAnormKEGGRes <- as.data.frame(RNAnormKEGGRes)[1:3,]
RNAsenKEGGRes <- as.data.frame(RNAsenKEGGRes)[1:3,]
ProtNormKEGGRes <- as.data.frame(ProtNormKEGGRes)[1:2,]
ProtSenKEGGRes <- as.data.frame(ProtSenKEGGRes)[1:3,]
MetNormKEGGRes <- as.data.frame(MetNormKEGGRes)[1:3,]
MetSenKEGGRes <- as.data.frame(MetSenKEGGRes)[1:3,]


RNAnormKEGGRes$Set <- "RNA_Control"
RNAsenKEGGRes$Set <- "RNA_OE"
ProtNormKEGGRes$Set <- "Protein_Control"
ProtSenKEGGRes$Set <- "Protein_OE"
MetNormKEGGRes$Set <- "Metabolite_Control"
MetSenKEGGRes$Set <- "Metabolite_OE"

RNAnormKEGGRes$MarkerType <- "RNA"
RNAsenKEGGRes$MarkerType <- "RNA"
ProtNormKEGGRes$MarkerType <- "Protein"
ProtSenKEGGRes$MarkerType <- "Protein"
MetNormKEGGRes$MarkerType <- "Metabolite"
MetSenKEGGRes$MarkerType <- "Metabolite"


RNAnormKEGGRes <- RNAnormKEGGRes[,c("Set","MarkerType","ID","Description","pvalue","Count")]
RNAsenKEGGRes <- RNAsenKEGGRes[,c("Set","MarkerType","ID","Description","pvalue","Count")]
ProtNormKEGGRes <- ProtNormKEGGRes[,c("Set","MarkerType","ID","Description","pvalue","Count")]
ProtSenKEGGRes <- ProtSenKEGGRes[,c("Set","MarkerType","ID","Description","pvalue","Count")]
MetNormKEGGRes <- MetNormKEGGRes[,c("Set","MarkerType","ID","Description","pvalue","Count")]
MetSenKEGGRes <- MetSenKEGGRes[,c("Set","MarkerType","ID","Description","pvalue","Count")]

KEGGRes <- rbind(RNAnormKEGGRes,RNAsenKEGGRes,ProtNormKEGGRes,ProtSenKEGGRes,MetNormKEGGRes,MetSenKEGGRes)

ggplot(KEGGRes,aes(x = Set,y = Description,color = -log10(pvalue),size = Count))+
    geom_point()+
    theme_bw(base_size = 12)+
    theme(aspect.ratio = 3/2)+
    theme(axis.title = element_blank())+
    MetBrewer::scale_color_met_c("Cassatt1",direction = -1)+
    theme(axis.text.x  = element_text(angle = 45, hjust = 1))+
    theme(legend.title = element_text(face = 'bold'))
ggsave("/home/data5/dingyanheng/COP/pub_branch/GDF11/sh_rm/DIABLO_KEGG.pdf",width = 6,height = 6)
