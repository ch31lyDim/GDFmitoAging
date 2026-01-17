library(MEGENA)
datExpr = read.table("/home/data5/dingyanheng/COP/GDF11/Transcriptome/data/TPM_name_mat",sep = '\t',header = T)
rownames(datExpr) = datExpr$Name
datExpr = datExpr[,-1]
datExpr_sd = apply(datExpr,MARGIN = 1,sd)
datExpr = datExpr[which(datExpr_sd!=0),]
datExpr = as.matrix(datExpr)
# input parameters
n.cores <- 62; # number of cores/threads to call for PCP
doPar <-TRUE; # do we want to parallelize?
method = "pearson" # method for correlation. either pearson or spearman. 
FDR.cutoff = 0.05 # FDR threshold to define significant correlations upon shuffling samples. 
module.pval = 0.05 # module significance p-value. Recommended is 0.05. 
hub.pval = 0.05 # connectivity significance p-value based random tetrahedral networks
cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs. 
hub.perm = 250; # number of permutations for calculating connectivity significance p-value. 

# annotation to be done on the downstream
annot.table=NULL
id.col = 1
symbol.col= 2
cor(t(datExpr[1:200,]))
run.par = doPar & (getDoParWorkers() == 1) 
if (run.par)
{
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  # check how many workers are there
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
}
ijw <- calculate.correlation(datExpr,FDR.cutoff = 0.05,doPerm = cor.perm,output.corTable = FALSE,output.permFDR = T,doPar = doPar,num.cores = n.cores)
#### register multiple cores if needed: note that set.parallel.backend() is deprecated. 

el <- calculate.PFN(ijw[,1:3],doPar = doPar,num.cores = n.cores,keep.track = FALSE)
g <- graph.data.frame(el,directed = FALSE)
##### perform MCA clustering.
MEGENA.output <- do.MEGENA(g,
 mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
 min.size = 10,max.size = vcount(g)/2,
 doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
 save.output = FALSE)

###### unregister cores as these are not needed anymore.
if (getDoParWorkers() > 1)
{
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
summary.output <- MEGENA.ModuleSummary(MEGENA.output,
    mod.pvalue = module.pval,hub.pvalue = hub.pval,
    min.size = 10,max.size = vcount(g)/2,
    annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
    output.sig = TRUE)

if (!is.null(annot.table))
{
  # update annotation to map to gene symbols
  V(g)$name <- paste(annot.table[[symbol.col]][match(V(g)$name,annot.table[[id.col]])],V(g)$name,sep = "|")
  summary.output <- output[c("mapped.modules","module.table")]
  names(summary.output)[1] <- "modules"
}
saveRDS(summary.output,"/home/data5/dingyanheng/COP/pub_branch/GDF11/module_local.rds")



library(clusterProfiler)
library(org.Cf.eg.db)
library(MetBrewer)
library(ggplot2)
library(ggpubr)
library(enrichplot)
mod = readRDS("/data1/dyh/GDF11/GDF11/module.rds")
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

GO_en_pairwise <- pairwise_termsim(GO_en)
emapplot(GO_en_pairwise,showCategory = 5,layout = "kk",circular = T,
          color = "-log10(p.adjust)",label_format = 50)

cluster_colors <- c(
  "c1_475" = "#E63946", # 亮红色(突出显示)
  "c1_476" = "#457B9D",
  "c1_477" = "#1D3557", 
  "c1_478" = "#A8DADC",
  "c1_479" = "#cbecc0",
  "c1_480" = "#2A9D8F"
)

p <- emapplot(GO_en_pairwise, showCategory = 5,
             layout = "fr",
             color = "-log10(p.adjust)",
             cex_category = 1.2,
             shadowtext = TRUE,
             label_format = 3)
# Apply additional styling
p_enhanced <- p + 
  scale_fill_manual(values = cluster_colors) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 20)),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )
ggsave("/data1/dyh/GDF11/MEGENA/module_enrich_GO_ema1.pdf",height = 8,width = 8)

en_475 <- enrichGO(gene = same_level$c1_475,
                       OrgDb = org.Cf.eg.db,
                       keyType = "SYMBOL",
                       ont = "ALL",
                       pAdjustMethod = "none",
                       pvalueCutoff = 0.1,
                       qvalueCutoff = 1)
en_475 <- pairwise_termsim(en_475)
res_475_1 <- as.data.frame(en_475)
res_475_1 <- arrange(res_475_1, FoldEnrichment)
#filter columns 'Description' with 'metabolic' in it
res_475 <- res_475_1[c(grepl("metabolic", res_475_1$Description, ignore.case = TRUE)),]
res_475 <- rbind(res_475_1[1:5,], res_475)
res_475 <- arrange(res_475, FoldEnrichment)
res_475 <- unique(res_475)
library(ggpubr)
res_475$labelx=rep(0,nrow(res_475))
res_475$labely=seq(nrow(res_475),1)
ggbarplot(res_475, y = "FoldEnrichment", x = 'Description', fill = "p.adjust") +
  # 添加文本到柱子内部
  geom_text(aes(label = Description, y = 5), 
            hjust = 0,  # 水平居中
            color = "black",
            size = 5,
            fontface = 'bold') +  # 调整文本大小
  theme_classic2(base_size = 18) +
  theme(aspect.ratio = 1/1,
        # 隐藏y轴上的描述文本和刻度线
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() +
  scale_fill_gradientn(colors = (met.brewer("Cassatt1")))+
  scale_y_continuous(expand = c(0,0))

ggsave("/data1/dyh/GDF11/plots/1F.pdf",height = 8,width = 8)
