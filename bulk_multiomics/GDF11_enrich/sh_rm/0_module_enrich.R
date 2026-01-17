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
                       pvalueCutoff = 0.05)
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
