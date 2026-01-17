library(MEGENA)
library(igraph)
library(ggplot2)
library(MetBrewer)
datExpr = read.table("/home/data3/dingyanheng/COP/data/TPM_name_mat",sep = '\t',header = T)
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
saveRDS(summary.output,"/home/data3/dingyanheng/COP/MEGENA/module.rds")

col = met.brewer('Monet',9)
GDF11_mode = lapply(summary.output[[1]],function(x){return("GDF11"%in%x)})
GDF11_mode = names(GDF11_mode)[which(GDF11_mode==TRUE)]

GDF11_mode

pnet.obj <- plot_module(output.summary = summary.output,PFN = g,subset.module = GDF11_mode[1],
    layout = "kamada.kawai",label.hubs.only = TRUE,
    # gene.set = list("GDF11"="GDF11"),color.code =  "#f41010",
    output.plot = FALSE,out.dir = "modulePlot",label.scaleFactor = 20,
    col.names = c(rep(NA,7),col[5]),
    hubLabel.col = "red",hubLabel.sizeProp = 2,show.topn.hubs = 26,show.legend = F,
    label.alpha = 1)

pnet.obj[[1]]+
  theme(aspect.ratio = 1/1)+
  theme_void()
ggsave("/home/data3/dingyanheng/COP/MEGENA/module/mod1.pdf",height = 15,width = 15)



pnet.obj <- plot_module(output.summary = summary.output,PFN = g,subset.module = GDF11_mode[2],
    layout = "kamada.kawai",label.hubs.only = F,
    gene.set = list("GDF11"="GDF11"),color.code =  "#f41010",
    output.plot = FALSE,out.dir = "modulePlot",label.scaleFactor = 20,
    # col.names = c(rep(NA,7),col[5]),
    hubLabel.col = "red",hubLabel.sizeProp = 2,show.topn.hubs = 26,show.legend = F,
    label.alpha = 1)
pnet.obj[[1]]+
  theme(aspect.ratio = 1/1)+
  theme_void()
ggsave("/home/data3/dingyanheng/COP/MEGENA/module/mod2.pdf",height = 15,width = 15)


pnet.obj <- plot_module(output.summary = summary.output,PFN = g,subset.module = GDF11_mode[3],
    layout = "kamada.kawai",label.hubs.only = F,
    gene.set = list("GDF11"="GDF11"),color.code =  "#f41010",
    output.plot = FALSE,out.dir = "modulePlot",label.scaleFactor = 20,
    # col.names = c(rep(NA,7),col[5]),
    hubLabel.col = "red",hubLabel.sizeProp = 2,show.topn.hubs = 26,show.legend = F,
    label.alpha = 1)
pnet.obj[[1]]+
  theme(aspect.ratio = 1/1)+
  theme_void()
ggsave("/home/data3/dingyanheng/COP/MEGENA/module/mod3.pdf",height = 15,width = 15)
