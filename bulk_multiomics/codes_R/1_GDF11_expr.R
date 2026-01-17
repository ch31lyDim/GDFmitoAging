library(data.table)
prot = fread("/data1/dyh/GDF11/all_sample.xls")
prot = prot[prot$Gene=="GDF11",]
Design = rep(c(rep("Normal",3),rep("H2O2",3)),3)
Genotype = c(rep("Control",6),rep("sh",6),rep("OverExp",6))

prot = prot[,-c(1,2,3)]
prot_df = data.frame(log2_Intensity = log2(t(prot)[,1]+1),
                     Design = Design,
                     Genotype = Genotype)
prot_df$Genotype = factor(prot_df$Genotype,levels=c("Control","OverExp","Hairpin"))
ggline(prot_df,x="Genotype",y="log2_Intensity",color="Design",
          position = position_dodge(0.6),add = c("mean_se","point"),width = 0.6)+
    stat_compare_means(ref.group = 'Control',label="p.signif",size=10)+
    theme_classic(base_size = 20)+
    scale_color_nejm()
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/GDF11_expr.pdf",height = 7,width = 9)
