library(ggvenn)
DEP = read.table("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP.txt",header = T)
rownames(DEP) = DEP$name
# DEP$C_H_vs_C_NC_p.val = p.adjust(DEP$C_H_vs_C_NC_p.val,'fdr')
# DEP$OE_H_vs_OE_NC_p.val = p.adjust(DEP$OE_H_vs_OE_NC_p.val,'fdr')
# DEP$SH_H_vs_SH_NC_p.val = p.adjust(DEP$SH_H_vs_SH_NC_p.val,'fdr')

Control_DEP = DEP$name[abs(DEP$C_H_vs_C_NC_ratio)>log2(2) & DEP$C_H_vs_C_NC_p.val<0.05]
OverExpr_DEP = DEP$name[abs(DEP$OE_H_vs_OE_NC_ratio)>log2(2) & DEP$OE_H_vs_OE_NC_p.val<0.05]
sh_DEP = DEP$name[abs(DEP$SH_H_vs_SH_NC_ratio)>log2(2) & DEP$SH_H_vs_SH_NC_p.val<0.05]

x = list(Control=Control_DEP,
         OverExpr=OverExpr_DEP,
         sh=sh_DEP)
ggvenn(x,set_name_size = 8,text_size = 6,stroke_linetype = 0,fill_alpha = 0.2)+
scale_fill_nejm()
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_venn.pdf")

Control_UP = DEP$name[DEP$C_H_vs_C_NC_ratio>log2(2) & DEP$C_H_vs_C_NC_p.val<0.05]
OverExpr_UP = DEP$name[DEP$OE_H_vs_OE_NC_ratio>log2(2) & DEP$OE_H_vs_OE_NC_p.val<0.05]
sh_UP = DEP$name[DEP$SH_H_vs_SH_NC_ratio>log2(2) & DEP$SH_H_vs_SH_NC_p.val<0.05]

x = list(Control_UP=Control_UP,
         OverExpr_UP=OverExpr_UP,
         sh_UP=sh_UP)
ggvenn(x,set_name_size = 8,text_size = 6,stroke_linetype = 0,fill_alpha = 0.2)+
scale_fill_nejm()
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_venn_UP.pdf")

Control_DN = DEP$name[DEP$C_H_vs_C_NC_ratio< -log2(2) & DEP$C_H_vs_C_NC_p.val<0.05]
OverExpr_DN = DEP$name[DEP$OE_H_vs_OE_NC_ratio< -log2(2) & DEP$OE_H_vs_OE_NC_p.val<0.05]
sh_DN = DEP$name[DEP$SH_H_vs_SH_NC_ratio< -log2(2) & DEP$SH_H_vs_SH_NC_p.val<0.05]

x = list(Control_DN=Control_DN,
         OverExpr_DN=OverExpr_DN,
         sh_DN=sh_DN)
ggvenn(x,set_name_size = 8,text_size = 6,stroke_linetype = 0,fill_alpha = 0.2)+
scale_fill_nejm()
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/DEP_venn_DN.pdf")


O_vs_C_DEP = intersect(OverExpr_DEP,Control_DEP)
O_vs_C_DEP_df = na.omit(data.frame(name = O_vs_C_DEP,
                           log2FC_O = DEP[O_vs_C_DEP,"OE_H_vs_OE_NC_ratio"],
                           log2FC_C = DEP[O_vs_C_DEP,"C_H_vs_C_NC_ratio"]))
O_vs_C_DEP_df$sig = ifelse(O_vs_C_DEP_df$log2FC_O*O_vs_C_DEP_df$log2FC_C>0,"Same","Diff")
O_vs_C_DEP_diff = O_vs_C_DEP_df[O_vs_C_DEP_df$sig=="Diff",'name']
O_vs_C_DEP_same = O_vs_C_DEP_df[O_vs_C_DEP_df$sig=="Same",'name']
write.table(O_vs_C_DEP_diff,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/O_vs_C_DEP_diff",
            col.names = F,row.names = F,quote=F,sep = '\t')
write.table(O_vs_C_DEP_same,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/O_vs_C_DEP_same",
            col.names = F,row.names = F,quote=F,sep = '\t')            

ggscatter(O_vs_C_DEP_df,x="log2FC_O",y="log2FC_C",add = 'reg.line',color = "sig")+
        geom_vline(xintercept = 0,col = "gray", linetype = 'dashed')+
        geom_hline(yintercept = 0,col = "gray", linetype = 'dashed')+
        scale_color_manual(values = c("Same"="#021f68","Diff"="#7d0d0d"))+
        geom_text_repel(mapping = aes(label=name,color = sig))+
        stat_cor(method = 'pearson', aes(x = log2FC_O, y = log2FC_C,color = sig),size = 4)
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/O_vs_C_DEP_scatter.pdf",width = 6,height = 6)


S_vs_C_DEP = intersect(sh_DEP,Control_DEP)
S_vs_C_DEP_df = na.omit(data.frame(name = S_vs_C_DEP,
                           log2FC_S = DEP[S_vs_C_DEP,"SH_H_vs_SH_NC_ratio"],
                           log2FC_C = DEP[S_vs_C_DEP,"C_H_vs_C_NC_ratio"]))
S_vs_C_DEP_df$sig = ifelse(S_vs_C_DEP_df$log2FC_S*S_vs_C_DEP_df$log2FC_C>0,"Same","Diff")
S_vs_C_DEP_diff = S_vs_C_DEP_df[S_vs_C_DEP_df$sig=="Diff",'name']
S_vs_C_DEP_same = S_vs_C_DEP_df[S_vs_C_DEP_df$sig=="Same",'name']
write.table(S_vs_C_DEP_diff,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/S_vs_C_DEP_diff",
            col.names = F,row.names = F,quote=F,sep = '\t')
write.table(S_vs_C_DEP_same,"/home/data5/dingyanheng/COP/GDF11/Proteome/res/ProtSets/S_vs_C_DEP_same",
            col.names = F,row.names = F,quote=F,sep = '\t')
ggscatter(S_vs_C_DEP_df,x="log2FC_S",y="log2FC_C",add = 'reg.line',color = "sig")+
        geom_vline(xintercept = 0,col = "gray", linetype = 'dashed')+
        geom_hline(yintercept = 0,col = "gray", linetype = 'dashed')+
        scale_color_manual(values = c("Same"="#021f68","Diff"="#7d0d0d"))+
        geom_text_repel(mapping = aes(label=name,color = sig))+
        stat_cor(method = 'pearson', aes(x = log2FC_S, y = log2FC_C,color = sig),size = 4)
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/S_vs_C_DEP_scatter.pdf",width = 6,height = 6)


O_vs_S_DEP = intersect(OverExpr_DEP,sh_DEP)
O_vs_S_DEP_df = na.omit(data.frame(name = O_vs_S_DEP,
                           log2FC_O = DEP[O_vs_S_DEP,"OE_H_vs_OE_NC_ratio"],
                           log2FC_S = DEP[O_vs_S_DEP,"SH_H_vs_SH_NC_ratio"]))
O_vs_S_DEP_df$sig = ifelse(O_vs_S_DEP_df$log2FC_O*O_vs_S_DEP_df$log2FC_S>0,"Same","Diff")
ggscatter(O_vs_S_DEP_df,x="log2FC_O",y="log2FC_S",add = 'reg.line',color = "sig")+
        geom_vline(xintercept = 0,col = "gray", linetype = 'dashed')+
        geom_hline(yintercept = 0,col = "gray", linetype = 'dashed')+
        scale_color_manual(values = c("Same"="#021f68","Diff"="#7d0d0d"))+
        geom_text_repel(mapping = aes(label=name,color = sig))+
        stat_cor(method = 'pearson', aes(x = log2FC_O, y = log2FC_S,color = sig),size = 4)
ggsave("/home/data5/dingyanheng/COP/GDF11/Proteome/res/O_vs_S_DEP_scatter.pdf",width = 6,height = 6)




