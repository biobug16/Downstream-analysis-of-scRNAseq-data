suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(cowplot)
  library(ggplot2)
  library(pheatmap)
  library(rafalib)
  # remotes::install_github("powellgenomicslab/scPred")
  library(scPred)
})

set.seed= 1623

test=readRDS("/data/Lab/mehanib2/primary_recurrent_paired/all_merged_seurat_with_subclass.RDS")


test2=subset(test, subset = Anno != 'Unknown')

pair=readRDS("/data/Lab/mehanib2/scRNA_Jan2023/GSE174554/Analysis/all_paired_tumor_merged_after_droplet_removal_merging_normalize_scale_pca_n_harmony_0.05res.RDS")

pair2=subset(pair, subset = IDH != 'IDH mutant')



transfer.anchors <- FindTransferAnchors(reference = test2, query = pair2,dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = test2$Anno,dims = 1:30)
pair2 <- AddMetaData(object = pair2, metadata = predictions)

pdf("paired_snRNASeq_glioma_with_label_transfer.pdf",width=8)
DimPlot(pair2, group.by = "predicted.id", label = T, repel = F,cols = c("purple", "yellow", "orange","red","blue","green","gray"),raster=FALSE)
dev.off()

saveRDS(pair2,"glioma_snRNA_pair_with_label_transfer_metadata.RDS")

pair2=readRDS("glioma_snRNA_pair_with_label_transfer_metadata.RDS")
pair2@meta.data=read.csv("glioma_snRNA_pair_with_label_transfer_metadata.csv",header=T,row.names=1)
myeloid=subset(pair2, subset = predicted.id == 'Myeloid')
malignant=subset(pair2, subset = predicted.id == 'Malignant')


sign=read.csv("comm_genes_wt_myeloid_sign_for_gsva.csv",header=T)
myeloid1=AddModuleScore(myeloid,sign)


ES.seurat <- enrichIt(obj = myeloid, 
                      gene.sets = sign, 
                      groups = 1000, cores = 12, 
                      min.size = 1,method = "UCell")

ES.seurat$mye_Diff=ES.seurat$sig_high-ES.seurat$sig_low
a=cbind(myeloid@meta.data, ES.seurat)
pdf("wt_myeloid_pair_type.pdf",width=20)
ggplot(a, aes(x=as.factor(Number), y= mye_Diff, fill=Type)) + geom_boxplot()+theme(text=element_text(size=21))+ stat_compare_means(aes(group = Type), label = "p.signif")

ES.seurat2 <- enrichIt(obj = malignant, 
                      gene.sets = sign2, 
                      groups = 1000, cores = 12, 
                      min.size = 1,method = "UCell")

ES.seurat2$malig_Diff=ES.seurat2$sig_high-ES.seurat2$sig_low
b=cbind(malignant@meta.data, ES.seurat2)



myeloid1=AddModuleScore(myeloid,sign)
colnames(myeloid1@meta.data)[26]="sig_high"
colnames(myeloid1@meta.data)[27]="sig_low"
myeloid1@meta.data$mye_Diff=myeloid1@meta.data[,26]-myeloid1@meta.data[,27]
a=myeloid1@meta.data
pdf("wt_myeloid_pair_type_addmod.pdf",width=20)
ggplot(a, aes(x=as.factor(Number), y= mye_Diff, fill=Type)) + geom_boxplot()+theme(text=element_text(size=21))+ stat_compare_means(aes(group = Type), label = "p.signif")
dev.off()

pdf("wt_myeloid_pair_tumor_type1.pdf",width=10)
ggplot(a, aes(x=Type, y= mye_Diff, fill=Type)) + geom_boxplot2(width = 0.8, width.errorbar = 0.5)+ stat_compare_means(aes(group = Type), label = "p.signif",label.x = 1.5, label.y = 0.65,size=15)+theme(text=element_text(size=25))
dev.off()

pdf("wt_myeloid_pair_type_poor.pdf",width=10)
ggplot(a, aes(x=Type, y= sig_high, fill=Type)) + geom_boxplot2(width = 0.8, width.errorbar = 0.5)+ stat_compare_means(aes(group = Type), label = "p.signif",label.x = 1.5, label.y = 0.15,size=15)+ ylab("mye_poor_surv")+theme(text=element_text(size=21))
dev.off()



pdf("wt_myeloid_pair_type_better.pdf",width=10)
ggplot(a, aes(x=Type, y= sig_low, fill=Type)) + geom_boxplot2(width = 0.8, width.errorbar = 0.5)+ stat_compare_means(aes(group = Type), label = "p.signif",label.x = 1.5, label.y = 0.4,size=15)+ ylab("mye_better_surv")+theme(text=element_text(size=21))
dev.off()


pdf("wt_myeloid_pair_type_poor.pdf",width=20)
ggplot(a, aes(x=as.factor(Number), y= sig_high, fill=Type)) + geom_boxplot()+theme(text=element_text(size=21))+ stat_compare_means(aes(group = Type), label = "p.signif")+ ylab("mye_poor_surv")
dev.off()




malignant1=AddModuleScore(malignant,sign2)
colnames(malignant1@meta.data)[26]="sig_high"
colnames(malignant1@meta.data)[27]="sig_low"
malignant1@meta.data$malig_Diff= malignant1@meta.data[,26]-malignant1@meta.data[,27]
b=malignant1@meta.data
pdf("wt_malignant_pair_type.pdf",width=20)
ggplot(b, aes(x=as.factor(Number), y= malig_Diff, fill=Type)) + geom_boxplot()+theme(text=element_text(size=21))+ stat_compare_means(aes(group = Type), label = "p.signif")
dev.off()

pdf("wt_malignant_pair_type_poor.pdf",width=10)
ggplot(b, aes(x=Type, y= sig_high, fill=Type)) + geom_boxplot2(width = 0.8, width.errorbar = 0.5)+ stat_compare_means(aes(group = Type), label = "p.signif",label.x = 1.5, label.y = 0.15,size=15)+ ylab("malig_poor_surv")+theme(text=element_text(size=21))
dev.off()


pdf("wt_malignant_pair_type_better.pdf",width=10)
ggplot(b, aes(x=Type, y= sig_low, fill=Type)) + geom_boxplot2(width = 0.8, width.errorbar = 0.5)+ stat_compare_means(aes(group = Type), label = "p.signif",label.x = 1.5, label.y = 0.15,size=15)+ ylab("malig_better_surv")+theme(text=element_text(size=21))
dev.off()


pdf("wt_malignant_pair_tumor_type1.pdf",width=10)
ggplot(b, aes(x=Type, y= malig_Diff, fill=Type)) + geom_boxplot2(width = 0.8, width.errorbar = 0.5)+ stat_compare_means(aes(group = Type), label = "p.signif",label.x = 1.5, label.y = 0.21,size=15)+theme(text=element_text(size=21))
dev.off()


