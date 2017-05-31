## ----available.data.types, eval=FALSE------------------------------------
#      available.data.types()

## ----process.multiple.studies 1, eval=FALSE------------------------------
#  process.multiple.studies(genes, cancers, high.throughput.data.type, data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value", "Median"), shorteded.cancer.names = TRUE, genelimit="none", resolution=600, RowCex=0.8, ColCex=0.8, heatmapMargines=c(24,17), cutoff=2, angle.for.heatmap.cancernames=45, heatmap.color = "RdBu", reverse.heatmap.color = TRUE, rewrite.output.list = TRUE, round=TRUE, top.genes = TRUE,  validate.genes = TRUE, Use.CancerCode.as.Name = FALSE, simplify.visulization=FALSE, simplifiction.cuttoff=FALSE)

## ----process.multiple.studies 2, eval=FALSE------------------------------
#  genes <- list(K.demethylases = c("KDM1A", "KDM1B", "KDM2A"), K.acetyltransferases = c("CLOCK", "CREBBP", "ELP3", "EP300"))
#  
#  cancernames <- c("Acute Myeloid Leukemia (TCGA, Provisional)", "Adrenocortical Carcinoma (TCGA, Provisional)", "Bladder Urothelial Carcinoma (TCGA, Provisional)", "Brain Lower Grade Glioma (TCGA, Provisional)", "Breast Invasive Carcinoma (TCGA, Provisional)")
#  
#  process.multiple.studies(genes, cancernames, "RNA-seq")

## ----process.multiple.studies 3, eval=FALSE------------------------------
#  cancernames <- matrix(c("Acute Myeloid Leukemia (TCGA, Provisional)", "acute myeloid leukemia", "Adrenocortical Carcinoma (TCGA, Provisional)", "adrenocortical carcinoma", "Bladder Urothelial Carcinoma (TCGA, Provisional)", "bladder urothelial carcinoma", "Brain Lower Grade Glioma (TCGA, Provisional)", "brain lower grade glioma", "Breast Invasive Carcinoma (TCGA, Provisional)",  "breast invasive carcinoma"), nrow = 5, ncol=2 , byrow = TRUE)

## ----process.multiple.studies 4, eval=FALSE------------------------------
#  process.multiple.studies(genes, cancernames, "RNA-seq", data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value"), shorteded.cancer.names = FALSE, resolution=300, RowCex=1, ColCex=1, heatmapMargines=c(15,5), cutoff=1.5, angle.for.heatmap.cancernames=30, heatmap.color = "redgreen")

## ----process.multiple.studies 5, eval=FALSE------------------------------
#  process.multiple.studies(genes, cancernames, "RNA-seq", heatmap.color = "redgreen", rewrite.output.list = FALSE)

## ----process.one.study 1, eval=FALSE-------------------------------------
#  process.one.study(genes, cancername, high.throughput.data.type, data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value", "Median"), transposedHeatmap=FALSE, desired.case.list="None", genelimit="none", resolution=600, RowCex=0.8, ColCex=0.8, heatmapMargines=c(10,10), cutoff="default", angle.for.heatmap.cancernames=45, heatmap.color = "RdBu", reverse.heatmap.color = TRUE, rewrite.output.list = TRUE, round=TRUE, top.genes = TRUE, validate.genes = TRUE, simplify.visulization=FALSE, simplifiction.cuttoff=FALSE)

## ----process.one.study 2, eval=FALSE-------------------------------------
#  genes <- list(K.demethylases = c("KDM1A", "KDM1B", "KDM2A"), K.acetyltransferases = c("CLOCK", "CREBBP", "ELP3", "EP300"))
#  
#  process.one.study(genes, "Breast Invasive Carcinoma (TCGA, Cell 2015)", "RNA-seq")

## ----process.one.study 3, eval=FALSE-------------------------------------
#  process.one.study(genes, cancernames, "RNA-seq", desired.case.list = c(3,4,5), data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value"), resolution=300, RowCex=1, ColCex=1, heatmapMargines=c(15,5), cutoff=1.5, angle.for.heatmap.cancernames=30, heatmap.color = "redgreen")

## ----process.one.study 4, eval=FALSE-------------------------------------
#  Enter the numeric index(es):                  2,3,4,5

## ----process.one.study 5, eval=FALSE-------------------------------------
#  process.one.study(genes, "Breast Invasive Carcinoma (TCGA, Cell 2015)", "RNA-seq", heatmap.color = "redgreen", rewrite.output.list = FALSE)

