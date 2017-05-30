## ----available.data.types, eval=FALSE------------------------------------
#      available.data.types()

## ----process.multiple.studies 1, eval=FALSE------------------------------
#  genes <- list(K.demethylases = c("KDM1A", "KDM1B", "KDM2A"), K.acetyltransferases = c("CLOCK", "CREBBP", "ELP3", "EP300"))
#  
#  cancernames <- c("Acute Myeloid Leukemia (TCGA, Provisional)", "Adrenocortical Carcinoma (TCGA, Provisional)", "Bladder Urothelial Carcinoma (TCGA, Provisional)", "Brain Lower Grade Glioma (TCGA, Provisional)", "Breast Invasive Carcinoma (TCGA, Provisional)")
#  
#  process.multiple.studies(genes, cancernames, "RNA-seq")

## ----process.multiple.studies 2, eval=FALSE------------------------------
#  cancernames <- matrix(c("Acute Myeloid Leukemia (TCGA, Provisional)", "acute myeloid leukemia", "Adrenocortical Carcinoma (TCGA, Provisional)", "adrenocortical carcinoma", "Bladder Urothelial Carcinoma (TCGA, Provisional)", "bladder urothelial carcinoma", "Brain Lower Grade Glioma (TCGA, Provisional)", "brain lower grade glioma", "Breast Invasive Carcinoma (TCGA, Provisional)",  "breast invasive carcinoma"), nrow = 5, ncol=2 , byrow = TRUE)

## ----process.multiple.studies 4, eval=FALSE------------------------------
#  process.multiple.studies(genes, cancernames, "RNA-seq", data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value"), shorteded.cancer.names = FALSE, resolution=300, RowCex=1, ColCex=1, heatmapMargines=c(15,5), cutoff=1.5, angle.for.heatmap.cancernames=30, heatmap.color = "redgreen")

## ----process.multiple.studies 3, eval=FALSE------------------------------
#  process.multiple.studies(genes, cancernames, "RNA-seq", heatmap.color = "redgreen", resetOldExpressionProfile = FALSE)

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

