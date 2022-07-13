#' @title Check Expression/methylation Profile for various cancer studies.
#'
#' @description This function Obtains the requested data for the given genes
#' across multiple cancer studie. It can check whether or not all genes are
#' included in cancer studies and and, if not, looks for the alternative gene
#' names. Then it calculates frequency percentage, frequency ratio, mean value
#' and median value of samples greather than specific value in the selected
#' cancer studies. Furthermore, it looks for the five genes that comprise the
#' highest values in each cancer study.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cbaf \cr
#' Type: \tab Package \cr
#' Version: \tab 1.19.5 \cr
#' Date: \tab 2022-07-14 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#'
#' @include cbaf-obtainMultipleStudies.R cbaf-automatedStatistics.R
#' cbaf-heatmapOutput.R cbaf-xlsxOutput.R
#'
#' @usage processMultipleStudies(genesList, submissionName, studiesNames,
#'   desiredTechnique, cancerCode = FALSE, validateGenes = TRUE, calculate =
#'   c("frequencyPercentage", "frequencyRatio", "meanValue"), cutoff=NULL,
#'   round=TRUE, topGenes = TRUE, shortenStudyNames = TRUE, geneLimit = 50,
#'   rankingMethod = "variation", heatmapFileFormat = "TIFF", resolution = 600,
#'   RowCex = "auto", ColCex = "auto", heatmapMargines = "auto",
#'   rowLabelsAngle = 0, columnLabelsAngle = 45, heatmapColor = "RdBu",
#'   reverseColor = TRUE, transposedHeatmap = FALSE, simplifyBy = FALSE,
#'   genesToDrop = FALSE, transposeResults = FALSE)
#'
#'
#'
#' @param genesList a list that contains at least one gene group
#'
#' @param submissionName a character string containing name of interest. It is
#' used for naming the process.
#'
#' @param studiesNames a character vector or a matrix that containes desired
#' cancer names. The character vector containes standard names of cancer studies
#'  that can be found on cbioportal.org, such as
#'  \code{"Acute Myeloid Leukemia (TCGA, NEJM 2013)"}. Alternatively, a matrix
#'  can be used if users prefer user-defined cancer names. In this case, the
#'  first column of matrix comprises the standard cancer names while the second
#'  column must contain the desired cancer names.
#'
#' @param desiredTechnique a character string that is one of the following
#' techniques: \code{"RNA-Seq"}, \code{"RNA-SeqRTN"}, \code{"microRNA-Seq"},
#' \code{"microarray.mRNA"}, \code{"microarray.microRNA"} or
#' \code{"methylation"}.
#'
#' @param cancerCode a logical value that tells the function to use cbioportal
#' abbreviated cancer names instead of complete cancer names, if set to be
#' \code{"TRUE"}. For example, \code{"laml_tcga_pub"} is the abbreviated name
#' for \code{"Acute Myeloid Leukemia (TCGA, NEJM 2013)"}.
#'
#' @param validateGenes a logical value that, if set to be \code{TRUE}, causes
#'  the function to check each cancer study to find whether or not each gene has
#'  a record. If a cancer doesn't have a record for specific gene, function
#'  looks for alternative gene names that cbioportal might use instead of the
#'  given gene name.
#'
#' @param calculate a character vector that containes the statistical procedures
#' users prefer the function to compute. The complete results can be obtained
#' by \code{c("frequencyPercentage", "frequencyRatio", "meanValue",
#' "medianValue")}. This will tell the function to compute the following:
#' \code{"frequencyPercentage"}, which is the percentge of samples having the
#' value greather than specific cutoff divided by the total sample size for
#' every study / study subgroup;
#' \code{"frequency ratio"}, which shows the number of selected samples divided
#' by the total number of samples that give the frequency percentage for every
#' study / study subgroup. It shows the selected and total sample sizes.;
#' \code{"Mean Value"}, that contains mean value of selected samples for each
#' study;
#' \code{"Median Value"}, which shows the median value of selected samples for
#' each study.
#' The default input is \code{calculate = c("frequencyPercentage",
#' "frequencyRatio", "meanValue")}.
#'
#' @param cutoff a number used to limit samples to those that are greather than
#' this number (cutoff). The default value for methylation data is \code{0.8}
#' while gene expression studies use default value of \code{2}. For methylation
#' studies, it is \code{average of relevant locations}, for the rest, it is
#' \code{"log z-score"}. To change the cutoff to any desired number, change the
#' option to \code{cutoff = desiredNumber} in which desiredNumber is the number
#' of interest.
#'
#' @param round a logical value that, if set to be \code{TRUE}, will force the
#' function to round all the calculated values to two decimal places. The
#' default value is \code{TRUE}.
#'
#' @param topGenes a logical value that, if set as \code{TRUE}, causes the
#' function to create three dataframes that contain the five top genes for each
#' cancer. To get all the three data.frames, \code{"frequencyPercentage"},
#' \code{"meanValue"} and \code{"medianValue"} must have been included for
#' \code{calculate}.
#'
#' @param shortenStudyNames a logical vector. If the value is set as TRUE,
#' function will try to remove the last part of the cancer names aiming to
#' shorten them. The removed segment usually contains the name of scientific
#' group that has conducted the experiment.
#'
#' @param geneLimit if large number of genes exist in at least one gene group,
#' this option can be used to limit the number of genes that are shown on
#' heatmap. For instance, \code{geneLimit=50} will limit the heatmap to 50 genes
#'  showing the most variation across multiple study / study subgroups. The
#'  default value is \code{50}.
#'
#' @param rankingMethod a character value that determines how genes will be
#' ranked prior to drawing heatmap. \code{"variation"} orders the genes based on
#' unique values in one or few cancer studies while \code{"highValue"} ranks the
#'  genes when they contain high values in multiple / many cancer studies. This
#'  option is useful when number of genes are too much so that user has to limit
#'  the number of genes on heatmap by \code{geneLimit}.
#'
#' @param heatmapFileFormat This option enables the user to select the desired
#' image file format of the heatmaps. The default value is \code{"TIFF"}. Other
#' supported formats include \code{"JPG"}, \code{"BMP"}, \code{"PNG"}, and
#' \code{"PDF"}.
#'
#' @param resolution a number. This option can be used to adjust the resolution
#' of the output heatmaps as 'dot per inch'. The defalut value is 600.
#'
#' @param RowCex a number that specifies letter size in heatmap row names,
#' which ranges from 0 to 2. If \code{RowCex = "auto"}, the function will
#' automatically determine the best RowCex.
#'
#' @param ColCex a number that specifies letter size in heatmap column names,
#' which ranges from 0 to 2. If \code{ColCex = "auto"}, the function will
#' automatically determine the best ColCex.
#'
#' @param heatmapMargines a numeric vector that is used to set heatmap margins.
#'  If \code{heatmapMargines = "auto"}, the function will automatically
#'  determine the best possible margines. Otherwise, enter the desired margine as
#'  e.g. c(10,10.)
#'
#' @param rowLabelsAngle a number that determines the angle with which the
#' gene names are shown in heatmaps. The default value is 0 degree.
#'
#' @param columnLabelsAngle a number that determines the angle with which the
#' studies/study subgroups names are shown in heatmaps. The default value is 45
#' degree.
#'
#' @param heatmapColor a character string that defines heatmap color. The
#' default value is \code{'RdBu'}. \code{'RdGr'} is also a popular color in
#' genomic studies. To see the rest of colors, please type
#' \code{library(RColorBrewer)} and then \code{display.brewer.all()}.
#'
#' @param reverseColor a logical value that reverses the color gradiant for
#' heatmap(s).
#'
#' @param transposedHeatmap a logical value that transposes heatmap rows to
#' columns and vice versa.
#'
#' @param simplifyBy a number that tells the function to change the values
#' smaller than that to zero. The purpose behind this option is to facilitate
#' recognizing candidate genes. Therefore, it is not suited for publications. It
#' has the same unit as \code{cutoff}.
#'
#' @param genesToDrop a character vector. Gene names within this vector will be
#' omitted from heatmap.The default value is \code{FALSE}.
#'
#' @param transposeResults a logical value that enables the function to replace
#' the columns and rows of data.
#'
#'
#'
#' @return a BiocFileCache object that containes some or all of the following
#' groups, based on what user has chosen: \code{obtainedData},
#' \code{validationResults}, \code{frequencyPercentage},
#' \code{Top.Genes.of.Frequency.Percentage}, \code{frequencyRatio},
#' \code{meanValue}, \code{Top.Genes.of.Mean.Value}, \code{medianValue},
#' \code{Top.Genes.of.Median.Value}. It also saves these results in one excel
#' files for convenience. Based on preference, three heatmaps for frequency
#' percentage, mean value and median can be generated. If more than one group of
#'  genes is entered, output for each group will be strored in a separate
#'  sub-directory.
#'
#' @examples
#' genes <- list(K.demethylases = c("KDM1A", "KDM1B", "KDM2A", "KDM2B", "KDM3A",
#'  "KDM3B", "JMJD1C", "KDM4A"), K.methyltransferases = c("SUV39H1", "SUV39H2",
#'  "EHMT1", "EHMT2", "SETDB1", "SETDB2", "KMT2A", "KMT2A"))
#'
#' studies <- c("Acute Myeloid Leukemia (TCGA, Provisional)",
#' "Adrenocortical Carcinoma (TCGA, Provisional)",
#' "Bladder Urothelial Carcinoma (TCGA, Provisional)",
#' "Brain Lower Grade Glioma (TCGA, Provisional)",
#' "Breast Invasive Carcinoma (TCGA, Provisional)")
#'
#' processMultipleStudies(genes, "test2", studies, "RNA-Seq",
#' calculate = c("frequencyPercentage", "frequencyRatio"), heatmapMargines =
#' c(16,10), RowCex = 1, ColCex = 1)
#'
#' @author Arman Shahrisa, \email{shahrisa.arman@hotmail.com} [maintainer,
#' copyright holder]
#' @author Maryam Tahmasebi Birgani, \email{tahmasebi-ma@ajums.ac.ir}
#'
#' @export



################################################################################
################################################################################
######### Evaluation of Frequency, Mean and Median for multiple Cancers ########
################################################################################
################################################################################

processMultipleStudies <- function(

  genesList,

  submissionName,

  studiesNames,

  desiredTechnique,

  cancerCode = FALSE,

  validateGenes = TRUE,

  calculate = c("frequencyPercentage", "frequencyRatio", "meanValue"),

  cutoff=NULL,

  round=TRUE,

  topGenes = TRUE,

  shortenStudyNames = TRUE,

  geneLimit = 50,

  rankingMethod = "variation",

  heatmapFileFormat = "TIFF",

  resolution = 600,

  RowCex = "auto",

  ColCex = "auto",

  heatmapMargines = "auto",

  rowLabelsAngle = 0,

  columnLabelsAngle = 45,

  heatmapColor = "RdBu",

  reverseColor = TRUE,

  transposedHeatmap = FALSE,

  simplifyBy = FALSE,

  genesToDrop = FALSE,

  transposeResults = FALSE

  ){

  ##############################################################################
  ### Obtaining data

  obtainMultipleStudies(

    genesList = genesList,

    submissionName = submissionName,

    studiesNames = studiesNames,

    desiredTechnique = desiredTechnique,

    cancerCode = cancerCode,

    validateGenes = validateGenes

    )

  message("")


  ##############################################################################
  ### Calculating statistics

  automatedStatistics(

    submissionName = submissionName,

    obtainedDataType = "multiple studies",

    calculate = calculate,

    cutoff = cutoff,

    round = round,

    topGenes = topGenes

    )

  message("")


  ##############################################################################
  ##############################################################################
  ### Create new directory for submission

  present.directory <- getwd()

  new.directory <- paste0(

    present.directory, "/", submissionName, " output for multiple studies"

  )


  dir.create(new.directory , showWarnings = FALSE)

  setwd(new.directory)


  ##############################################################################
  ### Preparing for heatmap output

  heatmapOutput(

    submissionName = submissionName,

    shortenStudyNames = shortenStudyNames,

    geneLimit = geneLimit,

    rankingMethod = rankingMethod,

    heatmapFileFormat = heatmapFileFormat,

    resolution = resolution,

    RowCex = RowCex,

    ColCex = ColCex,

    heatmapMargines = heatmapMargines,

    rowLabelsAngle = rowLabelsAngle,

    columnLabelsAngle = columnLabelsAngle,

    heatmapColor = heatmapColor,

    reverseColor = reverseColor,

    transposedHeatmap = transposedHeatmap,

    simplifyBy = simplifyBy,

    genesToDrop = genesToDrop

    )

  message("")


  ##############################################################################
  ### Preparing for excel output

  xlsxOutput(submissionName = submissionName,

             transposeResults = transposeResults)


  ##############################################################################
  ##############################################################################
  ### Change the directory to the first directory

  setwd(present.directory)

}
