#' @title Checking Expression/methylation Profile for various cancer studies.
#'
#' @description This function Obtains the requested data for the given genes across multiple cancer studies. It can check
#' whether or not all genes are included in cancer studies and, if not, looks for the alternative gene names. Tha main part
#' of function calculates frequency percentage, frequency ratio, mean value and median value of samples greather than specific
#' value in the selected cancer studies. Furthermore, it looks for the genes that comprise the highest values in each cancer study.
#' It can generate three meatmaps for
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cbaf \cr
#' Type: \tab Package \cr
#' Version: \tab 0.99.0 \cr
#' Date: \tab 2017-06-28 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#'
#' @include cbaf-obtainMultipleStudies.R cbaf-automatedStatistics.R cbaf-heatmapOutput.R cbaf-xlsxOutput.R
#'
#' @usage processMultipleStudies(genesList, submissionName, studiesNames, desiredTechnique,
#' cancerCode = FALSE, validateGenes = TRUE,
#' calculate = c("frequencyPercentage", "frequencyRatio", "meanValue", "medianValue"),
#' cutoff=NULL, round=TRUE, topGenes = TRUE, shortenStudyNames = TRUE, genelimit = "none",
#' resolution = 600, RowCex = 0.8, ColCex = 0.8, heatmapMargines = c(10,10), angleForYaxisNames = 45,
#' heatmapColor = "RdBu", reverseColor = TRUE, transposedHeatmap = FALSE, simplify = FALSE,
#' simplifictionCuttoff = FALSE)
#'
#'
#'
#' @param genesList a list that contains at least one gene group
#'
#' @param submissionName a character string containing name of interest. It is used for naming the process.
#'
#' @param studiesNames a character vector or a matrix that containes desired cancer names. The character vector containes standard
#' names of cancer studies that can be found on cbioportal.org, such as 'Acute Myeloid Leukemia (TCGA, NEJM 2013)'. Alternatively,
#' a matrix can be used if users prefer user-defined cancer names, in which the first column of matrix comprises the standard
#' cancer names while the second column must contain the desired cancer names.
#'
#' @param desiredTechnique a character string that is one of the following techniques: 'RNA-seq', 'microRNA-Seq', 'microarray.mRNA',
#' 'microarray.microRNA' or 'methylation'.
#'
#' @param cancerCode a logical value that tells the function to use cbioportal abbreviated cancer names instead of complete cancer
#' names, if set to be 'TRUE'. For example, 'laml_tcga_pub' is the abbreviated name for 'Acute Myeloid Leukemia (TCGA, NEJM 2013)'.
#'
#' @param validateGenes a logical value that, if set to be 'TRUE', function will check each cancer study to find whether
#' or not each gene has a record. If the given cancer doesn't have a record for specific gene, it checks for alternative gene
#' names that cbioportal might use instead of the given gene name.
#'
#' @param calculate a character vector that containes the statistical precedures users prefer the function to compute.
#' Default input is \code{c("frequencyPercentage", "frequencyRatio", "Mean.Value", "medianValue")}. This will tell the function to
#' compute the following:
#' 'frequencyPercentage', which is the number of samples having the value greather than specific cutoff divided by the total sample
#' size for every study / study subgroup;
#' 'frequency ratio', which shows the number of selected samples divided by the total number of samples that give the frequency
#' percentage for every study / study subgroup -to know selecected and total sample sizes only;
#' 'Mean Value', that contains mean value of selected samples for each study;
#' 'Median Value', which shows the median value of selected samples for each study.
#'
#' @param cutoff a number used to limit samples to those that are greather than specific number (cutoff). The default value for
#' methylation data is 0.6 while gene expression studies use default value of 2. For methylation studies, it is
#' \code{observed/expected ratio}, for the rest, it is 'z-score'. TO change the cutoff to any desired number, change the
#' option to \code{cutoff = desiredNumber} in which desiredNumber is the number of interest.
#'
#' @param round a logical value that, if set to be TRUE, will force the function to round all the calculated values
#' to two decimal places. The default value is TRUE.
#'
#' @param topGenes a logical value that, if set as TRUE, causes the function to create three data.frame that contain the five
#' top genes for each cancer. To get all the three data.frames, "Frequency.Percentage", "Mean.Value" and "Median" must have been
#' included for \code{calculate}.
#'
#' @param shortenStudyNames a logical vector. If the value is set as true, function will try to remove the end part of
#' cancer names aiming to shorten them. The removed segment usually contains the name of scientific group that has conducted
#' the experiment.
#'
#' @param genelimit if large number of genes exist in at least one gene group, this option can be used to limit the number of
#' genes that are shown on hitmap. For instance, \code{genelimit=50} will limit the heatmap to 50 genes showing the most variation
#' across multiple study / study subgroups.
#' The default value is \code{none}.
#'
#' @param resolution a number. This option can be used to adjust the resolution of the output heatmaps as 'dot per inch'.
#' The defalut value is 600.
#'
#' @param RowCex a number that specifies letter size in heatmap row names.
#'
#' @param ColCex a number that specifies letter size in heatmap column names.
#'
#' @param heatmapMargines a numeric vector that can be used to set heatmap margins. The default value is
#' \code{heatmapMargines=c(15,07)}.
#'
#' @param angleForYaxisNames a number that determines the angle with which the studies/study subgroups names are shown in heatmaps.
#' The default value is 45 degree.
#'
#' @param heatmapColor a character string that defines heatmao color. The default value is "RdBu". "redgreen" is also a popular
#' color in genomic studies. To see the rest of colors, please type \code{display.brewer.all()}.  Default value is 'TRUE'.
#'
#' @param reverseColor a logical value that reverses the color gradiant for heatmap(s).
#'
#' @param transposedHeatmap a logical value that changes row and colums of heatmap.
#'
#' @param simplify a logical value that tells the function whether or not to change values under
#' \code{simplifiction.cuttoff} to zero. The purpose behind this option is to facilitate seeing candidate genes. Therefore, it is
#' not suited for publications. Default value is 'FALSE'.
#'
#' @param simplifictionCuttoff a logical value that, if \code{simplify.visulization = TRUE}, needs to be set as a desired cuttoff
#' for \code{simplify.visulization}. It has the same unit as \code{cutoff}.
#'
#'
#'
#' @return a list that containes some or all of the following groups, based on what user has chosen: ValidationResults,
#' Frequency.Percentage, Top.Genes.of.Frequency.Percentage, Frequency.Ratio, mean.Value, Top.Genes.of.Mean.Value, median.Value,
#' Top.Genes.of.Median. It also saves these groups in one excel files for convenience. Based on preference, three heatmaps for
#' Frequency.Percentage, Mean.Value and Median can be generated. If more than one group of genes is entered, output for each group
#' will be strored in a separate sub-directory.
#'
#' @examples
#' genes <- list(K.acetyltransferases = c("CLOCK", "CREBBP", "ELP3", "EP300"))
#'
#' studies <- c("Acute Myeloid Leukemia (TCGA, Provisional)",
#' "Adrenocortical Carcinoma (TCGA, Provisional)")
#'
#' processMultipleStudies(genes, "test2", studies, "RNA-seq",
#' calculate = c("frequencyPercentage", "frequencyRatio", "meanValue",
#' "medianValue"), heatmapMargines = c(15,10))
#'
#' @author Arman Shahrisa, \email{shahrisa.arman@hotmail.com} [maintainer, copyright holder]
#' @author Maryam Tahmasebi Birgani, \email{tahmasebi-ma@ajums.ac.ir}
#'
#' @export



###################################################################################################
###################################################################################################
############ Evaluation of Median, Frequency and ExpressionMean for multiple Cancers ##############
###################################################################################################
###################################################################################################

processMultipleStudies <- function(genesList, submissionName, studiesNames, desiredTechnique, cancerCode = FALSE, validateGenes = TRUE,

                                   calculate = c("frequencyPercentage", "frequencyRatio", "meanValue", "medianValue"), cutoff=NULL, round=TRUE,

                                   topGenes = TRUE, shortenStudyNames = TRUE, genelimit = "none", resolution = 600, RowCex = 0.8, ColCex = 0.8,

                                   heatmapMargines = c(10,10), angleForYaxisNames = 45, heatmapColor = "RdBu", reverseColor = TRUE,

                                   transposedHeatmap = FALSE, simplify = FALSE, simplifictionCuttoff = FALSE){

  ##########################################################################
  ### Obtaining data

  obtainMultipleStudies(genesList = genesList, submissionName = submissionName, studiesNames = studiesNames, desiredTechnique = desiredTechnique,

                        cancerCode = cancerCode, validateGenes = validateGenes)



  ##########################################################################
  ### Calculating statistics

  automatedStatistics(submissionName = submissionName, obtainedDataType = "multiple studies", calculate = calculate, cutoff = cutoff,

                      round = round, topGenes = topGenes)



  ################################################################################
  ################################################################################
  ### Create new directory for submission

  present.directory <- getwd()

  dir.create(paste(present.directory, "/", submissionName, " output for multiple studies", sep = ""), showWarnings = FALSE)

  setwd(paste(present.directory, "/", submissionName, " output for multiple studies", sep = ""))



  ##########################################################################
  ### Preparing for heatmap output

  heatmapOutput(submissionName = submissionName, shortenStudyNames = shortenStudyNames, genelimit = genelimit, resolution = resolution,

                RowCex = RowCex, ColCex = ColCex, heatmapMargines = heatmapMargines, angleForYaxisNames = angleForYaxisNames,

                heatmapColor = heatmapColor, reverseColor = reverseColor, transposedHeatmap = transposedHeatmap, simplify = simplify,

                simplifictionCuttoff = simplifictionCuttoff)



  ##########################################################################
  ### Preparing for excel output

  xlsxOutput(submissionName = submissionName)



  ################################################################################
  ################################################################################
  ### Change the directory to the first directory

  setwd(present.directory)

}
