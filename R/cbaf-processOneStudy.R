#' @title Checking Expression/methylation Profile for various subgroups of a cancer study.
#'
#' @description This function Obtains the requested data for the given genes across multiple subgroups of a cancer. It can
#' check whether or not all genes are included in subgroups of a cancer study and, if not, looks for the alternative gene names.
#' Then it calculates frequency percentage, frequency ratio, mean value and median value of samples greather than
#' specific value in the selected subgroups of the cancer. Furthermore, it looks for the five genes that comprise the highest values
#' in each cancer subgroup.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cbaf \cr
#' Type: \tab Package \cr
#' Version: \tab 0.99.1 \cr
#' Date: \tab 2017-08-19 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#'
#' @include cbaf-obtainOneStudy.R cbaf-automatedStatistics.R cbaf-heatmapOutput.R cbaf-xlsxOutput.R
#'
#' @usage processOneStudy(genesList, submissionName, studyName, desiredTechnique,
#' desiredCaseList = FALSE, validateGenes = TRUE,
#' calculate = c("frequencyPercentage", "frequencyRatio", "meanValue", "medianValue"),
#' cutoff=NULL, round=TRUE, topGenes = TRUE, shortenStudyNames = TRUE,
#' genelimit = "none", resolution = 600, RowCex = 0.8, ColCex = 0.8,
#' heatmapMargines = c(10,10), angleForYaxisNames = 45, heatmapColor = "RdBu",
#' reverseColor = TRUE, transposedHeatmap = FALSE, simplify = FALSE,
#' simplifictionCuttoff = FALSE, genesToDrop = NULL)
#'
#'
#'
#' @param genesList a list that contains at least one gene group
#'
#' @param submissionName a character string containing name of interest. It is used for naming the process.
#'
#' @param studyName a character string showing the desired cancer name. It is an standard cancer study name that can be found on
#' cbioportal.org, such as 'Acute Myeloid Leukemia (TCGA, NEJM 2013)'.
#'
#' @param desiredTechnique a character string that is one of the following techniques: 'RNA-seq', 'microRNA-Seq', 'microarray.mRNA',
#' 'microarray.microRNA' or 'methylation'.
#'
#' @param desiredCaseList a numeric vector that contains the index of desired cancer subgroups, assuming the user knows index of
#' desired subgroups. If not, desiredCaseList has been set as 'none', function will show the available subgroups and ask the
#' user to enter the desired ones during the process. The default value is 'none'.
#'
#' @param validateGenes a logical value that, if set to be 'TRUE', causes the function to check each cancer study to find whether
#' or not each gene has a record. If a cancer doesn't have a record for specific gene, function looks for alternative gene
#' names that cbioportal might use instead of the given gene name.
#'
#' @param calculate a character vector that containes the statistical procedures users prefer the function to compute.
#' Default input is \code{c("frequencyPercentage", "frequencyRatio", "meanValue", "medianValue")}. This will tell the function to
#' compute the following:
#' 'frequencyPercentage', which is the percentge of samples having the value greather than specific cutoff divided by the total sample
#' size for every study / study subgroup;
#' 'frequency ratio', which shows the number of selected samples divided by the total number of samples that give the frequency
#' percentage for every study / study subgroup. It shows the selected and total sample sizes.;
#' 'Mean Value', that contains mean value of selected samples for each study;
#' 'Median Value', which shows the median value of selected samples for each study.
#'
#' @param cutoff a number used to limit samples to those that are greather than specific number (cutoff). The default value for
#' methylation data is 0.6 while gene expression studies use default value of 2. For methylation studies, it is
#' \code{observed/expected ratio}, for the rest, it is 'z-score'. To change the cutoff to any desired number, change the
#' option to \code{cutoff = desiredNumber}, in which desiredNumber is the number of interest.
#'
#' @param round a logical value that, if set to be TRUE, will force the function to round all the calculated values
#' to two decimal places. The default value is TRUE.
#'
#' @param topGenes a logical value that, if set as TRUE, causes the function to create three data.frame that contain the five
#' top genes for each cancer. To get all the three data.frames, "frequencyPercentage", "meanValue" and "medianValue" must have been
#' included for \code{calculate}.
#'
#' @param shortenStudyNames a logical vector. If the value is set as TRUE, function will try to remove the last part of the
#' cancer names aiming to shorten them. The removed segment usually contains the name of scientific group that has conducted
#' the experiment.
#'
#' @param genelimit if large number of genes exist in at least one gene group, this option can be used to limit the number of
#' genes that are shown on heatmap. For instance, \code{genelimit=50} will limit the heatmap to 50 genes showing the most variation
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
#' @param heatmapMargines a numeric vector that is used to set heatmap margins. The default value is
#' \code{heatmapMargines=c(15,07)}.
#'
#' @param angleForYaxisNames a number that determines the angle with which the studies/study subgroups names are shown in heatmaps.
#' The default value is 45 degree.
#'
#' @param heatmapColor a character string that defines heatmap color. The default value is "RdBu". "redgreen" is also a popular
#' color in genomic studies. To see the rest of colors, please type \code{library(RColorBrewer)} and then \code{display.brewer.all()}.
#'
#' @param reverseColor a logical value that reverses the color gradiant for heatmap(s).
#'
#' @param transposedHeatmap a logical value that transposes heatmap rows to columns and vice versa.
#'
#' @param simplify a logical value that tells the function whether or not to change values under
#' \code{simplifiction.cuttoff} to zero. The purpose behind this option is to facilitate seeing candidate genes. Therefore, it is
#' not suited for publications. Default value is 'FALSE'.
#'
#' @param simplifictionCuttoff a logical value that, if \code{simplify.visulization = TRUE}, needs to be set as a desired cuttoff
#' for \code{simplify.visulization}. It has the same unit as \code{cutoff}.
#'
#' @param genesToDrop a character vector. Gene names within this vector will be omitted from heatmap.
#'
#'
#'
#' @return a BiocFileCache object that containes some or all of the following groups, based on what user has chosen: Obtained.Data, Validation.Results,
#' Frequency.Percentage, Top.Genes.of.Frequency.Percentage, Frequency.Ratio, mean.Value, Top.Genes.of.Mean.Value, median.Value,
#' Top.Genes.of.Median.Value. It also saves these results in one excel files for convenience. Based on preference, three heatmaps for
#' frequency percentage, mean value and median can be generated. If more than one group of genes is entered, output for each group
#' will be strored in a separate sub-directory.
#'
#' @examples
#' # Sample BiocFileCache object which is created by the following code:
#' bfc_test <- BiocFileCache::BiocFileCache(system.file("extdata", "test", package = "cbaf"))
#'
#'
#' # Example of function usage:
#' genes <- list(K.demethylases = c("KDM1A", "KDM1B", "KDM2A", "KDM2B", "KDM3A",
#'  "KDM3B", "JMJD1C", "KDM4A"), K.methyltransferases = c("SUV39H1", "SUV39H2",
#'  "EHMT1", "EHMT2", "SETDB1", "SETDB2", "KMT2A", "KMT2A"))
#'
#' processOneStudy(genes, "test", "Breast Invasive Carcinoma (TCGA, Cell 2015)",
#' "RNA-seq", desiredCaseList = c(2,3,4,5), calculate = c("frequencyPercentage",
#' "frequencyRatio"), heatmapMargines = c(16, 10), RowCex = 1, ColCex = 1)
#'
#' @author Arman Shahrisa, \email{shahrisa.arman@hotmail.com} [maintainer, copyright holder]
#' @author Maryam Tahmasebi Birgani, \email{tahmasebi-ma@ajums.ac.ir}
#'
#' @export



###################################################################################################
###################################################################################################
########## Evaluation of Median, Frequency and ExpressionMean for Subtypes of a Cancer ############
###################################################################################################
###################################################################################################

processOneStudy <- function(genesList, submissionName, studyName, desiredTechnique, desiredCaseList = FALSE, validateGenes = TRUE,

                            calculate = c("frequencyPercentage", "frequencyRatio", "meanValue", "medianValue"), cutoff=NULL, round=TRUE,

                            topGenes = TRUE, shortenStudyNames = TRUE, genelimit = "none", resolution = 600, RowCex = 0.8, ColCex = 0.8,

                            heatmapMargines = c(10,10), angleForYaxisNames = 45, heatmapColor = "RdBu", reverseColor = TRUE,

                            transposedHeatmap = FALSE, simplify = FALSE, simplifictionCuttoff = FALSE, genesToDrop = NULL){

  ##########################################################################
  ### Obtaining data

  obtainOneStudy(genesList = genesList, submissionName = submissionName, studyName = studyName, desiredTechnique = desiredTechnique,

                 desiredCaseList = desiredCaseList, validateGenes = validateGenes)



  ##########################################################################
  ### Calculating statistics

  automatedStatistics(submissionName = submissionName, obtainedDataType = "single study", calculate = calculate, cutoff = cutoff,

                      round = round, topGenes = topGenes)



  ################################################################################
  ################################################################################
  ### Create new directory for submission

  present.directory <- getwd()

  dir.create(paste(present.directory, "/", submissionName, " output for a single study", sep = ""), showWarnings = FALSE)

  setwd(paste(present.directory, "/", submissionName, " output for a single study", sep = ""))



  ##########################################################################
  ### Preparing for heatmap output

  heatmapOutput(submissionName = submissionName, shortenStudyNames = shortenStudyNames, genelimit = genelimit, resolution = resolution,

                RowCex = RowCex, ColCex = ColCex, heatmapMargines = heatmapMargines, angleForYaxisNames = angleForYaxisNames,

                heatmapColor = heatmapColor, reverseColor = reverseColor, transposedHeatmap = transposedHeatmap, simplify = simplify,

                simplifictionCuttoff = simplifictionCuttoff, genesToDrop = genesToDrop)



  ##########################################################################
  ### Preparing for excel output

  xlsxOutput(submissionName = submissionName)



  ################################################################################
  ################################################################################
  ### Change the directory to the first directory

  setwd(present.directory)

}
