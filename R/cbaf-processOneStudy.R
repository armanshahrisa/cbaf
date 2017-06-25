#' @title Checking Expression/methylation Profile for various subgroups of a cancer study.
#'
#' @description This function Obtains the requested data for the given genes across multiple subgroups of a cancer. It can
#' check whether or not all genes are included in subgroups of a cancer study and, if not, looks for the alternative gene names.
#' The main part of function calculates frequency percentage, frequency ratio, mean expression and median of samples greather than
#' specific value in the selected subgroups of the cancer. Furthermore, it looks for the genes that comprise the highest values
#' in each cancer subgroup.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cBioAutomatedTools \cr
#' Type: \tab Package \cr
#' Version: \tab 0.99.0 \cr
#' Date: \tab 2017-06-22 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#' @import cgdsr xlsxjars xlsx gplots RColorBrewer rafalib Biobase genefilter
#'
#' @usage process.one.study(genes, cancername, high.throughput.data.type,
#' data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value", "Median"),
#' transposedHeatmap=FALSE, desired.case.list="None", genelimit="none",
#' resolution=600, RowCex=0.8, ColCex=0.8, heatmapMargines=c(10,10), cutoff=NULL,
#' angle.for.heatmap.cancernames=45, heatmap.color = "RdBu", reverse.heatmap.color = TRUE,
#' rewrite.output.list = TRUE, round=TRUE, top.genes = TRUE, validate.genes = TRUE,
#' simplify.visulization=FALSE, simplifiction.cuttoff=FALSE)
#'
#' @param genes a list that contains at least one gene group
#'
#' @param cancername a character string showing the desired cancer name. It is an standard cancer study name that can be found on
#' cbioportal.org, such as \code{Acute Myeloid Leukemia (TCGA, NEJM 2013)}.
#'
#' @param high.throughput.data.type a character string that is one of the following techniques: 'RNA-seq', 'microRNA-Seq',
#' 'microarray.mRNA', 'microarray.microRNA' or 'methylation'.
#'
#' @param data.presented.as a character vector that containes the statistical precedures users prefer the function to compute.
#' default is \code{c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value", "Median")}. This will tell the function to
#' compute the following:
#' frequency precentage, which is the number of samples having value greather than specific cutoff divided by the total sample
#' size for each cancer;
#' frequency ratio, which shows the number of selected samples divided by the total number of samples that give the frequency
#' percentage
#' for each cancer-to know selecected and total sample sizes only;
#' Mean Expression, that contains mean value of selected samples for each cancer;
#' Median, which shows the median value of selected samples for each cancer.
#'
#' @param transposedHeatmap a logical value that can transpose the rows and columns of heatmaps.
#'
#' @param desired.case.list a numeric vector that contains the index of rdesired cancer subgroups, if user knows index of desired
#' subgroups. If set to be \code{none}, function will ask the user to enter them during the process. The default value is
#' \code{none}.
#'
#' @param genelimit if large number of genes exists in at least one gene group, this option can be use to limit the number of
#' genes to be shown on hitmap. For instance, \code{genelimit=50} will limit the heatmap to 50 genes showing the most variation.
#' The default value is \code{none}.
#'
#' @param resolution a number. This option can be used to adjust the resolution of the output heatmaps as 'dot per inch'. The
#' defalut value is 600.
#'
#' @param RowCex a number that specifies letter size in heatmap row names.
#'
#' @param ColCex a number that specifies letter size in heatmap column names.
#'
#' @param heatmapMargines a numeric vectors that can be use to set heatmap margins. The default value is
#' \code{heatmapMargines=c(10,10)}.
#'
#' @param cutoff a number used to limit samples to those that are greather than specific number (cutoff). The default value for
#' methylation data is 0.6 while gene expression studies use default value of 2. For methylation studies, it is
#' \code{observed/expected ratio}, for the rest, it is \code{z-score} - Function can understand the correct unit. TO change the
#' cutoff to any desired number, change the option to \code{cutoff = desiredNumber} in which \code{desiredNumber} is the number
#' of interest.
#'
#' @param angle.for.heatmap.cancernames a number that determines the angle with which the cancer names are shown in heatmaps.
#' The default value is 45 degree.
#'
#' @param heatmap.color a character string matches standard color names. The default value is "RdBu". "redgreen" is also a popular
#' color in genomic studies. To see the rest of colors, please type \code{display.brewer.all()}.
#'
#' @param reverse.heatmap.color a logical value that reverses the color gradiant for heatmap.
#'
#' @param round a logical value that, if set to be \code{TRUE}, will force the function to round all the calculated values
#' to two decimal places. The default value is \code{TRUE}.
#'
#' @param top.genes a logical value that, if set as \code{TRUE}, cause the function to create three data.frame that contain the
#' gene names with the highest values for each cancer. To get all the three data.frames, `Frequency.Percentage`, `Mean.Value` and
#' `Median` must have been included for \code{data.presented.as}.
#'
#' @param validate.genes a logical value that, if set to be \code{TRUE}, function will checks each cancer study to finds whether
#' or not each gene has a record. If the given cancer doesn't have a record for specific gene, it checks for alternative gene
#' names that cbioportal might use instead of the given gene name.
#'
#' @param simplify.visulization a logical value that tells the function whether or not to change values under
#' \code{simplifiction.cuttoff} to zero. It only affects heatmaps to assist finding the candidate genes faster. Therefore, it is
#' not suited for publications.
#'
#' @param simplifiction.cuttoff a logical value that, if \code{simplify.visulization = TRUE}, needs to be set as a desired cuttoff
#' for \code{simplify.visulization}. It has the same unit as \code{cutoff}.
#'
#' @return a list that containes some or all of the following groups, based on what user has chosen: \code{Validation.Results},
#' \code{Frequency.Percentage}, \code{Top.Genes.of.Frequency.Percentage}, \code{Frequency.Ratio}, \code{Mean.Value},
#' \code{Top.Genes.of.Mean.Value}, \code{Median}, \code{Top.Genes.of.Median}. It also saves these groups in one excel
#' file for convenience. Based on preference, three heatmaps for \code{Frequency.Percentage}, \code{Mean.Value} and
#' \code{Median} can be generated. If more than one gene group is entered, output for each group will be strored in a separate sub-directory.
#'
#' @examples
#' # Creating a list that contains one gene group: 'K.demethylases'
#' genes <- list(K.demethylases = c("KDM1A", "KDM1B", "KDM2A"))
#'
#' # Running the function to obtain and process the selected data
#' process.one.study(genes, "Breast Invasive Carcinoma (TCGA, Cell 2015)", "RNA-seq", desired.case.list = c(3,4,5),
#' data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value"), heatmap.color = "redgreen")
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

                            transposedHeatmap = FALSE, simplify = FALSE, simplifictionCuttoff = FALSE){

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

  dir.create(paste(present.directory, "/", submissionName, " output for  single study", sep = ""), showWarnings = FALSE)

  setwd(paste(present.directory, "/", submissionName, " output for a single study", sep = ""))



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
