#' @title Generate heatmaps for various studies/subgroups of a study.
#'
#' @description This function can prepare heatmap for 'frequency percentage',
#' 'mean value' and 'median value' data provided by
#' automatedStatistics() function.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cbaf \cr
#' Type: \tab Package \cr
#' Version: \tab 1.1.1 \cr
#' Date: \tab 2017-11-11 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#'
#'
#' @importFrom genefilter rowVars
#'
#' @importFrom RColorBrewer brewer.pal
#'
#' @importFrom gplots heatmap.2 redgreen
#'
#' @importFrom BiocFileCache bfcnew bfcquery bfcpath
#'
#' @importFrom grDevices colorRampPalette dev.off tiff png bmp jpeg
#'
#' @importFrom utils head setTxtProgressBar txtProgressBar
#'
#'
#'
#' @include cbaf-obtainOneStudy.R cbaf-obtainMultipleStudies.R
#' cbaf-automatedStatistics.R
#'
#'
#'
#' @usage heatmapOutput(submissionName, shortenStudyNames = TRUE,
#'   genelimit = "none", rankingMethod = "variation", heatmapFileFormat = "TIFF"
#'   , resolution = 600, RowCex = 0.8, ColCex = 0.8, heatmapMargines = c(15,07),
#'   angleForYaxisNames = 45, heatmapColor = "RdBu", reverseColor = TRUE,
#'   transposedHeatmap = FALSE, simplify = FALSE, simplifictionCuttoff = FALSE,
#'   genesToDrop = NULL)
#'
#'
#'
#' @param submissionName a character string containing name of interest. It is
#' used for naming the process.
#'
#' @param shortenStudyNames a logical vector. If the value is set as TRUE,
#' function will try to remove the last part of the cancer names aiming to
#' shorten them. The removed segment usually contains the name of scientific
#' group that has conducted the experiment.
#'
#' @param genelimit if large number of genes exist in at least one gene group,
#' this option can be used to limit the number of genes that are shown on
#' heatmap. For instance, \code{genelimit=50} will limit the heatmap to 50 genes
#' that show the most variation across multiple study / study subgroups. The
#' default value is \code{none}.
#'
#' @param rankingMethod a character value that determines how genes will be
#' ranked prior to drawing heatmap. \code{"variation"} orders the genes based on
#' unique values in one or few cancer studies while \code{"highValue"} ranks the
#'  genes when they cotain high values in multiple / many cancer studies. This
#'  option is useful when number of genes are too much so that user has to limit
#'  the number of genes on heatmap by \code{genelimit}.
#'
#' @param heatmapFileFormat This option enables the user to select the desired
#' image file format of the heatmaps. The default value is \code{"TIFF"}. Other
#' suppoeted formats include \code{"PNG"}, \code{"BMP"}, and \code{"JPG"}.
#'
#' @param resolution a number. This option can be used to adjust the resolution
#' of the output heatmaps as 'dot per inch'. The defalut value is 600.
#'
#' @param RowCex a number that specifies letter size in heatmap row names.
#'
#' @param ColCex a number that specifies letter size in heatmap column names.
#'
#' @param heatmapMargines a numeric vector that is used to set heatmap margins.
#' The default value is \code{heatmapMargines=c(15,07)}.
#'
#' @param angleForYaxisNames a number that determines the angle with which the
#' studies/study subgroups names are shown in heatmaps. The default value is 45
#' degree.
#'
#' @param heatmapColor a character string that defines heatmap color. The
#' default value is "RdBu". "redgreen" is also a popular color in genomic
#' studies. To see the rest of colors, please type \code{library(RColorBrewer)}
#' and then \code{display.brewer.all()}.
#'
#' @param reverseColor a logical value that reverses the color gradiant for
#' heatmap(s).
#'
#' @param transposedHeatmap a logical value that transposes heatmap rows to
#' columns and vice versa.
#'
#' @param simplify a logical value that tells the function whether or not to
#' change values under \code{simplifictionCuttoff} to zero. The purpose behind
#' this option is to facilitate seeing candidate genes. Therefore, it is not
#' suited for publications.
#'
#' @param simplifictionCuttoff a logical value that, if
#' \code{simplify.visulization = TRUE}, needs to be set as a desired cuttoff for
#' \code{simplify.visulization}. It has the same unit as \code{cutoff}.
#'
#' @param genesToDrop a character vector. Gene names within this vector will be
#' omitted from heatmap.
#'
#'
#'
#' @return Based on preference, three heatmaps for \code{"Frequency.Percentage"}
#' , \code{"Mean.Value"} and \code{"Median.value"} can be generated. If more
#' than one group of genes are entered, output for each group will be strored in
#'  a separate sub-directory.
#'
#'
#'
#' @examples
#' genes <- list(K.demethylases = c("KDM1A", "KDM1B", "KDM2A", "KDM2B", "KDM3A",
#'  "KDM3B", "JMJD1C", "KDM4A"), K.methyltransferases = c("SUV39H1", "SUV39H2",
#'  "EHMT1", "EHMT2", "SETDB1", "SETDB2", "KMT2A", "KMT2A"))
#'
#' obtainOneStudy(genes, "test", "Breast Invasive Carcinoma (TCGA, Cell 2015)",
#' "RNA-seq", desiredCaseList = c(3,4))
#'
#' automatedStatistics("test", obtainedDataType = "single study", calculate =
#' c("frequencyPercentage", "frequencyRatio"))
#'
#' heatmapOutput(submissionName = "test", heatmapMargines = c(16, 10),
#' RowCex = 1, ColCex = 1)
#'
#'
#'
#' @author Arman Shahrisa, \email{shahrisa.arman@hotmail.com} [maintainer,
#' copyright holder]
#' @author Maryam Tahmasebi Birgani, \email{tahmasebi-ma@ajums.ac.ir}s
#'
#' @export



################################################################################
################################################################################
################## Generating heatmap for the processed data ###################
################################################################################
################################################################################

heatmapOutput <- function(

  submissionName,

  shortenStudyNames = TRUE,

  genelimit = "none",

  rankingMethod = "variation",

  heatmapFileFormat = "TIFF",

  resolution = 600,

  RowCex = 0.8,

  ColCex = 0.8,

  heatmapMargines = c(15,07),

  angleForYaxisNames = 45,

  heatmapColor = "RdBu",

  reverseColor = TRUE,

  transposedHeatmap = FALSE,

  simplify = FALSE,

  simplifictionCuttoff = FALSE,

  genesToDrop = NULL

  ){

  ##############################################################################
  ########## Prerequisites

  # Check submissionName

  if(!is.character(submissionName)){

    stop("'submissionName' must be entered as a character string for naming the process")

  }



  # Check heatmap image file format

  if(!(heatmapFileFormat %in% c("TIFF", "PNG", "JPG", "BMP"))){

    stop("'heatmapFileFormat' must be one of the supported image formats: 'TIFF', 'PNG', 'JPG', or 'BMP'")

  }



  ##############################################################################
  ########## Decide whether function should stops now!

  # Check wheather the requested data exists

  database <- system.file("extdata", submissionName, package="cbaf")

  if(!dir.exists(database)){

    stop("Please run one of the obtainSingleStudy() or obtainMultipleStudies() functions and then the automatedStatistics() function first")

  } else if(dir.exists(database)){

    bfc <- BiocFileCache(

      file.path(system.file("extdata", package = "cbaf"), submissionName)

      )

    if(!nrow(bfcquery(bfc, c("Parameters for automatedStatistics()"))) == 1){

      stop("Please run the automatedStatistics() function first")

    }

  }



  # obtain parameters for prevous function

  previousFunctionParam <-

    readRDS(

      bfcpath(bfc, bfcquery(bfc, c("Parameters for automatedStatistics()"))$rid)

      )




  # fetch an old parameter from the previous function

  desiredTechnique <- previousFunctionParam$desiredTechnique

  cutoff <- previousFunctionParam$cutoff



  # setting the value for cutoff

  if(desiredTechnique == "methylation"){

    cutoff.phrase <- "obs/exp cutoff"

  } else{

    cutoff.phrase <- "z-score cutoff"

  }





  # Store the new parameteres

  newParameters <-list()

  newParameters$submissionName <- submissionName

  newParameters$shortenStudyNames <- shortenStudyNames

  newParameters$genelimit <- genelimit

  newParameters$resolution <- resolution

  newParameters$RowCex <- RowCex

  newParameters$ColCex <- ColCex

  newParameters$heatmapMargines <- heatmapMargines

  newParameters$angleForYaxisNames <- angleForYaxisNames

  newParameters$heatmapColor <- heatmapColor

  newParameters$reverseColor <- reverseColor

  newParameters$transposedHeatmap <- transposedHeatmap

  newParameters$simplify <- simplify

  newParameters$simplifictionCuttoff <- simplifictionCuttoff

  newParameters$genesToDrop <- genesToDrop





  # Check wheather the requested data exists

  number.of.rows.parameters <-

    nrow(bfcquery(bfc, "Parameters for heatmapOutput()"))


  if(number.of.rows.parameters == 1){

    oldParameters <-

      readRDS(

        bfcpath(bfc, bfcquery(bfc, c("Parameters for heatmapOutput()"))$rid)

        )

    # Check whether the previous function is skipped

    if(previousFunctionParam$lastRunStatus == "skipped"){

      if(identical(oldParameters, newParameters)){

        continue <- FALSE

      } else{

        continue <- TRUE

      }

    } else{

      continue <- TRUE

    }

  } else{

    continue <- TRUE

  }





  # Getting the source data

  statisticsData <-

    readRDS(bfcpath(bfc, bfcquery(bfc, c("Calculated statistics"))$rid))

  if(!is.list(statisticsData)){

    stop("Input database must be a list.")

  }





  # get the working directory

  parent.directory <- getwd()



  ##############################################################################
  ########## Set the function ready to work

  # Report

  message("***", " Preparing the requested heatmaps for ", submissionName, " ***")

  if(simplify & is.numeric(simplifictionCuttoff) & !is.numeric(genelimit)){

    warning("--- Only significant results will be used to draw heatmaps ---")

  }



  # Count number of skipped heatmaps

  skipped <- 0



  # Create progressbar

  possible.subgroups <- c("Frequency.Percentage", "Mean.Value", "Median.Value")


  idx <- names(statisticsData[[1]]) %in% possible.subgroups

  total.number <- length(statisticsData)*length((statisticsData[[1]])[idx])


  heatmapOutputProgressBar <-

    txtProgressBar(min = 0, max = total.number, style = 3)

  ExtH <- 0





  ##############################################################################
  ########## Core segment

  # Save heatmaps in separate folder

  for(gr in seq_along(statisticsData)){

    # Subset data that can be presented as heatmap

    subset.name <- names(statisticsData)[gr]


    possible.subgroups.idx <-

      names(statisticsData[[gr]]) %in% possible.subgroups


    subset.data <- (statisticsData[[gr]])[possible.subgroups.idx]



    # Create a directory and set it as desired folder

    child.directory <- paste(

      gr, ". ", sub(x = subset.name, pattern = "\\.", replacement = "-"), sep=""

      )

    new.directory <- paste(parent.directory, child.directory, sep = "/")

    dir.create(new.directory, showWarnings = FALSE)

    setwd(new.directory)



    for(possible in seq_along(subset.data)){

      # subset statistics

      statistics.data <- subset.data[[possible]]



      # Remove desired genes

      if(!is.null(genesToDrop)){

        if(!is.character(genesToDrop)){

          stop("Please enter the desired genes to drop as a character vector.")

        } else{

          filtered.colnames <- !(colnames(statistics.data) %in% genesToDrop)

          statistics.data <- statistics.data[, filtered.colnames, drop = FALSE]

        }

      }



      name.statistics.data <- names(subset.data)[possible]


      # determine ourput file name

      output.file.name <- paste0(

        gsub(x = name.statistics.data, pattern = "\\.", replacement = "-"),

        ", ",

        gsub(x = subset.name, pattern = "\\.", replacement = " "),

        " (",

        cutoff.phrase,

        "=",

        cutoff,

        ")",

        if(heatmapFileFormat == "TIFF"){

          ".tiff"

        }else if(heatmapFileFormat == "PNG"){

          ".png"

        }else if(heatmapFileFormat == "JPG"){

          ".jpg"

        }else if(heatmapFileFormat == "BMP"){

          ".bmp"

        }

      )






      # Check continue permission

      if(continue | !continue & !file.exists(output.file.name)){




        # Check whether study names should be shorted

        if(shortenStudyNames){

          rownames(statistics.data) <-

            sapply(

              strsplit(

                as.character(rownames(statistics.data)),

                split=" (",

                fixed=TRUE

                ),

              function(x) (x[1])

              )

        }



        heatmap.data <- t(statistics.data)

        # Removing NA

        not.just.na <- apply(heatmap.data, 1, function(x) any(!is.na(x)==TRUE))

        heatmap.data <- heatmap.data[not.just.na,]

        # Removing NaN

        heatmap.data[is.nan(heatmap.data)] <- 0

        # Removing rows that contain only 0

        if(is.matrix(heatmap.data)){

          heatmap.data <- heatmap.data[rowSums(heatmap.data, na.rm = TRUE)!=0,]

        }

        if(is.matrix(heatmap.data)){

          if(nrow(heatmap.data) > 1 & ncol(heatmap.data) > 1){




            # !!! Simplifying heatmap for easy assessment !!!

            if(simplify &

               is.numeric(simplifictionCuttoff) &

               !is.numeric(genelimit)){

              heatmap.data[heatmap.data < simplifictionCuttoff] <- 0

            } else if(simplify &

                      !is.numeric(simplifictionCuttoff)){

              stop("Please set your desired cutoff for simplification in 'simplifictionCuttoff'")

            } else if(simplify &

                      is.numeric(simplifictionCuttoff) &

                      is.numeric(genelimit)){

              warning("There is no need to limit gene number because simplification option possibly limits the gene number even below the specified number")

            }




            # Limiting the number of genes in heatmap to get better resolution

            if(genelimit=="none" | genelimit > ncol(heatmap.data)){

              heatmap.data <- heatmap.data

            } else if(is.numeric(genelimit) & genelimit <= ncol(heatmap.data) &

                      rankingMethod == "variation"){

              ordering <- order(abs(rowVars(heatmap.data)), decreasing=TRUE)

              heatmap.data <- heatmap.data[ordering[seq_len(genelimit)],]

            } else if(is.numeric(genelimit) & genelimit <= ncol(heatmap.data) &

                      rankingMethod == "highValue"){

              ordering <- order(abs(rowSums(heatmap.data)), decreasing=TRUE)

              heatmap.data <- heatmap.data[ordering[seq_len(genelimit)],]

            } else{

              stop("Please type gene number limit or if whole genes are desired please type none")
            }





            # Heatmap color

            if(reverseColor){

              if(heatmapColor == "redgreen"){

                hmcol <- rev(redgreen(75))

              } else {

                hmcol <- rev(colorRampPalette(brewer.pal(9, heatmapColor))(100))

              }


            } else if (!reverseColor){

              if(heatmapColor == "redgreen"){

                hmcol <- redgreen(75)

              } else {

                hmcol <- colorRampPalette(brewer.pal(9, heatmapColor))(100)

              }

            }




            # Drawing heatmap

            if(heatmapFileFormat == "TIFF"){


              tiff(

                filename=paste(getwd(), output.file.name, sep="/"),

                width=9.5,

                height= 11,

                units = "in",

                res=resolution

              )


            }else if(heatmapFileFormat == "PNG"){


              png(

                filename=paste(getwd(), output.file.name, sep="/"),

                width=9.5,

                height= 11,

                units = "in",

                res=resolution

              )


            }else if(heatmapFileFormat == "BMP"){


              bmp(

                filename=paste(getwd(), output.file.name, sep="/"),

                width=9.5,

                height= 11,

                units = "in",

                res=resolution

              )


            }else if(heatmapFileFormat == "JPG"){


              jpeg(

                filename=paste(getwd(), output.file.name, sep="/"),

                width=9.5,

                height= 11,

                units = "in",

                res=resolution

              )


            }



            if(!transposedHeatmap){

              heatmap.2(

                heatmap.data,

                labCol=colnames(heatmap.data),

                na.color= "light gray",

                trace="none",

                symbreaks = TRUE,

                col= hmcol,

                cexRow = RowCex,

                cexCol= ColCex,

                margins = heatmapMargines,

                srtCol = angleForYaxisNames

                )

            } else if(transposedHeatmap){

              heatmap.2(

                t(heatmap.data),

                labCol=colnames(heatmap.data),

                na.color="light gray",

                trace="none",

                symbreaks = TRUE,

                col=hmcol,

                cexRow = RowCex,

                cexCol= ColCex,

                margins = heatmapMargines,

                srtCol = angleForYaxisNames)

            }

            dev.off()

          }

        }

      }else{

        skipped <- skipped + 1

      }

      # Update progressbar

      ExtH <- ExtH + 1

      setTxtProgressBar(heatmapOutputProgressBar, ExtH)

    }

  }

  # Close progressbar

  close(heatmapOutputProgressBar)



  # report number of skipped heatmaps

  if(skipped > 0 & skipped != 1){

    message("--- ", as.character(skipped), " out of ", as.character(total.number)," heatmaps were skipped: They already exist. ---")

  } else if(skipped > 0 & skipped == 1){

    message("--- ", as.character(skipped), " out of ", as.character(total.number)," heatmaps was skipped: It already exist. ---")

  }



  # Store the last parameter

  oldParamHeatmapOutput <- newParameters


  # Store the parameters for this run

  if(number.of.rows.parameters == 0){

    saveRDS(

      oldParamHeatmapOutput,

      file=bfcnew(bfc, "Parameters for heatmapOutput()", ext="RDS")

      )

  } else if(number.of.rows.parameters == 1){

    saveRDS(

      oldParamHeatmapOutput,

      file=bfc[[bfcquery(bfc, "Parameters for heatmapOutput()")$rid]]

      )

  }



  # change directory to parent directory

  setwd(parent.directory)

}
