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
#' Date: \tab 2017-06-01 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#' @import cgdsr Biobase
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
#' @param rewrite.output.list a logical value. This option can be used to modify heatmap, for instance change margin,
#' without obtaining data from internet again. If set to false, function will use the values previously stored in the
#' global environment to draw heatmaps and save excel file. The default value is \code{TRUE}.
#'
#' @param round a logical value that, if set to be \code{TRUE}, will force the function to round all the calculated values
#' to two decimal places. The default value is \code{TRUE}.
#'
#' @param top.genes a logical value that, if set as \code{TRUE}, cause the function to create three data.frame that contain the
#' gene names with the highest values for each cancer. To get all the three data.frames, `Frequency.Percentage`, `Mean.Value` and
#' `Median` must have been included for \code{data.presented.as}.
#'
#' more than 5 genes!
#'
#'
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



#########################################################################
#########################################################################
########### Automatically calculate statistical measurements ############
#########################################################################
#########################################################################

automatedStatistics<- function(obtainedDataType = "multiple studies", submissionName, desiredTechnique, calculate = c("frequencyPercentage", "frequencyRatio", "meanValue", "medianValue"),

                              cutoff=NULL, round=TRUE, topGenes = TRUE, validateGenes = TRUE){

  ##########################################################################
  ########## Prerequisites

  # Obtain the unprocessed data list

  if(obtainedDataType == "multiple studies"){

    databaseSymbol <- "obM"

    validationSymbol <- "vaM"

    haultType <- "ohaultM"

  } else if(obtainedDataType == "single study"){

    databaseSymbol <- "obS"

    validationSymbol <- "vaS"

    haultType <- "ohaultS"

  } else{

    stop("'obtainedDataType' must be entered as either 'multiple studies' or 'single study'.")

  }



  # Check submissionName

  if(!is.character(submissionName)){

    stop("'submissionName' must be entered as a character string for naming the process")

  }



  # high-throughput data type

  if(is.character(desiredTechnique)){

    if(!(desiredTechnique %in% c("RNA-seq", "microRNA-Seq", "Microarray.mRNA", "Microarray.microRNA", "methylation")) | length(desiredTechnique)!= 1){

      stop("'desiredTechnique' must contain one of the following techniques: 'RNA-seq', 'microRNA-Seq', 'Microarray.mRNA', 'Microarray.microRNA' or

           'methylation'")

    }

  } else {

    stop("'desiredTechnique' must be entered as a character string describing a technique name")

  }



  # setting the calue for cutoff

  if(desiredTechnique == "methylation"){

    cutoff.phrase <- "obs/exp cutoff"

    if(is.null(cutoff)){

      cutoff <- 0.6

    }

  } else{

    cutoff.phrase <- "zscore cutoff"

    if(is.null(cutoff)){

      cutoff <- 2

    }

  }



  # Check hault older

  if(exists(paste(haultType, ".", submissionName, sep = ""))){

    if(!is.null(paste(haultType, ".", submissionName, sep = ""))){



      if(exists(paste("cut", ".", submissionName, sep = ""))){

        oldCutoff <- get(paste(paste("cut", ".", submissionName, sep = "")))

        if(identical(oldCutoff, cutoff)){

          haultDecision1 <- NULL

        }

      }



      if(exists(paste("cal", ".", submissionName, sep = ""))){

        oldCalculate <- sort(get(paste(paste("cal", ".", submissionName, sep = ""))))

        newCalculate <- sort(calculate)


        if(identical(oldCalculate, newCalculate)){

          haultDecision2 <- NULL

        }

      }



      # Halt the function

      if(all(exists("haultDecision1"), exists("haultDecision2"))){

        return("--- Function 'automatedStatistics()' was skipped: The requested statistical parameters are already calculated ---")

      }

    }

  }





  # Check calculate

  if(is.character(calculate)){

    if(!any(calculate %in% c("frequencyPercentage", "frequencyRatio", "meanValue", "medianValue"))){

      stop("'calculate' must contain at least one of the following: 'frequencyPercentage', 'frequencyRatio', 'meanValue' and 'medianValue'.")

    }

  }



  # Getting the source data

  sourceDataList <- get(paste(databaseSymbol, ".", submissionName, sep = ""))

  if(!is.list(sourceDataList)){

    stop(paste("Input database must be a list.", sep = ""))

  }



  ##########################################################################
  ########## Set the function ready to work

  # creating output fortmat

  processedList <- list()

  options(stringsAsFactors = FALSE)



  # Create a progressbar

  automatedStatisticsProgressBar <- txtProgressBar(min = 0, max = length(sourceDataList[[1]])*length(sourceDataList), style = 3)



  ##########################################################################
  ########## Core segment

  # calculating the first 'for' loop for different gene groups

  for(gg in 1:length(sourceDataList)){

    geneNumber <- ncol(sourceDataList[[gg]][[1]])

    temList <- list()



    for(cs in 1:length(sourceDataList[[1]])){

      # start working on one study

      source.data.subset <- sourceDataList[[gg]][[cs]]

      source.data.subset.name <- names(sourceDataList[[gg]])[cs]

      genes.involved <- colnames(sourceDataList[[gg]][[cs]])





      # Creating and filling the empty matrix with frequency.percentage data

      if("frequencyPercentage" %in% calculate){

        # creating empty matrix

        frequency.percentage.for.a.subset <- matrix(, nrow = 1, ncol = geneNumber)

        dimnames(frequency.percentage.for.a.subset) <- list(source.data.subset.name, genes.involved)



        # calculate frequency percentage

          for(fp in 1:geneNumber){

            # Check all members are under cutoff

            if(!is.na(mean(as.vector(source.data.subset[,fp])[abs(source.data.subset[,fp])])) &

               is.nan(mean(as.vector(source.data.subset[,fp])[abs(source.data.subset[,fp]) >= cutoff], na.rm=TRUE))){

              frequency.percentage.for.a.subset[1, fp] <- 0

              # Check all members are NaN

            } else if (length((source.data.subset[,fp])[!is.nan(source.data.subset[,fp])]) == 0 & all(!is.finite(source.data.subset[,fp])) &

                       is.nan(mean(as.vector(source.data.subset[,fp])[abs(source.data.subset[,fp]) >= cutoff], na.rm=TRUE))){

              frequency.percentage.for.a.subset[1, fp] <- NaN

              # Check all members are NA

            } else if (length((source.data.subset[,fp])[!is.nan(source.data.subset[,fp])]) > 0 & all(!is.finite(source.data.subset[,fp])) &

                       is.nan(mean(as.vector(source.data.subset[,fp])[abs(source.data.subset[,fp]) >= cutoff], na.rm=TRUE))){

              frequency.percentage.for.a.subset[1, fp] <- NA

              # Mean is bigger than 0

            } else if (mean(as.vector(source.data.subset[,fp])[abs(source.data.subset[,fp]) >= cutoff], na.rm=TRUE) > 0 &

                       !is.nan(mean(as.vector(source.data.subset[,fp])[abs(source.data.subset[,fp]) >= cutoff], na.rm=TRUE))){

              frequency.percentage.for.a.subset[1, fp] <- 100*mean(as.vector(abs(source.data.subset[,fp]) >= cutoff), na.rm=TRUE)

              # Mean is smaller than 0

            } else if (mean(as.vector(source.data.subset[,fp])[abs(source.data.subset[,fp]) >= cutoff], na.rm=TRUE) < 0 &

                       !is.nan(mean(as.vector(source.data.subset[,fp])[abs(source.data.subset[,fp]) >= cutoff], na.rm=TRUE))){

              frequency.percentage.for.a.subset[1, fp] <- -100*(mean(as.vector(abs(source.data.subset[,fp]) >= cutoff), na.rm=TRUE))

            }

          }



          # Merging calculations

          if(cs==1){

            if(round==TRUE){

              temList$Frequency.Percentage <- round(frequency.percentage.for.a.subset, digits = 2)

            }else{

              temList$Frequency.Percentage <- frequency.percentage.for.a.subset

            }

          } else if(cs!=1){

            if(round==TRUE){

              temList$Frequency.Percentage <- rbind(temList$Frequency.Percentage, round(frequency.percentage.for.a.subset, digits = 2))

            }else{

              temList$Frequency.Percentage <- rbind(temList$Frequency.Percentage, frequency.percentage.for.a.subset)

            }

          }


          if(topGenes == TRUE){

            # Check if manual naming is requested

            pre.topGenes <- frequency.percentage.for.a.subset

            # Removing NaN and NA

            pre.topGenes[is.nan(pre.topGenes) | is.na(pre.topGenes)] <- 0

            # Creating data.frame for the top 5 genes

            post.topGenes <- data.frame(study = source.data.subset.name)

            # Correcting column name for the first column

            if(obtainedDataType == "multiple studies"){

              colnames(post.topGenes)[1] <- "Cancer Study"

            } else if(obtainedDataType == "single study"){

              colnames(post.topGenes)[1] <- "Cancer Study Subgroup"

            }

            # Finding the top 5 values

            topGenes.values <- head(unique(sort(pre.topGenes, decreasing = TRUE)), n = 5)

            for(topV in 1:length(topGenes.values)){

              topGene.name <- colnames(pre.topGenes)[pre.topGenes %in% topGenes.values[topV]]

              # check whether ttwo or more genes have the same rank

              if(length(topGene.name) > 1){

                topGene.name <- paste(topGene.name, collapse = ", ")

              }

              # rounding

              if(round==TRUE){

                complete.top <- data.frame(topGene = topGene.name, topValue = round(topGenes.values[topV], digits = 2))

              } else{

                complete.top <- data.frame(topGene = topGene.name, topValue = topGenes.values[topV])

              }

              # correcting column names

              colnames(complete.top) <- c(paste(topV, "th ", "Gene", sep=""), paste(topV, "th ", "Value", sep=""))

              # complete dataframe

              post.topGenes <- cbind(post.topGenes, complete.top)

            }


            # fixing the problem caused by more thank one gene with same rank

            if(length(topGenes.values) < 5){

              # Repeat unit

              fix.dataframe <- data.frame(topGene = "-", topValue = "-")

              # number of new units

              newUnits <- 5 - length(topGenes.values)

              # finding current number of units

              oldUnits <- length(topGenes.values)

              for(empty in 1:newUnits){

                colnames(fix.dataframe) <- c(paste(oldUnits + empty, "th ", "Gene", sep=""), paste(oldUnits + empty, "th ", "Value", sep=""))

                post.topGenes <- cbind(post.topGenes, fix.dataframe)

              }

            }

            # assigning the value to the second level list

            if(cs == 1){

              temList$Top.Genes.of.Frequency.Percentage <- post.topGenes

            } else if(cs!=1){

              temList$Top.Genes.of.Frequency.Percentage <- rbind(temList$Top.Genes.of.Frequency.Percentage, post.topGenes)

            }

          }

      }










      # Creating and filling the empty matrix with frequency.ratio data

      if("frequencyRatio" %in% calculate){

        # creating empty matrix

        frequency.ratio.for.a.subset <- matrix(, nrow = 1, ncol = geneNumber)

        dimnames(frequency.ratio.for.a.subset) <- list(source.data.subset.name, genes.involved)



        # calculate frequency ratio

        for(fr in 1:geneNumber){

          # Check all members are under cutoff

          if(!is.na(mean(as.vector(source.data.subset[,fr])[abs(source.data.subset[,fr])])) &

             is.nan(mean(as.vector(source.data.subset[,fr])[abs(source.data.subset[,fr]) >= cutoff], na.rm=TRUE))){

            frequency.ratio.for.a.subset[1, fr] <- paste(as.character(0), " out of ", as.character(length(as.vector(source.data.subset[,fr]))), sep="")

            # Check all members are NaN

          } else if (length((source.data.subset[,fr])[!is.nan(source.data.subset[,fr])]) == 0 & all(!is.finite(source.data.subset[,fr])) &

                     is.nan(mean(as.vector(source.data.subset[,fr])[abs(source.data.subset[,fr]) >= cutoff], na.rm=TRUE))){

            frequency.ratio.for.a.subset[1, fr] <- NaN

            # Check all members are NA

          } else if (length((source.data.subset[,fr])[!is.nan(source.data.subset[,fr])]) > 0 & all(!is.finite(source.data.subset[,fr])) &

                     is.nan(mean(as.vector(source.data.subset[,fr])[abs(source.data.subset[,fr]) >= cutoff], na.rm=TRUE))){

            frequency.ratio.for.a.subset[1, fr] <- NA

            # Mean is number

          } else if (!is.nan(mean(as.vector(source.data.subset[,fr])[abs(source.data.subset[,fr]) >= cutoff], na.rm=TRUE))){

            frequency.ratio.for.a.subset[1, fr] <- paste(as.character(length(na.omit(as.vector(source.data.subset[,fr])[abs(source.data.subset[,fr]) >= cutoff]))),

                                                                " out of ", as.character(length(as.vector(source.data.subset[,fr]))), sep="")

          }

        }



        # Merging calculations

        if(cs==1){

            temList$Frequency.Ratio <- frequency.ratio.for.a.subset


        } else if(cs!=1){

            temList$Frequency.Ratio <- rbind(temList$Frequency.Ratio, frequency.ratio.for.a.subset)

        }

      }










      # Creating and filling the empty matrix with mean.value data

      if("meanValue" %in% calculate){

        # creating empty matrix

        mean.value.for.a.subset <- matrix(, nrow = 1, ncol = geneNumber)

        dimnames(mean.value.for.a.subset) <- list(source.data.subset.name, genes.involved)



        # calculate frequency percentage

        for(mv in 1:geneNumber){

          # Check all members are under cutoff

          if(!is.na(mean(as.vector(source.data.subset[,mv])[abs(source.data.subset[,mv])])) &

             is.nan(mean(as.vector(source.data.subset[,mv])[abs(source.data.subset[,mv]) >= cutoff], na.rm=TRUE))){

            mean.value.for.a.subset[1, mv] <- 0

            # Check all members are NaN

          } else if (length((source.data.subset[,mv])[!is.nan(source.data.subset[,mv])]) == 0 & all(!is.finite(source.data.subset[,mv])) &

                     is.nan(mean(as.vector(source.data.subset[,mv])[abs(source.data.subset[,mv]) >= cutoff], na.rm=TRUE))){

            mean.value.for.a.subset[1, mv] <- NaN

            # Check all members are NA

          } else if (length((source.data.subset[,mv])[!is.nan(source.data.subset[,mv])]) > 0 & all(!is.finite(source.data.subset[,mv])) &

                     is.nan(mean(as.vector(source.data.subset[,mv])[abs(source.data.subset[,mv]) >= cutoff], na.rm=TRUE))){

            mean.value.for.a.subset[1, mv] <- NA

            # Mean is number

          } else if (!is.nan(mean(as.vector(source.data.subset[,mv])[abs(source.data.subset[,mv]) >= cutoff], na.rm=TRUE))){

            mean.value.for.a.subset[1, mv] <- mean(as.vector(source.data.subset[,mv])[abs(source.data.subset[,mv]) >= cutoff], na.rm=TRUE)

          }

        }



        # Merging calculations

        if(cs==1){

          if(round==TRUE){

            temList$Mean.Value <- round(mean.value.for.a.subset, digits = 2)

          }else{

            temList$Mean.Value <- mean.value.for.a.subset

          }

        } else if(cs!=1){

          if(round==TRUE){

            temList$Mean.Value <- rbind(temList$Mean.Value, round(mean.value.for.a.subset, digits = 2))

          }else{

            temList$Mean.Value <- rbind(temList$Mean.Value, mean.value.for.a.subset)

          }

        }


        if(topGenes == TRUE){

          # Check if manual naming is requested

          pre.topGenes <- mean.value.for.a.subset

          # Removing NaN and NA

          pre.topGenes[is.nan(pre.topGenes) | is.na(pre.topGenes)] <- 0

          # Creating data.frame for the top 5 genes

          post.topGenes <- data.frame(study = source.data.subset.name)

          # Correcting column name for the first column

          if(obtainedDataType == "multiple studies"){

            colnames(post.topGenes)[1] <- "Cancer Study"

          } else if(obtainedDataType == "single study"){

            colnames(post.topGenes)[1] <- "Cancer Study Subgroup"

          }

          # Finding the top 5 values

          topGenes.values <- head(unique(sort(pre.topGenes, decreasing = TRUE)), n = 5)

          for(topV in 1:length(topGenes.values)){

            topGene.name <- colnames(pre.topGenes)[pre.topGenes %in% topGenes.values[topV]]

            # check whether ttwo or more genes have the same rank

            if(length(topGene.name) > 1){

              topGene.name <- paste(topGene.name, collapse = ", ")

            }

            # rounding

            if(round==TRUE){

              complete.top <- data.frame(topGene = topGene.name, topValue = round(topGenes.values[topV], digits = 2))

            } else{

              complete.top <- data.frame(topGene = topGene.name, topValue = topGenes.values[topV])

            }

            # correcting column names

            colnames(complete.top) <- c(paste(topV, "th ", "Gene", sep=""), paste(topV, "th ", "Value", sep=""))

            # complete dataframe

            post.topGenes <- cbind(post.topGenes, complete.top)

          }


          # fixing the problem caused by more thank one gene with same rank

          if(length(topGenes.values) < 5){

            # Repeat unit

            fix.dataframe <- data.frame(topGene = "-", topValue = "-")

            # number of new units

            newUnits <- 5 - length(topGenes.values)

            # finding current number of units

            oldUnits <- length(topGenes.values)

            for(empty in 1:newUnits){

              colnames(fix.dataframe) <- c(paste(oldUnits + empty, "th ", "Gene", sep=""), paste(oldUnits + empty, "th ", "Value", sep=""))

              post.topGenes <- cbind(post.topGenes, fix.dataframe)

            }

          }

          # assigning the value to the second level list

          if(cs == 1){

            temList$Top.Genes.of.Mean.Value <- post.topGenes

          } else if(cs!=1){

            temList$Top.Genes.of.Mean.Value <- rbind(temList$Top.Genes.of.Mean.Value, post.topGenes)

          }

        }

      }










      # Creating and filling the empty matrix with median.value data

      if("medianValue" %in% calculate){

        # creating empty matrix

        median.value.for.a.subset <- matrix(, nrow = 1, ncol = geneNumber)

        dimnames(median.value.for.a.subset) <- list(source.data.subset.name, genes.involved)



        # calculate frequency percentage

        for(mdv in 1:geneNumber){

          # Check all members are under cutoff

          if(!is.na(mean(as.vector(source.data.subset[,mdv])[abs(source.data.subset[,mdv])])) &

             is.nan(mean(as.vector(source.data.subset[,mdv])[abs(source.data.subset[,mdv]) >= cutoff], na.rm=TRUE))){

            median.value.for.a.subset[1, mdv] <- 0

            # Check all members are NaN

          } else if (length((source.data.subset[,mdv])[!is.nan(source.data.subset[,mdv])]) == 0 & all(!is.finite(source.data.subset[,mdv])) &

                     is.nan(mean(as.vector(source.data.subset[,mdv])[abs(source.data.subset[,mdv]) >= cutoff], na.rm=TRUE))){

            median.value.for.a.subset[1, mdv] <- NaN

            # Check all members are NA

          } else if (length((source.data.subset[,mdv])[!is.nan(source.data.subset[,mdv])]) > 0 & all(!is.finite(source.data.subset[,mdv])) &

                     is.nan(mean(as.vector(source.data.subset[,mdv])[abs(source.data.subset[,mdv]) >= cutoff], na.rm=TRUE))){

            median.value.for.a.subset[1, mdv] <- NA

            # Mean is number

          } else if (!is.nan(mean(as.vector(source.data.subset[,mdv])[abs(source.data.subset[,mdv]) >= cutoff], na.rm=TRUE))){

            median.value.for.a.subset[1, mdv] <- median(as.vector(source.data.subset[,mdv])[abs(source.data.subset[,mdv]) >= cutoff], na.rm=TRUE)

          }

        }



        # Merging calculations

        if(cs==1){

          if(round==TRUE){

            temList$Median.Value <- round(median.value.for.a.subset, digits = 2)

          }else{

            temList$Median.Value <- median.value.for.a.subset

          }

        } else if(cs!=1){

          if(round==TRUE){

            temList$Median.Value <- rbind(temList$Median.Value, round(median.value.for.a.subset, digits = 2))

          }else{

            temList$Median.Value <- rbind(temList$Median.Value, median.value.for.a.subset)

          }

        }


        if(topGenes == TRUE){

          # Check if manual naming is requested

          pre.topGenes <- median.value.for.a.subset

          # Removing NaN and NA

          pre.topGenes[is.nan(pre.topGenes) | is.na(pre.topGenes)] <- 0

          # Creating data.frame for the top 5 genes

          post.topGenes <- data.frame(study = source.data.subset.name)

          # Correcting column name for the first column

          if(obtainedDataType == "multiple studies"){

            colnames(post.topGenes)[1] <- "Cancer Study"

          } else if(obtainedDataType == "single study"){

            colnames(post.topGenes)[1] <- "Cancer Study Subgroup"

          }

          # Finding the top 5 values

          topGenes.values <- head(unique(sort(pre.topGenes, decreasing = TRUE)), n = 5)

          for(topV in 1:length(topGenes.values)){

            topGene.name <- colnames(pre.topGenes)[pre.topGenes %in% topGenes.values[topV]]

            # check whether ttwo or more genes have the same rank

            if(length(topGene.name) > 1){

              topGene.name <- paste(topGene.name, collapse = ", ")

            }

            # rounding

            if(round==TRUE){

              complete.top <- data.frame(topGene = topGene.name, topValue = round(topGenes.values[topV], digits = 2))

            } else{

              complete.top <- data.frame(topGene = topGene.name, topValue = topGenes.values[topV])

            }

            # correcting column names

            colnames(complete.top) <- c(paste(topV, "th ", "Gene", sep=""), paste(topV, "th ", "Value", sep=""))

            # complete dataframe

            post.topGenes <- cbind(post.topGenes, complete.top)

          }


          # fixing the problem caused by more thank one gene with same rank

          if(length(topGenes.values) < 5){

            # Repeat unit

            fix.dataframe <- data.frame(topGene = "-", topValue = "-")

            # number of new units

            newUnits <- 5 - length(topGenes.values)

            # finding current number of units

            oldUnits <- length(topGenes.values)

            for(empty in 1:newUnits){

              colnames(fix.dataframe) <- c(paste(oldUnits + empty, "th ", "Gene", sep=""), paste(oldUnits + empty, "th ", "Value", sep=""))

              post.topGenes <- cbind(post.topGenes, fix.dataframe)

            }

          }

          # assigning the value to the second level list

          if(cs == 1){

            temList$Top.Genes.of.Median.Value <- post.topGenes

          } else if(cs!=1){

            temList$Top.Genes.of.Median.Value <- rbind(temList$Top.Genes.of.Median.Value, post.topGenes)

          }

        }

      }

      # Update progressbar

      setTxtProgressBar(automatedStatisticsProgressBar, gg*cs)

    }

    # assign the statistics list fot a subgroup of processedList

    processedList[[gg]] <- temList

    names(processedList)[gg] <- names(sourceDataList)[gg]

  }

  # close progressbar

  close(automatedStatisticsProgressBar)

  # Export the obtained data as list

  assign(paste("PrData", ":", submissionName, sep = ""), processedList, envir = globalenv())

  # store cutoff

  assign(paste(paste("cut", ".", submissionName, sep = "")), cutoff, envir = globalenv())

  # store calculate

  assign(paste(paste("cal", ".", submissionName, sep = "")), calculate, envir = globalenv())

}
