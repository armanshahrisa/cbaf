#' @title Checking Expression/methylation Profile for various subgroups of a cancer study.
#'
#' @description This function Obtaines the requested data for the given genes across multiple subgroups of a cancer. It can
#' check whether or not all genes are included in subgroups of a cancer study and, if not, looks for the alternative gene names.
#' Tha main part of function calculates frequency percentage, frequency ratio, mean expression and median of samples greather than
#' specific value in the selected subgroups of the cancer. Furthermore, it looks for the genes that comprise the highest values
#' in each cancer subgroup.
#'
#' @details
#' \tabular{lllll}{
#' Package: \tab cBioAutomatedTools \cr
#' Type: \tab Package \cr
#' Version: \tab 0.99.0 \cr
#' Date: \tab 2017-05-31 \cr
#' License: \tab Artistic-2.0 \cr
#' }
#'
#' @usage process.one.study(genes, cancername, high.throughput.data.type,
#' data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value", "Median"),
#' transposedHeatmap=FALSE, desired.case.list="None", genelimit="none",
#' resolution=600, RowCex=0.8, ColCex=0.8, heatmapMargines=c(10,10), cutoff="default",
#' angle.for.heatmap.cancernames=45, heatmap.color = "RdBu", reverse.heatmap.color = TRUE,
#' rewrite.output.list = TRUE, round=TRUE, top.genes = TRUE, validate.genes = TRUE,
#' Use.CancerCode.as.Name = FALSE, simplify.visulization=FALSE, simplifiction.cuttoff=FALSE)
#'
#' @param genes a list that contains at least one gene group
#'
#' @param cancername a character string showing the desired cancer name. It is an standard cancer name that can be found on
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
#' @param transposedHeatmap a logical value that can exchange the rows and columns of heatmaps.
#'
#' @param desired.case.list a numeric vector that contains the index of rdesired cancer subgroups. If set to be \code{FALSE},
#' function will ask the user to enter them during the process. The default value is \code{FALSE}.
#'
#' @param genelimit if large number of genes exists in at least one gene group, this option can be use to limit the number of
#' genes to be shown on hitmap. For instance, \code{genelimit=50} will limit the heatmap to 50 genes showing the most variation.
#' The default value is \code{FALSE}.
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
#' @param validate.genes a logical value that, if set to be \code{TRUE}, function will checks each cancer study to finds whether
#' or not each gene has a record. If the given cancer doesn't have a record for specific gene, it checks for alternative gene
#' names that cbioportal might use instead of the given gene name.
#'
#' @param Use.CancerCode.as.Name a logical value that tells the function to use abbreviated cancer names instead of complete
#' cancer names, if set to be \code{TRUE}. For example, \code{laml_tcga_pub} is the shortened name for
#' \code{Acute Myeloid Leukemia (TCGA, NEJM 2013)}.
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
#' genes <- list(K.demethylases = c("KDM1A", "KDM1B", "KDM2A"),
#' K.acetyltransferases = c("CLOCK", "CREBBP", "ELP3", "EP300"))
#'
#' cancername <- "Breast Invasive Carcinoma (TCGA, Cell 2015)"
#'
#' process.one.study(genes, cancername, "RNA-seq")
#'
#' process.one.study(genes, cancername, "RNA-seq", desired.case.list = c(3,4,5),
#' data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value"),
#' resolution=300, RowCex=1, ColCex=1, heatmapMargines=c(15,5),
#' cutoff=1.5, angle.for.heatmap.cancernames=30, heatmap.color = "redgreen")
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

process.one.study <- function(genes, cancername, high.throughput.data.type, data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Value", "Median"),

                                         transposedHeatmap=FALSE, desired.case.list="None", genelimit="none", resolution=600, RowCex=0.8, ColCex=0.8,

                                         heatmapMargines=c(10,10), cutoff="default", angle.for.heatmap.cancernames=45, heatmap.color = "RdBu", reverse.heatmap.color = TRUE,

                                         rewrite.output.list = TRUE, round=TRUE, top.genes = TRUE, validate.genes = TRUE, Use.CancerCode.as.Name = FALSE,

                                         simplify.visulization=FALSE, simplifiction.cuttoff=FALSE){


  ##########################################################################
  ### Checking Input data

  # Genes

  if(!is.list(genes)){

    stop("'genes' must be entered as a list containing at list one group of genes with descriptive group name for logistical purposes")

  }

  # cancers

  if(!is.character(cancername)){

    stop("'cancername' must be entered as a character string")

  }

  # high.throughput.data.type

  if(is.character(high.throughput.data.type)){

    if(!(high.throughput.data.type %in% c("RNA-seq", "microRNA-Seq", "Microarray.mRNA", "Microarray.microRNA"))){

      stop("'high.throughput.data.type' must contain one of the following techniques: 'RNA-seq', 'microRNA-Seq', 'Microarray.mRNA' or 'Microarray.microRNA'")

    }

  } else {

    stop("'high.throughput.data.type' must be entered in character string describing a technique name")

  }




  ##########################################################################
  ### Prerequisites

  # Get directory address

  dir.address <- getwd()




  # Choice of high-throughput data type

  if(high.throughput.data.type == "RNA-seq"){

    L1.characteristics <- c("Tumor Samples with mRNA data (RNA Seq V2)", "Tumors with mRNA data (RNA Seq V2)", "Tumor Samples with mRNA data (RNA Seq)", "Tumors with mRNA data (RNA Seq)")

    L2.characteristics <- c("mRNA Expression z-Scores (RNA Seq V2 RSEM)", "mRNA Expression z-Scores (RNA Seq RPKM)")

  } else if(high.throughput.data.type == "microRNA-Seq"){

    L1.characteristics <- c("Tumors with microRNA data (microRNA-Seq)")

    L2.characteristics <- c("microRNA expression Z-scores")

  } else if(high.throughput.data.type == "microarray.mRNA"){

    L1.characteristics <- c("Tumor Samples with mRNA data (Agilent microarray)", "Tumors with mRNA data (Agilent microarray)", "Tumor Samples with mRNA data (U133 microarray only)", "Tumors with mRNA data", "Tumors with mRNA")

    L2.characteristics <- c("mRNA Expression z-Scores (microarray)", "mRNA Expression z-Scores (U133 microarray only)", "mRNA expression z-scores (Illumina)", "mRNA expression Z-scores (all genes)", "mRNA Expression Z-Scores vs Normals", "mRNA Expression z-Scores (combined microarray)")

  } else if(high.throughput.data.type == "microarray.microRNA"){

    L1.characteristics <- c("Tumors with microRNA")

    L2.characteristics <- c("mRNA Expression Z-Scores vs Normals", "mRNA Expression z-Scores (combined microarray)")

  } else if(high.throughput.data.type == "methylation"){

    L1.characteristics <- c("Tumor Samples with methylation data (HM450)", "Tumors with methylation data (HM450)", "Tumor Samples with methylation data (HM27)", "Tumors with methylation data (HM27)", "Tumors with methylation data")

    L2.characteristics <- c("Methylation (HM450)", "Methylation (HM27)", "Methylation")

  } else{

    stop("High-throughput.data field can not be left empety. It should be chosen as 'RNA-seq', 'microRNA-Seq', 'microarray.mRNA', 'microarray.microRNA'or 'methylation'")

  }



  # Creating a vector for cancer names

  if(high.throughput.data.type == "methylation"){

    cutoff.phrase <- "obs/exp cutoff"

    if(is.character(cutoff)){

      cutoff <- 0.6

    }

  } else{

    cutoff.phrase <- "zscore cutoff"

    if(is.character(cutoff)){

      cutoff <- 2

    }

  }




  # Set cgdsr

  mycgds = CGDS("http://www.cbioportal.org/")




  ##########################################################################
  ### Core segment

  for(g in 1:length(names(genes))){

    # Reset directory

    setwd(dir.address)

    # Genes for each group

    c.genes <- as.character(unlist(genes[g]))






    if(rewrite.output.list == TRUE){




      print(paste("[[[", "Obtaining expression data for", names(genes)[g], "genes", "]]]", sep = " "))



      # First step of procedure outside the 'for' control structure as it is repeated in a single cancer

      mycancerstudy = getCancerStudies(mycgds)[which(getCancerStudies(mycgds)[,2]==cancername),1]






      # Finding the second characteristic of data in the cancer

      # The second characteristic is the same in all case list and therefore brought outside

      s.condition <- (getGeneticProfiles(mycgds,mycancerstudy)[,2])[getGeneticProfiles(mycgds,mycancerstudy)[,2] %in% L2.characteristics]

      mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[which(getGeneticProfiles(mycgds,mycancerstudy)[,2] == if(length(s.condition) >= 1){

        s.condition[1]

      } else if (length(s.condition) == 0){

        stop(paste(cancername, "doesn't have an appropriate 'level 2' condition for", high.throughput.data.type, "data!", sep=" "))

      }) ,1]






      # Chosing the desired case lists

      if(is.character(desired.case.list)){

        print(paste("Please enter the numeric index of desired case list(s) for", cancername, "separated by comma. For instance two indexes 1 and 2 should be enterd as 1, 2", sep=" "))

        Choices <- getCaseLists(mycgds,mycancerstudy)[,2]

        print(paste(1:length(Choices), Choices, sep=". "))

        inputCases <- readline(prompt = "Enter the numeric index(es):                  ")

        inputCases <- as.numeric(unlist(strsplit(inputCases, ",")))

      } else {

        inputCases <- desired.case.list

      }






      # Creating a vector which contains names of inputCases

      inputCases.names <- getCaseLists(mycgds,mycancerstudy)[inputCases ,2]

      # Shorten cancer name

      Shortened.cancername <- gsub(" ", ".", sapply(strsplit(as.character(cancername), split=" (", fixed=TRUE), function(x) (x[1])))













      # create progress bar

      pb <- txtProgressBar(min = 0, max = length(inputCases), style = 3)






      ## Getting the required gene expresssion profile ...

      # 'for' control structure for obtaining data and calculating the requested parameters

      for(i in 1:length(inputCases)){

        # Setting the first characteristics of data according to the desired case list

        ind <- getCaseLists(mycgds,mycancerstudy)[inputCases[i] ,2]

        mycaselist = getCaseLists(mycgds,mycancerstudy)[inputCases[i] ,1]







        # Obtaining Expression x-scores fore the requested genes

        ProfileData <- data.matrix(getProfileData(mycgds,c.genes[order(c.genes)],mygeneticprofile,mycaselist))

        # Get the order of genes in obtained expression data

        ordered.Genes <- colnames(ProfileData)

        Expected.Genes <- c.genes[order(c.genes)]






        # Create a matrix contained genes that are not included in required cancer

        if(validate.genes == TRUE){

          # Alter c.genes to be compatible with gene names in cBioPortal output

          Altered.sorted.gene.names <- sort(gsub("-", ".", c.genes))

          # Obtain name of genes that are absent in requested cancer

          Absent.genes <- Altered.sorted.gene.names[!Altered.sorted.gene.names %in% colnames(ProfileData)]

          # For loop for determining changed genes

          Alternative.gene.names <- vector("character", length = length(Absent.genes))

          for(ab in 1:length(Absent.genes)){

            Absent.gene.name.ProfileData <- colnames(data.matrix(getProfileData(mycgds, Absent.genes[ab], mygeneticprofile,mycaselist)))

            # Check wheter gene has an alternative name or missed in the database

            if(length(Absent.gene.name.ProfileData) == 1){

              Alternative.gene.names[ab] <- Absent.gene.name.ProfileData

            } else if(length(Absent.gene.name.ProfileData) != 1){

              Alternative.gene.names[ab] <- "-"

            }

          }

          # Naming Alternative.gene.names

          names(Alternative.gene.names) <- Absent.genes

          # Seperating genes with alternative names from those that are absent

          Gene.with.no.data <- Alternative.gene.names[Alternative.gene.names == "-"]

          Gene.with.data <- Alternative.gene.names[Alternative.gene.names != "-"]




          # Switch for Absent.genes

          if(length(Absent.genes) != 0 & any(Alternative.gene.names == "-")){



            ## Create profileData which containes absent genes

            # Find sample numbers in ProfileData

            numberOfRows <- nrow(ProfileData)

            # Create NA variable to be bind to Profile data with same sample size

            NotFound <- matrix(rep(NA, length(Gene.with.no.data)*numberOfRows), ncol = length(Gene.with.no.data), nrow = numberOfRows)

            # make sample names identical to sample names of ProfileData

            colnames(NotFound) <- names(Gene.with.no.data)

            # Merge two ProfileData

            CompletedProfileData <- cbind(ProfileData, NotFound)

            # Sort the merged ProfileData according to alphabetic name of genes (In order to mix to ProfileData)

            CompletedProfileData <- CompletedProfileData[,sort(colnames(CompletedProfileData))]

          } else if(length(Absent.genes) == 0 | (length(Absent.genes) != 0 & any(Alternative.gene.names != "-"))){

            CompletedProfileData <- ProfileData

          }

        }






        # Creating Empty matrixes to be filled with 'for' control structure

        if(validate.genes == TRUE){

          Complete.ordered.Genes <- colnames(CompletedProfileData)

          Serial.ProfileData.for.Validation <- vector("character", length=length(Expected.Genes))

          names(Serial.ProfileData.for.Validation) <- Complete.ordered.Genes

        }

        if("Frequency.Percentage" %in% data.presented.as){

          Serial.ProfileData.with.Frequency.Percentage <- vector("numeric", length=length(ordered.Genes))

        }

        if("Frequency.Ratio" %in% data.presented.as){

          Serial.ProfileData.with.Frequency.Ratio <- vector("character", length=length(ordered.Genes))

        }

        if("Mean.Value" %in% data.presented.as){

          Serial.ProfileData.with.Mean.Value <- vector("numeric", length=length(ordered.Genes))

        }

        if("Median" %in% data.presented.as){

          Serial.ProfileData.with.Median <- vector("numeric", length=length(ordered.Genes))

        }






        ### Calculating validity, frequency, Mean.Value and Median values for each requested gene

        if(validate.genes == TRUE){

          for(expect in 1:length(Expected.Genes)){

            ## Validating Genes

            # Correct those that are not found

            if(length((CompletedProfileData[,expect])[!is.nan(CompletedProfileData[,expect])]) > 0 & all(!is.finite(CompletedProfileData[,expect])) &

               is.nan(mean(as.vector(CompletedProfileData[,expect])[abs(CompletedProfileData[,expect])], na.rm=TRUE))){

              Serial.ProfileData.for.Validation[expect] <- "-"

            } else {

              Serial.ProfileData.for.Validation[expect] <- "Found"

            }

          }

          # Replacing name and alternative names

          for(re in 1:length(Gene.with.data)){

            replace.index <- which(names(Serial.ProfileData.for.Validation) == Gene.with.data[re])

            Serial.ProfileData.for.Validation[replace.index] <- paste("Alternative name: ", Gene.with.data[re], sep="")

            names(Serial.ProfileData.for.Validation)[replace.index] <- names(Gene.with.data)[re]

          }

          # Sorting genes alphabetically

          Serial.ProfileData.for.Validation <- Serial.ProfileData.for.Validation[sort(names(Serial.ProfileData.for.Validation))]

        }



        for(j in 1:length(ordered.Genes)){

          # Frequency.Percentage

          if("Frequency.Percentage" %in% data.presented.as){

            # Check all members are under cutoff

            if(!is.na(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j])])) &

               is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Frequency.Percentage[j] <- 0

              # Check all members are NaN

            } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) == 0 & all(!is.finite(ProfileData[,j])) &

                       is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Frequency.Percentage[j] <- NaN

              # Check all members are NA

            } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) > 0 & all(!is.finite(ProfileData[,j])) &

                       is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Frequency.Percentage[j] <- NA

              # Mean is bigger than 0

            } else if (mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE) > 0 &

                       !is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Frequency.Percentage[j] <- 100*mean(as.vector(abs(ProfileData[,j]) >= cutoff), na.rm=TRUE)

              # Mean is smaller than 0

            } else if (mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE) < 0 &

                       !is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Frequency.Percentage[j] <- -100*(mean(as.vector(abs(ProfileData[,j]) >= cutoff), na.rm=TRUE))

            }

          }


          # Frequency Ratio

          if("Frequency.Ratio" %in% data.presented.as){

            # Check all members are under cutoff

            if(!is.na(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j])])) &

               is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Frequency.Ratio[j] <- paste(as.character(0), "/", as.character(length(as.vector(ProfileData[,j]))), sep="")

              # Check all members are NaN

            } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) == 0 & all(!is.finite(ProfileData[,j])) &

                       is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Frequency.Ratio[j] <- NaN

              # Check all members are NA

            } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) > 0 & all(!is.finite(ProfileData[,j])) &

                       is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Frequency.Ratio[j] <- NA

              # Mean is number

            } else if (!is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Frequency.Ratio[j] <- paste(as.character(length(na.omit(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff]))),

                                                                  "/", as.character(length(as.vector(ProfileData[,j]))), sep="")

            }

          }


          # Mean.Value

          if("Mean.Value" %in% data.presented.as){

            # Check all members are under cutoff

            if(!is.na(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j])])) &

               is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Mean.Value[j] <- 0

              # Check all members are NaN

            } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) == 0 & all(!is.finite(ProfileData[,j])) &

                       is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Mean.Value[j] <- NaN

              # Check all members are NA

            } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) > 0 & all(!is.finite(ProfileData[,j])) &

                       is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Mean.Value[j] <- NA

              # Mean is number

            } else if (!is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Mean.Value[j] <- mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE)

            }

          }


          # Median

          if("Median" %in% data.presented.as){

            # Check all members are under cutoff

            if(!is.na(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j])])) &

               is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Median[j] <- 0

              # Check all members are NaN

            } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) == 0 & all(!is.finite(ProfileData[,j])) &

                       is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Median[j] <- NaN

              # Check all members are NA

            } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) > 0 & all(!is.finite(ProfileData[,j])) &

                       is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Median[j] <- NA

              # Mean is number

            } else if (!is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Median[j] <- median(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE)

            }

          }

        }


        # Converting each ProfileData to independant cartical matrix with one column

        if(validate.genes == TRUE){

          Original.validate.gene.names <- names(Serial.ProfileData.for.Validation)

          dim(Serial.ProfileData.for.Validation) <- c(length(Expected.Genes),1)

          rownames(Serial.ProfileData.for.Validation) <- Original.validate.gene.names

        }

        if("Frequency.Percentage" %in% data.presented.as){

          dim(Serial.ProfileData.with.Frequency.Percentage) <- c(length(ordered.Genes),1)

        }

        if("Frequency.Ratio" %in% data.presented.as){

          dim(Serial.ProfileData.with.Frequency.Ratio) <- c(length(ordered.Genes),1)

        }

        if("Mean.Value" %in% data.presented.as){

          dim(Serial.ProfileData.with.Mean.Value) <- c(length(ordered.Genes),1)

        }

        if("Median" %in% data.presented.as){

          dim(Serial.ProfileData.with.Median) <- c(length(ordered.Genes),1)

        }




        # Combining serries of columns to obtain four matrixes presenting 'Frequency.Percentage', 'Frequency in ratio', 'Mean.Value', 'Median' in a specific cancer

        if(i==1){

          if(validate.genes == TRUE){

            Validation.matrix <- Serial.ProfileData.for.Validation

          }

          if("Frequency.Percentage" %in% data.presented.as){

            Frequency.Percentage.Matrix <- Serial.ProfileData.with.Frequency.Percentage

          }

          if("Frequency.Ratio" %in% data.presented.as){

            Frequency.Ratio.Matrix <- Serial.ProfileData.with.Frequency.Ratio

          }

          if("Mean.Value" %in% data.presented.as){

            Mean.Value.Matrix <- Serial.ProfileData.with.Mean.Value

          }

          if("Median" %in% data.presented.as){

            Median.Matrix <- Serial.ProfileData.with.Median

          }

        }else if(i!=1){

          if(validate.genes == TRUE){

            Validation.matrix <- cbind(Validation.matrix, Serial.ProfileData.for.Validation)

          }

          if("Frequency.Percentage" %in% data.presented.as){

            Frequency.Percentage.Matrix <- cbind(Frequency.Percentage.Matrix, Serial.ProfileData.with.Frequency.Percentage)

          }

          if("Frequency.Ratio" %in% data.presented.as){

            Frequency.Ratio.Matrix <- cbind(Frequency.Ratio.Matrix, Serial.ProfileData.with.Frequency.Ratio)

          }

          if("Mean.Value" %in% data.presented.as){

            Mean.Value.Matrix <- cbind(Mean.Value.Matrix, Serial.ProfileData.with.Mean.Value)

          }

          if("Median" %in% data.presented.as){

            Median.Matrix <- cbind(Median.Matrix, Serial.ProfileData.with.Median)

          }

        }






        # update progress bar

        setTxtProgressBar(pb, i)

      }

      # Close progress bar

      close(pb)






      ### Adding dimention names to four matrixes

      print("Preparing a 'list' to fill with the requested data")






      ## Creating the required matrixes and Overal list

      # Create empty list

      Data.list <- list()

      # filling the list

      if(validate.genes == TRUE){

        colnames(Validation.matrix) <- inputCases.names

        Data.list$Validation.Results <- t(Validation.matrix)

      }

      if("Frequency.Percentage" %in% data.presented.as){

        dimnames(Frequency.Percentage.Matrix) <- list(ordered.Genes, inputCases.names)

        Data.list$Frequency.Percentage <- t(Frequency.Percentage.Matrix)

        if(top.genes==TRUE){

          # Check if manual naming is requested

          Exp.matrix <- t(Frequency.Percentage.Matrix)

          # Removing NaN

          Exp.matrix[is.nan(Exp.matrix)] <- 0

          # Removing NA

          Exp.matrix[is.na(Exp.matrix)] <- 0

          # Removing rows that contain only 0

          Exp.matrix <- Exp.matrix[rowSums(Exp.matrix)!=0,]

          # Empty data.frame to be filled with data

          Frequency.Percentage.top.genes.dataframe <- data.frame(Gene.with.highest.value=vector("character", length = nrow(Exp.matrix)), Frequency.Percentage=vector("numeric", length = nrow(Exp.matrix)), stringsAsFactors = FALSE)

          # Finding genes with highest value

          for (m in 1:nrow(Exp.matrix)){

            Frequency.Percentage.top.genes.dataframe[m,1] <- names(which.max(Exp.matrix[m,]))

          }

          # Finding values corresponding to these genes

          for (n in 1:nrow(Exp.matrix)){

            Frequency.Percentage.top.genes.dataframe[n,2] <- max(Exp.matrix[n,])

          }

          # Creating data.frame for desired profile

          rownames(Frequency.Percentage.top.genes.dataframe) <- rownames(Exp.matrix)

          Data.list$Top.Genes.of.Frequency.Percentage <- Frequency.Percentage.top.genes.dataframe

        }

      }

      if("Frequency.Ratio" %in% data.presented.as){

        dimnames(Frequency.Ratio.Matrix) <- list(ordered.Genes, inputCases.names)

        Data.list$Frequency.Ratio <- t(Frequency.Ratio.Matrix)

      }

      if("Mean.Value" %in% data.presented.as){

        dimnames(Mean.Value.Matrix) <- list(ordered.Genes, inputCases.names)

        Data.list$Mean.Value <- t(Mean.Value.Matrix)

        if(top.genes==TRUE){

          # Check if manual naming is requested

          Exp.matrix <- t(Mean.Value.Matrix)

          # Removing NaN

          Exp.matrix[is.nan(Exp.matrix)] <- 0

          # Removing NA

          Exp.matrix[is.na(Exp.matrix)] <- 0

          # Removing rows that contain only 0

          Exp.matrix <- Exp.matrix[rowSums(Exp.matrix)!=0,]

          # Empty data.frame to be filled with data

          Mean.Value.top.genes.dataframe <- data.frame(Gene.with.highest.value=vector("character", length = nrow(Exp.matrix)), Mean.Value=vector("numeric", length = nrow(Exp.matrix)), stringsAsFactors = FALSE)

          # Finding genes with highest value

          for (m in 1:nrow(Exp.matrix)){

            Mean.Value.top.genes.dataframe[m,1] <- names(which.max(Exp.matrix[m,]))

          }

          # Finding values corresponding to these genes

          for (n in 1:nrow(Exp.matrix)){

            Mean.Value.top.genes.dataframe[n,2] <- max(Exp.matrix[n,])

          }

          # Creating data.frame for desired profile

          rownames(Mean.Value.top.genes.dataframe) <- rownames(Exp.matrix)

          Data.list$Top.Genes.of.Mean.Value <- Mean.Value.top.genes.dataframe

        }

      }

      if("Median" %in% data.presented.as){

        dimnames(Median.Matrix) <- list(ordered.Genes, inputCases.names)

        Data.list$Median <- t(Median.Matrix)

        if(top.genes==TRUE){

          # Check if manual naming is requested

          Exp.matrix <- t(Median.Matrix)

          # Removing NaN

          Exp.matrix[is.nan(Exp.matrix)] <- 0

          # Removing NA

          Exp.matrix[is.na(Exp.matrix)] <- 0

          # Removing rows that contain only 0

          Exp.matrix <- Exp.matrix[rowSums(Exp.matrix)!=0,]

          # Empty data.frame to be filled with data

          Median.top.genes.dataframe <- data.frame(Gene.with.highest.value=vector("character", length = nrow(Exp.matrix)), Median=vector("numeric", length = nrow(Exp.matrix)), stringsAsFactors = FALSE)

          # Finding genes with highest value

          for (m in 1:nrow(Exp.matrix)){

            Median.top.genes.dataframe[m,1] <- names(which.max(Exp.matrix[m,]))

          }

          # Finding values corresponding to these genes

          for (n in 1:nrow(Exp.matrix)){

            Median.top.genes.dataframe[n,2] <- max(Exp.matrix[n,])

          }

          # Creating data.frame for desired profile

          rownames(Median.top.genes.dataframe) <- rownames(Exp.matrix)

          Data.list$Top.Genes.of.Median <- Median.top.genes.dataframe

        }

      }



      # Saving data list to enviroment with an informative name

      assign(paste(names(genes)[[g]], ".", "Data.List",  "_", "for", "_", Shortened.cancername, sep=""), Data.list, envir = globalenv())

    }












    ### Plotting Heatmap and preparing excel files

    print("Saving the output files")

    Shortened.cancername <- gsub(" ", ".", sapply(strsplit(as.character(cancername), split=" (", fixed=TRUE), function(x) (x[1])))

    List.to.Go <- get(paste(names(genes)[[g]], ".", "Data.List",  "_", "for", "_", Shortened.cancername, sep=""))

    # Create the appropriate directory and set it as desired folder

    newdir <- paste(g, ". ", sub(x = names(genes)[g], pattern = "\\.", replacement = "-"), sep="")

    dir.create(paste(dir.address, newdir, sep = "/"), showWarnings = FALSE)

    setwd(paste(dir.address, newdir, sep = "/"))

    # Remove any preexisting files in folder

    file.remove(list.files())

    # Export data

    if(simplify.visulization == TRUE & is.numeric(simplifiction.cuttoff) & !is.numeric(genelimit)){

      print("Only significant results will be used to draw heatmap")

    }

    # Create progress bar

    pb2 <- txtProgressBar(min = 0, max = length(names(List.to.Go)), style = 3)

    ## 'for' to prepare ouputs

    for(segment in 1:length(names(List.to.Go))){

      ## Subsetting matrixes from List by 'for' control structure

      Temporary.source <- List.to.Go[[segment]]

      Temporary.source.name <- names(List.to.Go)[[segment]]



      ## Creating Heatmap

      # Transpose Temporary.source

      heatmap.data <- t(Temporary.source)

      # Removing NA

      heatmap.data <- heatmap.data[(apply(heatmap.data, 1, function(x) any(!is.na(x)==TRUE))),]

      # Removing NaN

      heatmap.data[is.nan(heatmap.data)] <- 0

      # Removing rows that contain only 0

      if(Temporary.source.name %in% c("Frequency.Percentage","Mean.Value" ,"Median") & is.matrix(heatmap.data)){

        heatmap.data <- heatmap.data[rowSums(heatmap.data, na.rm = TRUE)!=0,]

      }

      if(Temporary.source.name %in% c("Frequency.Percentage","Mean.Value" ,"Median") & is.matrix(heatmap.data)){

        if(nrow(heatmap.data) > 1 & ncol(heatmap.data) > 1){




          # !!! Simplifying heatmap for easy assessment !!!

          if(simplify.visulization == TRUE & is.numeric(simplifiction.cuttoff) & !is.numeric(genelimit)){

            heatmap.data[heatmap.data < simplifiction.cuttoff] <- 0

          } else if(simplify.visulization == TRUE & !is.numeric(simplifiction.cuttoff)){

            stop("Please set your desired cut-off for simplification")

          } else if(simplify.visulization == TRUE & is.numeric(simplifiction.cuttoff) & is.numeric(genelimit)){

            warning("There is no need to limit gene number because simplification option possibly limits gene number even below the specified number")

          }




          # Limiting the number of genes in heatmap to get better resolution

          if(genelimit=="none"){

            heatmap.data <- heatmap.data

          } else if(is.numeric(genelimit)){

            ordering <- order(abs(rowVars(heatmap.data)), decreasing=TRUE)

            heatmap.data <- heatmap.data[ordering[1:genelimit],]

          } else{

            stop("Please type gene number limit or if whole genes are desired please type none")
          }



          # Determining that cancer names for heatmap should be obtained from vector or matrix

          tissue <- colnames(heatmap.data)




          # Heatmap color

          if(reverse.heatmap.color == TRUE){

            if(heatmap.color == "redgreen"){

              hmcol <- rev(redgreen(75))

            } else {

              hmcol <- rev(colorRampPalette(brewer.pal(9, heatmap.color))(100))

            }


          } else if (reverse.heatmap.color == FALSE){

            if(heatmap.color == "redgreen"){

              hmcol <- redgreen(75)

            } else {

              hmcol <- colorRampPalette(brewer.pal(9, heatmap.color))(100)

            }

          }




          # Drawing heatmap

          png(filename=paste(getwd(), paste(gsub(x = names(genes)[g], pattern = "\\.", replacement = "-"), ", ", gsub(x = Temporary.source.name, pattern = "\\.", replacement = " "), " (","Z-score", "=", cutoff, ")", ".PNG" , sep=""), sep="/"), width=9.5, height= 11, units = "in", res=resolution)



          if(transposedHeatmap==FALSE){

            heatmap.2(heatmap.data, labCol=tissue, na.color="light gray", trace="none", symbreaks = T, col=hmcol, cexRow = RowCex, cexCol= ColCex,

                      margins = heatmapMargines, srtCol = angle.for.heatmap.cancernames)

          } else if(transposedHeatmap==TRUE){

            heatmap.2(t(heatmap.data), labCol=tissue, na.color="light gray", trace="none", symbreaks = T, col=hmcol, cexRow = RowCex, cexCol= ColCex,

                      margins = heatmapMargines, srtCol = angle.for.heatmap.cancernames)

          }

          dev.off()

        }

      }




      ## Exporting the expression profile

      Exp.matrix2 <- Temporary.source

      # Rounding values

      if(round==TRUE & Temporary.source.name %in% c("Frequency.Percentage","Mean.Value" ,"Median")){

        Exp.matrix2 <- round(Exp.matrix2, digits=2)

      }

      # Copy data to a variable to modify

      Save.matrix <- Exp.matrix2

      if(!(Temporary.source.name %in% c("Top.Genes.of.Frequency.Percentage","Top.Genes.of.Mean.Value" ,"Top.Genes.of.Median"))){

        # Replace 'NaN' values with character string

        Save.matrix[is.nan(Save.matrix)] <- "NaN"

        # Replace 'NA' values with character string

        Save.matrix[is.na(Save.matrix)] <- "NA"

      }

      if(Temporary.source.name %in% c("Top.Genes.of.Frequency.Percentage","Top.Genes.of.Mean.Value" ,"Top.Genes.of.Median")){

        colnames(Save.matrix)[1:2] <- c("Gene with highest value", if(Temporary.source.name == "Top.Genes.of.Frequency.Percentage"){"Frequency Percentage"}

                                        else if(Temporary.source.name == "Top.Genes.of.Mean.Value"){"Mean Expression"}

                                        else if(Temporary.source.name == "Top.Genes.of.Median"){"Median"})

      }


      # Saving the expression profile

      write.xlsx(Save.matrix, file = paste(sub(x = names(genes)[g], pattern = "\\.", replacement = "-"), " (",cutoff.phrase, "=", cutoff, ")" , ".xlsx", sep = ""), sheetName = gsub(x = Temporary.source.name, pattern = "\\.", replacement = " "), append = TRUE)




      # Update progressbar

      setTxtProgressBar(pb2, segment)

    }

    # Close progress bar

    close(pb2)

    # Print information

    if(any(c("Frequency.Percentage", "Mean.Value", "Median") %in% names(List.to.Go))){

      print(paste("Requested excel and PNG files were saved in", getwd(), sep = " "))

    } else if(!any(c("Frequency.Percentage", "Mean.Value", "Median") %in% names(List.to.Go))){

      print(paste("Requested excel files were saved in", getwd(), sep = " "))

    }

  }

  setwd(dir.address)
}


