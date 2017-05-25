#' Checking Expression/methylation Profile for various cancers
#'
#' This function checks all the cancers from cbioportal.org to determine
#' datasets that they contain.
#'
#' @return A matrix that contain all cancers and their available datasets. It is
#' available in the global enviroment (user's workspace). For convenience, an excel
#' file will also be generated in the working directory.
#' @details
#' This function checks all the cancers that are registered in 'cbioportal.org' to
#' examine whether or not they contain RNA-seq, microRNA-seq, microarray(mRNA),
#' microarray(miRNA) and methylation datasets.
#' @usage Dataset.availability()
#' @export



###################################################################################################
###################################################################################################
############### Evaluation of Median, Frequency and ExpressionMean for All Cancers ################
###################################################################################################
###################################################################################################

Expression.Data.For.Cancers <- function(genes, cancers, High.throughput.data, data.presented.as = c("Frequency.Percentage", "Frequency.Ratio", "Mean.Expression", "Median"),

                                        shorteded.cancer.names = TRUE, genelimit="none", resolution=600, RowCex=0.8, ColCex=0.8, heatmapMargines=c(24,17), cutoff=2,

                                        angle.for.heatmap.cancernames=45, heatmap.color = "RdBu", reverse.heatmap.color = TRUE, resetOldExpressionProfile = TRUE, round=TRUE,

                                        top.genes = TRUE, validate.genes = TRUE, Use.CancerCode.as.Name = FALSE, simplify.visulization=FALSE, simplifiction.cuttoff=FALSE){


  ##########################################################################
  ### Checks whether the required packages are installed and installs if not

  # CRAN packages

  list.of.packages <- c("cgdsr", "gplots", "RColorBrewer", "rafalib", "xlsx")

  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

  if(length(new.packages)) install.packages(new.packages)




  # Bioconcuctor packages

  list.of.packages <- c("Biobase", "genefilter")

  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

  if(length(new.packages)){

    source("http://www.bioconductor.org/biocLite.R")

    biocLite(new.packages)




    # Configuring installed packages


    if (Sys.getenv("JAVA_HOME")!="")

      Sys.setenv(JAVA_HOME="")

    library(rJava)

    install.packages("xlsxjars", INSTALL_opts = "--no-multiarch")

  }




  # Notification for installing the required packages in LINUX

  if(!any((installed.packages()) %in% "xlsx")){

    print("If you are using Ubuntu please first install 'r-cran-rjava' and 'r-cran-xml' packages by typing 'sudo apt-get install r-cran-rjava' and 'sudo apt-get install r-cran-xml'")

  }




  ##########################################################################
  ### Checking Input data

  # Genes

  if(!is.list(genes)){

    stop("'genes' must be entered as a list containing at list one group of genes with descriptive group name for logistical purposes")

  }

  # cancers

  if(!is.character(cancers)){

    stop("'cancers' should contain at least one cancer name as a character string")

  }

  # High.throughput.data

  if(is.character(High.throughput.data)){

    if(!(High.throughput.data %in% c("RNA-seq", "microRNA-Seq", "Microarray.mRNA", "Microarray.microRNA"))){

      stop("'High.throughput.data' must contain one of the following techniques: 'RNA-seq', 'microRNA-Seq', 'Microarray.mRNA' or 'Microarray.microRNA'")

    }

  } else {

    stop("'High.throughput.data' must be entered in character string describing a technique name")

  }




  ##########################################################################
  ### Prerequisites

  # Get directory address

  dir.address <- getwd()




  # Choice of high-throughput data type

  if(High.throughput.data == "RNA-seq"){

    L1.characteristics <- c("Tumor Samples with mRNA data (RNA Seq V2)", "Tumors with mRNA data (RNA Seq V2)", "Tumor Samples with mRNA data (RNA Seq)", "Tumors with mRNA data (RNA Seq)")

    L2.characteristics <- c("mRNA Expression z-Scores (RNA Seq V2 RSEM)", "mRNA Expression z-Scores (RNA Seq RPKM)")

  } else if(High.throughput.data == "microRNA-Seq"){

    L1.characteristics <- c("Tumors with microRNA data (microRNA-Seq)")

    L2.characteristics <- c("microRNA expression Z-scores")

  } else if(High.throughput.data == "microarray.mRNA"){

    L1.characteristics <- c("Tumor Samples with mRNA data (Agilent microarray)", "Tumors with mRNA data (Agilent microarray)", "Tumor Samples with mRNA data (U133 microarray only)", "Tumors with mRNA data", "Tumors with mRNA")

    L2.characteristics <- c("mRNA Expression z-Scores (microarray)", "mRNA Expression z-Scores (U133 microarray only)", "mRNA expression z-scores (Illumina)", "mRNA expression Z-scores (all genes)", "mRNA Expression Z-Scores vs Normals", "mRNA Expression z-Scores (combined microarray)")

  } else if(High.throughput.data == "microarray.microRNA"){

    L1.characteristics <- c("Tumors with microRNA")

    L2.characteristics <- c("mRNA Expression Z-Scores vs Normals", "mRNA Expression z-Scores (combined microarray)")

  } else if(High.throughput.data == "methylation"){

    L1.characteristics <- c("Tumor Samples with methylation data (HM450)", "Tumors with methylation data (HM450)", "Tumor Samples with methylation data (HM27)", "Tumors with methylation data (HM27)", "Tumors with methylation data")

    L2.characteristics <- c("Methylation (HM450)", "Methylation (HM27)", "Methylation")

  } else{

    stop("High-throughput.data field can not be left empety. It should be chosen as 'RNA-seq', 'microRNA-Seq', 'microarray.mRNA', 'microarray.microRNA'or 'methylation'")

  }



  # Creating a vector for cancer names

  if(High.throughput.data == "methylation"){

    cutoff.phrase <- "obs/exp cutoff"

  } else{

    cutoff.phrase <- "zscore cutoff"

  }



  # Creating a vector for cancer names

  if(!is.matrix(cancers)){

    cancers <- cancers[order(cancers)]

    cancers.matrix <- cancers

  } else if(is.matrix(cancers)){

    cancers.matrix <- cancers

    cancers <- (cancers[,1])[order(cancers[,1])]

  }




  # Load the Required Packages

  library("cgdsr")

  mycgds = CGDS("http://www.cbioportal.org/")

  library(Biobase)

  library(abind)

  library(gplots)

  library(RColorBrewer)

  library(rafalib)

  library(genefilter)

  library(xlsx)




  ##########################################################################
  ### Core segment

  for(g in 1:length(names(genes))){

    # Reset directory

    setwd(dir.address)

    # Genes for each group

    c.genes <- as.character(unlist(genes[g]))




    if(resetOldExpressionProfile == TRUE){




      print(paste("[[[", "Obtaining expression data for", names(genes)[g], "genes", "]]]", sep = " "))



      # create progress bar

      pb <- txtProgressBar(min = 0, max = length(cancers), style = 3)




      # Empety vector equal to the number of requested cancers

      Requested.Cancers <- vector("character", length(cancers))




      ## Getting the required gene expresssion profile ...

      # 'for' control structure for obtaining data and calculating the requested parameters

      for(i in 1:length(cancers)){

        ##--## print(as.character(cancers[i]))

        mycancerstudy = getCancerStudies(mycgds)[which(getCancerStudies(mycgds)[,2]==as.character(cancers[i])),1]



        # Finding the first characteristics of data in the cancer

        f.condition <- (getCaseLists(mycgds,mycancerstudy)[,2])[getCaseLists(mycgds,mycancerstudy)[,2] %in% L1.characteristics]

        mycaselist = getCaseLists(mycgds,mycancerstudy)[which(getCaseLists(mycgds,mycancerstudy)[,2] == if(length(f.condition) >= 1){

          f.condition[1]

        } else if(length(f.condition) == 0){

          stop(paste(cancers[i], "doesn't contain any", High.throughput.data, "data!", sep=" "))

        }) ,1]



        # Finding the second characteristics of data in the cancer

        s.condition <- (getGeneticProfiles(mycgds,mycancerstudy)[,2])[getGeneticProfiles(mycgds,mycancerstudy)[,2] %in% L2.characteristics]

        mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[which(getGeneticProfiles(mycgds,mycancerstudy)[,2] == if(length(s.condition) >= 1){

          s.condition[1]

        } else if (length(s.condition) == 0){

          stop(paste(cancers[i], "doesn't have an appropriate 'level 2' condition for", High.throughput.data, "data!", sep=" "))

        }) ,1]



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

        if("Mean.Expression" %in% data.presented.as){

          Serial.ProfileData.with.Mean.Expression <- vector("numeric", length=length(ordered.Genes))

        }

        if("Median" %in% data.presented.as){

          Serial.ProfileData.with.Median <- vector("numeric", length=length(ordered.Genes))

        }






        ### Calculating validity, frequency, Mean.Expression and Median values for each requested gene

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


          # Mean.Expression

          if("Mean.Expression" %in% data.presented.as){

            # Check all members are under cutoff

            if(!is.na(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j])])) &

               is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Mean.Expression[j] <- 0

              # Check all members are NaN

            } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) == 0 & all(!is.finite(ProfileData[,j])) &

                       is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Mean.Expression[j] <- NaN

              # Check all members are NA

            } else if (length((ProfileData[,j])[!is.nan(ProfileData[,j])]) > 0 & all(!is.finite(ProfileData[,j])) &

                       is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Mean.Expression[j] <- NA

              # Mean is number

            } else if (!is.nan(mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE))){

              Serial.ProfileData.with.Mean.Expression[j] <- mean(as.vector(ProfileData[,j])[abs(ProfileData[,j]) >= cutoff], na.rm=TRUE)

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

        if("Mean.Expression" %in% data.presented.as){

          dim(Serial.ProfileData.with.Mean.Expression) <- c(length(ordered.Genes),1)

        }

        if("Median" %in% data.presented.as){

          dim(Serial.ProfileData.with.Median) <- c(length(ordered.Genes),1)

        }




        # Combining serries of columns to obtain four matrixes presenting 'Frequency.Percentage', 'Frequency in ratio', 'Mean.Expression', 'Median' in a specific cancer

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

          if("Mean.Expression" %in% data.presented.as){

            Mean.Expression.Matrix <- Serial.ProfileData.with.Mean.Expression

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

          if("Mean.Expression" %in% data.presented.as){

            Mean.Expression.Matrix <- cbind(Mean.Expression.Matrix, Serial.ProfileData.with.Mean.Expression)

          }

          if("Median" %in% data.presented.as){

            Median.Matrix <- cbind(Median.Matrix, Serial.ProfileData.with.Median)

          }

        }




        # Number of tumor samples for each cancer

        if(shorteded.cancer.names == TRUE){

          Requested.Cancers[i] <- sapply(strsplit(as.character(cancers[i]), split=" (", fixed=TRUE), function(x) (x[1]))

        } else if(shorteded.cancer.names == FALSE){

          Requested.Cancers[i] <- as.character(cancers[i])

        }

        # update progress bar

        setTxtProgressBar(pb, i)

      }

      # Close progress bar

      close(pb)






      ### Adding dimention names to four matrixes

      print("Preparing a 'list' to fill with the requested data")



      # Check how cancer names will be delivered to matrix

      if(Use.CancerCode.as.Name == FALSE){

        Cancer.Names <- Requested.Cancers

      } else if(Use.CancerCode.as.Name == TRUE){

        CancerCode <- vector("character", length = length(cancers))

        for(s in 1:length(cancers)){

          CancerCode[s] <- getCancerStudies(mycgds)[which(getCancerStudies(mycgds)[,2]==as.character(cancers[s])),1]

        }

        Cancer.Names <- CancerCode

      }



      ## Creating the required matrixes and Overal list

      # Create empty list

      Data.list <- list()

      # filling the list

      if(validate.genes == TRUE){

        colnames(Validation.matrix) <- Cancer.Names

        Data.list$Validation.Results <- t(Validation.matrix)

      }

      if("Frequency.Percentage" %in% data.presented.as){

        dimnames(Frequency.Percentage.Matrix) <- list(ordered.Genes, Cancer.Names)

        Data.list$Frequency.Percentage <- t(Frequency.Percentage.Matrix)

        if(top.genes==TRUE){

          # Check if manual naming is requested

          Exp.matrix <- t(Frequency.Percentage.Matrix)

          if(is.matrix(cancers.matrix)){

            rownames(Exp.matrix) <- cancers.matrix[,2]

          }

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

        dimnames(Frequency.Ratio.Matrix) <- list(ordered.Genes, Cancer.Names)

        Data.list$Frequency.Ratio <- t(Frequency.Ratio.Matrix)

      }

      if("Mean.Expression" %in% data.presented.as){

        dimnames(Mean.Expression.Matrix) <- list(ordered.Genes, Cancer.Names)

        Data.list$Mean.Expression <- t(Mean.Expression.Matrix)

        if(top.genes==TRUE){

          # Check if manual naming is requested

          Exp.matrix <- t(Mean.Expression.Matrix)

          if(is.matrix(cancers.matrix)){

            rownames(Exp.matrix) <- cancers.matrix[,2]

          }

          # Removing NaN

          Exp.matrix[is.nan(Exp.matrix)] <- 0

          # Removing NA

          Exp.matrix[is.na(Exp.matrix)] <- 0

          # Removing rows that contain only 0

          Exp.matrix <- Exp.matrix[rowSums(Exp.matrix)!=0,]

          # Empty data.frame to be filled with data

          Mean.Expression.top.genes.dataframe <- data.frame(Gene.with.highest.value=vector("character", length = nrow(Exp.matrix)), Mean.Expression=vector("numeric", length = nrow(Exp.matrix)), stringsAsFactors = FALSE)

          # Finding genes with highest value

          for (m in 1:nrow(Exp.matrix)){

            Mean.Expression.top.genes.dataframe[m,1] <- names(which.max(Exp.matrix[m,]))

          }

          # Finding values corresponding to these genes

          for (n in 1:nrow(Exp.matrix)){

            Mean.Expression.top.genes.dataframe[n,2] <- max(Exp.matrix[n,])

          }

          # Creating data.frame for desired profile

          rownames(Mean.Expression.top.genes.dataframe) <- rownames(Exp.matrix)

          Data.list$Top.Genes.of.Mean.Expression <- Mean.Expression.top.genes.dataframe

        }

      }

      if("Median" %in% data.presented.as){

        dimnames(Median.Matrix) <- list(ordered.Genes, Cancer.Names)

        Data.list$Median <- t(Median.Matrix)

        if(top.genes==TRUE){

          # Check if manual naming is requested

          Exp.matrix <- t(Median.Matrix)

          if(is.matrix(cancers.matrix)){

            rownames(Exp.matrix) <- cancers.matrix[,2]

          }

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

      assign(paste(names(genes)[[g]], "Data", "List", sep = "."), Data.list, envir = globalenv())

    }












    ### Plotting Heatmap and preparing excel files

    List.to.Go <- get(paste(names(genes)[[g]], "Data", "List", sep = "."))

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

      if(Temporary.source.name %in% c("Frequency.Percentage","Mean.Expression" ,"Median") & is.matrix(heatmap.data)){

        heatmap.data <- heatmap.data[rowSums(heatmap.data, na.rm = TRUE)!=0,]

      }

      if(Temporary.source.name %in% c("Frequency.Percentage","Mean.Expression" ,"Median") & is.matrix(heatmap.data)){

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

          if(!is.matrix(cancers.matrix)){

            tissue <- colnames(heatmap.data)

          } else if(is.matrix(cancers.matrix)){

            tissue <- cancers.matrix[,2]

            colnames(heatmap.data) <- cancers.matrix[,2]

          } else {

            stop("If you desire to enter your own cancer names, make sure that personal.cancer.names is a matrix of two columns,

                 fist column contains the real cancer names and the second column contains the desired cancer names")

          }




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

          png(filename=paste(getwd(), paste(gsub(x = names(genes)[g], pattern = "\\.", replacement = "-"), " ", gsub(x = Temporary.source.name, pattern = "\\.", replacement = " "), " (",cutoff.phrase, "=", cutoff, ")" , ".png", sep=""), sep="/"), width=9.5, height= 11, units = "in", res=resolution)

          heatmap.2(heatmap.data, labCol=tissue, na.color="light gray", trace="none", symbreaks = T, col=hmcol, cexRow = RowCex, cexCol= ColCex,

                    margins = heatmapMargines, srtCol = angle.for.heatmap.cancernames)

          dev.off()

        }

        }




      ## Exporting the expression profile

      Exp.matrix2 <- Temporary.source

      if(is.matrix(cancers.matrix)){

        rownames(Exp.matrix2) <- cancers.matrix[,2]

      }

      # Rounding values

      if(round==TRUE & Temporary.source.name %in% c("Frequency.Percentage","Mean.Expression" ,"Median")){

        Exp.matrix2 <- round(Exp.matrix2, digits=2)

      }

      # Copy data to a variable to modify

      Save.matrix <- Exp.matrix2

      if(!(Temporary.source.name %in% c("Top.Genes.of.Frequency.Percentage","Top.Genes.of.Mean.Expression" ,"Top.Genes.of.Median"))){

        # Replace 'NaN' values with character string

        Save.matrix[is.nan(Save.matrix)] <- "NaN"

        # Replace 'NA' values with character string

        Save.matrix[is.na(Save.matrix)] <- "NA"

      }

      if(Temporary.source.name %in% c("Top.Genes.of.Frequency.Percentage","Top.Genes.of.Mean.Expression" ,"Top.Genes.of.Median")){

        colnames(Save.matrix)[1:2] <- c("Gene with highest value", if(Temporary.source.name == "Top.Genes.of.Frequency.Percentage"){"Frequency Percentage"}

                                        else if(Temporary.source.name == "Top.Genes.of.Mean.Expression"){"Mean Expression"}

                                        else if(Temporary.source.name == "Top.Genes.of.Median"){"Median"})

      }



      # Saving the expression profile

      write.xlsx(Save.matrix, file = paste(gsub(x = names(genes)[g], pattern = "\\.", replacement = "-"), " (",cutoff.phrase, "=", cutoff, ")" , ".xlsx", sep = ""), sheetName = gsub(x = Temporary.source.name, pattern = "\\.", replacement = " "), append = TRUE)




      # Update progressbar

      setTxtProgressBar(pb2, segment)

      }

    # Close progress bar

    close(pb2)

    # Print information

    if(any(c("Frequency.Percentage", "Mean.Expression", "Median") %in% names(List.to.Go))){

      print(paste("Requested excel and PNG files were saved in", getwd(), sep = " "))

    } else if(!any(c("Frequency.Percentage", "Mean.Expression", "Median") %in% names(List.to.Go))){

      print(paste("Requested excel files were saved in", getwd(), sep = " "))

    }

  }

  setwd(dir.address)
}
