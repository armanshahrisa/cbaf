shortenStudiesNames = TRUE


if(shortenStudiesNames == TRUE){

  groupNames <- sapply(strsplit(as.character(studiesNames), split=" (", fixed=TRUE), function(x) (x[1]))

  }








aaaa <- list()
water <- list()

aaaa$ water <- list()


aaaa[[1]][[1]] <- c(1,2)
names(aaaa[[1]])[1] <- "Hello"




x <- -6

dre <- function(x){

  if (x < 0){
    return("skipped")
  }

  print("Hi!")

}

dre2 <- function(){
  print("Well done!")
}


dde <- function(x){
  dre(x=x)
  dre2()

}


dre(-6)


dde(3)

dde(-3)














genesList <- list(b=c("HOTAIR", "GAS5", "SNHG6"))

genesList <- list(b=c("HOTAIR", "GAS5"))

cancers <- c("Acute Myeloid Leukemia (TCGA, NEJM 2013)", "Acute Myeloid Leukemia (TCGA, Provisional)")


obtainMultipleStudies(genesList, "test", cancers, desiredTechnique = "RNA-seq")

