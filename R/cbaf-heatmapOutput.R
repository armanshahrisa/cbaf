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





matrix(, nrow = 15, ncol = 0)








genesList <- list(b=c("HOTAIR", "GAS5", "SNHG6"), hhc = c("PVT1", "GAS5", "DBBBGTZZXCC", "NSD2"))

genesList <- list(b=c("HOTAIR", "GAS5"))

cancers <- c("Acute Myeloid Leukemia (TCGA, NEJM 2013)", "Acute Myeloid Leukemia (TCGA, Provisional)")



obtainMultipleStudies(genesList, "test", cancers, desiredTechnique = "RNA-seq")








ObM.test$b$`Acute_Myeloid_Leukemia_(TCGA,_Provisional)`[,sort(colnames(ObM.test$b$`Acute_Myeloid_Leukemia_(TCGA,_Provisional)`))]




vbn <- matrix(c(1,2,3,4), ncol=4, nrow=1)
dimnames(vbn) <- list(c("can1"), c("a", "b", "c", "d"))

mkj <- matrix(c(5,6,7,8), ncol=4, nrow=1)
dimnames(mkj) <- list(c("can2"), c("a", "b", "c", "d"))

rbind(vbn, mkj)


