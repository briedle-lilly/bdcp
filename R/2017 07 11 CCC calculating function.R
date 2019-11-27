
##Writing a function which will give us the 7 desired CCC scores in each setting.

CCC.fxn <- function(arg1, arg2){
  n <- length(arg1)
  CCC <- 2*(n-1 / n)*cov(arg1,arg2) / ((n-1 / n)*var(arg1) + (n-1 / n)*var(arg2) + (mean(arg1) - mean(arg2))**2)
  return(round(CCC, digits = 5))
}

specific.CCC.fxn <- function(dataset){
  Wald.full <- CCC.fxn(dataset[[1]],dataset[[7]])
  
  Wald.partial <- CCC.fxn(dataset[[2]],dataset[[8]])
  
  score.full <- CCC.fxn(dataset[[3]],dataset[[9]])
  
  score.partial <- CCC.fxn(dataset[[4]],dataset[[10]])
  
  LRTKL.full <- CCC.fxn(dataset[[5]],dataset[[11]])
  
  LRTPIKL.partial <- CCC.fxn(dataset[[6]],dataset[[13]])
  
  LRTKL.partial <- CCC.fxn(dataset[[6]],dataset[[12]])
  
  return(c(Wald.full, Wald.partial, score.full, score.partial, LRTKL.full, LRTPIKL.partial, LRTKL.partial))
}

specific.CCC.fxn(dataset = linregtruenulln500)

load(file = "C:\\Users\\bnr\\Desktop\\Dissertation Stuff\\R programs\\saved pvalues and DCPs\\linreglightn100beta0.175sig60")
load(file = "C:\\Users\\bnr\\Desktop\\Dissertation Stuff\\R programs\\saved pvalues and DCPs\\linreglightn250beta0.1sig57")
load(file = "C:\\Users\\bnr\\Desktop\\Dissertation Stuff\\R programs\\saved pvalues and DCPs\\linreglightn1000beta0.05sig55")
     
load(file = "C:\\Users\\bnr\\Desktop\\Dissertation Stuff\\R programs\\saved pvalues and DCPs\\logreglightn100beta0.07")
load(file = "C:\\Users\\bnr\\Desktop\\Dissertation Stuff\\R programs\\saved pvalues and DCPs\\logreglightn250beta0.04")
load(file = "C:\\Users\\bnr\\Desktop\\Dissertation Stuff\\R programs\\saved pvalues and DCPs\\logreglightn1000beta0.01")

load(file = "C:\\Users\\bnr\\Desktop\\Dissertation Stuff\\R programs\\saved pvalues and DCPs\\linregtruenulln500")
load(file = "C:\\Users\\bnr\\Desktop\\Dissertation Stuff\\R programs\\saved pvalues and DCPs\\logregtruenulln500")


specific.CCC.fxn(dataset = linreglightn100beta0.175sig60 )
specific.CCC.fxn(dataset = linreglightn250beta0.1sig57)
specific.CCC.fxn(dataset = linreglightn1000beta0.05sig55)

specific.CCC.fxn(dataset = logreglightn100beta0.07)
specific.CCC.fxn(dataset = logreglightn250beta0.04)
specific.CCC.fxn(dataset = logreglightn1000beta0.01)


