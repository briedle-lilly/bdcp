
library(MASS)
###here are the linear regression functions which will produce the score vector and the information matrix.  When we only need portions of the score and info
#we will subset it within the larger function.  These functions will always produce the entire score and information.

lin.reg.obs.info.fxn <- function(beta.est, sig2.est, y, X){
  beta.dim <- length(beta.est)  #gives dimension of information matrix
  #the size of the matrix needs to account for the si2.est also
  info.matrix <- matrix(data = NA, nrow = beta.dim + 1, ncol = beta.dim + 1)
  
  #negative of second derivatives of beta twice
  info.matrix[c(1:beta.dim),c(1:beta.dim)] <- 1/sig2.est*(t(X) %*% X)
  
  #derivatives wrt beta sigma squared and then beta
  info.matrix[beta.dim +1, c(1:beta.dim)] <- 1 / sig2.est^2 * t(y - X %*% beta.est) %*% X
  
  #derivatives wrt to beta then sigma squared
  info.matrix[c(1:beta.dim),beta.dim + 1] <- 1 / sig2.est^2 * t(X) %*% (y - X %*% beta.est)
  
  ##derivative wrt to sig squared twice
  info.matrix[beta.dim + 1, beta.dim + 1] <- -length(y) / (2*sig2.est^2) + 1 / (sig2.est^3)*t(y - X %*% beta.est) %*% (y - X %*% beta.est)
  
  return(info.matrix)
}


##This is function which will produce the score vector. Beta.est and sig2.est allow for the score vector to be evaluated at different values.
#dl / dbeta = -1/sig2 *[-X^T*y + X^T*X*beta]
#dl / dsig2 = -n / (2*sig2) + 1 / (2*sig2^2)*RSS (model specific RSS)
lin.reg.score.fxn <- function(beta.est, sig2.est, y, X){
  score.vector <- vector(length = length(beta.est) + 1)
  
  #derivative wrt beta
  score.vector[1:length(beta.est)] <- - 1 / sig2.est *(-t(X) %*% y + t(X) %*% X %*% beta.est)
  
  #derivative wrt sigsquare
  score.vector[length(beta.est) + 1] <- - length(y) / (2*sig2.est) + 1 / (2*sig2.est^2)*t(y - X %*% beta.est) %*% (y - X %*% beta.est)
  
  return(score.vector)
}



#########SHOULD CREATE A SEPARATE BOOTSTRAPPING FUNCTION WHICH WE WILL CALL B TIMES TO CALCULATE B ESTIMATED DISCREPANCY COMPARISONS.  THIS SHOULD HELP ENSURE
#THAT WE ONLY STORE THE INFORMATION FROM ONE BOOTSTRAP SAMPLE AT A TIME. 
draw.boot.and.calc.disc.comparisons.fxn <- function(start.data, start.coefs.m123, start.sig2.hat.m123, start.sig2.hat.m1, obs.info.atbetahat.m123, inv.inv.info.atbetahat.m123, Wald.full.null.null.disc,
                                                    Wald.partial.null.null.disc, score.full.null.null.disc, score.partial.null.null.disc,  KL.full.null.mnull,
                                                    KL.partial.null.para.int.m1){
  
  
  
  n <- nrow(start.data)
  start.X <- matrix(cbind(as.vector(rep(1,n)),start.data[,1],start.data[,2],start.data[,3]), ncol = 4)
  start.y <- start.data[,4]
  
  
  #By including this in a repeat statement, we ensure there will be at least 2 ones and two zeroes in the bootstrap dataset.
  
    boot.indices <- sample(1:n, replace = TRUE)
    boot.data <- as.data.frame(rbind(start.data[boot.indices,]))
    y.b <- boot.data[,4]
  
  colnames(boot.data) <- c("x1.b" , "x2.b" , "x3.b", "y.b")
  boot.X <- as.matrix(cbind(as.vector(rep(1,n)),boot.data[,1:3]))
  
  
  #this is the estimated parameter vector for the nuisance parameters under the model m1 (i.e. the null model in the partial null setting).  We need to calculate
  #this when we are working with the overall KLD (not the parameter of interest overall KLD).  With the parameter of interest overall KLD, we plug in the pseudo-true
  #parameter values, which when we move to the bootstrap world to estimate the discrepancy, we would plug in the original sample's MLEs for the nuisance parameters.
  #Thus, we would not need bootstrap parameter estimates for the null model when we are using the PARAMETERS OF INTEREST overall KLD.  Note that we would also need
  #this bootstrap parameter vector for the partial null overall score discrepancy setting (i.e. not using the parameters of interest overall score discrepancy). 
  #However, we are not going to explore the standard, overall score discrepancy in the partial null setting.
  boot.beta.m1 <- c(lm(y.b ~ x1.b, data = boot.data)$coefficients,0,0)
  boot.sig2.hat.m1 <- as.numeric(t(y.b - boot.X %*% boot.beta.m1) %*% (y.b - boot.X %*% boot.beta.m1) / n)
  #This is the alternative model's bootstrap parameter vector.  We will be using this with the overall score and KL discrepancies (NOT THE PARAMETER OF INTEREST
  #OVERALL DISCREPANCIES)
  boot.beta.m123 <- as.vector(lm(y.b ~ x1.b + x2.b + x3.b, data = boot.data)$coefficients)
  boot.sig2.hat.m123 <- as.numeric(t(y.b - boot.X %*% boot.beta.m123) %*% (y.b - boot.X %*% boot.beta.m123) / n)
  #THIS IS THE ALTERNATIVE MODEL'S BOOTSTRAP PARAMETER VECTOR FOR THE NUISANCE PARAMETERS, CONDITIONED ON THE PARAMETERS OF INTEREST BEING EQUAL
  #TO THEIR MLES UNDER THE ORIGINAL SAMPLE. tHIS IS OUR THETA.HAT.STAR.C(1) 
  #In order to get this fit, for x.2 and x.3, I first need to define an offset that is equal to beta.hat.zero + x1*beta.hat.one (This is like the fixed intercept for
  #the model).  Thus, I need to not allow an intercept in this model fit. Note that boot.data[,1] corresponds to variable x1.  Note that we also have to constrain that 
  #sig2.est = sig2.hat.  However, since the mean is independent of the variance (which can be fairly easily seen by the form of the log-likelihood), then the beta which minimizes
  #is independent of the value of sigma squared anyway.  
  boot.and.orig.cond.beta.m123 <- c(start.coefs.m123[c(1,2)], 
                                    lm(y.b ~ x2.b + x3.b + offset(start.coefs.m123[1] + start.coefs.m123[2]*boot.X[,2]) - 1, data = boot.data)$coefficients)
  
  ###We now need to calculate the various bootstrap-estimated discrepancies.  Recall that for the Wald discrepancy the null model's bootstrap-estimated discrepancy is 
  #the wald statistic.  For the score discrepancy, the full null score discrepancy is the score stat and in the partial null setting we use the parameter of interest
  #overall score discrepancy, which yields a null model's estimated discrepancy equal to the score statistic.  For the KL discrepancy, the relationship is a little more
  #complicated.
  Wald.disc.full.null.m123 <- t(c(boot.beta.m123, boot.sig2.hat.m123) - c(start.coefs.m123, start.sig2.hat.m123)) %*% obs.info.atbetahat.m123 %*% (c(boot.beta.m123, boot.sig2.hat.m123) - c(start.coefs.m123, start.sig2.hat.m123))
  Wald.disc.partial.null.m123 <- t(boot.beta.m123[c(3,4)] - start.coefs.m123[c(3,4)]) %*% inv.inv.info.atbetahat.m123 %*% (boot.beta.m123[c(3,4)] - start.coefs.m123[c(3,4)])
  
  
  ###NOW WORKING ON THE SCORE DISCREPANCIES.  
  #In this case because we use the original data X and y, and the bootstrap parameter estimates.
  #For full null, we use the entire bootstrap parameter vector boot.beta.m123.  For the partial null setting, we use the nuisance parameters are set equal to
  #their MLEs under the original sample (start.coefs.m123[c(1,2)]) and their conditional bootstrap-estimated parameter of interest.  Thus, for the partial null we 
  #use boot.cond.nuis.para.beta.m123 for the parameter vector
  
  #Note that (start.data[,4]) is y.  The scores are based on the original data and the boot estimates, so we use X and y in our expressions, and the bootstrap parameter estimates 
  #(boot.beta and boot.sig2.hat).  For the bootstrap-estimated score discrepancy, in the full null setting we evaluate the score vector and observed information at the bootstrap MLEs.
  #We evaluate them at the (theta.hat.star.C1 and theta.hat.two) for the partial null setting.
  U.boot.full.null.m123 <- lin.reg.score.fxn(beta.est = boot.beta.m123, sig2.est = boot.sig2.hat.m123, y = start.y, X = start.X)
  #The partial null setting uses boot.and.orig.cond.beta.m123 for the beta estimates (which is a combo of the original MLEs for the nuisance parameters and the conditional MLEs
  #for the parameters of interest based on the bootstrap sample, where the conditioning is on the nuisance parameter MLEs under the original sample.)
  U.boot.partial.null.m123 <- lin.reg.score.fxn(beta.est = boot.and.orig.cond.beta.m123 , sig2.est = start.sig2.hat.m123, y = start.y, X = start.X)
  
  #The informations are based the original data and bootstrap parameter estimates so we use X and the p.hat (which incorporates the bootstrap parameter estimates)
  obs.info.atbetahatstar.m123 <- lin.reg.obs.info.fxn(beta.est = boot.beta.m123, sig2.est = boot.sig2.hat.m123, y = start.y, X = start.X)
  #The information for the partial null should only include be calculated for the parameters of interest.  We then use [I_11]^-1.  This will be evaluated at the
  #boot.cond.nuis.para.beta.m123, which is the MLES of the nuisance under the original sample and the paras of interest under the bootstrap, conditioned on nuisance MLES.  Note that
  #sig2 is a nuisance parameter so it is evaluated at the MLE under the original sample.
  obs.info.para.int.atbetahatstarC.m123 <- lin.reg.obs.info.fxn(beta.est = boot.and.orig.cond.beta.m123 , sig2.est = start.sig2.hat.m123, y = start.y, X = start.X)[c(3,4),c(3,4)]
  
  ##Need inverse infos
  inv.obs.info.atbetahatstar.m123 <- solve(obs.info.atbetahatstar.m123)
  inv.obs.info.para.int.atbetahatstarC.m123 <- solve(obs.info.para.int.atbetahatstarC.m123)
  
  ##Now calculate the alternative models' score discrepancies
  score.disc.full.null.m123 <- t(U.boot.full.null.m123) %*% inv.obs.info.atbetahatstar.m123 %*% U.boot.full.null.m123
  score.disc.partial.null.m123 <- t(U.boot.partial.null.m123)[3:4] %*% inv.obs.info.para.int.atbetahatstarC.m123 %*% U.boot.partial.null.m123[3:4]
  
  ##WORKING ON THE KL DISCREPANCY CALCULATIONS NOW
  ##We calculate the eta values to make the discrepancies easier to write.  The alternative model's full and partial null OVERALL KL discrepancies are equivalent.  
  #Under the PARAMETER OF INTEREST OVERALL KLD we have to use a different parameter vector (theta.hat.star.C1, theta.hat.2)
  eta.boot.and.orig.overall.KL.m123 <- start.X %*% boot.beta.m123
  eta.boot.and.orig.para.int.overall.KL.m123 <- start.X %*% boot.and.orig.cond.beta.m123
  
  KL.overall.disc.m123 <- n*log(2*pi) + n*log(boot.sig2.hat.m123) + 1 / boot.sig2.hat.m123 * t(start.y - eta.boot.and.orig.overall.KL.m123) %*% (start.y - eta.boot.and.orig.overall.KL.m123)
                          
  #Recall that this bootstrap-estimated discrepancy is going to use the nuisance parameter (i.e. beta0, beta1 and sig2) estimates under the original sample and the bootstrap
  #parameter estimates under the bootstrap data, conditioned on the nuisance parameters equaling their MLEs
  KL.para.int.overall.disc.m123 <- n*log(2*pi) + n*log(start.sig2.hat.m123)   + 1 / start.sig2.hat.m123 * t(start.y - eta.boot.and.orig.para.int.overall.KL.m123) %*% (start.y - eta.boot.and.orig.para.int.overall.KL.m123)
                                   
  
  #FOR THE PARTIAL NULL SETTING USING THE OVERALL KLD (AS OPPOSED TO THE PARA OF INT. OVERALL KLD), 
  #THE NULL MODEL (i.e. model m1) ACTUALLY DEPENDS ON THE BOOTSTRAP SAMPLE, SO WE CALCULATE THAT NOW. NOTE THE VARIANCE ESTIMATOR MUST ALSO BE ESTIMATED BY THE BOOTSTRAP.
  eta.boot.and.orig.overall.KL.m1 <- start.X %*% boot.beta.m1
  partial.null.KL.overall.disc.m1 <- n*log(2*pi) + n*log(boot.sig2.hat.m1)   + 1 / boot.sig2.hat.m1 * t(start.y - eta.boot.and.orig.overall.KL.m1) %*% (start.y - eta.boot.and.orig.overall.KL.m1)
  
  
  
  
  
  
  
  
  
  ##NOW CALCULATING THE DISCREPANCY COMPARISONS. Recall that for the Wald and score discrepancies the null models' discrepancy are actually the respective test statistics.
  ##THIS IS NOT TRUE FOR THE KL DISCREPANCIES.
  Wald.full.null.disc.comparison <- ifelse(Wald.full.null.null.disc < Wald.disc.full.null.m123, 1, 0)
  Wald.partial.null.disc.comparison <- ifelse(Wald.partial.null.null.disc < Wald.disc.partial.null.m123, 1, 0)
  
  score.full.null.disc.comparison <- ifelse(score.full.null.null.disc < score.disc.full.null.m123, 1, 0)
  score.partial.null.disc.comparison <- ifelse(score.partial.null.null.disc < score.disc.partial.null.m123, 1, 0)
  
  KL.full.null.disc.comparison <- ifelse(KL.full.null.mnull < KL.overall.disc.m123, 1, 0)
  KL.partial.null.overall.disc.comparison <- ifelse(partial.null.KL.overall.disc.m1 < KL.overall.disc.m123, 1, 0)
  KL.partial.null.para.int.overall.disc.comparison <- ifelse(KL.partial.null.para.int.m1 < KL.para.int.overall.disc.m123, 1, 0)
  
  return(c(Wald.full.null.disc.comparison,Wald.partial.null.disc.comparison,score.full.null.disc.comparison,score.partial.null.disc.comparison,
           KL.full.null.disc.comparison, KL.partial.null.overall.disc.comparison, KL.partial.null.para.int.overall.disc.comparison))
  
}





##This function will do all a completey cycle of calculating logistic regression p-values and bootstrap-estimated discrepancy comparison probabilities.  Our larger
#function will perform this process for iter iterations.

###GOING TO GO INSIDE THE FUNCTION FOR NOW SO WE WILL COMMENT OUT THE FUNCTION NOTATION
single.iter.lin.reg.fxn <- function(n =n, B = B, true.beta0 = true.beta0, true.beta1 = true.beta1, true.beta2 = true.beta2, true.beta3 = true.beta3
                                    , true.beta4 = true.beta4, true.sig2, mnull.sig2){
  
  
  #iter = 10
  #n = 100
  #B = 100
  #true.beta0 = 0
  #true.beta1 = 0
  #true.beta2 = 0
  #true.beta3 = 0
  #true.beta4 = 0
  #true.sig2 = 50
  #mnull.sig2 = 50

  
  ##Code if we want all predictors be independent of each other
  #x1.success.prob <- 0.6
  #x2.success.prob <- 0.4
  #x3.min <- -0.8
  #x3.max <- 0
  #x4.min <- 0
  #x4.max <- 2
  
  
  #x1 <- rbinom(n = n, size = 1, prob = x1.success.prob)
  #x2 <- rbinom(n = n, size = 1, prob = x2.success.prob)
  #x3 <- runif(n = n, min = x3.min, max = x3.max)
  #x4 <- runif(n = n, min = x4.min, max = x4.max)
  
  ###Code if we want to have the predictors have some correlation.
  corr <- 0.8
  var <- 100
  cov <- 100*0.8*0.8
  sig.mat <- matrix(c(rep(c(100,rep(cov,4)),3),100), nrow = 4, ncol = 4)
  X.includex4 <- mvrnorm(n = n, mu = c(2,2,-2,-2), Sigma = sig.mat)
  
  x1 <- X.includex4[,1]
  x2 <- X.includex4[,2]
  x3 <- X.includex4[,3]
  x4 <- X.includex4[,4]
  
  
  ##Calculating the original (sample from which we draw the bootstrap) eta values using the true betas and the observed x values
  orig.eta <- true.beta0 + x1*true.beta1 + x2*true.beta2 + x3*true.beta3 + x4* true.beta4
  
  
  
  ##Drawing a sample of y values for the ith iteration.  The sample will be distributed normally with mean X*beta (with beta including beta4).
    y <- orig.eta + rnorm(n = n, mean = 0, sd = sqrt(true.sig2))
   
  
  
  
  
  
  
  
  ###Storing the design matrix of the largest model (orig.X) and the the data x1,x2,x3,y (full.orig.data)
  full.orig.data <- matrix(cbind(x1,x2,x3,y), ncol = 4)
  orig.X <- matrix(cbind(1,x1,x2,x3), ncol = 4)
  
  
  ###Storing the coefficients of the fitted models for the original data
  orig.coef.mnull <- rep(0,4)
  orig.coef.m1 <- c(lm(y ~ x1)$coefficients, rep(0,2))
  orig.coef.m123 <- as.vector(lm(y ~ x1 + x2 + x3)$coefficients)
  
  ##calculating an orig.eta value for each model and for each of iter iterations.  This is the fitted value of eta for each of the models.  Note that
  #because the null model has 0 for all betas, then its fitted eta vector will simply be n zeroes.
  orig.eta.mnull <- rep(0,n)
  orig.eta.m1 <- orig.X %*% orig.coef.m1
  orig.eta.m123 <- orig.X %*% orig.coef.m123
  
  #sigma2 estimates
  sig2.hat.m1 <- as.numeric(t(y - orig.eta.m1) %*% (y - orig.eta.m1) / n)
  sig2.hat.m123 <- as.numeric(t(y - orig.eta.m123) %*% (y - orig.eta.m123) / n)
  
  #Calculating the LRT test stats and p-value##
  #We know that logL = -n/2log(2pi) - n/2log(sig2) - 1 / (2*sig2) * eta^T*eta and eta_i = x_i*beta.
  #For the null model where we also pre-specify sigma squared, then the likelihood ratio test does not simplify to simply the log of the ratio of the variance
  #estimates because for the null model, the "estimator" of sigma squared will not cancel with the sum of squares term in the likelihood.
  LRT.stat.m123vmnull <- n*log(mnull.sig2 / sig2.hat.m123) - n + 1 / mnull.sig2*sum(y^2)
  LRT.pval.m123vmnull <- 1 - pchisq(LRT.stat.m123vmnull, df = 5)
  #We also want the null model's discrepancy for the full null setting, where the null model pre-specifies all values to be zero.
 ##Recall that the KLD for the null model in the full null setting evaluated the -2loglik at the the pre-specified parameter values.  Also eta.mnull = 0,
  #which simplifies the last term to the sum of y^2.
  KLD.full.null.mnull <- -2*(-n/2 * log(2*pi) - n / 2 * log(mnull.sig2) - 1 / (2*mnull.sig2) * sum(y^2))
  
  LRT.stat.m123vm1 <- n*log(sig2.hat.m1 / sig2.hat.m123)
  LRT.pval.m123vm1 <- 1 - pchisq(LRT.stat.m123vm1, df = 2)
  ##We want the null model's discrepancy for the partial null setting where we are using the para of interest overall disc, so that the nuisance parameters (beta0,
  #beta1 and sig2) are estimated with the original sample, but the parameters of interest are set equal to zero.
  KLD.partial.null.para.int.m1 <- -2*(-n / 2 * log(2*pi) - n / 2 *log(sig2.hat.m1) - 1 / (2*sig2.hat.m1) * t(y - orig.eta.m1) %*% (y - orig.eta.m1))

  ###Working on Wald test p-values
  #For the Wald statistic in the full null setting we need the observed information evaluated at the MLE.  For the partial null setting we need to take the
  #inverse of the entire information matrix.  Then grab the (1,1) component of the inverse information (i.e. the parameters of interests component) and invert it.
  #Notice that all of this is being done on the larger model, which in our case for both the full and partial null settings is the model m123.
  obs.info.m123 <- lin.reg.obs.info.fxn(beta.est = orig.coef.m123, sig2.est = sig2.hat.m123, y = y, X = orig.X)
 
  ###THis is equivalent to I^{-1}(beta.hat).  This is what we will be using for the Wald statistic because it uses the observed info evaluated at the MLEs.  
  #The score statistic uses the information evaluated at the pre-defined values.  Thus for the full null model, we evaulate at the pre-specified beta values (i.e. 0).  For the
  #partial null model we evaluate info at beta2 = beta3 = 0 and let beta0 and beta1 equal their MLEs assuming beta2 and beta3 are zero.
  inv.info.m123 <- solve(obs.info.m123)
  
  
  ##These are giving the inverse of the inverse info where we only take the entries that apply to parameters being tested for the given hypothesis test.
  ##We only need to take the inverse of the elements which correspond to the parameters being tested.
  inv.inv.In11.m123vm1 <- solve(inv.info.m123[c(3:4),c(3:4)])
  
  #CALCULATING THE WALD TEST P-VALUES
  #Recall that for the full null setting (i.e. m123 vs mnull), the middle piece of the Wald stat is actually the information (or we could think of it as taking the
  #inverse information, then inverting the piece corresponding to the parameters of interest (i.e. all the parameters).  Thus, we would be inverting the 
  #entire matrix and then inverting the entire matrix again.
  Wald.stat.m123vmnull <- t(c(orig.coef.m123,sig2.hat.m123) - c(rep(0,4),mnull.sig2)) %*% obs.info.m123 %*% (c(orig.coef.m123,sig2.hat.m123) - c(rep(0,4),mnull.sig2))
  Wald.pval.m123vmnull<- 1 - pchisq(Wald.stat.m123vmnull, df = 5)
  
  Wald.stat.m123vm1 <- t(orig.coef.m123[c(3,4)]) %*% inv.inv.In11.m123vm1 %*% orig.coef.m123[c(3,4)]
  Wald.pval.m123vm1 <- 1 - pchisq(Wald.stat.m123vm1, df = 2)
  
  #######CALCULATING THE SCORE TEST STATISTIC AND ITS P-VALUE
  #we know that to calculate the score statistic we need to determine what the score and information are.  We then evaluate the score and information at the 
  #pre-specified values of beta (in the full null setting) or at the pre-specified values of the parameters of interest and the conditional MLEs for the 
  #nuisance parameters.
  #RECALL THAT FOR THE SCORE EVERYTHING NEEDS TO BE EVALUATED AT THE PRE-SPECIFIED VALUES THETA_0 NOT THETA_HAT!!!!!!
  
  #Defining the observed probabilities for a patient and for the pre-specified beta values.  Also, I defined each of the coefficients (i.e. beta hats) so that 
  #they have length 4, placing zeroes for the betas set equal to zero for the non-full models.  For the null model in the full null setting, all beta = 0.
  ###obs.info.m123 is equivalent to I^{-1}(beta.hat).  This is what we will be using for the Wald statistic because it uses the observed info evaluated at the MLEs.  
  #The score statistic uses the information evaluated at the pre-defined values.  Thus for the full null model, we evaulate at the pre-specified beta values (i.e. 0).  For the
  #partial null model we evaluate info at beta2 = beta3 = 0 and let beta0 and beta1 equal their MLEs assuming beta2 and beta3 are zero.
  U.m123vmnull <- lin.reg.score.fxn(beta.est = orig.coef.mnull, sig2.est = mnull.sig2, y = y, X = orig.X )
  U.m123vm1 <- lin.reg.score.fxn(beta.est = orig.coef.m1, sig2.est = sig2.hat.m1, y = y, X = orig.X) #We only want the score vector for para. of int.
  
   
    
  obs.info.atbeta0.mnull <- lin.reg.obs.info.fxn(beta.est = orig.coef.mnull, sig2.est = mnull.sig2, y = y, X = orig.X)
  obs.info.atbeta0.m1 <- lin.reg.obs.info.fxn(beta.est = orig.coef.m1, sig2.est = sig2.hat.m1, y = y, X = orig.X)
  
  ##For the full null, we only need the inverse information for the score statistic.  For the partial null, we need take the inverse of the information and then 
  #just grab the piece of the inverse information corresponding to the paramters of interest
  inv.obs.info.atbeta0.mnull <- solve(obs.info.atbeta0.mnull)
  inv.obs.info.atbeta0.m1 <- solve(obs.info.atbeta0.mnull)[c(3,4), c(3,4)]
  
  score.stat.m123vmnull <- t(U.m123vmnull) %*% inv.obs.info.atbeta0.mnull %*% U.m123vmnull
  score.stat.m123vm1   <- t(U.m123vm1[c(3,4)]) %*% inv.obs.info.atbeta0.m1 %*% U.m123vm1[c(3,4)]
  
  score.pval.m123vmnull <- 1 - pchisq(score.stat.m123vmnull, df = 5)
  score.pval.m123vm1 <- 1 - pchisq(score.stat.m123vm1, df = 2)
  
  
  
  ##THIS IS NOW ASKING THE FUNCTION WHICH DOES THE DRAWING OF THE BOOTSTRAP SAMPLE AND THE CALCULATING OF THE BOOTSTRAP DISCREPANCIES AND THEIR COMPARISON PROBS
  #TO BE DONE B TIMES.  BY CALLING THE FUNCTION B TIMES I BELIEVE IT WILL ONLY HAVE TO STORE THE BOOTSTRAP SAMPLE ON WHICH IT IS CURRENTLY WORKING.
  boot.disc.comparisons <- matrix(data = NA, nrow = B, ncol = 7) 
  for(b in 1:B){
    boot.disc.comparisons[b,] <- draw.boot.and.calc.disc.comparisons.fxn(
      start.data = full.orig.data, 
      start.coefs.m123 = orig.coef.m123, 
      start.sig2.hat.m123 = sig2.hat.m123,
      start.sig2.hat.m1 = sig2.hat.m1,
      obs.info.atbetahat.m123 = obs.info.m123,
      inv.inv.info.atbetahat.m123 = inv.inv.In11.m123vm1,
      Wald.full.null.null.disc = Wald.stat.m123vmnull, 
      Wald.partial.null.null.disc = Wald.stat.m123vm1,
      score.full.null.null.disc = score.stat.m123vmnull,
      score.partial.null.null.disc = score.stat.m123vm1,
      KL.full.null.mnull = KLD.full.null.mnull, 
      KL.partial.null.para.int.m1 = KLD.partial.null.para.int.m1)
  }
  
  
  
  boot.disc.comparison.probs <- apply(X = boot.disc.comparisons, 2, FUN = sum) / B
  
  
  
  return(c(Wald.pval.m123vmnull,
           Wald.pval.m123vm1, 
           score.pval.m123vmnull,
           score.pval.m123vm1,
           LRT.pval.m123vmnull,
           LRT.pval.m123vm1,
           boot.disc.comparison.probs[1],
           boot.disc.comparison.probs[2],
           boot.disc.comparison.probs[3],
           boot.disc.comparison.probs[4],
           boot.disc.comparison.probs[5],
           boot.disc.comparison.probs[6],
           boot.disc.comparison.probs[7]))
}






###This is the logistic regression function that does not work in parallel.  I could not figure out how to make this work in parallel so I will just have it run sequentially.
lin.reg.fxn <- function(iter, n, B, true.beta0, true.beta1, true.beta2, true.beta3, true.beta4, true.sig2, mnull.sig2){
  
  lin.reg.data <- matrix(data = NA, nrow = iter, ncol = 13)
  for(i in 1:iter){
    lin.reg.data[i,] <- single.iter.lin.reg.fxn(n = n, B = B, true.beta0 = true.beta0, true.beta1 = true.beta1, true.beta2 = true.beta2, true.beta3 = true.beta3
                                                , true.beta4 = true.beta4, true.sig2 = true.sig2, mnull.sig2 = mnull.sig2)
  } 
  
 iteration <- 1:iter 
  
  list.lin.reg.data <- data.frame("Wald.full.null.pval" = round(lin.reg.data[,1],5) ,
                            "Wald.partial.null.pval" = round(lin.reg.data[,2],5) , 
                            "score.full.null.pval" = round(lin.reg.data[,3],5),
                            "score.partial.null.pval" = round(lin.reg.data[,4],5),
                            "LRT.full.null.pval" = round(lin.reg.data[,5],5),
                            "LRT.partial.null.pval" = round(lin.reg.data[,6],5) ,
                            "Wald.full.null.disc.comparison.prob" = round(lin.reg.data[,7],5) ,
                            "Wald.partial.null.disc.comparison" = round(lin.reg.data[,8],5) ,
                            "score.full.null.disc.comparison.probs" = round(lin.reg.data[,9],5),
                            "score.partial.null.disc.comparison"= round(lin.reg.data[,10],5),
                            "KL.full.null.disc.comparison" = round(lin.reg.data[,11],5),
                            "KL.partial.null.overall.disc.comparison" = round(lin.reg.data[,12],5),
                            "KL.partial.null.para.int.overall.disc.comparison"= round(lin.reg.data[,13],5),
                            "iteration" = iteration)
  
  return(list.lin.reg.data)
  
}

linregtest <- lin.reg.fxn(iter = 10, 
                          n = 500, 
                          B = 100, 
                          true.beta0 = 0, 
                          true.beta1 = 0, 
                          true.beta2 = 0, 
                          true.beta3 = 0, 
                          true.beta4 = 0, 
                          true.sig2 = 50, 
                          mnull.sig2 = 50)




##previous working directory was: "C:/Users/bnr/Documents"
#I am saying this in case changing the working directory causes me some new unforeseen problems.
getwd()
setwd("C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs")
save(linregtest, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linregtruenull')

plotting.fxn <- function(dataset, maintitle){

par(mfrow = c(2,2), oma = c(0,0,3,0))

plot(x = dataset[[1]] ,y = dataset[[7]], main = "Wald EDCP v. Wald p-value in Full Null Setting" , ylab = "Overall Wald DCPs" 
     , xlab = "Wald p-values", xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1)

plot(x = dataset[[2]] ,y = dataset[[8]], main = "Wald EDCP v. Wald p-value in Partial Null Setting" , ylab = "Overall Wald DCPs" 
     , xlab = "Wald p-values", xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1)

plot(x = dataset[[3]] ,y = dataset[[9]], main = "Score EDCP v. Score p-value in Full Null Setting" , ylab = "Overall Score DCPs" 
     , xlab = "Score p-values", xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1)

plot(x = dataset[[4]] ,y = dataset[[10]], main = "Score EDCP v. Score p-value in Partial Null Setting" , ylab = "Para. of Int. Score DCPs" 
     , xlab = "Score p-values", xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1)

mtext(maintitle, outer = TRUE, cex = 1.5)

plot(x = dataset[[5]] ,y = dataset[[11]], main = "KL EDCP v. LRT p-value in Full Null Setting" , ylab = "Overall KL DCPs" 
     , xlab = "LRT p-values", xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1)

plot(x = dataset[[6]] ,y = dataset[[12]], main = "OKL EDCP v. LRT p-value in Partial Null Setting" , ylab = "Overall KL DCPs" 
     , xlab = "LRT p-values", xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1)

plot(x = dataset[[6]] ,y = dataset[[13]], main = "PIOKL EDCP v. LRT p-value in Partial Null Setting" , ylab = "Para. of Int. KL DCPs" 
     , xlab = "LRT p-values", xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1)

mtext(maintitle, outer = TRUE, cex = 1.5)

}

plotting.fxn(dataset = linregtest, maintitle = "Linear Regression Setting with True Null Hypothesis")

