

# log.reg.score.fxn <- function(X, y, beta){
#   #Calculating fitted probabilities based on model fit
#   p <- exp(X %*% beta) / (1 + exp(X %*% beta))
# 
#   #The score is then X^T(y - p)
#   U <- t(X) %*% (y - p)
# 
# #   return(U)
# # }
# 
# log.reg.obs.info.fxn <- function(X, y, beta){
#   p <- as.vector(exp(X %*% beta) / (1 + exp(X %*% beta)))
# 
#   W <- diag(p*(1-p), nrow = nrow(X), ncol = nrow(X))
# 
#   I <- t(X) %*% W %*% X
# 
#   return(I)
# }

# 
# #########SHOULD CREATE A SEPARATE BOOTSTRAPPING FUNCTION WHICH WE WILL CALL B TIMES TO CALCULATE B ESTIMATED DISCREPANCY COMPARISONS.  THIS SHOULD HELP ENSURE
# #THAT WE ONLY STORE THE INFORMATION FROM ONE BOOTSTRAP SAMPLE AT A TIME.
# log.reg.draw.boot.and.calc.disc.comparisons.fxn <- function(start.data, start.coefs.m123, obs.info.atbetahat.m123, inv.inv.info.atbetahat.m123, Wald.full.null.null.disc,
#                                                     Wald.partial.null.null.disc, score.full.null.null.disc, score.partial.null.null.disc,  KL.full.null.mnull,
#                                                     KL.partial.null.para.int.m1){
#   n <- nrow(start.data)
#   start.X <- matrix(cbind(as.vector(rep(1,n)),start.data[,1],start.data[,2],start.data[,3]), ncol = 4)
# 
# 
#   #By including this in a repeat statement, we ensure there will be at least 2 ones and two zeroes in the bootstrap dataset.
#   repeat{
#     boot.indices <- sample(1:n, replace = TRUE)
#     boot.data <- as.data.frame(rbind(start.data[boot.indices,]))
#     if(sum(boot.data[,4] == 1) > 1 & sum(boot.data[,4] == 0) > 1){
#       break
#     }
#   }
#   colnames(boot.data) <- c("x1.b" , "x2.b" , "x3.b", "y.b")
#   boot.X <- as.matrix(cbind(as.vector(rep(1,n)),boot.data[,1:3]))
#   start.y <- as.vector(start.data[,4])
# 
#   #this is the estimated parameter vector for the nuisance parameters under the model m1 (i.e. the null model in the partial null setting).  We need to calculate
#   #this when we are working with the overall KLD (not the parameter of interest overall KLD).  With the parameter of interest overall KLD, we plug in the pseudo-true
#   #parameter values, which when we move to the bootstrap world to estimate the discrepancy, we would plug in the original sample's MLEs for the nuisance parameters.
#   #Thus, we would not need bootstrap parameter estimates for the null model when we are using the PARAMETERS OF INTEREST overall KLD.  Note that we would also need
#   #this bootstrap parameter vector for the partial null overall score discrepancy setting (i.e. not using the parameters of interest overall score discrepancy).
#   #However, we are not going to explore the standard, overall score discrepancy in the partial null setting.
#   boot.beta.m1 <- c(glm(y.b ~ x1.b, family = binomial, data = boot.data)$coefficients,0,0)
#   #This is the alternative model's bootstrap parameter vector.  We will be using this with the overall score and KL discrepancies (NOT THE PARAMETER OF INTEREST
#   #OVERALL DISCREPANCIES)
#   boot.beta.m123 <- glm(y.b ~ x1.b + x2.b + x3.b, family = binomial, data = boot.data)$coefficients
#   #THIS IS THE ALTERNATIVE MODEL'S BOOTSTRAP PARAMETER VECTOR FOR THE NUISANCE PARAMETERS, CONDITIONED ON THE PARAMETERS OF INTEREST BEING EQUAL
#   #TO THEIR MLES UNDER THE ORIGINAL SAMPLE. tHIS IS OUR THETA.HAT.STAR.C(1)
#   #In order to get this fit, for x.2 and x.3, I first need to define an offset that is equal to beta.hat.zero + x1*beta.hat.one (This is like the fixed intercept for
#   #the model).  Thus, I need to not allow an intercept in this model fit. Note that boot.data[,1] corresponds to variable x1.
#   boot.and.orig.cond.beta.m123 <- c(start.coefs.m123[c(1,2)],
#                                     glm(y.b ~ x2.b + x3.b + offset(start.coefs.m123[1] + start.coefs.m123[2]*boot.X[,2]) - 1, family = binomial, data = boot.data)$coefficients)
# 
#   ###We now need to calculate the various bootstrap-estimated discrepancies.  Recall that for the Wald discrepancy the null model's bootstrap-estimated discrepancy is
#   #the wald statistic.  For the score discrepancy, the full null score discrepancy is the score stat and in the partial null setting we use the parameter of interest
#   #overall score discrepancy, which yields a null model's estimated discrepancy equal to the score statistic.  For the KL discrepancy, the relationship is a little more
#   #complicated.
#   Wald.disc.full.null.m123 <- t(boot.beta.m123 - start.coefs.m123) %*% obs.info.atbetahat.m123 %*% (boot.beta.m123 - start.coefs.m123)
#   Wald.disc.partial.null.m123 <- t(boot.beta.m123[c(3,4)] - start.coefs.m123[c(3,4)]) %*% inv.inv.info.atbetahat.m123 %*% (boot.beta.m123[c(3,4)] - start.coefs.m123[c(3,4)])


  # ###NOW WORKING ON THE SCORE DISCREPANCIES.  RECALL THAT U = X^t(y - p.hat) and I = X^T*W*X where W = diag(p.hat(1 - p.hat)).  p.hat = exp(X*beta) / (1 + exp(X*beta)).
  # #In this case because we use the original data X and y, and the bootstrap parameter estimates, then p.hat.boot = exp(X*beta.boot) / (1 + exp(X*beta.boot)).
  # #For full null, we use the entire bootstrap parameter vector boot.beta.m123.  For the partial null setting, we use the nuisance parameters are set equal to
  # #their MLEs under the original sample (start.coefs.m123[c(1,2)]) and their conditional bootstrap-estimated parameter of interest.  Thus, for the partial null we
  # #use boot.cond.nuis.para.beta.m123 for the parameter vector
  # #p.hat.boot.full.null.m123 <- as.vector(exp(start.X %*% boot.beta.m123) / (1 + exp(start.X %*% boot.beta.m123)))
  # #p.hat.boot.partial.null.m123 <- as.vector(exp(start.X %*% boot.and.orig.cond.beta.m123) / (1 + exp(start.X %*% boot.and.orig.cond.beta.m123)))
  # 
  # #Note that (start.data[,4]) is y.  The scores are based on the original data and the boot estimates, so we use X and y in our expressions, but p.hat incorporates
  # #the bootstrap parameter estimates.
  # U.full.null.m123 <- log.reg.score.fxn(X = start.X, y = start.y, beta = boot.beta.m123)
  # U.partial.null.m123 <- log.reg.score.fxn(X = start.X, y = start.y, beta =  boot.and.orig.cond.beta.m123)[c(3,4)]
  # 
  # #The informations are based the original data and bootstrap parameter estimates so we use X and the p.hat (which incorporates the bootstrap parameter estimates)
  # obs.info.atbetahatstar.m123 <- log.reg.obs.info.fxn(X = start.X, y = start.y, beta = boot.beta.m123)
  # #The information for the partial null should only include be calculated for the parameters of interest.  We then use [I_11]^-1.  This will be evaluated at the
  # #boot.cond.nuis.para.beta.m123, which is the MLES of the nuisance under the original sample and the paras of interest under the bootstrap, conditioned on nuisance MLES
  # obs.info.para.int.atbetahatstarC.m123 <- log.reg.obs.info.fxn(X = start.X, y = start.y, beta = boot.and.orig.cond.beta.m123)[c(3,4),c(3,4)]
  # 
  # ##Need inverse infos
  # inv.obs.info.atbetahatstar.m123 <- solve(obs.info.atbetahatstar.m123)
  # inv.obs.info.para.int.atbetahatstarC.m123 <- solve(obs.info.para.int.atbetahatstarC.m123)
  # 
  # ##Now calculate the alternative models' score discrepancies
  # score.disc.full.null.m123 <- t(U.full.null.m123) %*% inv.obs.info.atbetahatstar.m123 %*% U.full.null.m123
  # score.disc.partial.null.m123 <- t(U.partial.null.m123) %*% inv.obs.info.para.int.atbetahatstarC.m123 %*% U.partial.null.m123

  # ##WORKING ON THE KL DISCREPANCY CALCULATIONS NOW
  # ##We calculate the eta values to make the discrepancies easier to write.  The alternative model's full and partial null OVERALL KL discrepancies are equivalent.
  # #Under the PARAMETER OF INTEREST OVERALL KLD we have to use a different parameter vector (theta.hat.star.C1, theta.hat.2)
  # eta.boot.and.orig.overall.KL.m123 <- start.X %*% boot.beta.m123
  # eta.boot.and.orig.para.int.overall.KL.m123 <- start.X %*% boot.and.orig.cond.beta.m123
  # 
  # KL.overall.disc.m123 <- 2*(-sum(start.data[,4]*eta.boot.and.orig.overall.KL.m123) + sum(log(1 + exp(eta.boot.and.orig.overall.KL.m123))))
  # KL.para.int.overall.disc.m123 <- 2*(-sum(start.data[,4]*eta.boot.and.orig.para.int.overall.KL.m123) + sum(log(1 + exp(eta.boot.and.orig.para.int.overall.KL.m123))))
  # 
  # #FOR THE PARTIAL NULL SETTING USING THE OVERALL KLD, THE NULL MODEL ACTUALLY DEPENDS ON THE BOOTSTRAP SAMPLE, SO WE CALCULATE THAT NOW
  # eta.boot.and.orig.overall.KL.m1 <- start.X %*% boot.beta.m1
  # partial.null.KL.overall.disc.m1 <- 2*(-sum(start.data[,4]*eta.boot.and.orig.overall.KL.m1) + sum(log(1 + exp(eta.boot.and.orig.overall.KL.m1))))
  # 







  ##NOW CALCULATING THE DISCREPANCY COMPARISONS. Recall that for the Wald and score discrepancies the null models' discrepancy are actually the respective test statistics.
  ##THIS IS NOT TRUE FOR THE KL DISCREPANCIES.
#   Wald.full.null.disc.comparison <- ifelse(Wald.full.null.null.disc < Wald.disc.full.null.m123, 1, 0)
#   Wald.partial.null.disc.comparison <- ifelse(Wald.partial.null.null.disc < Wald.disc.partial.null.m123, 1, 0)
# 
#   score.full.null.disc.comparison <- ifelse(score.full.null.null.disc < score.disc.full.null.m123, 1, 0)
#   score.partial.null.disc.comparison <- ifelse(score.partial.null.null.disc < score.disc.partial.null.m123, 1, 0)
# 
#   KL.full.null.disc.comparison <- ifelse(KL.full.null.mnull < KL.overall.disc.m123, 1, 0)
#   KL.partial.null.overall.disc.comparison <- ifelse(partial.null.KL.overall.disc.m1 < KL.overall.disc.m123, 1, 0)
#   KL.partial.null.para.int.overall.disc.comparison <- ifelse(KL.partial.null.para.int.m1 < KL.para.int.overall.disc.m123, 1, 0)
# 
#   return(c(Wald.full.null.disc.comparison,Wald.partial.null.disc.comparison,score.full.null.disc.comparison,score.partial.null.disc.comparison,
#            KL.full.null.disc.comparison, KL.partial.null.overall.disc.comparison, KL.partial.null.para.int.overall.disc.comparison))
# 
# }
# 
# 
# 
# 
# 
# ##This function will do all a completey cycle of calculating logistic regression p-values and bootstrap-estimated discrepancy comparison probabilities.  Our larger
# #function will perform this process for iter iterations.
# 
# ###GOING TO GO INSIDE THE FUNCTION FOR NOW SO WE WILL COMMENT OUT THE FUNCTION NOTATION
# single.iter.log.reg.fxn <- function(n =n, B = B, true.beta0 = true.beta0, true.beta1 = true.beta1, true.beta2 = true.beta2, true.beta3 = true.beta3
#                                     , true.beta4 = true.beta4){
# 
# 
# ##################BEE SURE TO GO BACK AND DELETE THESE#####################################
#   #n = 100
#   ##B = 100
#   #true.beta0 = 0
#   #true.beta1 = 0
#   #true.beta2 = 0
#   #true.beta3 = 0
#   #true.beta4 = 0
# 
#   ##Code if we want all predictors be independent of each other
#     #x1.success.prob <- 0.6
#     #x2.success.prob <- 0.4
#     #x3.min <- -0.8
#     #x3.max <- 0
#     #x4.min <- 0
#     #x4.max <- 2
# 
# 
#     #x1 <- rbinom(n = n, size = 1, prob = x1.success.prob)
#     #x2 <- rbinom(n = n, size = 1, prob = x2.success.prob)
#     #x3 <- runif(n = n, min = x3.min, max = x3.max)
#     #x4 <- runif(n = n, min = x4.min, max = x4.max)
# 
#   ###Code if we want to have the predictors have some correlation.
#   corr <- 0.8
#   var <- 100
#   cov <- 100*0.8*0.8
#   sig.mat <- matrix(c(rep(c(100,rep(cov,4)),3),100), nrow = 4, ncol = 4)
#   X.includex4 <- MASS::mvrnorm(n = n, mu = c(2,2,-2,-2), Sigma = sig.mat)
# 
#   x1 <- X.includex4[,1]
#   x2 <- X.includex4[,2]
#   x3 <- X.includex4[,3]
#   x4 <- X.includex4[,4]
# 
#     ##Calculating the original (sample from which we draw the bootstrap) eta values using the true betas and the observed x values
#     orig.eta <- true.beta0 + x1*true.beta1 + x2*true.beta2 + x3*true.beta3 + x4* true.beta4
# 
# 
# 
#   ##Drawing a sample of y values for the ith iteration.  We are forcing there to be at least 5 zeroes and 5 ones in the original sample
#     #By including this in a repeat statement, we ensure there will be at least one 1 and one 0 in the bootstrap dataset.
#     repeat{
#       y <- rbinom(n = n, size = 1, prob = exp(orig.eta) / (1 + exp(orig.eta)) )
#       if(sum(y == 1) > 4 & sum(y == 0) > 4){
#         break
#       }
#     }
# 
# 
# 
# 
# 
# 
#     ###Storing all the data from the ith iteration
#     full.orig.data <- matrix(cbind(x1,x2,x3,y), ncol = 4)
# 
# 
# 
# 
#     x1 <- full.orig.data[,1]
#     x2 <- full.orig.data[,2]
#     x3 <- full.orig.data[,3]
#     y <- full.orig.data[,4]
# 
#     orig.X <- matrix(cbind(1,x1,x2,x3), ncol = 4)
# 
# 
#     ###Storing the coefficients of the fitted models for the original data
#     orig.coef.mnull <- rep(0,4)
#     orig.coef.m1 <- c(glm(y ~ x1, family = binomial)$coefficients, rep(0,2))
#     orig.coef.m123 <- glm(y ~ x1 + x2 + x3, family = binomial)$coefficients
# 
#     ##calculating an orig.eta value for each model and for each of iter iterations.  This is the fitted value of eta for each of the models.  Note that
#     #because the null model has 0 for all betas, then its fitted eta vector will simply be n zeroes.
#     orig.eta.mnull <- rep(0,n)
#     orig.eta.m1 <- orig.X %*% orig.coef.m1
#     orig.eta.m123 <- orig.X %*% orig.coef.m123
# 
#     #Calculating the LRT test stats and p-value##
# 
#     LRT.stat.m123vmnull <- 2 * ( sum(y*(orig.eta.m123 - orig.eta.mnull)) + sum(log((1 + exp(orig.eta.mnull)) / (1 + exp(orig.eta.m123)))))
#     LRT.pval.m123vmnull <- 1 - pchisq(LRT.stat.m123vmnull, df = 4 )
#     #We also want the null model's discrepancy for the full null setting, where the null model pre-specifies all values to be zero.
#     KLD.full.null.mnull <- 2*(-sum(y*orig.eta.mnull) + sum(log(1 + exp(orig.eta.mnull))))
# 
#     LRT.stat.m123vm1 <- 2 * ( sum(y*(orig.eta.m123 - orig.eta.m1)) + sum(log((1 + exp(orig.eta.m1)) / (1 + exp(orig.eta.m123)))))
#     LRT.pval.m123vm1 <- 1 - pchisq(LRT.stat.m123vm1, df = 2 )
#     ##We want the null model's discrepancy for the partial null setting where we are using the para of interest overall disc, so that the nuisance parameters (beta0
#     # and beta1) are estimated with the original sample, but the parameters of interest are set equal to zero.
#     KLD.partial.null.para.int.m1 <- 2*(-sum(y*orig.eta.m1) + sum(log(1 + exp(orig.eta.m1))))
# 
#     ###Working on Wald test p-values
#     ##This is giving the observed information matrix for each model that we will need.  For the Wald test we need the information matrix for whatever
#     ##the larger model is in the given hypothesis test.  we are always testing against model m123, so that is the only one we need.
# 
#     #This mlog.fourbetasv2 function is giving the NEGATIVE log likelihood of the logistic regression model, where we are sure to have the betas in the function.
#     #We know that logL = sum(y_i*eta_i) - sum(log(1 + exp(eta_i))) and eta_i = x_i*beta_i.  Also, we give it the negative of the log-likelihood because the nlm
#     #function MINIMIZES rather than maximizes.
# 
# 
# 
# 
# 
# 
#     #######################OLD WAY OF CALCULATING THE INFORMATION MATRIX
#     #mlog.fourbetasv2 <- function(beta, x, y){
#     #  -1*sum(y* (x %*% beta)) + sum(log(1 + exp(x %*% beta)))
#     #}
#     #This is the observd inforation evaluated at the MLEs beta.hat
#     #obs.info.m123.v2 <- nlm(f = mlog.fourbetasv2, p = orig.coef.m123,  x = orig.X, y = y , hessian = TRUE)$hessian
# 
#     ##This is giving the inverse of the information matrix
#     #inv.info.m123 <- solve(obs.info.m123)
# 
# 
#     #TRYING A NEW WAY OF CALCULATING THE INFORMATION MATRIX
#     p.m123 <- as.vector(exp(orig.X %*% orig.coef.m123) / (1 + exp(orig.X %*% orig.coef.m123)))
#     W.m123 <- diag(p.m123*(1 - p.m123), nrow = n, ncol =n)
#     obs.info.m123 <- t(orig.X) %*% W.m123 %*% orig.X
#     ###THis is equivalent to I^{-1}(beta.hat).  This is what we will be using for the Wald statistic because it uses the observed info evaluated at the MLEs.
#     #The score statistic uses the information evaluated at the pre-defined values.  Thus for the full null model, we evaulate at the pre-specified beta values (i.e. 0).  For the
#     #partial null model we evaluate info at beta2 = beta3 = 0 and let beta0 and beta1 equal their MLEs assuming beta2 and beta3 are zero.
#     inv.info.m123 <- solve(obs.info.m123)
# 
# 
# 
# 
# 
# 
# 
# 
#     ##These are giving the inverse of the inverse info where we only take the entries that apply to parameters being tested for the given hypothesis test.
#     ##We only need to take the inverse of the elements which correspond to the parameters being tested.
#     inv.inv.In11.m123vm1 <- solve(inv.info.m123[c(3:4),c(3:4)])
# 
#     #CALCULATING THE WALD TEST P-VALUES
#     #Recall that for the full null setting (i.e. m123 vs mnull), the middle piece of the Wald stat is actually the information (or we could think of it as taking the
#     #inverse information, then inverting the piece corresponding to the parameters of interest (i.e. all the parameters).  Thus, we would be inverting the
#     #entire matrix and then inverting the entire matrix again.
#     Wald.stat.m123vmnull <- t(orig.coef.m123) %*% obs.info.m123 %*% orig.coef.m123
#     Wald.pval.m123vmnull <- 1 - pchisq(Wald.stat.m123vmnull, df = 4)
# 
#     Wald.stat.m123vm1 <- t(orig.coef.m123[c(3,4)]) %*% inv.inv.In11.m123vm1 %*% orig.coef.m123[c(3,4)]
#     Wald.pval.m123vm1 <- 1 - pchisq(Wald.stat.m123vm1, df = 2)
# 
#     #######CALCULATING THE SCORE TEST STATISTIC AND ITS P-VALUE
#     #we know that to calculate the score statistic we need to determine what the score and information are.  We then evaluate the score and information at the
#     #pre-specified values of beta (in the full null setting) or at the pre-specified values of the parameters of interest and the conditional MLEs for the
#     #nuisance parameters.
#     #Let p(x_i,beta)_i be the probability of success given the observation's set of covariates and for a specific beta.  The beta we will use for the score
#     #is the pre-specified values.
#     #By doing some of my own math and looking up references,
#     #the score is equal to X^T(y - p), where p is the vector of estimated probabilities
#     #the observed information is X^T*W*X where W is an n x n diagonal matrix with p.hat*(1 - p.hat) along the diagonal.
#     #RECALL THAT FOR THE SCORE EVERYTHING NEEDS TO BE EVALUATED AT THE PRE-SPECIFIED VALUES THETA_0 NOT THETA_HAT!!!!!!
# 
#     #Defining the observed probabilities for a patient and for the pre-specified beta values.  Also, I defined each of the coefficients (i.e. beta hats) so that
#     #they have length 4, placing zeroes for the betas set equal to zero for the non-full models.
#     ###obs.info.m123 is equivalent to I^{-1}(beta.hat).  This is what we will be using for the Wald statistic because it uses the observed info evaluated at the MLEs.
#     #The score statistic uses the information evaluated at the pre-defined values.  Thus for the full null model, we evaulate at the pre-specified beta values (i.e. 0).  For the
#     #partial null model we evaluate info at beta2 = beta3 = 0 and let beta0 and beta1 equal their MLEs assuming beta2 and beta3 are zero.
# 
# 
#     U.m123vmnull <- log.reg.score.fxn(X = orig.X, y = y, beta = orig.coef.mnull)
#     U.m123vm1 <- log.reg.score.fxn(X = orig.X, y = y, beta = orig.coef.m1)
# 
#     obs.info.atbeta0.mnull <- log.reg.obs.info.fxn(X = orig.X, y = y, beta = orig.coef.mnull)
#     obs.info.atbeta0.m1 <- log.reg.obs.info.fxn(X = orig.X, y = y, beta = orig.coef.m1)
# 
#     ##For the full null, we only need the inverse information for the score statistic.  For the partial null, we need take the inverse of the information and then
#     #just grab the piece of the inverse information corresponding to the paramters of interest
#     inv.obs.info.atbeta0.mnull <- solve(obs.info.atbeta0.mnull)
#     inv.obs.info.atbeta0.m1 <- (solve(obs.info.atbeta0.m1))[c(3,4), c(3,4)]
# 
#     score.stat.m123vmnull <- t(U.m123vmnull) %*% inv.obs.info.atbeta0.mnull %*% U.m123vmnull
#     score.stat.m123vm1   <- t(U.m123vm1[c(3,4)]) %*% inv.obs.info.atbeta0.m1 %*% U.m123vm1[c(3,4)]
# 
#     score.pval.m123vmnull <- 1 - pchisq(score.stat.m123vmnull, df = 4)
#     score.pval.m123vm1 <- 1 - pchisq(score.stat.m123vm1, df = 2)
# 
#     #return(c(LRT.pval.m123vmnull, LRT.pval.m123vm1, Wald.pval.m123vmnull, Wald.pval.m123vm1, score.pval.m123vmnull, score.pval.m123vm1))
# 
#     ##THIS IS NOW ASKING THE FUNCTION WHICH DOES THE DRAWING OF THE BOOTSTRAP SAMPLE AND THE CALCULATING OF THE BOOTSTRAP DISCREPANCIES AND THEIR COMPARISON PROBS
#     #TO BE DONE B TIMES.  BY CALLING THE FUNCTION B TIMES I BELIEVE IT WILL ONLY HAVE TO STORE THE BOOTSTRAP SAMPLE ON WHICH IT IS CURRENTLY WORKING.
#     boot.disc.comparisons <- matrix(data = NA, nrow = B, ncol = 7)
#     for(b in 1:B){
#       boot.disc.comparisons[b,] <-
#         log.reg.draw.boot.and.calc.disc.comparisons.fxn(
#           start.data = full.orig.data,
#           start.coefs.m123 = orig.coef.m123,
#           obs.info.atbetahat.m123 = obs.info.m123,
#           inv.inv.info.atbetahat.m123 = inv.inv.In11.m123vm1,
#           Wald.full.null.null.disc = Wald.stat.m123vmnull,
#           Wald.partial.null.null.disc = Wald.stat.m123vm1,
#           score.full.null.null.disc = score.stat.m123vmnull,
#           score.partial.null.null.disc = score.stat.m123vm1,
#           KL.full.null.mnull = KLD.full.null.mnull,
#           KL.partial.null.para.int.m1 = KLD.partial.null.para.int.m1)
#     }
# 
# 
# 
#     boot.disc.comparison.probs <- apply(X = boot.disc.comparisons, 2, FUN = sum) / B
# 
# 
# 
# return(c(Wald.pval.m123vmnull,
#          Wald.pval.m123vm1,
#          score.pval.m123vmnull,
#          score.pval.m123vm1,
#          LRT.pval.m123vmnull,
#          LRT.pval.m123vm1,
#          boot.disc.comparison.probs[1],
#          boot.disc.comparison.probs[2],
#          boot.disc.comparison.probs[3],
#          boot.disc.comparison.probs[4],
#          boot.disc.comparison.probs[5],
#          boot.disc.comparison.probs[6],
#          boot.disc.comparison.probs[7]))
# }



# do_log_reg_sim <- function(iter,
#                            n,
#                            B,
#                            true_beta_vec){
# 
# true.beta0 <- true_beta_vec[[1]]
# true.beta1 <- true_beta_vec[[2]]
# true.beta2 <- true_beta_vec[[3]]
# true.beta3 <- true_beta_vec[[4]]
# true.beta4 <- true_beta_vec[[5]]
# 
#   log.reg.data <- matrix(data = NA, nrow = iter, ncol = 13)
#   for(i in 1:iter){
#     log.reg.data[i,] <- single.iter.log.reg.fxn(n = n, B = B, true.beta0 = true.beta0, true.beta1 = true.beta1, true.beta2 = true.beta2, true.beta3 = true.beta3
#                                                 , true.beta4 = true.beta4)
#   }
# 
#   iteration <- 1:iter
# 
# data.frame(
#     "Wald_full_null_pval" = round(log.reg.data[,1],5) ,
#     "Wald_partial_null_pval" = round(log.reg.data[,2],5) ,
#     "score_full_null_pval" = round(log.reg.data[,3],5),
#     "score_partial_null_pval" = round(log.reg.data[,4],5),
#     "LRT_full_null_pval" = round(log.reg.data[,5],5),
#     "LRT_partial_null_pval" = round(log.reg.data[,6],5) ,
#     "Wald_full_null_bdcp" = round(log.reg.data[,7],5) ,
#     "Wald_partial_null_bdcp" = round(log.reg.data[,8],5) ,
#     "score_full_null_bdcp" = round(log.reg.data[,9],5),
#     "score_partial_null_bdcp"= round(log.reg.data[,10],5),
#     "KL_full_null_bdcp" = round(log.reg.data[,11],5),
#     "KL_partial_null_overall_bdcp" = round(log.reg.data[,12],5),
#     "KL_partial_null_para_int_overall_bdcp"=
#       round(log.reg.data[,13],5),
#     "iteration" = iteration)
# 
# }
# 
# logregtest <- log.reg.fxn(iter = 25, n = 500, B = 500, true.beta0 = 0, true.beta1 = 0 , true.beta2 = 0, true.beta3 = 0, true.beta4 = 0)
# plotting.fxn(dataset = logregtest, maintitle = "log reg true null n500")
# 
# 
# logregtest2 <- log.reg.fxn(iter = 50, n = 1000, B = 200, true.beta0 = 0, true.beta1 = 0 , true.beta2 = 0, true.beta3 = 0, true.beta4 = 0)
# plotting.fxn(dataset = logregtest2, maintitle = "log reg bigger n true null")
