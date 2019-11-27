###Code to run when I leave 
# 
# setwd("C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs")
# 
# linregtruenulln50 <- lin.reg.fxn(iter = 50,
#                                  n = 50,
#                                  B = 10000,
#                                  true.beta0 = 0,
#                                  true.beta1 = 0, 
#                                  true.beta2 = 0,
#                                  true.beta3 = 0, 
#                                  true.beta4 = 0, 
#                                  true.sig2 = 50, 
#                                  mnull.sig2 = 50)
# save(linregtruenull, 
#      file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linregtruenulln50')
# plotting.fxn.0to1(dataset = linregtruenulln50, maintitle = "Linear Regression True Null n=50")
# 
# linregtruenulln100 <- lin.reg.fxn(iter = 50, n = 100, B = 10000, true.beta0 = 0, true.beta1 = 0, true.beta2 = 0, true.beta3 = 0, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 50)
# save(linregtruenulln100, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linregtruenulln100')
# 
# 
# linregtruenulln250 <- lin.reg.fxn(iter = 50, n = 250, B = 10000, true.beta0 = 0, true.beta1 = 0, true.beta2 = 0, true.beta3 = 0, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 50)
# save(linregtruenulln250, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linregtruenulln250')
# 
# 
# logregtruenulln100 <- log.reg.fxn(iter = 50, n = 100, B = 10000, true.beta0 = 0, true.beta1 = 0, true.beta2 = 0, true.beta3 = 0, true.beta4 = 0)
# save(logregtruenulln100, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logregtruenulln100')
# plotting.fxn.0to1(dataset = logregtruenulln100, maintitle = "Logistic Regression True Null n=100")
# 
# logregtruenulln250 <- log.reg.fxn(iter = 50, n = 250, B = 10000, true.beta0 = 0, true.beta1 = 0, true.beta2 = 0, true.beta3 = 0, true.beta4 = 0)
# save(logregtruenulln250, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logregtruenulln250')
# plotting.fxn.0to1(dataset = logregtruenulln250, maintitle = "Logistic Regression True Null n=250")
# 
# linregtruenulln500 <- lin.reg.fxn(iter = 50, n = 500, B = 10000, true.beta0 = 0, true.beta1 = 0, true.beta2 = 0, true.beta3 = 0, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 50)
# save(linregtruenulln500, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linregtruenulln500')
# #plotting.fxn.0to1(dataset = linregtruenulln250, maintitle = "Linear Regression True Null n=500")
# 
# logregtruenulln500 <- log.reg.fxn(iter = 50, n = 500, B = 10000, true.beta0 = 0, true.beta1 = 0, true.beta2 = 0, true.beta3 = 0, true.beta4 = 0)
# save(logregtruenulln500, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logregtruenulln500')
# #plotting.fxn.0to1(dataset = logregtruenulln250, maintitle = "Logistic Regression True Null n=500")
# logregtruenulln500
# 
# plotting.fxn(dataset = linregtruenulln100, maintitle = "hello")
# plotting.fxn(dataset = logregtruenulln250, maintitle = "hello")
# 
# ###############################LINEAR REGRESSION SLIGHYLY FALSE NULL SETTINGS############################################################
# linreglightn50beta0.2sig61 <- lin.reg.fxn(iter = 50, n = 50, B = 10000, true.beta0 = 0.20, true.beta1 = -0.20, true.beta2 = 0.20, true.beta3 = -0.20, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 61)
# save(linreglightn50beta0.2sig61, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linreglightn50beta0.2sig61')
# plotting.fxn(dataset = linreglightn50beta0.2sig61, maintitle = "n=50 beta2 = 0.2 sig2 = 61 Linear Regression")
# plotting.fxn.0to0.5(dataset = linreglightn50beta0.2sig61, maintitle = "n=50 beta2 = 0.2 sig2 = 61 Linear Regression")
# 
# linreglightn100beta0.175sig60 <- lin.reg.fxn(iter = 50, n = 100, B = 10000, true.beta0 = 0.175, true.beta1 = -0.175, true.beta2 = 0.175, true.beta3 = -0.175, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 60)
# save(linreglightn100beta0.175sig60, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linreglightn100beta0.175sig60')
# plotting.fxn(dataset = linreglightn100beta0.175sig60, maintitle = "n=100 beta2 = 0.175 sig2 = 60 Linear Regression")
# plotting.fxn.0to0.5(dataset = linreglightn100beta0.175sig60, maintitle = "n=100 beta2 = 0.175 sig2 = 60 Linear Regression")
# 
# linreglightn150beta0.15sig59 <- lin.reg.fxn(iter = 50, n = 150, B = 10000, true.beta0 = 0.15, true.beta1 = -0.15, true.beta2 = 0.15, true.beta3 = -0.15, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 59)
# save(linreglightn150beta0.15sig59, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linreglightn150beta0.15sig59')
# plotting.fxn(dataset = linreglightn150beta0.15sig59, maintitle = "n=150 beta2 = 0.15 sig2 = 59 Linear Regression")
# plotting.fxn.0to0.5(dataset = linreglightn150beta0.15sig59, maintitle = "n=150 beta2 = 0.15 sig2 = 59 Linear Regression")
# 
# linreglightn200beta0.125sig58 <- lin.reg.fxn(iter = 50, n = 200, B = 10000, true.beta0 = 0.125, true.beta1 = -0.125, true.beta2 = 0.125, true.beta3 = -0.125, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 58)
# save(linreglightn200beta0.125sig58, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linreglightn200beta0.125sig58')
# plotting.fxn(dataset = linreglightn200beta0.125sig58, maintitle = "n=200 beta2 = 0.125 sig2 = 58 Linear Regression")
# plotting.fxn.0to0.5(dataset = linreglightn200beta0.125sig58, maintitle = "n=200 beta2 = 0.125 sig2 = 58 Linear Regression")
# 
# 
# linreglightn250beta0.1sig57 <- lin.reg.fxn(iter = 50, n = 250, B = 10000, true.beta0 = 0.10, true.beta1 = -0.10, true.beta2 = 0.10, true.beta3 = -0.10, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 57)
# save(linreglightn250beta0.1sig57, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linreglightn250beta0.1sig57')
# plotting.fxn(dataset = linreglightn250beta0.1sig57, maintitle = "n=250 beta2 = 0.10 sig2 = 57 Linear Regression")
# plotting.fxn.0to0.5(dataset = linreglightn250beta0.1sig57, maintitle = "n=250 beta2 = 0.10 sig2 = 57 Linear Regression")
# 
# linreglightn500beta0.075sig56 <- lin.reg.fxn(iter = 50, n = 500, B = 10000, true.beta0 = 0.075, true.beta1 = -0.075, true.beta2 = 0.075, true.beta3 = -0.075, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 56)
# save(linreglightn500beta0.075sig56 , file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linreglightn500beta0.075sig56')
# plotting.fxn(dataset = linreglightn500beta0.075sig56, maintitle = "n=500 beta2 = 0.075 sig2 = 56 Linear Regression")
# plotting.fxn.0to0.5(dataset = linreglightn500beta0.075sig56, maintitle = "n=500 beta2 = 0.075 sig2 = 56 Linear Regression")
# 
# linreglightn1000beta0.05sig55 <- lin.reg.fxn(iter = 50, n = 1000, B = 10000, true.beta0 = 0.05, true.beta1 = -0.05, true.beta2 = 0.05, true.beta3 = -0.05, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 55)
# save(linreglightn1000beta0.05sig55 , file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linreglightn1000beta0.05sig55')
# plotting.fxn(dataset = linreglightn1000beta0.05sig55, maintitle = "n=1000 beta2 = 0.05 sig2 = 55 Linear Regression")
# plotting.fxn.0to0.5(dataset = linreglightn1000beta0.05sig55, maintitle = "n=1000 beta2 = 0.05 sig2 = 55 Linear Regression")
# 
# 
# ################################LOGISTIC REGRESSION SLIGHTLY FALSE NULL SETTING################################
# 
# 
# logreglightn100beta0.07 <- log.reg.fxn(iter = 50, n = 100, B = 10000, true.beta0 = 0.07, true.beta1 = -0.07, true.beta2 = 0.07, true.beta3 = -0.07, true.beta4 = 0)
# save(logreglightn100beta0.07 , file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logreglightn100beta0.07')
# plotting.fxn(dataset = logreglightn100beta0.07, maintitle = "n=100 beta2 = 0.07 Logistic Regression")
# plotting.fxn.0to0.5(dataset = logreglightn100beta0.07, maintitle = "n=100 beta2 = 0.07 Logistic Regression")
# 
# logreglightn150beta0.06 <- log.reg.fxn(iter = 50, n = 150, B = 10000, true.beta0 = 0.06, true.beta1 = -0.06, true.beta2 = 0.06, true.beta3 = -.06, true.beta4 = 0)
# save(logreglightn150beta0.06 , file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logreglightn150beta0.06')
# plotting.fxn(dataset = logreglightn150beta0.06, maintitle = "n=150 beta2 = 0.06 Logistic Regression")
# plotting.fxn.0to0.5(dataset = logreglightn150beta0.06, maintitle = "n=150 beta2 = 0.06 Logistic Regression")
# 
# logreglightn200beta0.05 <- log.reg.fxn(iter = 50, n = 200, B = 10000, true.beta0 = 0.05, true.beta1 = -0.05, true.beta2 = 0.05, true.beta3 = -.05, true.beta4 = 0)
# save(logreglightn200beta0.05 , file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logreglightn200beta0.05')
# plotting.fxn(dataset = logreglightn200beta0.05, maintitle = "n=200 beta2 = 0.05 Logistic Regression")
# plotting.fxn.0to0.5(dataset = logreglightn200beta0.05, maintitle = "n=200 beta2 = 0.05 Logistic Regression")
# 
# 
# logreglightn250beta0.04 <- log.reg.fxn(iter = 50, n = 250, B = 10000, true.beta0 = 0.04, true.beta1 = -0.04, true.beta2 = 0.04, true.beta3 = -.04, true.beta4 = 0)
# save(logreglightn250beta0.04 , file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logreglightn250beta0.04')
# plotting.fxn(dataset = logreglightn250beta0.04, maintitle = "n=250 beta2 = 0.04 Logistic Regression")
# plotting.fxn.0to0.5(dataset = logreglightn250beta0.04, maintitle = "n=250 beta2 = 0.04 Logistic Regression")
# 
# logreglightn500beta0.02 <- log.reg.fxn(iter = 50, n = 500, B = 10000, true.beta0 = 0.02, true.beta1 = -0.02, true.beta2 = 0.02, true.beta3 = -.02, true.beta4 = 0)
# save(logreglightn500beta0.02 , file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logreglightn500beta0.02')
# plotting.fxn(dataset = logreglightn500beta0.02, maintitle = "n=500 beta2 = 0.02 Logistic Regression")
# plotting.fxn.0to0.5(dataset = logreglightn500beta0.02, maintitle = "n=500 beta2 = 0.02 Logistic Regression")
# 
# logreglightn1000beta0.01 <- log.reg.fxn(iter = 50, n = 1000, B = 10000, true.beta0 = 0.01, true.beta1 = -0.01, true.beta2 = 0.01, true.beta3 = -.01, true.beta4 = 0)
# save(logreglightn1000beta0.01 , file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logreglightn1000beta0.01')
# plotting.fxn(dataset = logreglightn1000beta0.01, maintitle = "n=1000 beta2 = 0.01 Logistic Regression")
# plotting.fxn.0to0.5(dataset = logreglightn1000beta0.01, maintitle = "n=1000 beta2 = 0.01 Logistic Regression")
# 
# ###############################LINEAR REGRESSION STRONGLY FALSE FALSE NULL SETTINGS############################################################
# linregstrongn50beta0is0beta2is0.6sig52 <- lin.reg.fxn(iter = 50, n = 50, B = 20000, true.beta0 = 0, true.beta1 = -0, true.beta2 = 0.60, true.beta3 = -0.60, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 52)
# save(linregstrongn50beta0is0beta2is0.6sig52, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linregstrongn50beta0is0beta2is0.6sig52')
# plotting.fxn(dataset = linregstrongn50beta0is0beta2is0.6sig52, maintitle = "n=50 b0 = 0 b2 = 0.6 sig2 = 52 Linear Regression")
# plotting.fxn.0to0.1(dataset = linregstrongn50beta0is0beta2is0.6sig52, maintitle = "n=50 b0 = 0 b2 = 0.6 sig2 = 52 Linear Regression")
# 
# linregstrongn100beta0is0beta2is0.4sig52 <- lin.reg.fxn(iter = 50, n = 100, B = 20000, true.beta0 = 0, true.beta1 = -0, true.beta2 = 0.40, true.beta3 = -0.40, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 52)
# save(linregstrongn100beta0is0beta2is0.4sig52, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linregstrongn100beta0is0beta2is0.4sig52')
# plotting.fxn(dataset = linregstrongn100beta0is0beta2is0.4sig52, maintitle = "n=100 b0 = 0 b2 = 0.40 sig2 = 52 Linear Regression")
# plotting.fxn.0to0.1(dataset = linregstrongn100beta0is0beta2is0.4sig52, maintitle = "n=100 b0 = 0 b2 = 0.40 sig2 = 52 Linear Regression")
# 
# linregstrongn150beta0is0beta2is0.32sig52 <- lin.reg.fxn(iter = 50, n = 150, B = 20000, true.beta0 = 0, true.beta1 = -0, true.beta2 = 0.32, true.beta3 = -0.32, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 52)
# save(linregstrongn150beta0is0beta2is0.32sig52, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linregstrongn150beta0is0beta2is0.32sig52')
# plotting.fxn(dataset = linregstrongn150beta0is0beta2is0.32sig52, maintitle = "n=150 b0 = 0 b2 = 0.32 sig2 = 52 Linear Regression")
# plotting.fxn.0to0.1(dataset = linregstrongn150beta0is0beta2is0.32sig52, maintitle = "n=150 b0 = 0 b2 = 0.32 sig2 = 52 Linear Regression")
# 
# linregstrongn200beta0is0beta2is0.27sig52 <- lin.reg.fxn(iter = 50, n = 200, B = 20000, true.beta0 = 0, true.beta1 = -0, true.beta2 = 0.27, true.beta3 = -0.27, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 52)
# save(linregstrongn200beta0is0beta2is0.27sig52, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linregstrongn200beta0is0beta2is0.27sig52')
# plotting.fxn(dataset = linregstrongn200beta0is0beta2is0.27sig52, maintitle = "n=200 b0 = 0 b2 = 0.27 sig2 = 52 Linear Regression")
# plotting.fxn.0to0.1(dataset = linregstrongn200beta0is0beta2is0.27sig52, maintitle = "n=200 b0 = 0 b2 = 0.27 sig2 = 52 Linear Regression")
# 
# linregstrongn250beta0is0beta2is0.22sig52 <- lin.reg.fxn(iter = 50, n = 250, B = 20000, true.beta0 = 0, true.beta1 = -0, true.beta2 = 0.22, true.beta3 = -0.22, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 52)
# save(linregstrongn250beta0is0beta2is0.22sig52, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linregstrongn250beta0is0beta2is0.22sig52')
# plotting.fxn(dataset = linregstrongn250beta0is0beta2is0.22sig52, maintitle = "n=250 b0 = 0 b2 = 0.22 sig2 = 52 Linear Regression")
# plotting.fxn.0to0.1(dataset = linregstrongn250beta0is0beta2is0.22sig52, maintitle = "n=250 b0 = 0 b2 = 0.22 sig2 = 52 Linear Regression")
# 
# linregstrongn500beta0is0beta2is0.17sig52 <- lin.reg.fxn(iter = 50, n = 500, B = 20000, true.beta0 = 0, true.beta1 = -0, true.beta2 = 0.17, true.beta3 = -0.17, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 52)
# save(linregstrongn500beta0is0beta2is0.17sig52, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linregstrongn500beta0is0beta2is0.17sig52')
# plotting.fxn(dataset = linregstrongn500beta0is0beta2is0.17sig52, maintitle = "n=500 b0 = 0 b2 = 0.17 sig2 = 52 Linear Regression")
# plotting.fxn.0to0.1(dataset = linregstrongn500beta0is0beta2is0.17sig52, maintitle = "n=500 b0 = 0 b2 = 0.17 sig2 = 52 Linear Regression")
# 
# linregstrongn1000beta0is0beta2is0.12sig52 <- lin.reg.fxn(iter = 50, n = 1000, B = 20000, true.beta0 = 0, true.beta1 = -0, true.beta2 = 0.12, true.beta3 = -0.12, true.beta4 = 0, true.sig2 = 50, mnull.sig2 = 52)
# save(linregstrongn1000beta0is0beta2is0.12sig52, file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/linregstrongn1000beta0is0beta2is0.12sig52')
# plotting.fxn(dataset = linregstrongn1000beta0is0beta2is0.12sig52, maintitle = "n=1000 b0 = 0 b2 = 0.12 sig2 = 52 Linear Regression")
# plotting.fxn.0to0.1(dataset = linregstrongn1000beta0is0beta2is0.12sig52, maintitle = "n=1000 b0 = 0 b2 = 0.12 sig2 = 52 Linear Regression")
# 
# ################################LOGISTIC REGRESSION STRONGLY FALSE NULL SETTING################################
# 
# logregstrongn100beta0.11 <- log.reg.fxn(iter = 50, n = 100, B = 20000, true.beta0 = 0, true.beta1 = -0, true.beta2 = 0.11, true.beta3 = -.11, true.beta4 = 0)
# save(logregstrongn100beta0.11 , file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logregstrongn100beta0.11')
# plotting.fxn(dataset = logregstrongn100beta0.11, maintitle = "n=100 b0=0 b2=0.11 Logistic Regression")
# plotting.fxn.0to0.1(dataset = logregstrongn100beta0.11, maintitle = "n=100 b0=0 b2=0.11 Logistic Regression")
# 
# logregstrongn150beta0.09 <- log.reg.fxn(iter = 50, n = 150, B = 20000, true.beta0 = 0, true.beta1 = -0, true.beta2 = 0.09, true.beta3 = -.09, true.beta4 = 0)
# save(logregstrongn150beta0.09 , file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logregstrongn150beta0.09')
# plotting.fxn(dataset = logregstrongn150beta0.09, maintitle = "n=150 b0=0 b2=0.09 Logistic Regression")
# plotting.fxn.0to0.1(dataset = logregstrongn150beta0.09, maintitle = "n=150 b0=0 b2=0.09 Logistic Regression")
# 
# logregstrongn200beta0.07 <- log.reg.fxn(iter = 50, n = 200, B = 20000, true.beta0 = 0, true.beta1 = -0, true.beta2 = 0.07, true.beta3 = -.07, true.beta4 = 0)
# save(logregstrongn200beta0.07 , file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logregstrongn200beta0.07')
# plotting.fxn(dataset = logregstrongn200beta0.07, maintitle = "n=200 b0=0 b2=0.07 Logistic Regression")
# plotting.fxn.0to0.1(dataset = logregstrongn200beta0.07, maintitle = "n=200 b0=0 b2=0.07 Logistic Regression")
# 
# logregstrongn250beta0.06 <- log.reg.fxn(iter = 50, n = 250, B = 20000, true.beta0 = 0, true.beta1 = -0, true.beta2 = 0.06, true.beta3 = -.06, true.beta4 = 0)
# save(logregstrongn250beta0.06 , file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logregstrongn250beta0.06')
# plotting.fxn(dataset = logregstrongn250beta0.06, maintitle = "n=250 b0=0 b2=0.06 Logistic Regression")
# plotting.fxn.0to0.1(dataset = logregstrongn250beta0.06, maintitle = "n=250 b0=0 b2=0.06 Logistic Regression")
# 
# 
# logregstrongn500beta0.04 <- log.reg.fxn(iter = 50, n = 500, B = 20000, true.beta0 = 0, true.beta1 = -0, true.beta2 = 0.04, true.beta3 = -.04, true.beta4 = 0)
# save(logregstrongn500beta0.04 , file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logregstrongn500beta0.04')
# plotting.fxn(dataset = logregstrongn500beta0.04, maintitle = "n=500 b0=0 b2=0.04 Logistic Regression")
# plotting.fxn.0to0.1(dataset = logregstrongn500beta0.04, maintitle = "n=500 b0=0 b2=0.04 Logistic Regression")
# 
# 
# 
# logregstrongn1000beta0.0275 <- log.reg.fxn(iter = 50, n = 1000, B = 20000, true.beta0 = 0, true.beta1 = -0, true.beta2 = 0.0275, true.beta3 = -.0275, true.beta4 = 0)
# save(logregstrongn1000beta0.0275 , file = 'C:/Users/bnr/Desktop/Dissertation Stuff/R programs/saved pvalues and DCPs/logregstrongn1000beta0.0275')
# plotting.fxn(dataset = logregstrongn1000beta0.0275, maintitle = "n=1000 b0=0 b2=0.0275 Logistic Regression")
# plotting.fxn.0to0.1(dataset = logregstrongn1000beta0.0275, maintitle = "n=1000 b0=0 b2=0.0275 Logistic Regression")
