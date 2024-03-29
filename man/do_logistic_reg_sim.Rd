% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do_logistic_reg_sim.R
\name{do_logistic_reg_sim}
\alias{do_logistic_reg_sim}
\title{Title The function which produced the BDCP logistic regression simulation study results.}
\usage{
do_logistic_reg_sim(iter, n, B, true_beta_vec)
}
\arguments{
\item{iter}{The number of iterations of the simuation to run.  For instance, if iter = 100, then the function will draw 100 separate samples, calculating the respective p-values and BDCPs for each of the 100 iterations.  iter controls the number of rows in the return data.frame.}

\item{n}{The sample size for each of the iterations.}

\item{B}{The number of bootstrap samples used to estimate each of the BDCPs.  Note that B bootstrap samples will be drawn for each of the iter iterations, so a total of iter * B total bootstrap samples will be drawn.}

\item{true_beta_vec}{The vector of length 5 that controls the underlying data generating mechanism.  The true data-generating mechanism is of the form
y_i | x_i ~ bernoulli(\theta_i), where log(\theta_i / (1 - \theta_i)) = \beta_0 + \beta_1x_1 + \beta_2x_2 + \beta_3x_3 + \beta_4x_4.  Importantly, \beta_4 is not included as a predictor in any of the specified models.  Thus, set \beta_4 to a nonzero value in order to test an underspecified null setting.  Specifically, the full null hypothesis tests H_0: \beta_0 = \beta_1 = \beta_2 = \beta_3 = 0 versus the general alternative.  THe partial null hypothesis tests H_0: \beta_2 = \beta_3 = 0, letting \beta_0 and \beta_1 be maximized.}
}
\value{
An iter x 14 data.frame.  Each row of the data.frame represents the results from a single iteration.
}
\description{
Title The function which produced the BDCP logistic regression simulation study results.
}
\examples{
#' Model containing parameters /beta_0, /beta_1, /beta_2 and /beta_3 is over-specified because the true value of /beta_3 is zero.
#' over_specified <- do_logistic_reg_sim(
  iter = 10, 
  n = 500, 
  B = 50,
  true_beta_vec = c(0.2, -0.25, 0.30, 0, 0))
Model containing parameters /beta_0, /beta_1, /beta_2 and /beta_3 is properly-specified.
properly_specified <-
  do_logistic_reg_sim(iter = 10,
                      n = 500,
                      B = 50,
                      true_beta_vec = c(0.2, -0.25, 0.30, -0.30, 0))
Model containing parameters /beta_0, /beta_1, /beta_2 and /beta_3 is under-specified because true value of /beta_4 is nonzero.  
under_specified <- 
  do_logistic_reg_sim(iter = 10,
                      n = 500, 
                      B = 50,
                      true_beta_vec = c(0.2, -0.25, 0.30, -0.30, -0.20))
}
