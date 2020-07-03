# MLEB
Machine Learning Empirical Bayes

R code for an empirical bayes estimation for respondent level utilities ~ Multivariate Normal(mu, sigma).  Unlike simple empirical bayes, mu and sigma are estimated via an EM loop.  The initial starting point of mu is the aggregate solution across all respondents.  Sigma is a diagonal matrix of the standard errors for each variable, coupled with the Peter Lenk adjustment for number of levels of categorical variables.  Respondent level utilities are estimated, and we in turn update mu and sigma.  We use Bayesian updating for mu and sigma, jointly using a normal-inverse-wishart conjugate prior.  This gives simple closed form updating for mu and sigma.   

Beginning July 2020, MLEB was recoded to use parallel threads in R.  MLEB loops across respondents in several cases, and these now use foreach loops, requiring the doParallel package.
