% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SuperAdap.R
\name{SuperAdap}
\alias{SuperAdap}
\title{One-sided test in a sensitivity analysis in matched observational studies
via the two-stage programming method}
\usage{
SuperAdap(Q, Z, index, alpha, Gamma, alternative)
}
\arguments{
\item{Q}{An N by 2 matrix of scores of the two component sum test statistics,
where N is the total number of units in the study. That is, the n-th
entry of the first column of \code{Q} is the score of the first component test
statistic of unit n and the n-th entry of the second column of \code{Q} is the
score of the second component test statistic of unit n.}

\item{Z}{An N-length vector of the treatment indicators of all N units in the
study: 1 if treated and 0 if not. The n-th entry of \code{Z} is the treatment
indicator of unit n. Note that in each matched set there is one and only
one treated unit.}

\item{index}{An N-length vector of the matched set indexes of all N units in
the study, taking value from 1 to I, where I is the total number of
matched sets. That is, the n-th entry of \code{index} is the matched
set index of unit n.}

\item{alpha}{The level of the one-sided test.}

\item{Gamma}{The sensitivity parameter in the Rosenbaum bounds sensitivity analysis, which
is a prespecified number that is greater than or equal to 1.}

\item{alternative}{The direction of the alternative in a one-sided test, can be either
"greater", i.e., greater than, or "less", i.e., less than.}
}
\value{
An indicator of the null hypothesis to be rejected or not: "reject" or "failed to reject".
}
\description{
\code{SuperAdap} is the main function for performing a one-sided test
   under the Rosenbaum bounds sensitivity analysis framework via the two-stage
   programming method, which is an adaptive approach with a nice asymptotic
   property called "super-adaptivity". For details of the two-stage programming
   method, see the paper "Increasing Power for Observational Studies of Aberrant
   Response: An Adaptive Approach" by Heng, Kang, Small, and Fogarty.
}
\examples{
#We randomly generate a dataset with I matched sets along with
#the scores of the two component test statistics
I=200
n<-rep(0, I)
for (i in 1:I){
  n[i]=sample(c(2:6), size = 1)
}
N=sum(n)
index_1<-rep(0, N)
Z_1<-rep(0, N)
Q_1<-matrix(0, nrow = N, ncol = 2)
S<-rep(0, I)
for (i in 1:I){
  S[i]=sum(n[1:i])
}
index_1[1:S[1]]=1
Z_1[1]=1
Z_1[1:S[1]]=sample(Z_1[1:S[1]], size = S[1])
if (I>1){
  for (i in 2:I){
    index_1[(S[i-1]+1):S[i]]=i
    Z_1[(S[i-1]+1)]=1
    Z_1[(S[i-1]+1):S[i]]=sample(Z_1[(S[i-1]+1):S[i]], size = n[i])
  }
}
for (s in 1:N){
  if (Z_1[s]==1){
   Q_1[s, 1]=rnorm(1, mean = 0.5, sd=1)
   Q_1[s, 2]=rnorm(1, mean = 0.6, sd=1)
 }
 else {
   Q_1[s, 1]=rnorm(1, mean = 0.1, sd=1)
   Q_1[s, 2]=rnorm(1, mean = 0, sd=1)
 }
}
#We then run the SuperAdap function with the generated dataset
result_1=SuperAdap(Q = Q_1, Z = Z_1, index = index_1,
                      alpha = 0.05, Gamma = 1.5, alternative = "greater")
result_2=SuperAdap(Q = Q_1, Z = Z_1, index = index_1,
                      alpha = 0.05, Gamma = 5, alternative = "greater")
}
