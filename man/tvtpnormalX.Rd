\name{tvtpnormalX}
\alias{tvtpnormalX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dynamic Gaussian Copula Model with exogenous variable
}
\description{
This function estimates the parameter(s) of a bivariate dynamic gaussian copula with exogenous variable using  maximum likelihood estimation.
}
\usage{
dynamicnormalX(data,z, plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ The matrix T x 2 margins of the two random varibles}
  \item{data}{ plot the time varying dependence parameter (plot=TRUE, default)}
  \item{z}{The vector T x exogenous variable}
  }
\details{
This model can be used as an alternative way of connecting the marginal distributions to restore the joint distribution. The ARMA(1,10) process of Patton(2006) is extened to ARMAX(1,10), thus the additional exogenous variable is added in the time varying copula equaiton.}
\value{

\item{result }{Estimated parameter, Standard error, t-stat, p-value}
\item{AIC }{Akaiki Information Criteria}
\item{BIC }{Bayesain Information Criteria}
\item{Loglikelihood }{Log Likeihood function'}
\item{tvtpdep }{the vector of time varying dependence parameter}
}
\references{
Patton, A. J. (2006). Modelling asymmetric exchange rate dependence. International economic review, 47(2), 527-556.

Maneejuk, P., & Yamaka, W. (2019). Predicting Contagion from the US Financial Crisis to International Stock Markets Using Dynamic Copula with Google Trends. Mathematics, 7(11), 1032.
}
\author{
Woraphon Yamaka
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(VineCopula)
library("tseries")
library("quantmod")
library("PerformanceAnalytics")

tckk <- c("^N225", "CL=F","^DJI") # ticker names defined
numtk <- length(tckk);
ustart <- "2010-12-30";
uend <- "2020-2-29" # start and end date
all_dat <- list(); # empty list to fill in the data
for(i in 1:numtk)
{
  all_dat[[i]] <- xxx <- get.hist.quote(instrument = tckk[i], start=ustart, end=uend, quote = c("Close"), provider = "yahoo", compression = "m")
}
OIL=all_dat[[2]]+0.000001
N225=all_dat[[1]]
DOWJ=all_dat[[3]]
rOIL=diff(log(OIL))     # exogenous variable
rN225=diff(log(N225))
rDOWJ=diff(log(DOWJ))
# normal margins
u=pnorm(rDOWJ/sd(rDOWJ))
v=pnorm(rN225/sd(rN225))

# Correlation rho
cor(u,v)
# maximum likelihood estimates for comparison
BiCopEst(u, v, family = 1, method = "mle")
# Dynamic Gaussian Copula
data=cbind(u,v)
model=dynamicnormalx(data, z=rOIL, plot=TRUE)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
