\name{CorTransform}
\alias{CorTransform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transforming the time varying copula parameters
}
\description{
This function computes the theoretical tail dependence coefficients,  Kendall’s tau value, and l Blomqvist’s beta value of a bivariate copula for given parameter values.
}
\usage{
CorTransform(dep, family, transform="Tau")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dep}{the vector of T x 1 dependencey parameter}
  \item{family}{Copula function (1=Gaussian, 2=Student-t)}
  \item{transform}{ The transfomation methods : "Tau" , "Beta", and "Tail"}
}

\details{
No
}
\value{

\item{out}{time varying correlation measures}


}
\references{
Schepsmeier, U., Stoeber, J., Brechmann, E. C., Graeler, B., Nagler, T., Erhardt, T., ... & Killiches, M. (2015). Package ‘VineCopula’. R package version, 2(5).
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
# Example
# Gaussian family =1
model=dynamicnormal(data, plot=TRUE)
out=CorTransform(dep=model$tvtpdep, family=1, transform="Tau")

# Stundet-t family =2
model1=dynamicT(data, plot=TRUE)
out2=CorTransform(model1$tvtpdep, family=2, transform="Beta")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
