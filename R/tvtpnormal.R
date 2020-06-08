
# This is an  function named 'Dynamic Copula Model'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chiang Mai University
#
#
## Main Function
###############################==============

normal_pdf=function(u,v,k1){
cop <- BiCop(family = 1, par = k1)
pdf =BiCopPDF(u, v, cop)
return(pdf)
}


normtvtp=function(theta,data,rhobar,Call="optim"){

x = qnorm(data[,1])
y = qnorm(data[,2])
T = nrow(data)
u = data[,1]
v = data[,2]
w=theta[1]
a=theta[2]
b=theta[3]
d=theta[4]
kappa =rep(-0.99,T)
kappa[1] = rhobar;			# this is the MLE of kappa in the time-invariant version of this model
psi1=rep(0,T)

for (jj in  2:T){
   if (jj<=10){
      psi1[jj] = w+ a*kappa[jj-1] + b*(mean(abs(x[1:jj-1]*y[1:jj-1])))
   }else{
      psi1[jj] = w+ (a*kappa[jj-1]) + b*(mean(abs(x[(jj-10):(jj-1)]*y[(jj-10):(jj-1)])))
   }
   kappa[jj] = 1.998/(1+exp(-psi1[jj]))-0.999;		# a modified logistic transformation
}


rhohat = kappa;  # time-path of conditional copula parameter


# like normal cop
#CL = ((1-kappa^2)^(-0.5))*(exp((2*(1-kappa^2))^(-1)*((x^2)+(y^2)-2*kappa*x*y)))*(exp(0.5*((x^2)+(y^2))))




## log-like normal
CL = -1*(2*(1-kappa^2))^(-1)*((x^2)+(y^2)-2*kappa*x*y);
CL = CL + 0.5*((x^2)+(y^2));
CL = CL - 0.5*log(1-(kappa^2));

#CL=rep(0,T)
#for ( i in 1:T){
#CL[i]=normal_pdf(u[i],v[i],rhohat[i])
#}

CL = sum((CL));

n=T
if (is.infinite(CL))  # control for optimization
    CL<--n*100
if (is.nan(CL))  # control for optimization
    CL<--n*100
cat("Sum of log Likelihood for normdynamic ->",sprintf("%4.4f",c(CL,kappa[n])),"\n")

  if (Call=="optim")
    return(CL)

  if (Call=="filtering")
    {
    specOut<-list(LL=CL ,beta=c(theta),TVTP=rhohat)

    return(specOut)
    }
}


##---------------
dynamicnormal=function(data, plot){
lower =rep(-5,3)
upper =rep(5,3)
theta=c(w=0.01,beta1=0.002,beta2=-0.002)

model <- optim(theta,normtvtp,data=data,rhobar=0.5, Call="optim",
          control = list(maxit=100000,fnscale=-1),method="L-BFGS-B",
           lower =lower,upper =upper, hessian=TRUE)

# table of results
n=nrow(data)
coef<- model$par
model$se <- sqrt(-diag(solve(model$hessian)))
S.E.= model$se
(paramsWithTs = cbind (model$par , model$par/model$se ) )
stat=coef/S.E.
pvalue <- 2*(1 - pnorm(abs(stat)))
confident=matrix(0,length(coef),2)
for (i in 1:length(coef)){
confident[i,]=coef[i]+c(-1,1)*S.E.[i]*qt(0.975,60)}

result <- cbind(coef,S.E.,stat,pvalue)
BIC= -2*model$value+ (log(n)*length(coef))
AIC = -2*model$value + 2*length(coef)

rhot=normtvtp(coef,data=data,rhobar=0.5,Call="filtering")$TVTP
if ( plot==TRUE){
plot(ts( rhot[5:n]))}

output=list(
  result=result,
  AIC=AIC,
  BIC=BIC,
  Loglikelihood=model$value,
  tvtpdep=rhot
)
output
}

# Example
#library(VineCopula)
#library("tseries")
#library("quantmod")
#library("PerformanceAnalytics")

#tckk <- c("^N225", "CL=F") # ticker names defined
#numtk <- length(tckk);
#ustart <- "2010-12-30";
#uend <- "2020-2-29" # start and end date
#all_dat <- list(); # empty list to fill in the data
#for(i in 1:numtk)
#{
#  all_dat[[i]] <- xxx <- get.hist.quote(instrument = tckk[i], start=ustart, end=uend, quote = c("Close"), provider = "yahoo", compression = "m")
#}
#OIL=all_dat[[2]]+0.000001
#N225=all_dat[[1]]
#rOIL=diff(log(OIL))
#rN225=diff(log(N225))
# normal margins
#u=pnorm(rOIL/sd(rOIL))
#v=pnorm(rN225/sd(rN225))

# Correlation rho
#cor(u,v)
# maximum likelihood estimates for comparison
#BiCopEst(u, v, family = 1, method = "mle")
# Dynamic Gaussian Copula
#data=cbind(u,v)
#model=dynamicnormal(data, plot=TRUE)
