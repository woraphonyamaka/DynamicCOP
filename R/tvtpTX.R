# This is an  function named 'Dynamic Copula Model'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chiang Mai University
#
#
## Main Function
###############################==============

T_pdf=function(u,v,k1,k2){
cop <- BiCop(family = 2, par = k1,par2 = k2)
pdf =BiCopPDF(u, v, cop)
return(pdf)
}


Ttvtpx=function(theta,data,z,Call="optim"){

x = qstd(data[,1])
y = qstd(data[,2])
n=length(x)
T = nrow(data)
u = data[,1]
v = data[,2]
w=theta[1]
a=theta[2]
b=theta[3]
d=theta[4]
w1=theta[5]
a1=theta[6]
b1=theta[7]

kappa =cbind(rep(-0.99,n),rep(5,n))

kappa[,1] = 0.5;			# this is the MLE of kappa in the time-invariant version of this model
kappa[,2] = 5;               # degree of freedom
psi1=matrix(0,n,2)

for (jj in  2:n){
 if (jj<=10){
      psi1[jj,1] = w+ a*kappa[jj-1,1] + b*(mean(abs(x[1:jj-1]*y[1:jj-1])))+d*z[jj]
      psi1[jj,2] = w1+ a1*kappa[jj-1,2] + b1*(mean(abs(x[1:jj-1]*y[1:jj-1])))

kappa[jj,1] = 1.998/(1+exp(-psi1[jj]))-0.999;		          # a modified logistic transformation
kappa[jj,2] = (exp( psi1[jj,2])/(1+exp(psi1[jj,2])))*98 + 2.01   # a modified logistic transformation df

   }else{
 psi1[jj,1] = w+  a*kappa[jj-1,1]  + b* (mean(abs(x[(jj-10):(jj-1)]*y[(jj-10):(jj-1)])))+d*z[jj]
 psi1[jj,2] = w1+ a1*kappa[jj-1,2] + b1*(mean(abs(x[(jj-10):(jj-1)]*y[(jj-10):(jj-1)])))

kappa[jj,1] = 1.998/(1+exp(-psi1[jj]))-0.999;		      # a modified logistic transformation
kappa[jj,2] = (exp( psi1[jj,2])/(1+exp(psi1[jj,2])))*98 + 2.01   # a modified logistic transformation df

}

}



RHO = kappa[,1];  # time-path of conditional T copula parameter
NU= kappa[,2]        # time-path of conditional T copula parameter


#CL=rep(0,T)
#for ( tt in 1:T){
#   CL[tt] = gamma((NU[tt]+2)/2)/gamma(NU[tt]/2)*((NU[tt]*pi)^(-1))/sqrt(1-RHO[tt]^2)*((1 + (x[tt]^2+y[tt]^2-2*RHO[tt]*x[tt]*y[tt])/(NU[tt]*(1-RHO[tt]^2)))^(-(NU[tt]+2)/2))
#   CL[tt]= CL[i]/(dstd(x[tt],0,1,NU[tt])*dstd(y[tt],0,1,NU[tt]))+0.000000001;
#}


CL=rep(0,T)
for ( i in 1:T){
CL[i]=T_pdf(u[i],v[i],RHO[i],NU[i])
}

CL = sum(log(CL));

n=T
if (is.infinite(CL))  # control for optimization
    CL<--n*100
if (is.nan(CL))  # control for optimization
    CL<--n*100
cat("Sum of log Likelihood for TdynamicX ->",sprintf("%4.4f",c(CL,kappa[n])),"\n")

  if (Call=="optim")
    return(CL)

  if (Call=="filtering")
    {
    specOut<-list(LL=CL ,beta=c(theta),TVTP=RHO)

    return(specOut)
    }
}


##---------------
dynamicTx=function(data, z, plot){
lower =rep(-10,7)
upper =rep(10,7)
theta=c(w=0.01,beta1=0.002,beta2=-0.002,exo=0.5,
w2=0.01,beta3=0.002,beta4=-0.002)


model <- optim(theta,Ttvtpx,data=data, z=z, Call="optim",
          control = list(maxit=100000,fnscale=-1),method="L-BFGS-B",
           lower =lower,upper =upper, hessian=TRUE )
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

rhot=Ttvtpx(coef,data=data,z=z,Call="filtering")$TVTP
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

#tckk <- c("^N225", "CL=F","^DJI") # ticker names defined
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
#DOWJ=all_dat[[3]]
#rOIL=diff(log(OIL))     # exogenous variable
#rN225=diff(log(N225))
#rDOWJ=diff(log(DOWJ))
# normal margins
#u=pnorm(rDOWJ/sd(rDOWJ))
#v=pnorm(rN225/sd(rN225))

# Correlation rho
#cor(u,v)
# maximum likelihood estimates for comparison
#BiCopEst(u, v, family = 2, method = "mle")
# Dynamic Student-t  Copula
#data=cbind(u,v)
#model=dynamicTx(data, z=rOIL, plot=TRUE)


