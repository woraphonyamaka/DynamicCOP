#BiCopPar2TailDep Tail Dependence Coefficients of a Bivariate Copula

CorTransform=function(dep, family, transform){
  kappa=dep
  COR=c()
  n=length(dep)
  #### Beta dependence
  if (transform=="Beta"){
  for ( i in 1:n) {
    COR[i]=BiCopPar2Beta(family = family, par = kappa[i])
  }
  plot(ts(COR[5:n]), main="Time varying Beta", ylab="Beta values")
  COR
  }
  if (transform=="Tail"){
  ### Tail dependence
  for ( i in 1:n) {
    COR[i]=BiCopPar2TailDep(family, kappa[i])$upper
  }
  plot(ts(COR[5:n]), main="Time varying tail dependence", ylab="Tail values")
  COR
  }
  if (transform=="Tau"){
  ## Kendal Tau
  for ( i in 1:n) {
    COR[i]=BiCopPar2Tau(family = family, par = kappa[i])
  }
  plot(ts(COR[5:n]), main="Time varying Tau", ylab="Tau values")
  COR
  }
}

# Example
# Gaussian family =1
#model=dynamicnormal(data, plot=TRUE)
#out=CorTransform(dep=model$tvtpdep, family=1, transform="Tau")

# Stundet-t family =2
#model1=dynamicT(data, plot=TRUE)
#out2=CorTransform(model1$tvtpdep, family=2, transform="Beta")
