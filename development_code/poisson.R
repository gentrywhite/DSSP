library(sp)
library(gstat)
library(DSSP)
library(tidyverse)
library(mvtnorm)

##  Generate Data

data(meuse.all)
coordinates(meuse.all) <- ~ x + y

M<-make.M(scale(coordinates(meuse.all)))
COV<-solve(diag(1,nrow(meuse.all))+M$M)

set.seed(14061972)

la<-rmvnorm(1,mean=rep(0,nrow(meuse.all)),sigma=COV)

y<-rpois(length(la),la%>%exp()%>%t())

df<-tibble(Y=y,x = meuse.all$x, y = meuse.all$y,la = t(la))
coordinates(df)<-~x+y

##  Generate Model Fit

TMP<-DSSP(formula = Y~1,data = df,pars = c(0.001,0.001), N=1000)

TMP2<-DSSP(formula = la~1,data = df,pars = c(0.001,0.001), N=1000)

df <- tibble(nu = t(TMP$nu),delta = TMP$delta, eta = TMP$eta)
df2 <- tibble(nu = t(TMP2$nu),delta = TMP2$delta, eta = TMP2$eta)

gaussian_post<-function(pars)
{
  n<- length(y)
  nu <- pars[1:n]
  delta <- pars[n+1]%>%as.numeric()
  eta <- as.numeric(pars[n+2])
  dmvnorm(y,mean=nu,sigma = diag(delta,n),log = TRUE)
}
  
poisson_post<-function(pars)  
{
  n<- length(y)
  nu <- pars[1:n]
  delta <- pars[n+1]%>%as.numeric()
  eta <- as.numeric(pars[n+2])
  dpois(y,exp(nu),log = TRUE)%>%sum()
  
}

##  Now test against the PSIS

gp <- apply(df,1,gaussian_post)

pp <- apply(df,1,poisson_post)

w<-pp-gp

w<-exp(w)

w<-w/sum(w)

library(loo)

log_ratios <- -1 * array(pp-gp,dim=c(length(gp),1,1))
r_eff <- relative_eff(exp(-log_ratios))
psis_result <- psis(log_ratios, r_eff = r_eff)
str(psis_result)
plot(psis_result)

##  Meh

##  Now try to just use the eta and delta samples to draw from the "true" posterior using 
##  RU, ARS, SLICE etc.  Can't be worse than INLA. 

eta <- TMP$eta

delta <-TMP$delta

par <- 0.5*eta/delta

##  pars <- eta/2*delta

poisson_post <- function(x,y,pars,m)
{
  
  dpois(y,exp(x),log = TRUE)-x*m*x*pars-0.5*sum(log(M$M.eigen$values[1:161]))+(0.5*n)*log(pars)
}

library(HI)

m<-diag(M$M)
n<-length(m)
pars <-rep(par[1],n)
mat <- cbind(y,pars,m)

sample_pp<-function(X)
{
arms(0, function(x) poisson_post(x,X[1],X[2],X[3]), function(x) (x<10)*(x>-10), n.sample = 1)
} 

nu<-apply(mat,1,sample_pp)




