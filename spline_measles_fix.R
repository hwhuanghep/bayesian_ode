##this code give direct ode parameter estimation method for measles model
##based on spline bayesian mcmc method including
##derivative term
rm(list=ls())
graphics.off(); #close all graphics windows
require(deSolve)  #loads ODE solver package
require(mvtnorm)
require(nloptr) #for fitting
require(GenSA)
require(msm)

rk<-function(x,z){ 
  x=x/scale
  z=z/scale 
  ((z-0.5)^2-1/12)*((x-0.5)^2-1/12)/4-((abs(x-z)-0.5)^4-(abs(x-z)-0.5)^2/2+7/240)/24
}

rk1<-function(x,z){ 
  x=x/scale
  z=z/scale 
  out=((z-0.5)^2-1/12)*(x-0.5)/2-(4*(abs(x-z)-0.5)^3-(abs(x-z)-0.5))*sign(x-z)/24
  out/scale
}

# set up the penalized regression spline penalty matrix, given knot sequence xk
spl.S<-function(xk){ 
  q<-length(xk)+2;
  S<-matrix(0,q,q) # initialize matrix to 0
  S[3:q,3:q]<-outer(xk,xk,FUN=rk) # fill in non-zero part
  S
}

# set up the design matrix, given knot sequence xk
spl.X<-function(x,xk){
  n=length(x) 
  q<-length(xk)+2;
  S=matrix(1,n,q)
  S[,2]=x/scale
  S[,3:q]<-outer(x,xk,FUN=rk) 
  S
}

spl.X1<-function(x,xk){
  n=length(x) 
  q<-length(xk)+2;
  S=matrix(0,n,q)
  S[,2]=1/scale
  S[,3:q]<-outer(x,xk,FUN=rk1) 
  S
}

odeequations = function(t,y,x){
  beta = x[1];
  delta = x[2];
  ds = -beta*y[1]*y[2];
  dw = beta*y[1]*y[2] - delta*y[2];
  dydt = c(ds,dw)
  return(list(dydt))
}

odeformulas = function(x,theta,bmatrix){
  beta = theta[1];
  delta = theta[2];
  y = bmatrix%*%x
  ds = -beta*y[,1]*y[,2];
  dw = beta*y[,1]*y[,2] - delta*y[,2];
  dydt = cbind(ds,dw)
  return(dydt)
}

odematrix = function(y){
  m = matrix(0,2,2)
  m[1,] = c(-y[1]*y[2],0)
  m[2,] = c(y[1]*y[2],-y[2]) 
  return(m)
}

odevar1 <- function(x,theta,gamma,bmatrix,bmatrix1){
  beta = theta[1];
  delta = theta[2];
  s=bmatrix%*%x
  ds=bmatrix1%*%x
  hmatrix <- t(bmatrix1+bmatrix*beta*s[,2])%*%(bmatrix1+bmatrix*beta*s[,2])/gamma[1]
  hmatrix=hmatrix+beta^2*t(bmatrix*s[,2])%*%(bmatrix*s[,2])/gamma[2]
  mu <- beta*t(bmatrix*s[,2])%*%(ds[,2]+delta*s[,2])/gamma[2]
  return(list(matrix=hmatrix,mu=mu))
}

odevar2 = function(x,theta,gamma,bmatrix,bmatrix1){
  beta = theta[1];
  delta = theta[2];
  s=bmatrix%*%x
  ds=bmatrix1%*%x
  hmatrix=t(bmatrix1+bmatrix*(delta-beta*s[,1]))%*%(bmatrix1+bmatrix*(delta-beta*s[,1]))/gamma[2]
  hmatrix <- hmatrix+beta^2*t(bmatrix*s[,1])%*%(bmatrix*s[,1])/gamma[1]
  mu <- -beta*t(bmatrix*s[,1])%*%ds[,1]/gamma[1]
  return(list(matrix=hmatrix,mu=mu))
}

gpode = function(tobs,yobs,knots,missing,conditioning){
  n=nrow(yobs)
  k=ncol(yobs)
  q = length(knots)+2;
  #initial values for x 
  x=matrix(1,q,k)
  bmatrix=spl.X(tobs,knots)
  bmatrix1=spl.X1(tobs,knots)
  #make orthogonal matrix from bmatrix
  mothg=diag(q)
  a0=sqrt(sum(bmatrix[1,]^2))
  mothg[,1]=bmatrix[1,]/a0
  tmatrix=matrix(rnorm(q*q),q,q)
  for(i in 2:q){
    mothg[,i]=tmatrix[,i]
    for(j in 1:(i-1)){
      mothg[,i]=mothg[,i]-sum(mothg[,i]*mothg[,j])*mothg[,j] 
    }
    mothg[,i]=mothg[,i]/sqrt(sum(mothg[,i]^2))
  }
  ##initial value for all parameters
  sigma = rep(1,k)   #control goodness of fit
  gamma = rep(1,k)   #control ODE
  tau = 1
  theta = c(1e-4,7/7)   #four estimated parameters 
  #saved state  
  x.save = matrix(0,q,k)
  #saved parameters
  theta.save = rep(0,2)
  icount = 0
  nrep=10000
  nburn=5000
  nskip=10
  for(ii in 1:nrep){

    ##obtain ds given x
    ds=bmatrix1%*%x 

    #sample sigma
    squaresum = apply((yobs-bmatrix%*%x)^2,2,sum)
    if(!missing[1]){
      sigma[1] = 1/rgamma(1,shape=n/2+1,rate=2/squaresum[1])
    }
    if(!missing[2]){
      sigma[2] = 1/rgamma(1,shape=n/2+1,rate=2/squaresum[2])
    }
    sigma=rep(1e-9,2)

    #sample tau
    D=spl.S(knots)
    #D=diag(q)
    squareterm = sum(diag(t(x)%*%D%*%x))
    tau = 1/rgamma(1,shape=q*k/2+1,scale=2/squareterm)
    #print(tau)

    ##sample gamma
    muode = odeformulas(x,theta,bmatrix)
    squareterm = diag((t(ds-muode)%*%(ds-muode)))
    gamma[1] = 1/rgamma(1,shape=n/2+1,scale=2/squareterm[1]) + 1e-5    
    gamma[2] = 1/rgamma(1,shape=n/2+1,scale=2/squareterm[2]) + 1e-5    
    gamma=rep(1e-3,2)

    ##sample x[,1]
    ##get results from ODE
    odeout <- odevar1(x,theta,gamma,bmatrix,bmatrix1)
    hmatrix1 = odeout$matrix
    mu1 = odeout$mu
    ##get results from penalty on dx
    hmatrix2 = D/tau
    ##get results from observations
    hmatrix3 = t(bmatrix)%*%bmatrix/sigma[1]
    mu3 = t(bmatrix)%*%yobs[,1]/sigma[1]
    if(missing[1]){
    	sigmax = solve(hmatrix1+hmatrix2)
      mux = sigmax%*%mu1
    }else{
      sigmax = solve(hmatrix1+hmatrix2+hmatrix3)
      mux = sigmax%*%(mu1+mu3)
    }
    x[,1] = mux+chol(sigmax)%*%rnorm(q)
    if(missing[1]&conditioning[1]){
      mux=t(mothg)%*%mux
      sigmax=t(mothg)%*%sigmax%*%mothg
      mux1=mux[2:q]+sigmax[2:q,1]/sigmax[1,1]*(u0/a0-mux[1])
      sigma12=matrix(sigmax[2:q,1],ncol=1)
      sigmax1=sigmax[2:q,2:q]-sigma12%*%t(sigma12)/sigmax[1,1]
      xt=rep(0,q)
      xt[1]=u0/a0
      xt[2:q] = mux1+chol(sigmax1)%*%rnorm(q-1)
      x[,1]=mothg%*%xt      
    }
    
    ##sample x[,2]
    ##get results from ODE
    odeout <- odevar2(x,theta,gamma,bmatrix,bmatrix1)
    hmatrix1 = odeout$matrix
    mu1 = odeout$mu
    ##get results from penalty on dx
    hmatrix2 = D/tau
    ##get results from observations      
    hmatrix3 = t(bmatrix)%*%bmatrix/sigma[2]
    mu3 = t(bmatrix)%*%yobs[,2]/sigma[2]
    if(missing[2]){
      sigmax = solve(hmatrix1+hmatrix2)
      mux = sigmax%*%mu1
    }else{
      sigmax = solve(hmatrix1+hmatrix2+hmatrix3)
      mux = sigmax%*%(mu1+mu3)
    }
    x[,2] = mux+chol(sigmax)%*%rnorm(q)
    if(missing[2]&conditioning[2]){
      mux=t(mothg)%*%mux
      sigmax=t(mothg)%*%sigmax%*%mothg
      mux1=mux[2:q]+sigmax[2:q,1]/sigmax[1,1]*(i0/a0-mux[1])
      sigma12=matrix(sigmax[2:q,1],ncol=1)
      sigmax1=sigmax[2:q,2:q]-sigma12%*%t(sigma12)/sigmax[1,1]
      xt=rep(0,q)
      xt[1]=i0/a0
      xt[2:q] = mux1+chol(sigmax1)%*%rnorm(q-1)      
      x[,2]=mothg%*%xt      
    }    

    ##sample theta
    thetavar = 0
    thetamu = 0
    s=bmatrix%*%x
    ds=bmatrix1%*%x
    for(i in 1:n){ 
      thetavar = thetavar + (s[i,1]*s[i,2])^2*(1/gamma[1]+1/gamma[2])
      thetamu = thetamu - s[i,1]*s[i,2]*ds[i,1]/gamma[1]
      thetamu = thetamu + s[i,1]*s[i,2]*(ds[i,2]+theta[2]*s[i,2])/gamma[2]
    }
    thetamu = thetamu/thetavar
    thetavar = 1/thetavar
    theta[1] = rtnorm(1,mean=thetamu,sd=sqrt(thetavar),lower=0)

    #save the results
    if(ii>nburn&((ii-nburn)/nskip==as.integer((ii-nburn)/nskip))){
      print(theta)
      x.save = x.save + x
      theta.save = theta.save + theta
      icount = icount + 1
    }
  }
  return(list(x=x.save/icount,theta=theta.save/icount))
}

atolv=1e-7; 
rtolv=1e-7;

real=1
if(real){
setwd("H:/spring2017/ode")
data=read.table("measles.txt")
Y0 = cbind(c(0:52),data[37:89,2]);
tobs = Y0[,1]
n = nrow(Y0)
yobs = cbind(rep(1,n),Y0[,2]) 
scale=max(tobs)
k = ncol(yobs)
knots=1:9/10*scale
odesol0=cbind(tobs,yobs)
}

simu=1-real
if(simu){
pars = c(3.87e-7,1.4);
u0=3951600
i0=2011
Y0 = c(u0,i0);
tvec = 0:52;
odesol0=try(lsoda(Y0,tvec,odeequations,pars,atol=atolv,rtol=rtolv));

tobs = odesol0[,1]
yobs = odesol0[,-1] 
scale=max(tobs)
k = ncol(yobs)
n = nrow(yobs)
knots=1:9/10*scale
yobs = yobs + matrix(rnorm(k*n,0,0.15),n,k)
}

missing=c(1,0)
conditioning=c(0,0)
gpout = gpode(tobs,yobs,knots,missing,conditioning)
bmatrix=spl.X(tobs,knots)
bmatrix1=spl.X1(tobs,knots)
xhat=bmatrix%*%gpout$x
	
Y0 = xhat[1,]
pars = gpout$theta
odesol=try(lsoda(Y0,tobs,odeequations,pars,atol=atolv,rtol=rtolv));

#postscript("real.eps")
par(mfrow=c(2,1))
par(cex.lab=2)
par(cex.main=2)
par(cex.axis=2)
par(mar=c(2,4.3,2.3,1)+0.3)
if(simu){
plot(odesol0[,1],odesol0[,2],xlab="t",ylab="S")
lines(odesol[,1],odesol[,2])
lines(odesol[,1],xhat[,1],col="red")
}else{
plot(odesol[,1],odesol[,2],axes=F,xlab="",ylab="S",ylim=c(3.2e6,4e6),type="l",lwd=2)
axis(side=1,seq(0,60,by=10),tick=T,labels=F)
axis(side=2,seq(3.2e6,4e6,by=2e5),tick=T,labels=c("3.2","3.4","3.6","3.8","4"))
axis(side=3,0,labels=paste("x",expression(10^6),sep=""),tick=F)
box()
}
plot(odesol0[,1],odesol0[,3],xlab="",ylab="I",axes=F,pch=3,lwd=3)
lines(odesol[,1],odesol[,3],lwd=2)
axis(side=1,seq(0,60,by=10),tick=T,labels=F)
axis(side=2,seq(0,20000,by=5000),tick=T,labels=c("0","0.5","1","1.5","2"))
axis(side=3,0,labels=paste("x","10^4",sep=""),tick=F)
mtext(side=1,text="time",line=1,cex=2)
box()
#lines(odesol[,1],xhat[,2],col="red")
#graphics.off()



