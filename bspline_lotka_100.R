##this code give direct ode parameter estimation method for lotka-volterra model
##based on bspline bayesian mcmc method for 100 samples including
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
  alpha = x[1];
  beta = x[2];
  gamma=x[3]
  delta = x[4];
  ds = y[1]*(alpha - beta*y[2]);
  dw = -y[2]*(gamma - delta*y[1]);
  dydt = c(ds,dw)
  return(list(dydt))
}

odeformulas = function(x,theta,bmatrix){
  alpha = theta[1];
  beta = theta[2];
  gamma = theta[3];
  delta = theta[4];
  y=bmatrix%*%x
  ds = y[,1]*(alpha - beta*y[,2]);
  dw = -y[,2]*(gamma - delta*y[,1]);
  dydt = cbind(ds,dw)
  return(dydt)
}

odematrix = function(y){
  m = matrix(0,2,4)
  m[1,] = c(y[1],-y[1]*y[2],0,0)
  m[2,] = c(0,0,-y[2],y[1]*y[2]) 
  return(m)
}

odevar1 <- function(x,theta,gamma,bmatrix,bmatrix1){
  alpha = theta[1];
  beta = theta[2];
  gam = theta[3];
  delta = theta[4];
  s=bmatrix%*%x
  ds=bmatrix1%*%x
  hmatrix <- t(bmatrix1-bmatrix*(alpha-beta*s[,2]))%*%(bmatrix1-bmatrix*(alpha-beta*s[,2]))/gamma[1]
  hmatrix=hmatrix+delta^2*t(bmatrix*s[,2])%*%(bmatrix*s[,2])/gamma[2]
  mu <- delta*t(bmatrix*s[,2])%*%(ds[,2]+gam*s[,2])/gamma[2]
  return(list(matrix=hmatrix,mu=mu))
}

odevar2 = function(x,theta,gamma,bmatrix,bmatrix1){
  alpha = theta[1];
  beta = theta[2];
  gam = theta[3];
  delta = theta[4];
  s=bmatrix%*%x
  ds=bmatrix1%*%x
  hmatrix <- t(bmatrix1+(gam-delta*s[,1])*bmatrix)%*%(bmatrix1+(gam-delta*s[,1])*bmatrix)
  hmatrix=hmatrix/gamma[2]+beta^2*t(bmatrix*s[,1])%*%(bmatrix*s[,1])/gamma[1]
  mu <- beta*t(bmatrix*s[,1])%*%(alpha*s[,1]-ds[,1])/gamma[1]
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
  theta = rep(1,4)   #four estimated parameters 
  #saved state  
  x.save = matrix(0,q,k)
  #saved parameters
  theta.save = rep(0,4)
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
    if(missing[2]){ 
      gamma=rep(1e-7,2)
    }

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
    thetavar = matrix(0,4,4)
    thetamu = matrix(0,4,1)
    s=bmatrix%*%x
    ds=bmatrix1%*%x
    for(i in 1:n){ 
      mm = odematrix(s[i,])
      thetavar = thetavar + t(mm)%*%diag(1/gamma)%*%mm
      thetamu = thetamu + t(mm)%*%diag(1/gamma)%*%ds[i,]
    }
    thetamu = solve(thetavar+diag(1e-5,4))%*%thetamu
    thetavar = solve(thetavar+diag(1e-5,4))
    do.sep <- 1
    if(do.sep){
      for(i in 1:4){
        mutemp <- thetamu[i,1] + thetavar[i,-i]%*%solve(thetavar[-i,-i])%*%(theta[-i]-thetamu[-i,])
        vartemp <- thetavar[i,i] - thetavar[i,-i]%*%solve(thetavar[-i,-i])%*%thetavar[-i,i]
        theta[i] = rtnorm(1,mean=mutemp,sd=sqrt(vartemp),lower=0)
      }
    }else{
      theta = thetamu + chol(thetavar)%*%matrix(rnorm(4),4,1)
    }

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

#setwd("F:/spring2017/ode")

pars = c(3,2,3,1);
u0=5
i0=3
Y0 = c(u0,i0);
tvec = seq(0,4,by=0.1);
odesol1=try(lsoda(Y0,tvec,odeequations,pars,atol=atolv,rtol=rtolv));

tobs = odesol1[,1]
yobs = odesol1[,-1] 
scale=max(tobs)
k = ncol(yobs)
n = nrow(yobs)
knots=1:9/10*scale

nrep=2
thetas=NULL
for(i in 1:nrep){
  yobs = yobs + matrix(rnorm(k*n,0,0.2),n,k)
  missing=c(0,0)
  conditioning=c(0,0)
  gpout = gpode(tobs,yobs,knots,missing,conditioning)
  thetas = rbind(thetas,gpout$theta)
}
print(apply(thetas,2,mean))
print(apply(thetas,2,sd))
write.table(thetas,file="lotka_theta.txt",sep="\t",quote=F,row.names=F,col.names=F)

bmatrix=spl.X(tobs,knots)
xhat=bmatrix%*%gpout$x
Y0 = xhat[1,]
pars = apply(thetas,2,mean)
odesol=try(lsoda(Y0,tobs,odeequations,pars,atol=atolv,rtol=rtolv));


par(mfrow=c(2,1))
par(cex.lab=1.5)
par(cex.main=1.5)
par(cex.axis=1.5)
par(mar=c(4,4.3,0.3,1)+0.3)
plot(tobs,yobs[,1],xlab="t",ylab="S",ylim=range(yobs[,1]),pch=3,lwd=3,main="")
lines(odesol1[,1],odesol1[,2],lty=2)
lines(odesol[,1],odesol[,2],lty=2,col="red")
plot(tobs,yobs[,2],xlab="t",ylab="W",ylim=range(yobs[,2]),pch=3,lwd=3,main="")
lines(odesol1[,1],odesol1[,3],lty=2)
lines(odesol[,1],odesol[,3],lty=2,col="red")

