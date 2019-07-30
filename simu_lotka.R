##this code give direct ode parameter estimation method
##based on gaussian process bayesian mcmc method including
##derivative term
rm(list=ls())
graphics.off(); #close all graphics windows
require(deSolve)  #loads ODE solver package
require(mvtnorm)
require(nloptr) #for fitting
require(GenSA)
require(msm)

odeequations = function(t,y,x){
  alpha = x[1];
  beta = x[2];
  gamma = x[3];
  delta = x[4];
  ds = y[1]*(alpha - beta*y[2]);
  dw = -y[2]*(gamma - delta*y[1]);
  dydt = c(ds,dw)
  return(list(dydt))
}

odeformulas = function(y,x){
  alpha = x[1];
  beta = x[2];
  gamma = x[3];
  delta = x[4];
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

odevar1 <- function(dx,y,x,gamma){
  hmatrix <- diag((x[1]-x[2]*y[,2])^2+(x[4]*y[,2])^2)/gamma
  mu <- ((x[1]-x[2]*y[,2])*dx[,1]+(x[3]*y[,2]+dx[,2])*x[4]*y[,2])/((x[1]-x[2]*y[,2])^2+(x[4]*y[,2])^2)
  return(list(matrix=hmatrix,mu=mu))
}

odevar2 <- function(dx,y,x,gamma){
  hmatrix <- diag((x[2]*y[,1])^2+(x[4]*y[,1]-x[3])^2)/gamma
  mu <- ((x[4]*y[,1]-x[3])*dx[,2]-x[2]*y[,1]*(dx[,1]-x[1]*y[,1]))/((x[2]*y[,1])^2+(x[4]*y[,1]-x[3])^2)
  return(list(matrix=hmatrix,mu=mu))
}

gpode = function(tobs,yobs,k,k0){
  D = 2
  n = length(tobs);
  ##initial value
  sigma = rep(1,k)
  gamma = 1
  tau = 1
  theta = rep(1,4)
  if(k0!=k){
    missing=1
    x = cbind(rep(0,n),rep(0,n),yobs)  
    dx <- matrix(0,n,k)
    for(i in 1:(n-1)){
      dx[i,1] <- (x[i+1,1]-x[i,1])/(tobs[i+1]-tobs[i])
    }
    dx[n,1] <- dx[n-1,1]
    x[,2] <- (theta[1]*x[,1]-dx[,1])/theta[2]/x[,1]
    yobs=x
  }else{
    x=yobs
    missing=0
  }
  ##define GP covariance
  gpsigma.0 = diag(n)
  gpsigmad1.0 = diag(n)
  gpsigmad2.0 = diag(n)
  gpsigmadd.0 = diag(n)
  for(i in 1:n){
    for(j in 1:n){
      gpsigma.0[i,j] = exp(-(tobs[i] - tobs[j])^2/D^2)
      gpsigmad1.0[i,j] = - 2*(tobs[i] - tobs[j])*gpsigma.0[i,j]/D^2
      gpsigmad2.0[i,j] = - 2*(tobs[j] - tobs[i])*gpsigma.0[i,j]/D^2
      gpsigmadd.0[i,j] = (2 - 4*(tobs[i] - tobs[j])^2/D^2)*gpsigma.0[i,j]/D^2
    }
  }
  gpsigma.0 <- gpsigma.0 + diag(1e-3,n)
  gpsigmad1.0  <- gpsigmad1.0 
  gpsigmad2.0  <- gpsigmad2.0 
  gpsigmadd.0  <- gpsigmadd.0 + diag(1e-3,n)
  gpsigma = tau*gpsigma.0 
  gpsigmad1 = tau*gpsigmad1.0 
  gpsigmad2 = tau*gpsigmad2.0 
  gpsigmadd = tau*gpsigmadd.0 
  i.count = 0
  x.save = matrix(0,n,k)
  theta.save = rep(0,4)
  for(ii in 1:11000){
    ##sample dx given x
    mudxx = gpsigmad1%*%solve(gpsigma)%*%x
    sigmadxx = gpsigmadd - gpsigmad1%*%solve(gpsigma)%*%gpsigmad2
    lambda = gamma*diag(n)
    ##ode results
    muode = odeformulas(x,theta)
    dxvar = sigmadxx%*%solve(sigmadxx+lambda)%*%lambda
    dxmu = dxvar%*%(solve(sigmadxx)%*%mudxx+solve(lambda)%*%muode)
    dx = dxmu+chol(dxvar)%*%matrix(rnorm(n*k),n,k)
    ##sample tau
    mudxx.0 = gpsigmad1.0%*%solve(gpsigma.0)%*%x
    sigmadxx.0 = gpsigmadd.0 - gpsigmad1.0%*%solve(gpsigma.0)%*%gpsigmad2.0 
    squareterm = sum(diag(t(x)%*%solve(gpsigma.0)%*%x)) 
    + sum(diag(t(dx-mudxx.0)%*%solve(sigmadxx.0)%*%(dx-mudxx.0)))
    tau = 1/rgamma(1,shape=k*n+1,scale=2/squareterm)
    prob=-0.5*(squareterm/tau+log(det(tau*gpsigma.0))+log(det(tau*sigmadxx.0)))

    gpsigma = tau*gpsigma.0
    gpsigmad1 = tau*gpsigmad1.0
    gpsigmad2 = tau*gpsigmad2.0
    gpsigmadd = tau*gpsigmadd.0
    ##sample gamma
    squareterm = sum(diag((t(dx-muode)%*%(dx-muode))))
    gamma = 1/rgamma(1,shape=n*k/2+1,scale=2/squareterm)
    ##sample theta
    thetavar = matrix(0,4,4)
    thetamu = matrix(0,4,1)
    for(i in 1:n){ 
      mm = odematrix(x[i,])
      thetavar = thetavar + t(mm)%*%mm
      thetamu = thetamu + t(mm)%*%dx[i,]
    }
    thetamu = solve(thetavar+gamma*diag(1e-3,4))%*%thetamu
    thetavar = solve(thetavar+gamma*diag(1e-3,4))*gamma
    do.sep <- 0
    if(do.sep){
      for(i in 1:4){
        mutemp <- thetamu[i,1] + thetavar[i,-i]%*%solve(thetavar[-i,-i])%*%(theta[-i]-thetamu[-i,])
        vartemp <- thetavar[i,i] - thetavar[i,-i]%*%solve(thetavar[-i,-i])%*%thetavar[-i,i]
        theta[i] = rtnorm(1,mean=mutemp,sd=sqrt(vartemp),lower=0)
      }
    }else{
      theta = thetamu + chol(thetavar)%*%matrix(rnorm(4),4,1)
    }
    ##sample x[,1]
    ##get results from ODE
    odeout <- odevar1(dx,x,theta,gamma)
    hmatrix1 = odeout$matrix
    mu1 = odeout$mu
    ##get results from Gaussian conditional on dx
    mu2 = as.vector(gpsigmad2%*%solve(gpsigmadd)%*%dx[,1])
    hmatrix2 = gpsigma - gpsigmad2%*%solve(gpsigmadd)%*%gpsigmad1
    ##get results from observations
    hmatrix3 = diag(sigma[1],n)
    mu3 = as.vector(yobs[,1])
    sigmax = solve(hmatrix1+solve(hmatrix2)+solve(hmatrix3))
    mux = sigmax%*%(hmatrix1%*%mu1+solve(hmatrix2)%*%mu2+solve(hmatrix3)%*%mu3)
    x[,1] = mux+chol(sigmax)%*%rnorm(n)
    
    ##sample x[,2]
    ##get results from ODE
    odeout <- odevar2(dx,x,theta,gamma)
    hmatrix1 = odeout$matrix
    mu1 = odeout$mu
    ##get results from Gaussian conditional on dx
    mu2 = as.vector(gpsigmad2%*%solve(gpsigmadd)%*%dx[,2])
    hmatrix2 = gpsigma - gpsigmad2%*%solve(gpsigmadd)%*%gpsigmad1

    if(missing){
      sigmax = solve(hmatrix1+solve(hmatrix2))
      mux = sigmax%*%(hmatrix1%*%mu1+solve(hmatrix2)%*%mu2)
    }else{
    ##get results from observations      
      hmatrix3 = diag(sigma[2],n)
      mu3 = as.vector(yobs[,2])
      sigmax = solve(hmatrix1+solve(hmatrix2)+solve(hmatrix3))
      mux = sigmax%*%(hmatrix1%*%mu1+solve(hmatrix2)%*%mu2+solve(hmatrix3)%*%mu3)
    }
    x[,2] = mux+chol(sigmax)%*%rnorm(n)
    
    ##sample sigma
    for(i in 1:k){
      squaresum = sum((yobs[,i]-x[,i])^2)
      sigma[i] = 1/rgamma(1,shape=n/2+1,rate=2/squaresum)
    }

    ##sample D
    sampleD=1
    if(sampleD){
      Dnew=runif(1,D-delta,D+delta)

      if(Dnew<=0){
        Dnew=delta
      }
      if(Dnew>5){
        Dnew=5
      }
   
      gpsigma.new = diag(n)
      gpsigmad1.new = diag(n)
      gpsigmad2.new = diag(n)
      gpsigmadd.new = diag(n)
      for(i in 1:n){
        for(j in 1:n){
          gpsigma.new[i,j] = exp(-(tobs[i] - tobs[j])^2/Dnew^2)
          gpsigmad1.new[i,j] = - 2*(tobs[i] - tobs[j])*gpsigma.new[i,j]/Dnew^2
          gpsigmad2.new[i,j] = - 2*(tobs[j] - tobs[i])*gpsigma.new[i,j]/Dnew^2
          gpsigmadd.new[i,j] = (2 - 4*(tobs[i] - tobs[j])^2/Dnew^2)*gpsigma.new[i,j]/Dnew^2
        }
      }
      gpsigma.new <- gpsigma.new + diag(1e-3,n)
      gpsigmadd.new <- gpsigmadd.new + diag(1e-3,n)
      mudxx.new = gpsigmad1.new%*%solve(gpsigma.new)%*%x
      sigmadxx.new = gpsigmadd.new - gpsigmad1.new%*%solve(gpsigma.new)%*%gpsigmad2.new 
      squareterm = sum(diag(t(x)%*%solve(gpsigma.new)%*%x)) 
      + sum(diag(t(dx-mudxx.new)%*%solve(sigmadxx.new)%*%(dx-mudxx.new)))
      probnew=-0.5*(squareterm/tau+log(det(tau*gpsigma.new))+log(det(tau*sigmadxx.new)))
      if(exp(probnew-prob)>runif(1)){
        D=Dnew
        gpsigma.0 = gpsigma.new
        gpsigmad1.0 = gpsigmad1.new
        gpsigmad2.0 = gpsigmad2.new
        gpsigmadd.0 = gpsigmadd.new
        gpsigma = tau*gpsigma.0 
        gpsigmad1 = tau*gpsigmad1.0 
        gpsigmad2 = tau*gpsigmad2.0 
        gpsigmadd = tau*gpsigmadd.0 
      }
    }


    if(ii>1000&ii/10==as.integer(ii/10)){
 	print(theta)
      print(sigma)
      i.count = i.count + 1
      x.save = x.save + x
      theta.save = theta.save + theta
    }
  }
  return(list(x=x.save/i.count,theta=theta.save/i.count))
}

atolv=1e-7; 
rtolv=1e-7;
pars = c(2,1,4,1);
Y0 = c(5,3);
tvec = seq(0,4,by=0.2);
odesol=try(lsoda(Y0,tvec,odeequations,pars,atol=atolv,rtol=rtolv));
par(mfrow=c(2,1))
plot(odesol[,1],odesol[,2],type="l")
plot(odesol[,1],odesol[,3],type="l")

tobs = odesol[,1]
yobs = odesol[,-1] 
k = ncol(yobs)
n = nrow(yobs)
yobs = yobs + matrix(rnorm(k*n,0,0.15),n,k)

gpout = gpode(tobs,yobs[1,],2,1)
Y0 = gpout$x[1,]
pars = gpout$theta
odesol=try(lsoda(Y0,tvec,odeequations,pars,atol=atolv,rtol=rtolv));

par(mfrow=c(2,1))
par(cex.lab=1.5)
par(cex.main=1.5)
par(cex.axis=1.5)
par(mar=c(4,4.3,0.3,1)+0.3)

plot(tobs,yobs[,1],xlab="t",ylab="S")
lines(odesol[,1],odesol[,2])
plot(tobs,yobs[,2],xlab="t",ylab="W")
lines(odesol[,1],odesol[,3])


