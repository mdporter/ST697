################################################################################
# B-splines and penalized b-splines (PB-splines)
# ST 697 Fall 2015
# MD Porter (University of Alabama)
################################################################################


library(fda)
library(splines)

#-----------------------------------------------------------------------------
#-- Generate Data
#-----------------------------------------------------------------------------
## simulate some data - from mgcv::magic
set.seed(1)
n = 200
x = sort(runif(n))          # sort the x values so plots of splines is simple
f <- function(x) 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
y = f(x) + rnorm(n, 0, sd = 2)

plot(x,y)
curve(f,add=TRUE,lwd=2)


#-----------------------------------------------------------------------------
#- Polynomial basis functions
#-----------------------------------------------------------------------------

P2 = poly(x,2)
matplot(x,P2,typ='p',pch=19)
abline(h=0,col="lightgray");abline(v=mean(x),col="lightgray")

P3 = poly(x,3)
matplot(x,P3,typ='p',pch=19)
abline(h=0,col="lightgray");abline(v=mean(x),col="lightgray")


#--------------------
#- Polynomial fits
#--------------------

#- polynomial models
m0 = lm(y~1)               # intercept only
m1 = lm(y~x)               # + linear
m2 = lm(y~poly(x,2))       # + linear + quadractic
m3 = lm(y~poly(x,3))       # + linear + quadractic + cubic

#- plots of fits
plot(x,y)
curve(f,add=TRUE,col="grey",lwd=2)
tt = seq(min(x)-.1,max(x)+.1,length=100)
lines(tt, predict(m0, data.frame(x = tt)), col="black")
lines(tt, predict(m1, data.frame(x = tt)), col="red")
lines(tt, predict(m2, data.frame(x = tt)), col="green")
lines(tt, predict(m3, data.frame(x = tt)), col="blue")
title("Polynomial Fits")

m10 = lm(y~poly(x,10))
lines(tt, predict(m10, data.frame(x = tt)), lw=3,col="orange")

#-- Show components
beta = coef(m3)                # coefficients

plot(x,y,type='n')             # set-up plot (but don't plot anything type="n)
P3 = cbind(1,poly(x,3))        # create polynomial basis function
matlines(x,P3[,1] * beta[1],col="orange")  # constant
matlines(x,P3[,2] * beta[2],col="red")     # linear
matlines(x,P3[,3] * beta[3],col="green")   # quadratic
matlines(x,P3[,4] * beta[4],col="blue")    # cubic

matlines(x, P3 %*% beta, col="black",lwd=2) # fitted curve





#-----------------------------------------------------------------------------
#-- B-spline Basis Functions
#-----------------------------------------------------------------------------
require(fda) # used for bspline function
?create.bspline.basis


nkts = 5                              # number of knots
rng = c(0,1)                          # support
kts = seq(rng[1],rng[2],length=nkts)  # knots

#- Degree = 0 (piecewise constant - regressogram)
deg = 0                                      # degree                      
spl = create.bspline.basis(kts,nord=deg+1)   # make spline object
plot(spl)

#- Degree = 1
deg = 1                               
spl = create.bspline.basis(kts,nord=deg+1)   # make spline object
plot(spl)

#- Degree = 2
deg = 2                              
spl = create.bspline.basis(kts,nord=deg+1)   # make spline object
plot(spl)

#- Degree = 3 (cubic splines)
deg = 3
spl = create.bspline.basis(kts,nord=deg+1) 
plot(spl)


#-- Plot B-splines of several degrees
par(mfrow=c(4,1))
nkts = 11
kts = seq(rng[1],rng[2],length=nkts) 
for(deg in 0:3){
  plot(create.bspline.basis(kts,nord=deg+1)) 
  title(paste("B-splines of degree",deg))
}
par(mfrow=c(1,1))



#--------------------
#- B-Spline fits
#--------------------

deg = 3
nkts = 5
rng = range(x)
kts = seq(rng[1],rng[2],length=nkts)  

spl = create.bspline.basis(kts,nord=deg+1)   # make spline object
plot(spl)  
B = eval.basis(x,spl)
theta = solve(crossprod(B)) %*% crossprod(B,y)
      # solve(crossprod(B),crossprod(B,y))     # use solve() directly

plot(x,y); curve(f,add=TRUE,col="grey",lwd=2)
lines(x,B %*% theta)



#-- bs() method. A little different at the right side
# you can use the bs() function to implement b-splines directly into lm() and glm(), etc.
# it is slightly different than use the create.bspline.basis (by one column).
# we are going to add penalties next, so it is better to use the create.bsplines.basis,
#  for what we will do next

library(splines)  # for bs() function
B2 = bs(x,knots=kts[-c(1,nkts)],degree=deg,Boundary.knots=kts[c(1,nkts)])
matplot(x,B2,typ='l',lty=1)

B3 = bs(x,df=5)               # specify 5 degrees of freedom (parameters)
matplot(x,B3,typ='l',lty=1)

model.bs = lm(y ~ bs(x,df=5)-1)  # no intercept

plot(x,y); curve(f,add=TRUE,col="grey",lwd=2)
lines(x,predict(model.bs),col="red")




#-----------------------------------------------------------------------------
#-- Bspline Functions
#-----------------------------------------------------------------------------
library(fda)

######  bspline(x,y,nkts=50,deg=3,rng=range(x))
# Create p.spline object
#  INPUTS
#   x is 1-D data, doesn't need to be sorted
#   y is response data
#   nkts is number of knots. Will not work if nkts>length(x) (or too makes B singluar)
#   deg is the degree of B-spline. Order=deg+1. deg=3 gives cubic splines
#   rng is range of inputs
#  OUTPUTS
#   spl - B-spline basis function
#   beta - the fitted coefficients for the basis functions

bspline <- function(x,y,nkts=50,deg=3,rng=range(x)){
  kts = seq(rng[1],rng[2],length=nkts)                 # Consider extending a bit (i.e. kts*1.1)
  spl = create.bspline.basis(kts,nord=deg+1)           # plot(spl) to see basis functions
  B = eval.basis(x,spl)                                # Basis Matrix
  theta = solve(crossprod(B)) %*% crossprod(B,y)       # estimated coefficients
  object = list(spl=spl,B=B,theta=theta,kts=kts,call=match.call())
  return(object)
}


######  predict.bspline(bspl,x0)
# Predict from a bspline
#  INPUTS
#   bspl - a bspline object 
#   x0 - a vector of locations where you want a prediction
#  OUTPUTS
#   yhat - The estimated responses
predict.bspline <- function(bspl,x0){  
  if(length(x0)>5000){          # Approximate to minimize large matrices
    xx = seq(min(x0),max(x0),length=5000)
    B0 = eval.basis(xx,bspl$spl) 
    yhat = B0 %*% bspl$theta
    yhat = approx(xx,yhat,xout=x0)$y        
  } else{ 
    B0 = eval.basis(x0,bspl$spl) 
    yhat = B0 %*% bspl$theta 
  }    
  return(yhat)
}




#-----------------------------------------------------------------------------
#-- Compare 3 B-splines 
#-----------------------------------------------------------------------------
bspl.0 = bspline(x,y,nkts=5,deg=0)
bspl.5 = bspline(x,y,nkts=5,deg=3)
bspl.50 = bspline(x,y,nkts=50,deg=3)

par(mfrow=c(1,1))
tt = seq(min(x),max(x),length=2000)
plot(x,y); curve(f,add=TRUE,col="grey",lwd=2)
abline(v=bspl.0$kts,lty=3)
lines(tt,predict.bspline(bspl.0,tt),lwd=2)
lines(tt,predict.bspline(bspl.5,tt),col=2,lwd=2)
lines(tt,predict.bspline(bspl.50,tt),col=3,lwd=2)
legend('topright',c("5 kts:  deg=0","5 kts:  deg=3","50 kts: deg=3"),
       col=1:3,lty=1,lwd=2,bty='n',inset=.05)
title("B-spline fits")




#-- Details of basis functions and parameters
bspl = bspline(x,y,nkts=20,deg=3)
BB = with(bspl, sweep(B,2,theta,'*'))

plot(x,y)
matlines(x,BB,col="orange",lty=3)                    # basis functions
lines(x,bspl$B %*% bspl$theta,col="orange",lwd=2)    # fitted line
with(bspl,                                           # add coefficents   
 points(c(kts[1],kts,kts[length(kts)]),theta,pch=19,col=2,cex=1.25)
)







#-----------------------------------------------------------------------------
#-- Penalized B-Splines (PB-splines)
#-----------------------------------------------------------------------------

#-- Difference Matrix

#- pord = 1
D = diag(8)        # 8 basis functions
D = diff(D)        # difference matrix
D
crossprod(D)       # difference penalty

#- pord = 2
D = diag(8)       
for(k in 1:2) D = diff(D)
D
crossprod(D)


#- pord = 3
D = diag(8)       
for(k in 1:3) D = diff(D) 
D
crossprod(D)

#- pord = 0 (ridge penalty)
D = diag(8)       
crossprod(D)



#-- Function for penalized b-spline

p.bspline <- function(x,y,lambda=0,nkts=50,deg=3,pord=1,rng=range(x)){
  kts = seq(rng[1],rng[2],length=nkts)                 # Consider extending a bit (i.e. kts*1.1)
  spl = create.bspline.basis(kts,nord=deg+1)           # plot(spl) to see basis functions
  B = eval.basis(x,spl)                 # Basis Matrix
  D = diag(ncol(B))                     # Matix for penalty
  if(pord>0) for(k in 1:pord) D=diff(D) # difference matrix
  crossB = crossprod(B)
  crossD = crossprod(D)  
  Q = solve(crossB + lambda * crossD) 
  theta = Q %*% crossprod(B,y)          # estimate parameters
  df = sum(diag(Q %*% crossB))          # effective degree of freedom
  object = list(spl=spl,B=B,theta=theta,kts=kts,lambda=lambda,pord=pord,df=df,
                call=match.call())
  return(object)
}



#-- Prediction
pord = 2
deg = 3
ntks = 50
pspl.0 = p.bspline(x, y, lambda=0, pord=pord, nkts=ntks, deg=deg)
pspl.5 = p.bspline(x, y, lambda=5, pord=pord, nkts=ntks, deg=deg)
pspl.500 = p.bspline(x, y, lambda=500, pord=pord, nkts=ntks, deg=deg)


par(mfrow=c(1,1))
tt = seq(min(x),max(x),length=2000)
plot(x,y); curve(f,add=TRUE,col="grey",lwd=2)
lines(tt,predict.bspline(pspl.0,tt),lwd=2)
lines(tt,predict.bspline(pspl.5,tt),col="red",lwd=2)
lines(tt,predict.bspline(pspl.500,tt),col="green",lwd=2)
title("PB-spline fits (pord=2)")
legend("topright",paste("df =",round(c(pspl.0$df,pspl.5$df,pspl.500$df),1)),col=1:3,lwd=2,lty=1)




#-- Details of basis functions and parameters
#   Change lambda and pord and see what happens
pspl = p.bspline(x,y, lambda=1e6, pord=3, nkts=50, deg=3)
BB = with(pspl, sweep(B,2,theta,'*'))

plot(x,y)
matlines(x,BB,col="orange",lty=3)                    # basis functions
lines(x,pspl$B %*% pspl$theta,col="orange",lwd=2)    # fitted line
with(pspl,                                           # add coefficents   
 points(c(kts[1],kts,kts[length(kts)]),theta,pch=19,col=2,cex=1.25)
)





#-- Maximum penalty
#   When the penalty (lambda) gets large, the fit approaches a polynomial
#    of degree pord - 1. (spline degree sufficiently large)

deg = 3         # b-spline degree
ntks = 50       # number of knots
poly.degree = 1 # polynomial degree 

pspl = p.bspline(x,y, lambda=1e8, pord=poly.degree+1, nkts=ntks, deg=deg)
m3 = lm(y~poly(x,poly.degree))       # polynomial fit

tt = seq(min(x),max(x),length=2000)
plot(x,y); curve(f,add=TRUE,col="grey",lwd=2)
lines(tt,predict.bspline(pspl,tt),lwd=2)
lines(tt,predict(m3,newdata=data.frame(x=tt)),col="red",lty=2)






################################################################################
#-- GAM Example
################################################################################

#-- Load Data
url = "http://www.stat.cmu.edu/~larry/all-of-nonpar/=data/rock.dat"
data = read.table(url,header=TRUE)
X.train = data[,1:3]
Y.train = log(data[,4])

#X = X.train
X = scale(X.train)
Y = Y.train - mean(Y.train)
Xnames = colnames(X)

#-- OLS
beta.ls = coef(lm(Y~X))

#-- Ridge Regression (see ridge.R for details on setting lambda)
library(MASS)
beta.ridge = coef(lm.ridge(Y~X,lambda=12.3))

#-- GAM
library(mgcv)
fmla = as.formula(Y~paste("s",Xnames))
terms = paste(paste0("s(",Xnames,",sp=.001)"),collapse='+')
fmla = as.formula(paste("Y~",terms))
fit = gam(fmla,data=data.frame(Y,X))



#-- Component Plots
par(mfrow=c(3,1),mar=c(4,4.5,1,1),oma=c(0,0,2,0))
for(j in 1:ncol(X)){
  ord = order(X[,j])
  xs = X[ord,j]
  yhat.ls =  xs * beta.ls[j+1] #+ beta.ls[1] 
  yhat.ridge = xs * beta.ridge[j+1] #+ beta.ridge[1]
  yrng = c(-4,4)
  plot(range(xs),yrng,typ='n',
       xlab=paste(Xnames[j],"(scaled)"),ylab=bquote(f[.(j)]),las=1)
  abline(h=0,col='grey85')
  par(new=TRUE);plot(fit,select=j,ylim=yrng,yaxt='n',xaxt='n',xlab='',ylab='',
                     shade=TRUE,lwd=2,col=1)  
  points(xs,Y[ord])
  lines(xs,yhat.ls,col=4,lwd=2)
  #lines(xs,yhat.ridge,col=2,lwd=2)
  # rug(xs)
}
legend(.5,29.5, legend=c("OLS","GAM"), col=c(4,1), lwd=2,lty=1, horiz=TRUE,
       xpd=NA,bty='n',cex=1.25,xjust=.5)
