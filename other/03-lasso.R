library(glmnet)
library(MASS)


#-------------------------------------------------------------------------------
#-- Prostate Data - Comparing Lasso and Ridge
#-------------------------------------------------------------------------------
url = "https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data"
data = read.table(url,header=TRUE)

data.train = subset(data,train==TRUE)
X.train = as.matrix(data.train[,1:8])
Y.train = as.matrix(data.train[,"lpsa"])


#- Center and Scale X
X = scale(X.train)            # Center and scale predictors
Y = Y.train


#-- Fit LS model
fit.ls = lm(Y~X)


#-- Fit Lasso Model
library(glmnet)
fit.lasso = glmnet(X,Y,alpha=1)
B.lasso = t(coef(fit.lasso))   # coefficients
fit.lasso$penalty = rowSums(abs(B.lasso[,-1]))



#-- Fit Ridge Model
library(MASS)
lamseq = c(0,exp(seq(log(.01),log(2000),length=1000))) # exponential sequence
fit.ridge = lm.ridge(Y~X, lambda=lamseq )
B.ridge = coef(fit.ridge)             # coefficients
fit.ridge$penalty = rowSums(B.ridge[,-1]^2)  
XTX = crossprod(X)
fit.ridge$df = sapply(fit.ridge$lambda, 
                      function(l) sum(diag(solve(XTX + diag(l,ncol(X))) %*% XTX)))


#-- Plot Both Models
# pdf('lectures/figs/lasso-ridge.pdf', width=2*4, height=2*3)
par(mfrow=c(1,2))

#- Lasso Fit
xrng = range(fit.lasso$penalty) + c(0,.22)
yrng = range(B.lasso[,-1])
plot(xrng, yrng,typ='n', lty=1, las=1,
        xlab='L1 norm: sum of absolute betas', ylab='beta (standardized)')
abline(h=0, col='grey80')
abline(v=axTicks(1), col='grey80')
matlines(fit.lasso$penalty, B.lasso[,-1], lty=1)
text(max(fit.lasso$penalty), tail(B.lasso[,-1],1), labels=colnames(X),pos=4,cex=.8)
axis(3,axTicks(1), round(approx(fit.lasso$pen, fit.lasso$df, axTicks(1))$y,1))
mtext('Lasso',side=3,line=2.5)

#- Ridge Fit
xrng = range(fit.ridge$penalty) + c(0,.09)
yrng = range(B.ridge[,-1])
plot(xrng, yrng, typ='n', lty=1, las=1,
        xlab='L2 norm: sum of squared betas', ylab='beta (standardized)')
abline(h=0,col='grey80')
abline(v=axTicks(1),col='grey80')
matlines(fit.ridge$pen, B.ridge[,-1],lty=1)
text(max(fit.ridge$pen), B.ridge[1,-1],labels=colnames(X),pos=4,cex=.8)
axis(3,axTicks(1),round(approx(fit.ridge$pen,fit.ridge$df,axTicks(1))$y,1))
mtext('Ridge',side=3,line=2.5)
# dev.off()

#-- MSE vs. EDF
# pdf('lectures/figs/lasso-ridge-mse.pdf', width=2*4, height=2*3)
par(mfrow=c(1, 1))

yhat.lasso = predict(fit.lasso, newx = X)
r = -sweep(yhat.lasso, 1, Y, '-')
mse.lasso = colMeans(r^2)
plot(fit.lasso$df, mse.lasso, xlab="edf", ylab="mse (training)", las=1,
     main='MSE vs. EDF (not including intercept)')

yhat.ridge = cbind(1,X) %*% t(B.ridge)
r = -sweep(yhat.ridge, 1, Y, '-')
mse.ridge = colMeans(r^2)
lines(fit.ridge$df, mse.ridge, col="red")
legend("topright",c("lasso","ridge"),col=c(1,2),lty=c(NA,1),pch=c(1,NA))
# dev.off()


#-------------------------------------------------------------------------------
#-- Elastic Net
#-------------------------------------------------------------------------------

#-- Elastic Net (alpha=0 is ridge, alpha=1 is lasso)
fit.enet = glmnet(X, Y, alpha=0.5)

# pdf('lectures/figs/lasso-enet.pdf', width=2*4, height=2*3)
par(mfrow=c(1,1))
plot(fit.lasso, col="lightgray")
plot(fit.enet, add=TRUE)
legend("topleft",c("lasso"),col=c("gray"),lty=1)
# dev.off()



#-------------------------------------------------------------------------------
#-- Soft Threshold
#-------------------------------------------------------------------------------

soft <- function(b,thres) sign(b)*( abs(b) - thres)*( abs(b) > thres)

b = seq(-3,4,length=50)
plot(b,soft(b,thres=1),main="soft-threshold function lambda=1")
abline(h=0,col="lightgray")







