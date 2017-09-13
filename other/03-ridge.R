#### Ridge Regression ####

library(MASS)
library(glmnet)


#-- Correlated Data Example
library(MASS) # for lm.ridge function
set.seed(10)
n = 25
x1 = rnorm(n)
x2 = rnorm(n,mean=x1,sd=.01)
cor(x1,x2)

y = rnorm(n, mean=1+x1+x2,sd=2)  # m(x) = 1 + 1x_1 + 1x_2

#- OLS
lm(y~x1+x2)$coef

#- ridge  
lm.ridge(y~x1+x2,lambda=.5)   # set lambda=0.5



#-- Fit ridge regression model for sequence of lambdas

lam.seq = exp(seq(log(100),log(1e-5),length=500))
m = lm.ridge(y~x1+x2, lambda=lam.seq)   # ridge regression model

beta = coef(m)                   # matrix of estimated coefficients (for each lambda)
penalty = rowSums(beta[,-1]^2)   # total penalty P(\beta) (sum of squared betas)

ridge.stuff = data.frame(lam = lam.seq, intercept=beta[,1], x1=beta[,2],
                         x2=beta[,3], penalty=penalty, GCV=m$GCV, 
                         row.names=NULL)
head(ridge.stuff)



#-- ridge path for sequence of lambdas

matplot(log(m$lambda),coef(m)[,2:3],typ='l',lty=1,xlab="log(lambda)",ylab="beta",las=1)
legend("topright",c("beta_1","beta_2"),col=1:2,lty=1)
abline(h=1,lty=1,col="lightgray")

#-- penalty for sequence of lambdas
matplot(log(m$lambda),penalty,typ='l',lty=1,xlab="log(lambda)",ylab="Pen(beta)",las=1)
abline(h=1,lty=1,col="lightgray")

#-- GCV for sequence of lambdas
matplot(log(m$lambda),m$GCV,typ='l',lty=1,xlab="log(lambda)",ylab="GCV",las=1)


#-- ridge with glmnet. Notice: no formula interface. 

library(glmnet)
m2 = glmnet(cbind(x1, x2), y, alpha=0, lambda=lam.seq)

#- ridge path vs. lambda 
plot(m2, "lambda", las=1); abline(h=1,lty=1,col="lightgray")

#- ridge path vs. L1 norm of penalty
plot(m2, "norm", las=1); abline(h=1,lty=1,col="lightgray")
