####  03-subset selection  ####
library(tidyverse)


#-- read in prostate data
library(readr)
data.url = 'https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data'
prostate = readr::read_tsv(data.url)    # tab separated
train = prostate$train
test = !train

#-- predictors
preds = c("lcavol", "lweight", "age", "lbph", "svi", "lcp", "gleason", "pgg45")
response = "lpsa"
p = length(preds)


#-- fit best subset and eval on train and test sets
out = NULL
for(k in 1:p){
  var.mat = combn(preds, k)
  nk = ncol(var.mat)
  for(i in 1:nk){
    vars = c(response, var.mat[,i])
    m = lm(lpsa~., data=prostate[train, vars])
    MSE.train = mean(m$residuals^2)
    yhat.test = predict(m, newdata=prostate[test,  ])
    MSE.test = mean((prostate[test, response] - yhat.test)^2)
    out = rbind(out, data.frame(k, MSE.train, MSE.test))
  }
}



#-- Plot training MSE
library(dplyr)
yrng = range(c(out$MSE.train, out$MSE.test))
best.k = out %>% group_by(k) %>% summarize(MSE = min(MSE.train))
plot(out$k, out$MSE.train, 
     pch=19, col="gray", 
     xlab='k', ylab='MSE.train', 
     las=1, ylim=yrng)
lines(best.k$k, best.k$MSE, col='red', type='b')


#-- Plot test MSE
best.k = out %>% group_by(k) %>% summarize(MSE = min(MSE.test))
plot(out$k, out$MSE.test, 
     pch=19, col="gray", 
     xlab='k', ylab='MSE.test', 
     las=1, ylim=yrng)
lines(best.k$k, best.k$MSE, col='red', type='b')
