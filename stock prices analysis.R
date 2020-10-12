

lassoISTA <- function(X, Y, lambda){
  THRESH = 1e-15
  p = ncol(X)
  n = nrow(X)
  L = 2*eigen(t(X) %*% X)$values[1]
  
  XtX = t(X) %*% X
  XtY = t(X) %*% Y
  
  beta = rep(0, p)
  while (TRUE){
    beta_new = softThresh( beta - (2/L)*(XtX %*% beta - XtY), lambda/L)
    if (sqrt(sum((beta_new - beta)^2)) < THRESH)
      break
    else
      beta = beta_new
  }
  return(beta_new)
}


softThresh <- function(u, lambda){
  u[abs(u) <= lambda] = 0
  u[u > lambda] = u[u > lambda] - lambda
  u[u < -lambda] = u[u < -lambda] + lambda
  return(u)
}

###################
## 

library(glmnet)

stocks_df = read.csv("sp500_sm.csv")
names(stocks_df)

dep_comp = "TGT"
indep_comps = setdiff(names(stocks_df), c("date", dep_comp))

n = nrow(stocks_df)
X = cbind(as.matrix(stocks_df[, indep_comps]))
X = scale(X)
Y = as.vector(stocks_df[, dep_comp])
Y = Y - mean(Y)

lambda = 0.005*n

beta_lasso = lassoISTA(X, Y, lambda)
beta_glmnet = glmnet(X, Y, lambda = lambda/(2*n))$beta

glmnet_deviation = sum((beta_lasso - beta_glmnet)^2) / 
  min(sum(beta_lasso^2), sum(beta_glmnet^2))

lasso_selected = indep_comps[which(abs(beta_lasso) > 1e-10)]

beta_ols = solve((t(X) %*% X), t(X) %*% Y)
ols_top5 = indep_comps[order(abs(beta_ols), decreasing=TRUE)[1:5]]


print(sprintf("Relative deviation between glmnet soln and ISTA soln: %f", 
              glmnet_deviation))
print(paste("Lasso selected predictors: ", toString(lasso_selected)))
print(paste("OLS top 5 predictors: ", toString(ols_top5)))
