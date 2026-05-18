load(file = "introduction_to_gwas/model_extensions/multi_trait/multi-trait-fit.RData")

fit = to_save$fit
K = to_save$K
X = to_save$X

## genetic vars
genetic_var <- diag(fit$sigma$`u:id`)

# Residual variances
residual_var <- diag(fit$sigma$`u:units`)

## genetic covariance
cov_SL_SW = fit$sigma$`u:id`["SL","SW"]

## residual covariance
cov_res <- fit$sigma$`u:units`["SL","SW"]

cor_SL_SW <- cov_SL_SW / sqrt(genetic_var["SL"] * genetic_var["SW"])


res_cor <- cov_res / sqrt(residual_var["SL"] * residual_var["SW"])


## getting the X'(XX')^-1 matrix
Kinv <- solve(K)
Kinv <- solve(K + diag(rnorm(ncol(K))/1e+3, ncol(K), ncol(K)))

Xc = scale(X, center = TRUE, scale = FALSE) 

XKinv <- t(Xc)%*%Kinv

## individual genetic effects
g_sl <- fit$U$`u:id`$SL

## 1) get SNP effects --> β_SNP = X * K^-1 * g
beta_hat_sl <- XKinv %*% g_sl ## marker effects

## 2) get the standard errors: Var(β_SNP)=X(K^−1)Var(g)(K−1)t(X)
varg_sl = fit$VarU$`u:id`$SL
cov_beta_sl <- t(X) %*% Kinv %*% varg_sl %*% t(Kinv) %*% X
var_beta_sl = diag(cov_beta_sl)

# Variance of SNP effects
# PEV_u <- fit$PevU$`u:id`$SL
# Var_beta <- t(X) %*% Kinv %*% PEV_u %*% Kinv %*% X

# Standard errors
se_beta <- sqrt(var_beta_sl)

## 3) t-statistc = β_SNP  / SE(β_SNP)  
tstat = beta_hat_sl/se_beta

## 4) p-values (two-sided)
p_values <- 2 * (1 - pnorm(abs(tstat)))
hist(p_values)
