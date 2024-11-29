# if(TRUE){args=(commandArgs(TRUE))
# for(i in 1:length(args)){eval(parse(text=args[[i]]))}}else{
#   jobID=1
# }

### libraries required
require(MASS)
require(rstanarm)
require(glmnet)
require(caret)
require(monomvn)
require(caret)
require(BoomSpikeSlab)
require(SSLASSO)
require(bayesreg)
require(horseshoe)
require(susieR)
require(randomForest)
require(ppcor)
require(doParallel)
library(ranger)
library(Rfast)
###used functions
MSEY<-function(Y,X,B){print(sum((Y-X%*%B)^2)/nrow(Y))}
MSEB<-function(beta,B){print(sum((beta-B)^2))}
confusion_measure<-function(beta,B){
  B<-relevel(as.factor(ifelse(round(B,2) == 0, 0, 1)),ref="1")
  beta_<-ifelse(round(beta,2) == 0, 0, 1)
  if(sum(beta_)==0|sum(beta_)==length(beta_)){
    TP<-NA
    FP<-NA
    TN<-NA
    FN<-NA
    Accu<-NA
    prec<-NA
    recal<-NA
    Fscore<-NA
  }else{
    beta<-relevel(as.factor(ifelse(beta == 0, 0, 1)),ref="1")
    TP<-confusionMatrix(beta,B)$table[1,1]
    FP<-confusionMatrix(beta,B)$table[1,2]
    TN<-confusionMatrix(beta,B)$table[2,2]
    FN<-confusionMatrix(beta,B)$table[2,1]
    Accu<-(TP+TN)/(TP+TN+FP+FN)
    prec<-(TP)/(TP+FP)
    recal<-(TP)/(TP+FN)
    Fscore<-(2*prec*recal)/(prec + recal)
  }
  return(list(TP=TP,FP=FP,TN=TN,FN=FN,Accu=Accu,prec=prec,recal=recal,Fscore=Fscore))
}


# Function to calculate metrics (TP, TN, FP, FN, FDR, FNR, F1 Score)
calculate_metrics <- function(beta_hat, beta_true) {
  # True Positives: Predicted non-zero and true non-zero
  true_positives <- sum((beta_hat != 0) & (beta_true != 0))
  
  # False Positives: Predicted non-zero and true zero
  false_positives <- sum((beta_hat != 0) & (beta_true == 0))
  
  # False Negatives: Predicted zero and true non-zero
  false_negatives <- sum((beta_hat == 0) & (beta_true != 0))
  
  # True Negatives: Predicted zero and true zero
  true_negatives <- sum((beta_hat == 0) & (beta_true == 0))
  
  # FDR: False Discovery Rate
  fdr <- false_positives / (true_positives + false_positives)
  fdr <- ifelse(is.nan(fdr), 0, fdr)  # Handle NaN cases
  
  # FNR: False Negative Rate
  fnr <- false_negatives / (true_positives + false_negatives)
  fnr <- ifelse(is.nan(fnr), 0, fnr)  # Handle NaN cases
  
  # Precision: True Positives / (True Positives + False Positives)
  precision <- true_positives / (true_positives + false_positives)
  precision <- ifelse(is.nan(precision), 0, precision)  # Handle division by zero
  
  # Recall: True Positives / (True Positives + False Negatives)
  recall <- true_positives / (true_positives + false_negatives)
  recall <- ifelse(is.nan(recall), 0, recall)  # Handle division by zero
  
  # F1 Score: 2 * (Precision * Recall) / (Precision + Recall)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  f1_score <- ifelse(is.nan(f1_score), 0, f1_score)  # Handle NaN cases
  
  return(list(TP = true_positives, FP = false_positives, TN = true_negatives, FN = false_negatives,
              FDR = fdr, FNR = fnr, Precision = precision, Recall = recall, F1 = f1_score))
}
#set.seed(123)

#fn_iteration<-function(Sim){

#  for(iterr in 1:Sim){
#  set.seed(iterr)
# Function to create block-exchangeable covariance matrix
generate_covariance_matrix <- function(p,q, alpha) {
  alpha1 <- alpha[1]
  alpha2 <- alpha[2]
  alpha3 <- alpha[3]
  
  # Function to create R11 or R22 block matrices
  create_block_matrix <- function(size, diagonal_value, off_diagonal_value) {
    mat <- matrix(off_diagonal_value, nrow = size, ncol = size)
    diag(mat) <- diagonal_value
    return(mat)
  }
  
  # Create R11 and R22 (block matrices)
  R11 <- create_block_matrix(q, diagonal_value = 1, off_diagonal_value = alpha1)
  R22 <- create_block_matrix(p-q, diagonal_value = 1, off_diagonal_value = alpha3)
  
  # Create R12 matrix (all elements equal to alpha2)
  R12 <- matrix(alpha2, nrow = q, ncol = p-q)
  
  # Construct the full covariance matrix
  Covariance_matrix <- rbind(
    cbind(R11, R12),
    cbind(t(R12), R22)
  )
  
  return(Covariance_matrix)
}


#Allmodels<-function(Sim){ 
  ###data generation
  p <- 5000   # number of variables
  spars_list<-list(c(0.03),c(0.1),c(0.3))
  spars<-spars_list[[1]]
  q <- p*spars  # number of nonzero coefficients
  n <- 200    # number of observations
  outn=50
  ps<-q
  #spars<-ps/p
  #beta_s <- 0.4   # non-zero coefficient value
  beta<-numeric(p)
  beta[1:ps] <-c(rep(1,(ps/2)),rep(-1,(ps/2)))
  alpha_list <- list(c(0.3, 0.3, 0.3), c(0.5, 0.5, 0.5), c(0.7, 0.7, 0.7), c(0.6, 0.7, 0.9),c(0.8,0.5,0.9),c(0.3,0.3,0.3))
  alpha<-alpha_list[[1]]
  
  cov_matrix <- generate_covariance_matrix (p,q, alpha)
  cholemy<-chol(cov_matrix)
  
  for (Sim in 1:100) {
  iterr=Sim
  cat("simulation--:", iterr, "\n")
  set.seed(iterr)
  
  Z <- matrix(rnorm(n * p), ncol = p)  # Standard normal variables
  X <- Z%*%t(cholemy) #sweep(Z %*% t(cholemy), 2, rep(0,p), "+")  # Transform Z to have mean mu and covariance A
  #X <- mvrnorm(n, mu = rep(0, p), Sigma = cov_matrix)
  # beta <- c(rep(beta_s, q), rep(0, p - q))
  Y <- X %*% beta + rnorm(n)
  
  X_scale=X
  Y_scale=Y-mean(Y)
  #out of sample generation
  outX<- matrix(rnorm(outn * p), nrow = outn)
  outY <- outX %*% beta + rnorm(outn)
  outX_scale=outX
  outY_scale=outY-mean(outY)
  
  #lasso
  print("lasso")
  p=ncol(X_scale)
  n=nrow(X_scale)
  start_time <- Sys.time()
  #lasso_model <- glmnet(X_scale, Y_scale, alpha = 1)
  cv_model <- cv.glmnet(X_scale, Y_scale, alpha = 1)
  lambda_optimal <- cv_model$lambda.min
  lasso_model_optimal <- glmnet(X_scale, Y_scale, alpha = 1, lambda = lambda_optimal)
  end_time <- Sys.time()
  runtime1=difftime(end_time, start_time, units="secs")
  runtime1
  SSE1<-sum((beta-lasso_model_optimal$beta)^2)/p
  MSPE1<-sum((X_scale%*%lasso_model_optimal$beta-Y_scale)^2)/(n)
  MSPE1out<-sum((outX_scale%*%lasso_model_optimal$beta-outY_scale)^2)/(outn)
  
  out1<-list(runtime=runtime1,calculate_metrics(round(lasso_model_optimal$beta,2),round(beta,2)),SSE=SSE1,MSPE=MSPE1,MSPEout=MSPE1out)
  
  ##elasticnet
  print("EL")
  start_time <- Sys.time()
  eln_model <- glmnet(X_scale, Y_scale, alpha = 0.5)
  cv_model1 <- cv.glmnet(X_scale, Y_scale, alpha = 0.5)
  lambda_optimal <- cv_model1$lambda.min
  eln_model_optimal <- glmnet(X_scale, Y_scale, alpha = 0.5, lambda = lambda_optimal)
  end_time <- Sys.time()
  runtime2=difftime(end_time, start_time, units="secs")
  runtime2
  SSE2<-sum((beta-eln_model_optimal$beta)^2)/p
  MSPE2<-sum((X_scale%*%eln_model_optimal$beta-Y_scale)^2)/(n)
  MSPE2out<-sum((outX_scale%*%eln_model_optimal$beta-outY_scale)^2)/(outn)
  out2<-list(runtime=runtime2,calculate_metrics(round(eln_model_optimal$beta,2),round(beta,2)),SSE=SSE2,MSPE=MSPE2,MSPEout=MSPE2out)
  
  
  ###susie
  print("su")
  Y1<-drop(Y_scale)
  start_time <- Sys.time()
  res=susie(X_scale,Y1,L=10,intercept=FALSE,max_iter=5000)
  end_time <- Sys.time()
  runtime10=difftime(end_time, start_time, units="secs")
  runtime10
  SSE10<-sum((beta-susie_get_posterior_mean(res))^2)/p
  MSPE10<-sum((X_scale%*%susie_get_posterior_mean(res)-Y1)^2)/n
  MSPE10out<-sum((outX_scale%*%susie_get_posterior_mean(res)-outY_scale)^2)/(outn)
  out10<-list(runtime=runtime10,calculate_metrics(round(susie_get_posterior_mean(res),2),round(beta,2)),SSE=SSE10,MSPE=MSPE10,MSPEout=MSPE10out)
  
  #### random forest
  ## Two-step procedure: randomforest to screen a set of features
  ## and then apply the stepwise variable selection technique
  ## using the BIC criteria
  ## Step one
  
  print("TWIST")
  start_time <- Sys.time()
  mydata=data.frame(Y_scale, X_scale)
  #outrf=randomForest(Y_scale~., ntree=1000,mtry=500, data=mydata)
  outrf <- ranger(
    formula = Y_scale ~ .,               # Response variable and predictors
    data = mydata,                       # Your dataset
    mtry = 100,                          # Number of variables tried at each split
    num.trees = 500,                     # You can adjust the number of trees as needed
    min.node.size = 5,                   # Minimum size of terminal nodes
    splitrule =   "variance",              # Splitting rule for regression
    importance = "impurity"              # Optionally, you can specify how to calculate feature importance
  )
  impoutrf=importance(outrf)
  impvariables=sort(impoutrf, index.return=T, decreasing=TRUE)$ix
  mythreshold=min((n-2), 0.1*p) # to screen out important variables
  selectedvariablerf=impvariables[1:mythreshold]
  
  ## Step two
  mydata.rf=data.frame(X_scale[,selectedvariablerf], Y_scale)
  names(mydata.rf)<-c(selectedvariablerf, "response")
  full.rf=lm(response~., data=mydata.rf)
  final.rf <- step(full.rf, k=log(n), trace=0)# using the BIC criteria
  end_time <- Sys.time()
  runtime12=difftime(end_time, start_time, units="secs")
  summary(final.rf)
  sig.genes.rf=noquote(noquote(names(final.rf$coefficients[-1])))
  sig.genes.rf=as.numeric(gsub("`", "",sig.genes.rf))
  final_rf_coef<-numeric(p)
  final_rf_coef[sig.genes.rf] <- as.numeric(final.rf$coefficients[-1])
  # sig.genes.rf=noquote(noquote(names(full.rf$coefficients[-1])))
  # sig.genes.rf=as.numeric(gsub("`", "",sig.genes.rf))
  # final_rf_coef<-numeric(p)
  # final_rf_coef[sig.genes.rf] <- as.numeric(full.rf$coefficients[-1])
  
  SSE12<-sum((beta-final_rf_coef)^2)/p
  MSPE12<-sum((X_scale%*%final_rf_coef-Y_scale)^2)/(n)
  MSPE12out<-sum((outX_scale%*%final_rf_coef-outY_scale)^2)/(outn)
  #out12<-list(runtime=runtime12,confusion_measure(final_rf_coef,beta),SSE=SSE12,MSPE=MSPE12,MSPEout=MSPE12out)
  out12<-list(runtime=runtime12,calculate_metrics(round(final_rf_coef,2),round(beta,2)),SSE=SSE12,MSPE=MSPE12,MSPEout=MSPE12out)
  ###combining results for one simulation
  #setwd("~/2024/spanew") 
  df<-cbind(unlist(out1),unlist(out2),unlist(out10),unlist(out12))
  colnames(df)<-c("lasso","EL","SUSIE","RF")
  write.csv(df,paste("Allsap_","_Sim_",Sim,"_p_",p,"_spars_",spars,"_cov_",alpha[1],".csv",sep=""))
  beta_estimates<-as.data.frame(cbind(c(beta),c(susie_get_posterior_mean(res)),
                                      c(final_rf_coef),c(lasso_model_optimal$beta[,1]),c(eln_model_optimal$beta[,1])))
  colnames(beta_estimates)<-c("T","SUSIE","RF","lasso","EL")
  write.csv(beta_estimates,paste("beta_estimatessap_",Sim,"_p_",p,"_spars_",spars,"_cov_",alpha[1],".csv",sep=""))
  }
  
#}

#ans = fn_iteration(1)
#ans



# Sim=60
# registerDoParallel(cores=8)
# results<-foreach(i=1:Sim)%dopar%{
#   set.seed(i)
#   Allmodels(Sim=i)
# }



