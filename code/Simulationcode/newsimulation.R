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
###used functions
MSEY<-function(Y,X,B){print(sum((Y-X%*%B)^2)/nrow(Y))}
MSEB<-function(beta,B){print(sum((beta-B)^2))}
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

  
Allmodels<-function(p,spars,response,SNR,Sim){ 
  ###data generation
  iterr=Sim
  print(iterr)
  outn=50
  n=200
  p=1000;spars=0.03;SNR=5;response="C"
  ps<-ceiling(spars*p) #significant beta
  beta<-numeric(p)
  beta[1:ps] <-c(rep(1,(ps/2)),rep(-1,(ps/2)))
  set.seed(1)
  X <- matrix(rnorm(10000 * p), nrow = 10000)
  s1<-sd(X %*% beta)
  cat("sd:", s1, "\n")
  
  set.seed(iterr)
  
  if(response=="B"){
    X <- matrix(rbinom(n * p,1,prob=0.7)*(-2)+1, nrow = n)}else{
      X <- matrix(rnorm(n * p), nrow = n)
    }
  
  Y <- X %*% beta + rnorm(n,sd=s1/sqrt(SNR))
  
  X_scale=X
  Y_scale=Y-mean(Y)
  
  #out of sample generation
  outX<- matrix(rnorm(outn * p), nrow = outn)
  outY <- outX %*% beta + rnorm(outn,sd=s1/sqrt(SNR))
  outX_scale=outX
  outY_scale=outY-mean(outY)
  
  #lasso
  p=ncol(X_scale)
  n=nrow(X_scale)
  start_time <- Sys.time()
  #lasso_model <- glmnet(X_scale, Y_scale, alpha = 1)
  cv_model <- cv.glmnet(X_scale, Y_scale, alpha = 1)
  lambda_optimal <- cv_model$lambda.min
  lasso_model_optimal <- glmnet(X_scale, Y_scale, alpha = 1, lambda = lambda_optimal)
  end_time <- Sys.time()
  runtime1=difftime(end_time, start_time, units="secs")
  SSE1<-sum((beta-lasso_model_optimal$beta)^2)/p
  MSPE1<-sum((X_scale%*%lasso_model_optimal$beta-Y_scale)^2)/(n)
  MSPE1out<-sum((outX_scale%*%lasso_model_optimal$beta-outY_scale)^2)/(outn)
  out1<-list(runtime=runtime1,calculate_metrics(lasso_model_optimal$beta,beta),SSE=SSE1,MSPE=MSPE1,MSPEout=MSPE1out)

  ##elasticnet
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
  out2<-list(runtime=runtime2,calculate_metrics(eln_model_optimal$beta,beta),SSE=SSE2,MSPE=MSPE2,MSPEout=MSPE2out)

  
  ##alasso
  start_time <- Sys.time()
  # Step 1: Initial Lasso to get coefficients
  init_ridge <- glmnet(X_scale, Y_scale, alpha = 0) #ridge
  beta_init <- as.vector(init_ridge$beta[, which.min(init_ridge$lambda)])
  # Step 2: Define adaptive weights
  # Avoid division by zero (or near-zero) by setting a small epsilon
  epsilon <- 1e-6
  weights <- 1 / (abs(beta_init) + epsilon) 
  # Step 3: Fit Adaptive Lasso using weights
  cv_alasso <- cv.glmnet(X_scale, Y_scale, alpha = 1, penalty.factor = weights)
  lambda_alasso <- cv_alasso$lambda.min
  alasso_fit <- glmnet(X_scale, Y_scale, alpha = 1, lambda = lambda_alasso, penalty.factor = weights)
  end_time <- Sys.time()
  runtime3 <- difftime(end_time, start_time, units = "secs")
  # Step 4: Performance evaluation
  p <- ncol(X_scale)
  n <- nrow(X_scale)
  outn <- nrow(outX_scale)
  
  SSE3 <- sum((beta - alasso_fit$beta)^2) / p
  MSPE3 <- sum((X_scale %*% alasso_fit$beta - Y_scale)^2) / n
  MSPE3out <- sum((outX_scale %*% alasso_fit$beta - outY_scale)^2) / outn
  out3 <- list(runtime = runtime3, calculate_metrics(alasso_fit$beta, beta),
               SSE = SSE3,MSPE = MSPE3,MSPEout = MSPE3out)
  
  ##sslasso
  start_time <- Sys.time()
  ssL1<-SSLASSO(X_scale,Y_scale,variance="unknown",max.iter=5000)
  lambda0_start <- which(ssL1$iter < 100)[1]
  lambda0 <- seq(lambda0_start, 100, by = 1)
  ssL1<- SSLASSO(X_scale, Y_scale, variance = "unknown", lambda1 = 1, lambda0 = lambda0)
  end_time <- Sys.time()
  runtime4=difftime(end_time, start_time, units="secs")
  runtime4
  SSE4<-sum((beta-ssL1$beta[,1])^2)/p
  MSPE4<-sum((X_scale%*%ssL1$beta[,1]-Y_scale)^2)/n
  MSPE4out<-sum((outX_scale%*%ssL1$beta[,1]-outY_scale)^2)/(outn)
  out4<-list(runtime=runtime3,calculate_metrics(ssL1$beta[,1],beta),SSE=SSE4,MSPE=MSPE4,MSPEout=MSPE4out)

  ##horseshoe (HS) ####
  start_time <- Sys.time()
  model1<-horseshoe::horseshoe(Y_scale,X_scale,method.tau = "halfCauchy",method.sigma = "Jeffreys",burn=1000,nmc=5000)
  end_time <- Sys.time()
  runtime5=difftime(end_time, start_time, units="secs")
  runtime5
  SSE5<-sum((beta-model1$BetaHat)^2)/p
  MSPE5<-sum((X_scale%*%model1$BetaHat-Y_scale)^2)/n
  MSPE5out<-sum((outX_scale%*%model1$BetaHat-outY_scale)^2)/(outn)
  out5<-list(runtime=runtime5,calculate_metrics(model1$BetaHat,beta),SSE=SSE5,MSPE=MSPE5,MSPEout=MSPE5out)

  ##horseshoe+ (HSP)
  start_time <- Sys.time()
  hsp<-bayesreg(Y_scale~X_scale,as.data.frame(X_scale,Y_scale),model = "gaussian",prior="hs+",n.samples=5000,burnin=1000)
  end_time <- Sys.time()
  runtime6=difftime(end_time, start_time, units="secs")
  runtime6
  SSE6<-sum((beta-hsp$mu.beta)^2)/p
  MSPE6<-sum((X_scale%*%hsp$mu.beta-Y_scale)^2)/n
  MSPE6out<-sum((outX_scale%*%hsp$mu.beta-outY_scale)^2)/(outn)
  out6<-list(runtime=runtime6,calculate_metrics(hsp$mu.beta,beta),SSE=SSE6,MSPE=MSPE6,MSPEout=MSPE6out)

  
  ##Spike-and-Normal (SN)
  start_time <- Sys.time()
  prior_ssN1<-SpikeSlabPrior(x=X_scale,y=Y_scale, expected.model.size=ps,
                             diagonal.shrinkage = 0.5,   # use Zellner's prior
                             prior.df = spars)#shrink to zero
  ssN1 <- lm.spike(Y_scale ~ X_scale - 1, niter = 10000,prior=prior_ssN1,
                   error.distribution = "gaussian",ping=500,seed=iterr)
  end_time <- Sys.time()
  runtime7=difftime(end_time, start_time, units="secs")
  runtime7
  s <- summary(ssN1, order = FALSE)$coefficients
  pip_sn    <- s[, "inc.prob"]
  cond_mean <- s[, "mean.inc"]
  ssbeta<-numeric(p)
  ssbeta[pip_sn >= 0.5] <- cond_mean[pip_sn >= 0.5]
  SSE7<-sum((beta-ssbeta)^2)/p
  MSPE7<-sum((X_scale%*% ssbeta-Y_scale)^2)/n
  MSPE7out<-sum((outX_scale%*% ssbeta-outY_scale)^2)/(outn)
  out7<-list(runtime=runtime7,calculate_metrics(ssbeta,beta),SSE=SSE7,MSPE=MSPE7,MSPEout=MSPE7out)

  ##Bayesian Lasso (BL)
  start_time <- Sys.time()
  model_bl<-blasso(X_scale,Y_scale,T=200,case = "default",icept = FALSE)
  end_time <- Sys.time()
  runtime8=difftime(end_time, start_time, units="secs")
  runtime8
  SSE8<-sum((beta-colMeans(model_bl$beta))^2)/p
  MSPE8<-sum((X_scale%*%colMeans(model_bl$beta)-Y_scale)^2)/n
  MSPE8out<-sum((outX_scale%*%colMeans(model_bl$beta)-outY_scale)^2)/(outn)
  out8<-list(runtime=runtime8,calculate_metrics(colMeans(model_bl$beta),beta),SSE=SSE8,MSPE=MSPE8,MSPEout=MSPE8out)

  ##Regularized HS (RHS: rstanarm-hs)
  D <- ncol (X_scale)
  n <- nrow (X_scale)
  p0 <- ps # prior guess for the number of relevant variables
  tau0 <- p0 /(D-p0) / sqrt (n) # rstanarm will scale this by sigma automatically
  prior4<-hs(global_scale = tau0,slab_scale=2.5, slab_df=4)
  options (mc.cores = parallel :: detectCores ())
  start_time <- Sys.time()
  fit <- stan_glm (Y_scale ~-1+X_scale, family = gaussian (), data = data.frame (I(X_scale),Y_scale),
                   prior = prior4,
                   chains = 4,iter = 500,sparse = TRUE,seed=iterr)
  end_time <- Sys.time()
  runtime9=difftime(end_time, start_time, units="secs")
  SSE9<-sum((beta-fit$coefficients)^2)/p
  MSPE9<-sum((X_scale%*%fit$coefficients-Y_scale)^2)/n
  MSPE9out<-sum((outX_scale%*%fit$coefficients-outY_scale)^2)/(outn)
  out9<-list(runtime=runtime9,calculate_metrics(fit$coefficients,beta),SSE=SSE9,MSPE=MSPE9,MSPEout=MSPE9out)
 
   ###susie
   Y1<-drop(Y_scale)
   start_time <- Sys.time()
   res=susie(X_scale,Y1,L=ps,intercept=FALSE,max_iter=5000)
   end_time <- Sys.time()
   runtime10=difftime(end_time, start_time, units="secs") 
   runtime10
   pip_susie <- as.numeric(susie_get_pip(res))
   beta_susie<-numeric(p)
   beta_susie[pip_susie>0.5]<-susie_get_posterior_mean(res)[pip_susie>0.5]
   SSE10<-sum((beta-beta_susie)^2)/p
   MSPE10<-sum((X_scale%*%beta_susie-Y1)^2)/n
   MSPE10out<-sum((outX_scale%*%beta_susie-outY_scale)^2)/(outn)
   out10<-list(runtime=runtime10,calculate_metrics(beta_susie,beta),SSE=SSE10,MSPE=MSPE10,MSPEout=MSPE10out)
   
   
   ###SIS+lasso ####
   start_time <- Sys.time()
   data <- as.data.frame(X_scale)
   data$y <- Y_scale
   # Step 1: Screening - calculate marginal correlation
   cor_vals <- sapply(as.data.frame(X_scale), function(x) abs(cor(x, data$y)))
   ranked_vars <- names(sort(abs(cor_vals), decreasing = TRUE))
   
   # Select top d variables
   d <- ceiling(n / log(n))  # Around 38 variables
   sis_indices <- ranked_vars[1:d]
   
   # Step 2: Prepare data for regsubsets
   mydata_sis <- data.frame(Y_scale, data[, sis_indices , drop = FALSE])
   X_sis<-data[, sis_indices , drop = FALSE]
   colnames(mydata_sis)[1] <- "response"
   
   # Step 3: Apply forward stepwise selection (based on BIC)
   # Null model (intercept only)
   cvfit <- cv.glmnet(as.matrix(X_sis), Y_scale, alpha = 1, standardize = TRUE)
   lasso_coef <- coef(cvfit, s = "lambda.min")[-1]  # remove intercept
   selected_local <- which(lasso_coef != 0)
   selected_global <- sis_indices[selected_local]
   selected_global <- as.integer(gsub("^V", "", selected_global ))
   
   sis_coef <- numeric(p)
   sis_coef[selected_global] <- lasso_coef
  
   # Step 7: Evaluate
   SSE11 <- sum((beta - sis_coef)^2) / p
   MSPE11 <- sum((X_scale %*% sis_coef - Y_scale)^2) / n
   MSPE11out <- sum((outX_scale %*% sis_coef  - outY_scale)^2) / outn
   
   # Step 8: Wrap results
   runtime11 <- difftime(Sys.time(), start_time, units = "secs")
   out11 <- list(runtime = runtime11,
                calculate_metrics(sis_coef , beta),
                SSE = SSE11,
                MSPE = MSPE11,
                MSPEout = MSPE11out) 
   
   ###RandomForest
   ##stepone
   
   start_time <- Sys.time()
   
   mydata=data.frame(Y_scale, X_scale)
   fit_rf <- rfsrc(Y_scale ~ ., data = mydata, ntree = 500,forest = TRUE,importance = TRUE)
   
   # Calculate minimal depth
   varsel_result <- var.select(fit_rf)
   
   # Extract the minimal depth matrix and variable names
   depth_matrix <- varsel_result$md.obj$order  # usually a matrix with depth info
   var_names <- varsel_result$md.obj$names     # names corresponding to rows in 'order'
   
   # Create a data frame with variable names and their minimal depths
   depth_df <- data.frame(
     variable =paste0("X", 1:p),
     depth = depth_matrix[, 1],  # column 1 is minimal depth
     stringsAsFactors = FALSE
   )
   depth_df$variable<-toupper(depth_df$variable)
   # Filter for only the variables in topvars
   topvars <- varsel_result$topvars
   depth_df_topvars <- depth_df[depth_df$variable %in% topvars, ]
   
   # Sort by lowest depth (i.e., most important)
   depth_df_topvars <- depth_df_topvars[order(depth_df_topvars$depth), ]
   
   # Select the top 10% based on depth
   #myd= ceiling(n/sqrt(log(n))) #ceiling(n / log(n))
   myd= ceiling(n / log(n))
   n_top <- min(myd, nrow(depth_df_topvars)) #ceiling(0.10 * nrow(depth_df_topvars))
   top10pct_topvars <- depth_df_topvars$variable[1:n_top]
   
   # Result
   print(top10pct_topvars)
   
   #slc_feature=X[top10pct_topvars]
   selected_rf_vars <- as.numeric(gsub("X", "", top10pct_topvars))
   
   # Step 4: Use selected variables for prediction
   X_rf_train <- X_scale[, selected_rf_vars, drop = FALSE]
   X_rf_test <- outX_scale[, selected_rf_vars, drop = FALSE]
   
   
   colnames(X_rf_test) <- colnames(X_rf_train)
   
   # 7. Fit a linear model on selected variables (can also do RF again if desired)
   ft_rfsrc <- rfsrc(Y_scale ~ ., data = data.frame(Y_scale, X_rf_train), ntree = 500)
   #lm(Y_scale ~ ., data = data.frame(Y_scale = Y_scale, X_rf_train))
   # 8. Predict on training and test data
   train_pred <- predict(ft_rfsrc, newdata = data.frame(X_rf_train))
   test_pred <- predict(ft_rfsrc, newdata = data.frame(X_rf_test))
   
   # 9. End timer
   end_time <- Sys.time()
   runtime12 <- difftime(end_time, start_time, units = "secs")
   beta_rf <- rep(0, p)
   names(beta_rf) <- paste0("X", 1:p)
   beta_rf[selected_rf_vars] <- 1
   
   # 9. Metrics

   SSE12 <- sum((beta - beta_rf)^2) / p
   MSPE12 <- mean((train_pred$predicted - c(Y_scale))^2)
   MSPE12out <- mean((test_pred$predicted - c(outY_scale))^2)
   
   # 11. Store output
   out12 <- list(
     runtime = runtime12,
     calculate_metrics(beta_rf, beta),
     SSE = SSE12,
     MSPE = MSPE12,
     MSPEout = MSPE12out
   )
   
   ## Step two
   mydata.rf=data.frame(X_scale[,selected_rf_vars], Y_scale)
   names(mydata.rf)<-c(selected_rf_vars, "response")
   null.rf=lm(response~1, data=mydata.rf)
   full.rf=lm(response~., data=mydata.rf)
   final.rf <- step(null.rf,  scope = list(lower = null.rf, upper = full.rf), 
                    direction = "forward", trace = TRUE, k=log(n)) #     full.rf, k=log(n), trace=0)# using the BIC criteria
   end_time <- Sys.time()
   runtime13=difftime(end_time, start_time, units="secs")
   summary(final.rf)
   
   # 1. Remove backticks and extract variable numbers from names like "X420"
   sig.genes.rf <- gsub("`", "", names(final.rf$coefficients[-1]))  # remove backticks
   sig.genes.rf <- gsub("X", "", sig.genes.rf)                      # remove "X"
   sig.genes.rf <- as.numeric(sig.genes.rf)                         # convert to numeric
   
   
   #sig.genes.rf=noquote(noquote(names(final.rf$coefficients[-1])))
   #sig.genes.rf=as.numeric(gsub("`", "",sig.genes.rf))
   beta_rfsfs<-numeric(p)
   beta_rfsfs[sig.genes.rf] <- as.numeric(final.rf$coefficients[-1])
   SSE13<-sum((beta-beta_rfsfs)^2)/p
   MSPE13<-sum((X_scale%*%beta_rfsfs-Y_scale)^2)/(n)
   MSPE13out<-sum((outX_scale%*%beta_rfsfs-outY_scale)^2)/(outn)
   out13<-list(runtime=runtime13,calculate_metrics(beta_rfsfs,beta),SSE=SSE13,MSPE=MSPE13,MSPEout=MSPE13out)
   
   
   

   
   
   
  df<-cbind(unlist(out1),unlist(out2),unlist(out3),unlist(out4),unlist(out5),unlist(out6),unlist(out7),unlist(out8),unlist(out9),unlist(out10),unlist(out11),unlist(out12),unlist(out13))
  colnames(df)<-c("lasso","EL","alasso","SL","HS","HSP","SN","BL","RHS","SUSIE","SISL","RF","RFSFS")
  write.csv(df,paste("All_","_SNR_",SNR,"_Sim_",Sim,"_p_",p,"_spars_",spars,".csv",sep=""))
  beta_estimates<-as.data.frame( cbind(c(lasso_model_optimal$beta),c(eln_model_optimal$beta),c(alasso_fit$beta),c(ssL1$beta[,1]),
                                       c(model1$BetaHat),c(hsp$mu.beta),c(ssbeta),c(colMeans(model_bl$beta)),c(fit$coefficients),c(beta_susie),
                                       c(sis_coef),c(beta_rf),c(beta_rfsfs))) 
  colnames(beta_estimates)<-c("lasso","EL","alasso","SL","HS","HSP","SN","BL","RHS","SUSIE","SISL","RF","RFSFS")
  #write.csv(beta_estimates,paste("beta_estimates_",Sim,"_SNR_",SNR,"_p_",p,"_spars_",spars,".csv",sep=""))
  
  }
#}

#ans = fn_iteration(1)
#ans



Sim=100
registerDoParallel(cores=8)
results<-foreach(i=1:Sim)%dopar%{
set.seed(i)
Allmodels(p=1000,spars=0.1,"C",SNR=5,Sim=i)
}
