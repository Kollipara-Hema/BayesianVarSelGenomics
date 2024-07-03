### libraries required
require("MASS")
require("rstanarm")
require("glmnet")
require("caret")
require("monomvn")
require("caret")
require("BoomSpikeSlab")
require("SSLASSO")
require("bayesreg")
require("horseshoe")
require("susieR")

###used functions
MSEY<-function(Y,X,B){print(sum((Y-X%*%B)^2)/nrow(Y))}
MSEB<-function(beta,B){print(sum((beta-B)^2))}
confusion_measure<-function(beta,B){
  beta<-relevel(as.factor(ifelse(beta == 0, 0, 1)),ref="1")
  B_<-ifelse(round(B,2) == 0, 0, 1)
  if(sum(B_)==0|sum(B_)==length(B_)){
    TP<-NA
    FP<-NA
    TN<-NA
    FN<-NA
    Accu<-NA
    prec<-NA
    recal<-NA
    Fscore<-NA
  }else{
    B<-relevel(as.factor(ifelse(round(B,2) == 0, 0, 1)),ref="1")
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
#set.seed(123)

###Simulation starting

fn_iteration<-function(iterr){
  set.seed(iterr)
  Allmodels(200,1000,0.1,"B",5,Sim=iterr) 
}

Allmodels<-function(n,p,spars,response,SNR,Sim){
  ###data generation
  #Sim=iterr
  #n=200
  #p=1000
  #spars=0.1
  if(response=="B"){
    X <- matrix(rbinom(n * p,1,prob=0.7)*(-2)+1, nrow = n)}else{
      X <- matrix(rnorm(n * p), nrow = n)}
  
  ps<-ceiling(spars*p) #significant beta
  beta <- sample(c(runif(ps,-3,3), rep(0, p - ps))) 
  
  #SNR=5
  Y <- X %*% beta + rnorm(n,sd=sd(X %*% beta)/sqrt(SNR))
  
  X_scale=X
  Y_scale=Y-mean(Y)
  
  ## LASSO
  p=ncol(X_scale)
  n=nrow(X_scale)
  start_time <- Sys.time()
  lasso_model <- glmnet(X_scale, Y_scale, alpha = 1)
  cv_model <- cv.glmnet(X_scale, Y_scale, alpha = 1)
  lambda_optimal <- cv_model$lambda.min
  lasso_model_optimal <- glmnet(X_scale, Y_scale, alpha = 1, lambda = lambda_optimal)
  end_time <- Sys.time()
  runtime1=difftime(end_time, start_time, units="secs")
  runtime1
  SSE1<-sum((beta-lasso_model_optimal$beta)^2)/p
  MSPE1<-sum((X_scale%*%lasso_model_optimal$beta-Y_scale)^2)/(n)
  out1<-list(runtime=runtime1,confusion_measure(beta,lasso_model_optimal$beta),SSE=SSE1,MSPE=MSPE1)
  
  ## ELASTICNET  (EL)
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
  out2<-list(runtime=runtime2,confusion_measure(beta,eln_model_optimal$beta),SSE=SSE2,MSPE=MSPE2)
  
  ## SSLASSO (SL)
  start_time <- Sys.time()
  ssL1<-SSLASSO(X_scale,Y_scale,variance="unknown",max.iter=5000)
  lambda0_start <- which(ssL1$iter < 100)[1]
  lambda0 <- seq(lambda0_start, 100, by = 1)
  ssL1<- SSLASSO(X_scale, Y_scale, variance = "unknown", lambda1 = 1, lambda0 = lambda0)
  end_time <- Sys.time()
  runtime3=difftime(end_time, start_time, units="secs") 
  runtime3
  #SSE3<-sum((beta-ssL1$beta[,ncol(ssL1$beta)])^2)/p
  #MSPE3<-sum((X_scale%*%ssL1$beta[,ncol(ssL1$beta)]-Y_scale)^2)/n
  #out3<-list(runtime=runtime3,confusion_measure(beta,ssL1$beta[,ncol(ssL1$beta)]),SSE=SSE3,MSPE=MSPE3)
  
  SSE4<-sum((beta-ssL1$beta[,1])^2)/p
  MSPE4<-sum((X_scale%*%ssL1$beta[,1]-Y_scale)^2)/n
  out4<-list(runtime=runtime3,confusion_measure(beta,ssL1$beta[,1]),SSE=SSE4,MSPE=MSPE4)
  
  ## Horseshoe (HS)
  start_time <- Sys.time()
  model1<-horseshoe::horseshoe(Y_scale,X_scale,method.tau = "halfCauchy",method.sigma = "Jeffreys",burn=1000,nmc=5000)
  end_time <- Sys.time()
  runtime5=difftime(end_time, start_time, units="secs") 
  runtime5
  SSE5<-sum((beta-model1$BetaHat)^2)/p
  MSPE5<-sum((X_scale%*%model1$BetaHat-Y_scale)^2)/n
  out5<-list(runtime=runtime5,confusion_measure(beta,model1$BetaHat),SSE=SSE5,MSPE=MSPE5)
  
  ## Horseshoe+ (HSP)
  start_time <- Sys.time()
  hsp<-bayesreg(Y_scale~X_scale,as.data.frame(X_scale,Y_scale),model = "gaussian",prior="hs+",n.samples=5000,burnin=1000)
  end_time <- Sys.time()
  runtime6=difftime(end_time, start_time, units="secs") 
  runtime6
  SSE6<-sum((beta-hsp$mu.beta)^2)/p
  MSPE6<-sum((X_scale%*%hsp$mu.beta-Y_scale)^2)/n
  out6<-list(runtime=runtime6,confusion_measure(beta,hsp$mu.beta),SSE=SSE6,MSPE=MSPE6)
  
  ## SN
  start_time <- Sys.time()
  prior_ssN1<-SpikeSlabPrior(x=X_scale,y=Y_scale, expected.model.size=ps,          
                             diagonal.shrinkage = 0.5,   # use Zellner's prior
                             prior.df = spars)#shrink to zero
  ssN1 <- lm.spike(Y_scale ~ X_scale - 1, niter = 10000,prior=prior_ssN1,
                   error.distribution = "gaussian",ping=500,seed=iterr)
  end_time <- Sys.time()
  runtime7=difftime(end_time, start_time, units="secs") 
  runtime7
  SSE7<-sum((beta-colMeans(ssN1$beta))^2)/p
  MSPE7<-sum((X_scale%*%colMeans(ssN1$beta)-Y_scale)^2)/n
  out7<-list(runtime=runtime7,confusion_measure(beta,colMeans(ssN1$beta)),SSE=SSE7,MSPE=MSPE7)
  
  ## BLASSO (BL)
  start_time <- Sys.time()
  model_bl<-blasso(X_scale,Y_scale,T=5000,case = "default",icept = FALSE)
  end_time <- Sys.time()
  runtime8=difftime(end_time, start_time, units="secs")
  runtime8
  SSE8<-sum((beta-colMeans(model_bl$beta))^2)/p
  MSPE8<-sum((X_scale%*%colMeans(model_bl$beta)-Y_scale)^2)/n
  out8<-list(runtime=runtime8,confusion_measure(beta,colMeans(model_bl$beta)),SSE=SSE8,MSPE=MSPE8)
  
  ## Regularised Horseshoe (RHS)
  D <- ncol (X_scale)
  n <- nrow (X_scale)
  p0 <- ps # prior guess for the number of relevant variables
  tau0 <- p0 /(D-p0) / sqrt (n) # rstanarm will scale this by sigma automatically
  prior4<-hs(global_scale = tau0,slab_scale=2.5, slab_df=4)
  options (mc.cores = parallel :: detectCores ())
  start_time <- Sys.time()
  fit <- stan_glm (Y_scale ~-1+X_scale, family = gaussian (), data = data.frame (I(X_scale),Y_scale),
                   prior = prior4,
                   chains = 4,iter = 3000,sparse = TRUE,seed=iterr)
  end_time <- Sys.time()
  runtime9=difftime(end_time, start_time, units="secs") 
  SSE9<-sum((beta-fit$coefficients)^2)/p
  MSPE9<-sum((X_scale%*%fit$coefficients-Y_scale)^2)/n
  out9<-list(runtime=runtime9,confusion_measure(beta,fit$coefficients),SSE=SSE9,MSPE=MSPE9)
  
  ## SUSIE
  Y1<-drop(Y_scale)
  start_time <- Sys.time()
  res=susie(X_scale,Y1,L=ps,intercept=FALSE,max_iter=5000)
  end_time <- Sys.time()
  runtime10=difftime(end_time, start_time, units="secs") 
  runtime10
  SSE10<-sum((beta-susie_get_posterior_mean(res))^2)/p
  MSPE10<-sum((X_scale%*%susie_get_posterior_mean(res)-Y1)^2)/n
  out10<-list(runtime=runtime10,confusion_measure(beta,susie_get_posterior_mean(res)),SSE=SSE10,MSPE=MSPE10)
  
  
  ###combining results for one simulation
  
  df<-cbind(unlist(out1),unlist(out2),unlist(out4),unlist(out5),unlist(out6),unlist(out7),
            unlist(out8),unlist(out9),unlist(out10))
  colnames(df)<-c("lasso","EL","SL","HS","HSP","SN","blasso","RHS","SUSIE")
  write.csv(df,paste("All_","_SNR_",SNR,"_Sim_",Sim,"_p_",p,"_spars_",spars,".csv",sep=""))
  beta_estimates<-as.data.frame(cbind(c(ssL1$beta[,1]),c(model1$BetaHat),c(hsp$mu.beta),c(colMeans(ssN1$beta)),c(colMeans(model_bl$beta)),c(fit$coefficients),c(susie_get_posterior_mean(res)),c(lasso_model_optimal$beta[,1]),c(eln_model_optimal$beta[,1])))
  colnames(beta_estimates)<-c("SL","HS","HSP","SN","blasso","RHS","SUSIE","lasso","EL")
  write.csv(beta_estimates,paste("beta_estimates_",Sim,"_SNR_",SNR,"_p_",p,"_spars_",spars,".csv",sep=""))
}

sim_no=1 #no of simulations

ans = fn_iteration(sim_no)
ans