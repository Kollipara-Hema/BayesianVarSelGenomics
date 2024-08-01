if(TRUE){args=(commandArgs(TRUE))
for(i in 1:length(args)){eval(parse(text=args[[i]]))}}else{
  jobID=1
}

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

fn_iteration<-function(iterr){
  set.seed(iterr)
  
  
  ###data generation
  Sim=iterr
  outn=50
  n=200
  p=250;spars=0.03;SNR=1;response="B"
  ps<-ceiling(spars*p) #significant beta
  beta <- c(rep(1,(ps/2)),rep(-1,(ps/2)), rep(0, p - ps))
  
  
  if(response=="B"){
    X <- matrix(rbinom(n * p,1,prob=0.7)*(-2)+1, nrow = n)}else{
      X <- matrix(rnorm(n * p), nrow = n)
    }
  #X <- matrix(rbinom(n * p,1,prob=0.7)*(-2)+1, nrow = n)
  #X <- matrix(rnorm(n * p), nrow = n)
  
  
  #beta <- sample(c(runif(ps,-3,3), rep(0, p - ps))) 
  
  
  Y <- X %*% beta + rnorm(n,sd=2.47)
  
  X_scale=X
  Y_scale=Y-mean(Y)
  #out of sample generation
  outX<- matrix(rnorm(outn * p), nrow = outn)
  outY <- outX %*% beta + rnorm(outn,sd=2.47)
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
  runtime1
  SSE1<-sum((beta-lasso_model_optimal$beta)^2)/p
  MSPE1<-sum((X_scale%*%lasso_model_optimal$beta-Y_scale)^2)/(n)
  MSPE1out<-sum((outX_scale%*%lasso_model_optimal$beta-outY_scale)^2)/(outn)
  
  out1<-list(runtime=runtime1,confusion_measure(lasso_model_optimal$beta,beta),SSE=SSE1,MSPE=MSPE1,MSPEout=MSPE1out)
  
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
  out2<-list(runtime=runtime2,confusion_measure(eln_model_optimal$beta,beta),SSE=SSE2,MSPE=MSPE2,MSPEout=MSPE2out)
  
  ##sslasso
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
  MSPE4out<-sum((outX_scale%*%ssL1$beta[,1]-outY_scale)^2)/(outn)
  out4<-list(runtime=runtime3,confusion_measure(ssL1$beta[,1],beta),SSE=SSE4,MSPE=MSPE4,MSPEout=MSPE4out)
  
  ##horseshoe####
  start_time <- Sys.time()
  model1<-horseshoe::horseshoe(Y_scale,X_scale,method.tau = "halfCauchy",method.sigma = "Jeffreys",burn=1000,nmc=5000)
  end_time <- Sys.time()
  runtime5=difftime(end_time, start_time, units="secs") 
  runtime5
  SSE5<-sum((beta-model1$BetaHat)^2)/p
  MSPE5<-sum((X_scale%*%model1$BetaHat-Y_scale)^2)/n
  MSPE5out<-sum((outX_scale%*%model1$BetaHat-outY_scale)^2)/(outn)
  out5<-list(runtime=runtime5,confusion_measure(model1$BetaHat,beta),SSE=SSE5,MSPE=MSPE5,MSPEout=MSPE5out)
  
  ##hsp
  start_time <- Sys.time()
  hsp<-bayesreg(Y_scale~X_scale,as.data.frame(X_scale,Y_scale),model = "gaussian",prior="hs+",n.samples=5000,burnin=1000)
  end_time <- Sys.time()
  runtime6=difftime(end_time, start_time, units="secs") 
  runtime6
  SSE6<-sum((beta-hsp$mu.beta)^2)/p
  MSPE6<-sum((X_scale%*%hsp$mu.beta-Y_scale)^2)/n
  MSPE6out<-sum((outX_scale%*%hsp$mu.beta-outY_scale)^2)/(outn)
  out6<-list(runtime=runtime6,confusion_measure(hsp$mu.beta,beta),SSE=SSE6,MSPE=MSPE6,MSPEout=MSPE6out)
  
  ##SN
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
  MSPE7out<-sum((outX_scale%*%colMeans(ssN1$beta)-outY_scale)^2)/(outn)
  out7<-list(runtime=runtime7,confusion_measure(colMeans(ssN1$beta),beta),SSE=SSE7,MSPE=MSPE7,MSPEout=MSPE7out)
  
  ##blasso
  start_time <- Sys.time()
  model_bl<-blasso(X_scale,Y_scale,T=5000,case = "default",icept = FALSE)
  end_time <- Sys.time()
  runtime8=difftime(end_time, start_time, units="secs")
  runtime8
  SSE8<-sum((beta-colMeans(model_bl$beta))^2)/p
  MSPE8<-sum((X_scale%*%colMeans(model_bl$beta)-Y_scale)^2)/n
  MSPE8out<-sum((outX_scale%*%colMeans(model_bl$beta)-outY_scale)^2)/(outn)
  out8<-list(runtime=runtime8,confusion_measure(colMeans(model_bl$beta),beta),SSE=SSE8,MSPE=MSPE8,MSPEout=MSPE8out)
  
  ##rhs
  D <- ncol (X_scale)
  n <- nrow (X_scale)
  p0 <- ps # prior guess for the number of relevant variables
  tau0 <- p0 /(D-p0) / sqrt (n) # rstanarm will scale this by sigma automatically
  prior4<-hs(global_scale = tau0,slab_scale=2.5, slab_df=4)
  options (mc.cores = parallel :: detectCores ())
  start_time <- Sys.time()
  fit <- stan_glm (Y_scale ~-1+X_scale, family = gaussian (), data = data.frame (I(X_scale),Y_scale),
                   prior = prior4,
                   chains = 4,iter = 5000,sparse = TRUE,seed=iterr)
  end_time <- Sys.time()
  runtime9=difftime(end_time, start_time, units="secs") 
  SSE9<-sum((beta-fit$coefficients)^2)/p
  MSPE9<-sum((X_scale%*%fit$coefficients-Y_scale)^2)/n
  MSPE9out<-sum((outX_scale%*%fit$coefficients-outY_scale)^2)/(outn)
  out9<-list(runtime=runtime9,confusion_measure(fit$coefficients,beta),SSE=SSE9,MSPE=MSPE9,MSPEout=MSPE9out)
  
  ###susie
  Y1<-drop(Y_scale)
  start_time <- Sys.time()
  res=susie(X_scale,Y1,L=ps,intercept=FALSE,max_iter=5000)
  end_time <- Sys.time()
  runtime10=difftime(end_time, start_time, units="secs") 
  runtime10
  SSE10<-sum((beta-susie_get_posterior_mean(res))^2)/p
  MSPE10<-sum((X_scale%*%susie_get_posterior_mean(res)-Y1)^2)/n
  MSPE10out<-sum((outX_scale%*%susie_get_posterior_mean(res)-outY_scale)^2)/(outn)
  out10<-list(runtime=runtime10,confusion_measure(susie_get_posterior_mean(res),beta),SSE=SSE10,MSPE=MSPE10,MSPEout=MSPE10out)
  
  
  ### step 
  est.lasso=as.numeric(coef(lasso_model_optimal))[-1]
  mydata.lasso=data.frame(X_scale[,which(est.lasso!=0)], Y_scale)
  names(mydata.lasso)<-c(which(est.lasso!=0), "response")
  start_time <- Sys.time()
  #null.lasso=lm(response~1, data=mydata.lasso)
  full.lasso=lm(response~., data=mydata.lasso)
  #final.lasso <- step(null.lasso, score=list(lower=null.lasso, upper=full.lasso),direction="forward", k=log(n), trace=0)
  final.lasso <- step(full.lasso, k=log(n), trace=0)# using the BIC criteria
  end_time <- Sys.time()
  runtime11=difftime(end_time, start_time, units="secs")
  sig.genes.l=noquote(noquote(names(final.lasso$coefficients[-1])))
  sig.genes.l=as.numeric(gsub("`", "",sig.genes.l))
  final_lasso_coef<-numeric(p)
  final_lasso_coef[sig.genes.l] <- as.numeric(final.lasso$coefficients[-1])
  SSE11<-sum((beta-final_lasso_coef)^2)/p
  MSPE11<-sum((X_scale%*%final_lasso_coef-Y_scale)^2)/(n)
  MSPE11out<-sum((outX_scale%*%final_lasso_coef-outY_scale)^2)/(outn)
  out11<-list(runtime=runtime11,confusion_measure(final_lasso_coef,beta),SSE=SSE11,MSPE=MSPE11,,MSPEout=MSPE11out)
  
  
  #### random forest
  ## Two-step procedure: randomforest to screen a set of features
  ## and then apply the stepwise variable selection technique 
  ## using the BIC criteria
  ## Step one
  start_time <- Sys.time()
  mydata=data.frame(Y_scale, X_scale)
  outrf=randomForest(Y_scale~., ntree=500, data=mydata)
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
  SSE12<-sum((beta-final_rf_coef)^2)/p
  MSPE12<-sum((X_scale%*%final_rf_coef-Y_scale)^2)/(n)
  MSPE12out<-sum((outX_scale%*%final_rf_coef-outY_scale)^2)/(outn)
  out12<-list(runtime=runtime12,confusion_measure(final_rf_coef,beta),SSE=SSE12,MSPE=MSPE12,,MSPEout=MSPE12out)
  
  ####selection by partial correlation
  start_time <- Sys.time()
  pcorr<-pcor(mydata)
  partial_corrs <- c(pcorr$estimate[1,2:(p+1)])
  impvariables.pc=sort(abs(partial_corrs), index.return=T, decreasing=TRUE)$ix
  mythreshold=min((n-2), 0.1*p) # to screen out important variables 
  selectedvariablerf.pc=impvariables.pc[1:mythreshold]
  
  ## Step two
  mydata.rf.pc=data.frame(X_scale[,selectedvariablerf.pc], Y_scale)
  names(mydata.rf.pc)<-c(selectedvariablerf.pc, "response")
  full.rf.pc=lm(response~., data=mydata.rf.pc)
  final.rf.pc <- step(full.rf.pc, k=log(n), trace=0)# using the BIC criteria
  end_time <- Sys.time()
  runtime13=difftime(end_time, start_time, units="secs")
  summary(final.rf.pc)
  sig.genes.rf.pc=noquote(noquote(names(final.rf.pc$coefficients[-1])))
  sig.genes.rf.pc=as.numeric(gsub("`", "",sig.genes.rf.pc))
  final_rf_coef.pc<-numeric(p)
  final_rf_coef.pc[sig.genes.rf.pc] <- as.numeric(final.rf.pc$coefficients[-1])
  SSE13<-sum((beta-final_rf_coef.pc)^2)/p
  MSPE13<-sum((X_scale%*%final_rf_coef.pc-Y_scale)^2)/(n)
  MSPE13out<-sum((outX_scale%*%final_rf_coef.pc-outY_scale)^2)/(outn)
  out13<-list(runtime=runtime13,confusion_measure(final_rf_coef.pc,beta),SSE=SSE13,MSPE=MSPE13,MSPEout=MSPE13out)
  ###combining results for one simulation
  
  df<-cbind(unlist(out1),unlist(out2),unlist(out4),unlist(out5),unlist(out6),unlist(out7),
            unlist(out8),unlist(out9),unlist(out10),unlist(out11),unlist(out12),unlist(out13))
  colnames(df)<-c("lasso","EL","SL","HS","HSP","SN","blasso","RHS","SUSIE","STEP","RF","RFPC")
  write.csv(df,paste("All_","_SNR_",SNR,"_Sim_",Sim,"_p_",p,"_spars_",spars,".csv",sep=""))
  beta_estimates<-as.data.frame(cbind(c(ssL1$beta[,1]),c(model1$BetaHat),c(hsp$mu.beta),c(colMeans(ssN1$beta)),c(colMeans(model_bl$beta)),c(fit$coefficients),c(susie_get_posterior_mean(res)),
                                      c(final_lasso_coef),c(final_rf_coef),c(final_rf_coef.pc),c(lasso_model_optimal$beta[,1]),c(eln_model_optimal$beta[,1])))
  colnames(beta_estimates)<-c("SL","HS","HSP","SN","blasso","RHS","SUSIE","STEP","RF","RFPC","lasso","EL")
  write.csv(beta_estimates,paste("beta_estimates_",Sim,"_SNR_",SNR,"_p_",p,"_spars_",spars,".csv",sep=""))
}

ans = fn_iteration(jobID)
ans