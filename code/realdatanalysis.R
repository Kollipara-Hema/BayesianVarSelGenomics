library(MASS)
library(rstanarm)
library(glmnet)
library(caret)
library(monomvn)
library(caret)
library(BoomSpikeSlab)
library(SSLASSO)
library(UCSCXenaTools)
library(bayesreg)
library(horseshoe)
require(randomForest)
MSEY<-function(Y,X,B){print(sum((Y-X%*%B)^2)/nrow(Y))}
MSEB<-function(beta,B){print(sum((beta-B)^2))}

set.seed(321)
###lasso
set.seed(143)
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
MSEY(Y_scale,X_scale,lasso_model_optimal$beta)
sum((X_scale%*%lasso_model_optimal$beta-Y_scale)^2)/n
sum(lasso_model_optimal$beta!=0)
sum(round(lasso_model_optimal$beta,2)!=0)

###elasitc net
set.seed(143)
start_time <- Sys.time()
eln_model <- glmnet(X_scale, Y_scale, alpha = 0.5)
cv_model1 <- cv.glmnet(X_scale, Y_scale, alpha = 0.5)
lambda_optimal <- cv_model1$lambda.min
eln_model_optimal <- glmnet(X_scale, Y_scale, alpha = 0.5, lambda = lambda_optimal)
end_time <- Sys.time()
runtime2=difftime(end_time, start_time, units="secs")
runtime2
MSEY(Y_scale,X_scale,eln_model_optimal$beta)
sum((X_scale%*%eln_model_optimal$beta-Y_scale)^2)/n
sum(eln_model_optimal$beta!=0)
sum(round(eln_model_optimal$beta,2)!=0)

###blasso

### T=5000(MCMC samples)
###RJ to do feature selection
###M=maximum number of allowed covariates
### beta(initial beta)
###lambda2>0
###case="default"
###
set.seed(143)
start_time <- Sys.time()
model_bl1<-blasso(X_scale,Y_scale,T=5000,case = "default",icept = FALSE)
end_time <- Sys.time()
runtime3=difftime(end_time, start_time, units="secs")
runtime3
MSEY(Y_scale,X_scale,colMeans(model_bl1$beta))
sum((X_scale%*%colMeans(model_bl1$beta)-Y_scale)^2)/n
sum(colMeans(model_bl1$beta)!=0)
sum(round(colMeans(model_bl1$beta),2)!=0)

###SSLASSO
set.seed(143)
start_time <- Sys.time()
ssL1<-SSLASSO(X_scale,Y_scale,variance="unknown",max.iter=5000)
lambda0_start <- which(ssL1$iter < 100)[1]
lambda0 <- seq(lambda0_start, 100, by = 1)
ssL1<- SSLASSO(X_scale, Y_scale, variance = "unknown", lambda1 = 1, lambda0 = lambda0)
end_time <- Sys.time()
runtime4=difftime(end_time, start_time, units="secs") 
runtime4
MSEY(Y_scale,X_scale,ssL1$beta[,ncol(ssL1$beta)])
sum((X_scale%*%ssL1$beta[,ncol(ssL1$beta)]-Y_scale)^2)/n
sum(ssL1$beta[,ncol(ssL1$beta)]!=0)
sum(round(ssL1$beta[,ncol(ssL1$beta)],2)!=0)
sum((X_scale%*%ssL1$beta[,1]-Y_scale)^2)/n
sum(ssL1$beta[,1]!=0)
sum(round(ssL1$beta[,1],2)!=0)
MSEY(Y_scale,X_scale,ssL1$beta[,1])

##SUSIE

# install.packages("remotes")
#remotes::install_github("stephenslab/susieR")
library(susieR)
set.seed(143)
Y_scale<-drop(Y_scale)
start_time <- Sys.time()
res1=susie(X_scale,Y_scale,L=10,intercept=FALSE,max_iter=5000)
end_time <- Sys.time()
runtime5=difftime(end_time, start_time, units="secs") 
runtime5
susue_beta_est<-susie_get_posterior_mean(res1)

#susie_get_cs(res1)#eX_scaletractcrediblesetsfromfit 

#plot(Y_scale,predict(res1)) 
sum((Y_scale-predict(res1))^2)/n
sum(susie_get_posterior_mean(res1)!=0)
sum(round(susie_get_posterior_mean(res1),2)!=0)
sum(res1$pip>0.5)
susie_plot(res1,"PIP")

####horseshoe####
set.seed(143)
start_time <- Sys.time()
model1<-horseshoe::horseshoe(Y_scale,X_scale,method.tau = "halfCauchy",method.sigma = "Jeffreys",burn=1000,nmc=5000)
end_time <- Sys.time()
runtime6=difftime(end_time, start_time, units="secs") 
runtime6
MSEY(Y_scale,X_scale,model1$BetaHat)
sum((X_scale%*%model1$BetaHat-Y_scale)^2)/n
sum(model1$BetaHat!=0)
sum(round(model1$BetaHat,2)!=0)

#####horseshoe+#####
set.seed(143)
start_time <- Sys.time()
hsp<-bayesreg(Y_scale~X_scale,merged_data,model = "gaussian",prior="hs+",n.samples=5000,burnin=1000)

#hsp<-bayesreg(Y_scale~X_scale,merged_data,model = "gaussian",prior="hs+",n.samples=5000,burnin=1000)
end_time <- Sys.time()
runtime7=difftime(end_time, start_time, units="secs") 
runtime7
MSEY(Y_scale,X_scale,model1$BetaHat)
sum((X_scale%*%hsp$mu.beta-Y_scale)^2)/n
sum(hsp$mu.beta!=0)
sum(round(hsp$mu.beta,2)!=0)

#spike and slab
set.seed(143)
start_time <- Sys.time()
prior_ssN1<-SpikeSlabPrior(x=X_scale,y=Y_scale, expected.model.size=100,   # eX_scalepect 100 nonzero predictors          
                           diagonal.shrinkage = 0.5,   # use Zellner's prior
                           prior.df = spars)# shrink to zero)
#)
## default prior
## 
ssN1 <- lm.spike(Y_scale ~ X_scale - 1, niter = 10000,#prior=prior_ssN1,
                 error.distribution = "gaussian",ping=500,seed=123)

end_time <- Sys.time()
runtime8=difftime(end_time, start_time, units="secs") 
runtime8
MSEY(Y_scale,X_scale,ssN1$beta[which.min(ssN1$sse),])
MSEY(Y_scale,X_scale,colMeans(ssN1$beta))
sum((X_scale%*%colMeans(ssN1$beta)-Y_scale)^2)/n
sum(colMeans(ssN1$beta)!=0)
sum(round(colMeans(ssN1$beta),2)!=0)

### regularized horshoe #####
set.seed(321)
options (mc.cores = parallel :: detectCores ())
# set up the prior , use hY_scaleperprior tau~half - CauchY_scale (0, tau0 ^2)
D <- ncol (X_scale)
n <- nrow (X_scale)
p0 <- 10 # prior guess for the number of relevant variables
tau0 <- p0 /(D-p0) / sqrt (n)
ps=p0
prior4<-hs_plus(global_scale = tau0,
                slab_df = 4,slab_scale = 2.5)
start_time <- Sys.time()
fit4 <- stan_glm (Y_scale ~-1+X_scale, family = gaussian (), data = data.frame(I(X_scale),Y_scale), prior = prior4, warmup = 1000, adapt_delta = 0.95, chains = 4,iter = 5000,sparse = TRUE,seed= 10) #seed = 12345
end_time <- Sys.time()
runtime9=difftime(end_time, start_time, units="secs") 
runtime9
sum((X_scale%*%fit4$coefficients-Y_scale)^2)/(n)
sum(fit4$coefficients!=0)
sum(round(fit4$coefficients,2)!=0)
sigma(fit4)

RHS1<-round(fit4$coefficients,2)
RHS1[which(RHS1!=0)]

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
runtime11
sig.genes.l=noquote(noquote(names(final.lasso$coefficients[-1])))
sig.genes.l=as.numeric(gsub("`", "",sig.genes.l))
final_lasso_coef<-numeric(p)
final_lasso_coef[sig.genes.l] <- as.numeric(final.lasso$coefficients[-1])
sum((final_lasso_coef)!=0)
sum(round(final_lasso_coef,2)!=0)
sum((X_scale%*%final_lasso_coef-Y_scale)^2)/(n)


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
mythreshold=min((n-2), ceiling(0.1*p)) # to screen out important variables
selectedvariablerf=impvariables[1:mythreshold]

## Step two
mydata.rf=data.frame(X_scale[,selectedvariablerf], Y_scale)
names(mydata.rf)<-c(selectedvariablerf, "response")
full.rf=lm(response~., data=mydata.rf)
final.rf <- step(full.rf, k=log(n), trace=0)# using the BIC criteria
end_time <- Sys.time()
runtime12=difftime(end_time, start_time, units="secs")
runtime12
summary(final.rf)
sig.genes.rf=noquote(noquote(names(final.rf$coefficients[-1])))
sig.genes.rf=as.numeric(gsub("`", "",sig.genes.rf))
final_rf_coef<-numeric(p)
final_rf_coef[sig.genes.rf] <- as.numeric(final.rf$coefficients[-1])
sum((final_rf_coef)!=0)
sum(round(final_rf_coef,2)!=0)
sum((X_scale%*%final_rf_coef-Y_scale)^2)/(n)



beta_estimates<-data.frame(c(colMeans(model_bl1$beta)),c(ssL1$beta[,ncol(ssL1$beta)]),c(ssL1$beta[,1]),
                           c(susie_get_posterior_mean(res1)),c(model1$BetaHat),c(hsp$mu.beta),c(colMeans(ssN1$beta)),
                           c(fit4$coefficients),lasso_model_optimal$beta[,1],eln_model_optimal$beta[,1],c(final_lasso_coef),c(final_rf_coef))
colnames(beta_estimates)<-c("BL","SSL","SSL1","Susie","hs","hsp","SSN","rhs","lasso","el","STEP","RF")
row.names(beta_estimates)<-colnames(X_scale)
write.csv(beta_estimates,"MET_tumor_female.csv")
MF<-beta_estimates

