data_all <- list.files(path = "files path",  # Identify all CSV files
                       pattern = "*.csv", full.names = TRUE)


read.tcsv = function(file, header=TRUE, sep=",", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  ## empty strings are converted to NA
  out = read.csv(text=x, sep=sep, header=header, na.strings = "", ...)
  return(out)
  
}

SNR_250_0.1s<-matrix(NA,ncol=19)
colnames(SNR_250_0.1s)<-c("Method" , "runtime", "TP",  "FP","TN","FN","Accu","prec","recal","Fscore","MSBE","MSPE","FDR","FNR","SNR","p","sparse", "Sim", "name")
for(i in 1:(length(data_all))){
  print(i)
  d<-read.tcsv(data_all[i])
  head(d)
  colnames(d)<-c("Method" , "runtime" ,"TP","FP" ,"TN", "FN","Accu","prec", "recal" ,"Fscore" ,"MSBE" ,"MSPE" ) 
  d[, -1] <- lapply(d[, -1], function(x) as.numeric(gsub("NA", NA, x)))
  d$FDR<-ifelse(is.na(d$FP)==TRUE,NA,d$FP/(d$FP+d$TP))
  d$FNR<-ifelse(is.na(d$FN)==TRUE,NA,d$FN/(d$FN+d$TP))
  d$SNR<-substring(data_all[i],first = regexpr("_SNR_",data_all[i])+5, last = regexpr("_Sim_", data_all[i])-1)
  d$p<-substring(data_all[i],first = regexpr("_p_",data_all[i])+3, last = regexpr("_spars_", data_all[i])-1)
  d$sparse<-substring(data_all[i],first = regexpr("_spars_",data_all[i])+7, last = regexpr(".csv", data_all[i])-1)
  d$Sim<-substring(data_all[i],first = regexpr("_Sim_",data_all[i])+5, last = regexpr("_p_", data_all[i])-1)
  d$name<-substring(data_all[i],first = regexpr("All/",data_all[i])+9, last = regexpr(".csv", data_all[i])-1)
  SNR_250_0.1s<-rbind(d,SNR_250_0.1s)
  
}
library(tm)
df1<- SNR_250_0.1s[-which( SNR_250_0.1s$Method=="SL"),]
df<-na.omit(df1)

colnames(df)<-c("Model",colnames(df)[-1])
df$Model<-removeNumbers(toupper(df$Model))
df$spars<-paste0("Sparsity=",df$spars)
df$p<-paste0("P=",df$p)
df$SNR<-paste0("SNR=",df$SNR)
df$Model[which(df$Model=="BLASSO")]="BL"
df$Model[which(df$Model=="SLF")]="SL"
df$Model<-as.factor(df$Model)
df$spars<-as.factor(df$spars)
df$p<-as.factor(df$p)
df$SNR<-as.factor(df$SNR)

halfdata<-df
anotherhalf<-df
df<-rbind(halfdata,anotherhalf)
p.labs=c("P=1000","P=250")
names(p.labs)=c("1000","250")
SNR.labs=c("SNR=1","SNR=5")
names(SNR.labs)=c("1","5")
spars.labs=c("Sparsity=3%","Sparsity=10%")
names(spars.labs)=c("0.03","0.1")
df$new_p=factor(df$p,levels=c('P=250','P=1000'))
df$Model<-factor(df$Model,levels=c("LASSO","EL","BL","HS","HSP","RHS","SN","SL","SUSIE"))
#df$sp<-paste0(df$Metric," , ",df$spars)


p1<-ggplot(df, aes(x = Model, y = FDR ,fill=Model)) +
  geom_boxplot()+ # geom_violin(trim = FALSE) +
  #stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red")+
  stat_summary(fun.data = function(y) data.frame(y = max(y) + 0.05, label = round(mean(y),2)), 
               geom = "text", 
               position = position_dodge(width = 0.75), color = "red", size = 2.5)+
  facet_grid(spars+SNR+new_p~.,scales = "free") +theme(axis.text.x=element_text(size=8, angle=45))+  
  theme(strip.text = element_text(face="bold", size=9,lineheight=5.0), axis.text.x=element_text(angle=45,face="bold"))+labs(fill = "Method",y=NULL,x=NULL)+ggtitle("")
#theme_minimal()


p2<-ggplot(df, aes(x = Model, y = FNR ,fill=Model)) +
  geom_boxplot()+ # geom_violin(trim = FALSE) +
  #stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red")+
  stat_summary(fun.data = function(y) data.frame(y = max(y) + 0.05, label = round(mean(y),2)), 
               geom = "text", 
               position = position_dodge(width = 0.75), color = "red", size = 2.5)+
  facet_grid(spars+SNR+new_p~.,scales = "free") +theme(axis.text.x=element_text(size=8, angle=45))+  
  theme(strip.text = element_text(face="bold", size=9,lineheight=5.0), axis.text.x=element_text(angle=45,face="bold"))+labs(fill = "Method",y=NULL,x=NULL)+ggtitle("")
#theme_minimal()

p3<-ggplot(df, aes(x = Model, y = MSBE ,fill=Model)) +
  geom_boxplot()+ # geom_violin(trim = FALSE) +
  #stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red")+
  stat_summary(fun.data = function(y) data.frame(y = max(y) + 0.02, label = round(mean(y),2)), 
               geom = "text", 
               position = position_dodge(width = 0.75), color = "red", size = 2.5)+
  facet_grid(spars+SNR+new_p~.,scales = "free") +theme(axis.text.x=element_text(size=8, angle=45))+  
  theme(strip.text = element_text(face="bold", size=9,lineheight=5.0), axis.text.x=element_text(angle=45,face="bold"))+labs(fill = "Method",y=NULL,x=NULL)+ggtitle("")
#theme_minimal()facet_grid(SNR+new_p~spars,scales = "free")

p4<-ggplot(df, aes(x = Model, y = MSPE ,fill=Model)) +
  geom_boxplot()+ # geom_violin(trim = FALSE) +
  #stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red")+
  stat_summary(fun.data = function(y) data.frame(y = max(y) + 20, label = round(mean(y),2)), 
               geom = "text", 
               position = position_dodge(width = 0.75), color = "red", size = 2.5)+
  facet_grid(spars+SNR+new_p~.,scales = "free") +theme(axis.text.x=element_text(size=8, angle=45))+  
  theme(strip.text = element_text(face="bold", size=9,lineheight=5.0), axis.text.x=element_text(angle=45,face="bold"))+labs(fill = "Method",y=NULL,x=NULL)+ggtitle("")
#theme_minimal()

p5<-ggplot(df, aes(x = Model, y = Fscore ,fill=Model)) +
  geom_boxplot()+ # geom_violin(trim = FALSE) +
  #stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red")+
  #stat_summary(fun.data = function(y) data.frame(y = max(y) + 0.05, label = length(y)), 
  #            geom = "text", 
  #           position = position_dodge(width = 0.75))+
  stat_summary(fun.data = function(y) data.frame(y = max(y) + 0.05, label = round(mean(y),2)), 
               geom = "text", 
               position = position_dodge(width = 0.75), color = "red", size = 2.5)+
  facet_grid(spars+SNR+new_p~.,scales = "free") +theme(axis.text.x=element_text(size=8, angle=45))+  
  theme(strip.text = element_text(face="bold", size=9,lineheight=5.0), axis.text.x=element_text(angle=45,face="bold"))+labs(fill = "Method",y=NULL,x=NULL)+ggtitle("")
#theme_minimal()

p1$layers <- p1$layers[!sapply(p1$layers, function(layer) inherits(layer$geom, "GeomText"))]  
p2$layers <- p2$layers[!sapply(p2$layers, function(layer) inherits(layer$geom, "GeomText"))]  
p3$layers <- p3$layers[!sapply(p3$layers, function(layer) inherits(layer$geom, "GeomText"))]  
p4$layers <- p4$layers[!sapply(p4$layers, function(layer) inherits(layer$geom, "GeomText"))] 
p5$layers <- p5$layers[!sapply(p5$layers, function(layer) inherits(layer$geom, "GeomText"))]

ggpubr::ggarrange(p1+theme(strip.text = element_blank()),p2,ncol=2,
                  nrow=1,common.legend=TRUE,legend="right",widths= c(1,1), labels = c("FDR" ,"FNR"))

ggpubr::ggarrange(p3+theme(strip.text = element_blank()),p4,ncol=2,
                  nrow=1,common.legend=TRUE,legend="right",widths= c(1,1), labels = c("MSB","MSP"))
