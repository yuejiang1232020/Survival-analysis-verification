
library("glmnet")
library("survival")
library("survminer")
asType=c("AA","AD","AP","AT","ES","RI","ME","")##¡°¡±is ALL type  

setwd("C:\\Users\\17475\\Desktop\\Survival analysis verification")##Set to your own working path, double slash required for win10 systems
rt=read.table("uniSigExp.txt",header=T,sep="\t",row.names=1,check.names=F)##       
rt$futime[rt$futime<=0]=1
rt$futime=rt$futime/365
genes=colnames(rt)

for (i in 1:length(asType)) {
  gene=grep(paste0("\\|",asType[i]),genes,value=T)
  geneLength=ifelse(length(gene)>20,20,length(gene))
  rt1=rt[,c("futime","fustat",gene[1:geneLength])]
  
  x=as.matrix(rt1[,c(3:ncol(rt1))])
  y=data.matrix(Surv(rt1$futime,rt1$fustat))
  
  ###lasso model
  fit <- glmnet(x, y, family = "cox", maxit = 1000)
  cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
  coef <- coef(fit, s = cvfit$lambda.min)
  index <- which(coef != 0)
  actCoef <- coef[index]
  lassoGene=row.names(coef)[index]
  lassoGene=c("futime","fustat",lassoGene)
  lassoSigExp=rt[,lassoGene]
  lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
  write.table(lassoSigExp,file=paste0("lassoSigExp",".",asType[i],".txt"),sep="\t",row.names=F,quote=F)
  
  
  rt2=read.table(paste0("lassoSigExp",".",asType[i],".txt"),header=T,sep="\t",check.names=F,row.names=1)
  #COX model construction
  multiCox=coxph(Surv(futime, fustat) ~ ., data = rt2)
  multiCox=step(multiCox,direction = "both")
  multiCoxSum=summary(multiCox)
  
  #Output model parameters
  outTab=data.frame()
  outTab=cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  outTab=cbind(id=row.names(outTab),outTab)
  outTab=gsub("`","",outTab)
  write.table(outTab,file=paste0("multiCox",".",asType[i],".xls"),sep="\t",row.names=F,quote=F)
  
  #Output patient risk values
  riskScore=predict(multiCox,type="risk",newdata=rt2)
  coxGene=rownames(multiCoxSum$coefficients)
  coxGene=gsub("`","",coxGene)
  outCol=c("futime","fustat",coxGene)
  risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
  write.table(cbind(id=rownames(cbind(rt2[,outCol],riskScore,risk)),cbind(rt2[,outCol],riskScore,risk)),
              file=paste0("risk",".",asType[i],".txt"),
              sep="\t",
              quote=F,
              row.names=F)

  dt=read.table(paste0("risk",".",asType[i],".txt"),header=T,sep="\t")
  diff=survdiff(Surv(futime, fustat) ~risk,data = dt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  
  fit <- survfit(Surv(futime, fustat) ~ risk, data = dt)
  #Plotting survival curves
  surPlot=ggsurvplot(fit, 
             data=dt,
             conf.int=TRUE,
             pval=paste0("p=",pValue),
             pval.size=4,
             risk.table=TRUE,
             legend.labs=c("High risk", "Low risk"),
             legend.title="Risk",
             xlab="Time(years)",
             break.time.by = 1,
             risk.table.title="",
             palette=c("red", "blue"),
             risk.table.height=.25)
  pdf(file=paste0("survival",".",asType[i],".pdf"),onefile = FALSE,width = 6.5,height =5.5)
  print(surPlot)
  dev.off()

}
