rm(list=ls())

library(glmnet)
library(dplyr)
require(doMC)
#registerDoMC(cores=16)

registerDoMC(cores=8)


set.seed(12)
x.test<-sample_n(lncProtexprs,435)
y.test<-factor(x.test$Subgroup)
x.train<-lncProtexprs[-as.numeric(row.names(x.test)),]
y.train<-factor(x.train$Subgroup)
###############################################
### Test elastic regression
### May perform better when n > p and correlated covariates
###############################################
# generate a set of different levels of alpha to test in elastic net
# set up a list for the models of variable alpha

fit.Lasso <- glmnet(x.train, y.train, family="multinomial", alpha=1)
fit.Ridge <- glmnet(x.train, y.train, family="multinomial", alpha=0)
fit.Elnet.5 <- glmnet(x.train, y.train, family="multinomial", alpha=.5)
for (i in 0:10) {
  assign(paste("elasticFit", i, sep=""), cv.glmnet(x.train, y.train, type.measure="mse", 
                                                   alpha=i/10,family="multinomial"))
}




lncProtexprs<-read.table(file ="lncProtexprs.042218.txt",header=T)
y<-factor(lncProtexprs$`exprs$Subgroup`)

#set seed for reproducibility, split the test and training sets to 1:2 ratio.
set.seed(12)
test<-sample_n(lncProtexprs, 435)
#test<-sample_n(lncProtexprs, 330)
dim(test)
train<-lncProtexprs[-as.numeric(row.names(test)),]
dim(train)
train[1:4,1:4]
ytest<-factor(test$`exprs.Subgroup`)
head(ytest)
ytrain<-factor(train$`exprs.Subgroup`)
head(ytrain)

# Use cross validation to for parameter tuning to build various models. 
# Set the foldid parameter so that the aplhas can be compared.

foldid=sample(1:10,size=length(ytrain),replace=TRUE)
 
#cvlasso<-cv.glmnet(data.matrix(train[,3:493]), ytrain,foldid=foldid, family="multinomial", alpha=1,type.multinomial = "grouped",type.measure="class",parallel = TRUE)

#Inspect the coefficients for the model when alpha=1, lasso.
#coef(cvlasso,s = "lambda.min")$`1`

#cvridge<-cv.glmnet(data.matrix(train[,3:493]), ytrain,foldid=foldid, family="multinomial", alpha=0,type.multinomial = "grouped",type.measure="class",parallel = TRUE)

#Inspect the coefficients for the model when alpha=0, ridge.
#coef(cvridge,s="lambda.min")$`1`

# With alpha=0.5, if predictors are correlated within groups, the elastic net selects the groups as in or out together. 
#cvelastic0.5<-cv.glmnet(data.matrix(train[,3:493]), ytrain,foldid=foldid, family="multinomial", alpha=0.5,type.multinomial = "grouped",type.measure="class",parallel = TRUE)

#Inspect the coefficients for the model when alpha=0.5
#myCoefs<-coef(cvelastic0.5,s = "lambda.min")
#myResults <- data.frame(
#  features = myCoefs@Dimnames[[1]][which(myCoefs != 0 )], #intercept included
#  coefs    = myCoefs[which(myCoefs != 0 ) ]  #intercept included
#)
#write.table(myResults,file="cvelasticnet.coefs.txt",quote=F,sep="\t",row.names=F,col.names=T)

# With alpha = 0.8, elastic net "performs much like the lasso, but removes degeneracies and wild behaviour caused by extreme correlations in the predictors."  

#cvelastic0.8<-cv.glmnet(data.matrix(train[,3:493]), ytrain,foldid=foldid,family="multinomial", alpha=0.8,type.multinomial = "grouped",type.measure="class",parallel = TRUE)

#Inspect the non-zero coefficients for the model when alpha=0.8
#coef(cvelastic0.8,s = "lambda.min")$`1`

#pdf("/users/cnodzak/BoG_2018/glmnet.models.plots.042318.pdf")
#pdf("/users/cnodzak/BoG_2018/glmnet.models.plots.042318.n330.pdf")
#par(mfrow=c(3,2))
#plot(cvridge);plot(cvelastic0.5);plot(cvelastic0.8);plot(cvlasso)
#plot(log(cvlasso$lambda),cvlasso$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cvlasso$name)
#points(log(cvelastic0.5$lambda),cvelastic0.5$cvm,pch=19,col="green")
#points(log(cvelastic0.8$lambda),cvelastic0.8$cvm,pch=19,col="pink")
#points(log(cvridge$lambda),cvridge$cvm,pch=19,col="blue")
#legend("topleft",legend=c("alpha = 1","alpha = 0.8","alpha = 0.5","alpha = 0"),pch=19,col=c("red","pin#k","green","blue"))
#dev.off()

#GenesCVe0.5<-names(coef(cvelastic0.5,s="lambda.min")$`1`[coef(cvelastic0.5,s = "lambda.min")$`1`[,1] != 0,])
#GenesCVe0.8<-names(coef(cvelastic0.8,s="lambda.min")$`1`[coef(cvelastic0.8,s = "lambda.min")$`1`[,1] != 0,])
#GenesCVlasso<-names(coef(cvlasso,s="lambda.min")$`1`[coef(cvlasso,s = "lambda.min")$`1`[,1] != 0,])
#GenesCVridge<-rownames(coef(cvridge,s="lambda.min")$`1`[coef(cvridge,s = "lambda.min")$`1`[,1] != 0,])
#head(GenesCVridge)
#n<-max(length(GenesCVe0.5),length(GenesCVe0.8),length(GenesCVlasso),length(GenesCVridge))
#n
#length(GenesCVe0.5)
#length(GenesCVe0.5)<-n
#head(GenesCVe0.5)
#length(GenesCVe0.8)<-n
#head(GenesCVe0.8)
#length(GenesCVlasso)<-n
#head(GenesCVlasso)
#length(GenesCVridge)<-n
#head(GenesCVridge)
#genes<-cbind(GenesCVe0.5,GenesCVe0.8,GenesCVlasso,GenesCVridge)
#colnames(genes)<-c("elastic0.5","elastic0.8","lasso","ridge")

#write.table(genes,file="/users/cnodzak/BoG_2018/NonZeroGenes.glmnet.042318.n330.txt",sep="\t",quote=F,row.names=F,col.names=T)


#predlasso<-predict.cv.glmnet(cvlasso, newx = data.matrix(test[,3:493]), s = "lambda.min", type="class")
#prede0.8<-predict.cv.glmnet(cvelastic0.8, newx = data.matrix(test[,3:493]), s = "lambda.min", type="class")
#prede0.5<-predict.cv.glmnet(cvelastic0.5,newx = data.matrix(test[,3:493]), s = "lambda.min", type="class")
#predridge<-predict.cv.glmnet(cvridge,newx = data.matrix(test[,3:493]), s = "lambda.min", type="class")

#predictions<-cbind.data.frame(test$exprs.SampleID,test$exprs.Subgroup,predlasso,prede0.8,prede0.5,predridge)
#colnames(predictions)<-c("sampleID","Subgroup","Lasso","alpha=0.8","alpha=0.5","Ridge")
#write.table(predictions,file="/users/cnodzak/BoG_2018/glmnetModelPredictions.042318.n330.txt",row.names=F,col.names=T,quote=F,sep="\t")



################################
### Implement RandomForest model
################################
library(randomForest)
library(ggplot2)
ranf<-randomForest(as.factor(exprs.Subgroup) ~ . ,data=train[,2:493])
error_df<-data.frame(error_rate=ranf$err.rate[,'OOB'],num_trees=1:ranf$ntree)
ranf
test.pred <- predict(ranf, test[,2:493],type="response")
table(observed = ytest, predicted = test.pred)



#pdf("/users/cnodzak/BoG_2018/randomForest.errorPlot.042318.n330.pdf")
#print(ggplot(error_df,aes(x=num_trees,y=error_rate))+geom_line()+ggtitle("Accuracy of random forest with addition of more trees"))
#dev.off()

#ranfALL<-randomForest(as.factor(exprs.Subgroup) ~ . ,data=data.matrix(lncProtexprs[2:493]))
#error_dfALL<-data.frame(error_rate=ranfALL$err.rate[,'OOB'],num_trees=1:ranfALL$ntree)
#ranfALL

#ranfALL<-randomForest(as.factor(exprs.Subgroup) ~ . ,data=data.matrix(lncProtexprs[2:493]),ntree=175)
#error_dfALL<-data.frame(error_rate=ranfALL$err.rate[,'OOB'],num_trees=1:ranfALL$ntree)
#ranfALL


#pdf("/users/cnodzak/BoG_2018/randomForest.ALL.ntree175.errPlot.050218.pdf")
#print(ggplot(error_dfALL,aes(x=num_trees,y=error_rate))+geom_line()+ggtitle("Accuracy of random forest with addition of more trees"))
#dev.off()


# Assess variable importance in randomForest Models

#ranfALL.IMP<-randomForest(as.factor(exprs.Subgroup) ~ . ,data=data.matrix(lncProtexprs[2:493]),importance=TRUE,ntree=175)
#ranf.IMP<-randomForest(as.factor(exprs.Subgroup) ~ . ,data=data.matrix(train[2:493]),importance=TRUE)

#pdf("/users/cnodzak/BoG_2018/varIMPplot.train.Acc.042318.n330.pdf")
#varImpPlot(ranf.IMP,type=1,n.var=50)
#dev.off()

#pdf("/users/cnodzak/BoG_2018/varIMPplot.train.Gini.042318.n330.pdf")
#varImpPlot(ranf.IMP,type=2,n.var=50)
#dev.off()

pdf("/users/cnodzak/BoG_2018/varIMPplot.ALL.ntree175.Acc.050218.pdf")
varImpPlot(ranfALL.IMP,type=1,n.var=30)
dev.off()

pdf("/users/cnodzak/BoG_2018/varIMPplot.ALL.ntree175.Gini.050218.pdf")
varImpPlot(ranfALL.IMP,type=2,n.var=30)
dev.off()

