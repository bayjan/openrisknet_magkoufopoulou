args = commandArgs(trailingOnly = TRUE) # Provide arguments only using --args
if(length(args)!=3){
  print(args)
  print("Please specify locations of: work_dir, training_file and test_file with --args option of R");
  print("There must be at least 3 arguments specified. Quitting!");
  quit(save = "no", status = -1, runLast = TRUE);
}

print("Provided args:");
print(args);
work_dir = args[1];
training_file = args[2];
test_file = args[3];

# Allow to work with many variables, default is 5000 and I changed it to 50000
options(expressions=500000);

#setwd("/home/jbayjanov/projects/tgx/dixa_classification/data/caret_run")
setwd(work_dir)


library(caret)
# library(doMC)
# registerDoMC(cores = 24)


set.seed(998)
## read input

#training<-read.table("all_caret_train_data.tsv",header = TRUE, sep = "\t", row.names = 1)
training1 <-read.table(file=training_file, header = TRUE, sep = "\t", row.names = 1)
training = data.frame(t(training1))

dim(training)


table(training$Class)

### 10 fold 3 repeats cross-validation

control <- trainControl(method="repeatedcv", number=10, repeats=3,## Estimate class probabilities
                        classProbs = TRUE,
                        ## Evaluate performance using 
                        ## the following function
                        summaryFunction = twoClassSummary)

# train the Pam
set.seed(7)
modelPam <- train(Class~., data=training, method="pam", trControl=control,k=0)
models_list <- list(PAM=modelPam)

results <- resamples(models_list)
# summarize the distributions
summary(results)

png(file = "bwplot.png", bg = "transparent")
# boxplots of results
bwplot(results)
dev.off()


png(file = "dotplot.png", bg = "transparent")
# dot plots of results
dotplot(results)
dev.off()


## Give prediction results for all models
#allModels <- list(PAM=modelPam, KNN=modelKnn, SVM=modelSvm,RF=modelRf,svmL=modelsvmLinear,pls=modelpls)

allModels <- models_list; 
pred <- predict(allModels,type="prob")
write.table(pred,"training_set_prediction.txt",sep = "\t", col.names=NA)

### Variable importance for all models
for (i in allModels){
  ImpMeasure<-data.frame(varImp(i)$importance)
  write.table(ImpMeasure,paste(i$method,"_variable_importance.txt", sep = ""),sep="\t", col.names=NA)
}



## Confusion matrix and parameters choosen for each model

for (i in allModels){
  testPred <- predict(i, training)
  cf<- confusionMatrix(testPred, training$Class)
  
  write.table(as.matrix(cf, what = "xtabs"),paste(i$method,"_performance.txt", sep = ""),sep="\t", append = TRUE)
  write.table(as.matrix(cf, what = "overall"),paste(i$method,"_performance.txt", sep = ""),sep="\t", append = TRUE)
  write.table(as.matrix(cf, what = "classes"),paste(i$method,"_performance.txt", sep = ""),sep="\t", append = TRUE)
  
  write(i$method, file = "parameters.txt",append = TRUE, sep = "\t")
  write(rbind(colnames(as.matrix(i$bestTune)),as.matrix(i$bestTune)), file = "parameters.txt",append = TRUE, sep = "\t")
  write("      ", file = "parameters.txt",append = TRUE, sep = "\t")
}


### predict the validationset for each model
#testing<-read.table("TK6_test_transposed.txt",header = TRUE, sep = "\t", row.names = 1)
testing1<-read.table(file=test_file, header = TRUE, sep = "\t", row.names = 1)
testing = data.frame(t(testing1))

pred_val <- predict(allModels,newdata=testing,type="prob")


write.table(pred_val,"validation_set_prediction.txt",sep = "\t", col.names=NA)


for (i in allModels){
  valPred <- predict(i, testing)
  cf<- confusionMatrix(valPred, testing$Class)
  
  write.table(as.matrix(cf, what = "xtabs"),paste(i$method,"_validation_performance.txt", sep = ""),sep="\t", append = TRUE)
  write.table(as.matrix(cf, what = "overall"),paste(i$method,"_validation_performance.txt", sep = ""),sep="\t", append = TRUE)
  write.table(as.matrix(cf, what = "classes"),paste(i$method,"_validation_performance.txt", sep = ""),sep="\t", append = TRUE)
}
