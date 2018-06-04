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


library("pamr")
setwd(work_dir)
training_data = pamr.from.excel(training_file, 104,sample.labels=TRUE)
test_data = pamr.from.excel(test_file, 86,sample.labels=TRUE)
fit = pamr.train(training_data,threshold=0)
train_prediction = pamr.predict(fit, training_data$x, threshold=0)
test_prediction = pamr.predict(fit, test_data$x, threshold=0)
training_actual_predicted=cbind(training_data$samplelabels, training_data$y,as.character(train_prediction))
colnames(training_actual_predicted)=c("array_ids","actual","predicted")
test_actual_predicted=cbind(test_data$samplelabels, test_data$y,as.character(test_prediction))
colnames(test_actual_predicted)=c("array_ids","actual","predicted")
write.table(training_actual_predicted,file=paste(work_dir,"/training_actual_predicted_pam.tsv",sep=""), quote=F, sep="\t", row.names=F)
write.table(test_actual_predicted,file=paste(work_dir,"/test_actual_predicted_pam.tsv",sep=""), quote=F, sep="\t", row.names=F)
