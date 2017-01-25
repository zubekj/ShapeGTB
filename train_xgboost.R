library(xgboost)

data <- read.csv("./data/ShapeGTB_train_data_sequence.csv", header=TRUE)
label <- data$class
data <- data[,4:(ncol(data)-1)]

clf <- xgboost(data=as.matrix(data), label=label, max.depth=8, nrounds=300, eta=0.1, objective="binary:logistic", save_name="ShapeGTB.model")
xgb.save(clf, "ShapeGTB.model")
