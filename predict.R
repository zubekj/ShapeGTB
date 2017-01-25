library(xgboost)
source("feature_encoding.R")

# Example usage
args <- commandArgs(trailingOnly=TRUE)
snps <- read.csv(args[1], stringsAsFactors=FALSE, header=FALSE)
data <- encode_features(snps)

clf <- xgb.load("ShapeGTB.model")
results <- as.data.frame(predict(clf, as.matrix(data)))
colnames(results) <- "ShapeGTB_score"

write.csv(cbind(data, results), paste0(args[1], ".out"), row.names=FALSE)
