# ShapeGTB
Machine learning tool for prioritization of functional single-nucleotide polymorphisms in human promoters

Supplementary code for publication:

**The role of local DNA shape in prioritization of functional variants in human promoters**, Maja Malkowska, Julian Zubek, Dariusz Plewczynski, Lucjan Wyrwicz

## Requirements

* xgboost
* [DNAshapeR](http://tsupeichiu.github.io/DNAshapeR/)

## Files

* *data/* -- original data used to train and test the classifier
* *ShapeGTB.model* -- trained xgboost model
* *train_xgboost.R* -- script for retraining the model on the original data
* *feature_encoding.R* -- functions used to encode features used by the classifier
* *predict.R* -- an example of model usage

## Usage

Input data format is given by the file *test.csv*. Each line contains a sequence of 9 nucleotides (wild type) and the central nucleotide occurring in the mutated variant:

```
GCCAAGTGA,C
GAATGACTG,A
ACTGAATTC,T
...
```

You can calculate all the features and run the classifier using *predict.R* script:

> Rscript predict.R test.csv

Output will be saved as *test.csv.out*. Column *ShapeGTB_score* in the output table can be interpreted as a likelihood that the variant is functional (range 0-1).
