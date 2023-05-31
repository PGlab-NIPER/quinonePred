# quinonePred
QuinonePred is a machine learning based tool. This tool can predict whether a small drug-like molecule can produce reactive metabolites of the quinone type during metabolism.

# Contents
* ``data`` This folder contains the dump of classification models and the file required for calculating the descriptors
* ``predict.R`` This is the R script for making predictions on the test set

# System dependencies
* Java run time environment (64 bit)
* R language (64 bit) (make sure to add the path of R language to the sytem environment variables)

# Required R libraries
* Use following command in the command prompt to install required R libraries 
R -e "install.packages(c('kernlab','randomForest','caret','rcdk','xgboost'),repos='https://cloud.r-project.org', dependencies=TRUE")
