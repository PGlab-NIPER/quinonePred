# quinonePred
QuinonePred is a machine learning based tool. This tool can predict whether a small drug-like molecule can produce reactive metabolites of the quinone type during metabolism.

# Contents
* ``data`` This folder contains the dump of classification models and the file required for calculating the descriptors
* ``predict.R`` This is the R script for making predictions on the test set

# System dependencies
* Java run time environment (64 bit)
* R language (64 bit) 
* make sure to add the path of R language to the sytem environment variables.
For adding path of R language to system environment variable use following command:                        
``pathman /au C:\Program Files\R\R-4.1.1\bin\x64\``                                                                               
Above command is for R version 4.1.1, change it according to the R version installed               
For example for R version 4.2.0:                                                       
``pathman /au C:\Program Files\R\R-4.2.0\bin\x64\``

# Required R libraries
* Use following command in the command prompt to install required R libraries:                                                         
R -e "install.packages(c('kernlab','randomForest','caret','rcdk','xgboost'),repos='https://cloud.r-project.org', dependencies=TRUE")
