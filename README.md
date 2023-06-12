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
```bash
pathman /au C:\Program Files\R\R-4.1.1\bin\x64\
```                                                                               
Above command is for R version 4.1.1, change it according to the R version installed               
For example for R version 4.2.0:                                                       
```bash
pathman /au C:\Program Files\R\R-4.2.0\bin\x64\
```

# Required R libraries
* Use following command in the command prompt to install required R libraries:                                                         
```bash
R -e "install.packages(c('kernlab','randomForest','caret','rcdk','xgboost'),repos='https://cloud.r-project.org', dependencies=TRUE)"
```

# Usage
To make prediction using quinonePred first save the structures of query molecules in .sdf format

1. Download the repository and save all the files and folders of repository in one folder.
2. Then Use the follwing command by changing current working directory to the folder where repository predict.R is saved.
```bash
Rscript predict.R
```
3. User wil be prompted to select ``.sdf`` file of query molecules                                                                             
Once the .sdf file is selected and the process is completed, results will be saved in ``quinone_predictions.csv`` file.
