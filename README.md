# EVE <img src="/images/EVE.png" width="20%"> 


[![GitHub repo version](https://img.shields.io/badge/release-v0.3-blue.svg)](https://github.roche.com/ED-EB/EVE)


# Getting started


### Clone this repo to your HOME

```console

ssh rosalind.gene.com

cd ~
git clone --branch dev-models https://github.roche.com/ED-EB/EVE.git
```
_current setup assumes users always clone this repo to your unix HOME_. The branch dev-models contains the latest update.


### Install Required Packages


Install necessary R libraries to R3.6.1 on Rosalind, as actual model fitting uses R3.6.1. R3.6.1 is not accessible from [the Rstudio server](http://gred-rstudio-p01.sc1.roche.com:8080/) and has to be accessed from the commandline:

ml R/R_3.6.1_Bioc_3.10/R-3.6.1-Bioc-3.10-prd-20200217
R

Then verify that the following libraries are available in R3.6.1. You will also need the nextdoor package for lasso. R4.0.0 is available from Rosalind Rstudio server; we should be able to submit eve jobs (run the *.r files in the example fold) and generate reports (run the *.rmd files) there.


```r
pkgs <- c("prodlim", "survival", "survminer", "glmnet", "pec", "pheatmap",
"tidyverse", "randomForestSRC", "broom", "caret", "ggrepel", "pROC")

devtools::install_github("zouw2/NextDoor/nextdoor")
```


### Algorithms

1. **glmnet (lasso.r)** for binary, multinomial, gaussian and survival outcomes
2. **survival random forest (rfeSRCC.r)** for survival outcome based on _randomforestSRC_
3. **XGBoost (XGBoost.r)** for binary, multiclass, regression, and survival outcomes 

### Execute


Go to [`examples`](https://github.roche.com/ED-EB/EVE/tree/master/examples) folder, copy a *.R to your [project home](https://docs.google.com/spreadsheets/d/1OAmZDae7MF9NXBBwR6YpHjLxbUVFbFw7_y6HzFJegHY/edit#gid=0&range=A4). This R script 

1. collects input parameters, including an engine script which may or may not be implemented in R, 
2. prepares input data file if necessary
3. splits repeated CV or LOOCV by repeated calling the engine script on different CPU

For more definition of user inputs in a *.R, please see [Google Doc](https://docs.google.com/spreadsheets/d/1OAmZDae7MF9NXBBwR6YpHjLxbUVFbFw7_y6HzFJegHY/edit#gid=0).

After you modify *.R appropriately, then click **Source** in RStudio, 
it will start running the models on different CPUs.

Note: once the job is submitted you will see `log` folder under your project folder. 
When job is done, you will see `results` folder containing all outputs files. 
**Please do NOT Push these two folders to this repo because those files can be huge**

### Harvest Results

  - go to [`examples/`](https://github.roche.com/ED-EB/EVE/tree/dev-models/examples) (these .Rmd files are reporting template).
  - choose corresponding `.Rmd` file. 
  (for example, if `eve_binaryclass_lasso.R` was submitted, 
  then use `eve_binaryclass_lasso.Rmd` to gather the output)
  - modify the project path and name accordingly
  - knit to pdf and enjoy it

[Check out what the report looks like.](https://github.roche.com/ED-EB/EVE/blob/master/examples/reports/test_eval_multiclass_xgboost.pdf)


### Contribute to EVE

We highly recommend users add new models (e.g. add SVM, Neural Network, etc..), 
please make sure the model outputs prediction and feature importance (and grid search parameters if any) 
that follow [this data format](https://github.roche.com/ED-EB/EVE/tree/master/examples/model_output_format). 
With consistent output format, user can use the same reporting template(s) to harvest and compare the result across different models.
