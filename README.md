# EVE <img src="/images/EVE.png" width="20%"> 


[![GitHub repo version](https://img.shields.io/badge/release-v0.3-blue.svg)](https://github.roche.com/ED-EB/EVE)


# Getting started


### Clone this repo to your HOME

```console
ssh rescomp5003.gene.com
cd ~
git clone https://github.roche.com/ED-EB/EVE.git
```
_current setup assumes users always clone this repo to your unix HOME_


### Install Required Packages

Install some Python libraries

```console
cd ~/EVE
module load apps/python3
pip3 install --user -r requirements.txt
```

Install some R libraries

Open [RStudio on rescomp](http://rescomp5105.gene.com:8080), 
make sure you are using **R 3.5.1** or above, then:

```r
pkgs <- c("prodlim", "survival", "survminer", "glmnet", "pec", "pheatmap",
"tidyverse", "randomForestSRC", "broom", "caret", "ggrepel", "pROC")

install.packages(pkgs)
```


### Algorithms

1. **glmnet (lasso.r)** for binary, multinomial, gaussian and survival outcomes
2. **survival random forest (rfeSRCC.r)** for survival outcome based on _randomforestSRC_
3. **XGBoost (XGBoost_xxx.py)** for binary, multiclass, regression, and survival outcomes 

### Execute

Go to [`examples`](https://github.roche.com/ED-EB/EVE/tree/master/examples) folder, copy sbatch_xxxx.R to your [project home](https://docs.google.com/spreadsheets/d/1OAmZDae7MF9NXBBwR6YpHjLxbUVFbFw7_y6HzFJegHY/edit#gid=0&range=A4). This R script 
1. collects input parameters, including an engine script which may or may not be implemented in R, 
2. prepares input data file if necessary
3. splits repeated CV or LOOCV by repeated calling the engine script on different CPU

For more definition of user inputs in sbatch_xxx.R, please see [here](https://docs.google.com/spreadsheets/d/1OAmZDae7MF9NXBBwR6YpHjLxbUVFbFw7_y6HzFJegHY/edit#gid=0).

After you modify sbatch_xxx.R appropriately, then click **Source** in RStudio, 
it will start running the models on different CPUs.

Note: once the job is submitted you will see `log` folder under your project folder. 
When job is done, you will see `results` folder containing all outputs files. 
**Please do NOT Push these two folders to this repo because those files can be huge**

### Harvest Results

  - go to [`examples/reports/`](https://github.roche.com/ED-EB/EVE/tree/master/examples/reports) (these .Rmd files are reporting template).
  - choose corresponding `.Rmd` file. 
  (for example, if `sbatch_xgb_regression.R` was used, 
  then use `test_eval_regression_xgboost.Rmd`)
  - modify the project path and name accordingly
  - knit to pdf and enjoy it

[Check out what the report looks like.](https://github.roche.com/ED-EB/EVE/blob/master/examples/reports/test_eval_multiclass_xgboost.pdf)


### Contribute to EVE

We highly recommend users add new models (e.g. add SVM, Neural Network, etc..), 
please make sure the model outputs prediction and feature importance (and grid search parameters if any) 
that follow [this data format](https://github.roche.com/ED-EB/EVE/tree/master/examples/model_output_format). 
With consistent output format, user can use the same reporting template(s) to harvest and compare the result across different models.
