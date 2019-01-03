# EVE


# Getting started


### Clone this repo to your HOME

```diff
- __currently the paths used in EVE assume users always clone this repo to your unix HOME__
```

```console
ssh rescomp5003.gene.com
cd ~
git clone https://github.roche.com/ED-EB/EVE.git
```


### Install Required Packages

If you will use **Python**:

```console
cd ~/EVE

module load apps/python

pip install --user -r requirements.txt
```

If you will use **R**:

Open [RStudio on rescomp](http://rescomp5105.gene.com:8080), then:

```r
pkgs <- c("prodlim", "survival", "survminer", "glmnet", "pec", 
"tidyverse", "randomForestSRC", "broom", "caret", "ggrepel", "pROC")

install.packages(pkgs)
```


### Run jobs

#### (tests)

ToDo: Need to implement proper testing


  - for testing whether each functionality runs correctly
  - for testing performance of each algorithm (assume each piece of the code works correctly, 
  how do we quickly compare the performance of the model?)
  

### Currently available analysis engines

1. glmnet for binary, multinomial, gaussian and survival outcomes (lasso.r), based on glmnet
2. survival random forest (rfeSRCC.r), based on randomforestSRC
3. xgboost for categorical, gaussian and survival outcomes (XGBoost_xxx.py)

### Execute

Go to `tests` folder, copy sbatch_xxxx.R to your project folder. This r script 
1. collects input parameters, including an engine script which may or may not be implemented in R, 
2. prepares input data file if necessary
3. splits repeated CV or LOOCV by repeated calling the engine script on different CPU

For more definition of user inputs in sbatch_xxx.R, please see [here](https://docs.google.com/spreadsheets/d/1OAmZDae7MF9NXBBwR6YpHjLxbUVFbFw7_y6HzFJegHY/edit#gid=0).

Open sbatch_xxx.R, modify to fit your project and then click **Source** in RStudio. 
It can now work for both R and Python engines. 

Note: once the job is submitted you will see `log` folder under tests. 
When job is done, you will see `results` folder containing all outputs files. 
**Please do NOT Push these two folders to this repo because those files can be huge**

### Harvest Results

  - go to `tests/reports/` (these .Rmd files are reporting template).
  - choose corresponding `.Rmd` file. 
  (for example, if `sbatch_xgb_regression.R` was used, 
  then use `test_eval_regression_xgboost.Rmd`)
