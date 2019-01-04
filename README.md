# EVE

[![GitHub repo version](https://img.shields.io/github/release/qubyte/rubidium.svg)](https://github.roche.com/ED-EB/EVE)

# Getting started


### Clone this repo to your HOME

(currently the paths used in the code only allow users to clone this repo to HOME)

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
  

### Execute

Go to `tests` folder, open sbatch_xxxx.R

If it is python engine (e.g. xgboost), remember to use command line (for example):

```console
module load apps/python
Rscript ~/EVE/tests/sbatch_xgb_regression.R
```

If it is R engines, just open sbatch_xxx.R, and then click **Source** in RStudio.

Note: once the job is submitted you will see `log` folder under tests. 
When job is done, you will see `results` folder containing all outputs files. 
**Please do NOT Push these two folders to this repo because those files can be huge**

### Harvest Results

  - go to `tests/reports/` (these .Rmd files are reporting template).
  - choose corresponding `.Rmd` file. 
  (for example, if `sbatch_xgb_regression.R` was used, 
  then use `test_eval_regression_xgboost.Rmd`)
