+++
title = "Sample size calculation with PS software"
date = 2019-07-04T14:38:20+05:30
draft = false

# Tags and categories
# For example, use `tags = []` for no tags, or the form `tags = ["A Tag", "Another Tag"]` for one or more tags.
tags = ["sample size", "power", "clinical research"]
categories = ["Statistics"]

# Featured image
# Place your image in the `static/img/` folder and reference its filename below, e.g. `image = "example.jpg"`.
# Use `caption` to display an image caption.
#   Markdown linking is allowed, e.g. `caption = "[Image credit](http://example.org)"`.
# Set `preview` to `false` to disable the thumbnail in listings.
[header]
image = ""
caption = ""
preview = true

+++


## Introduction

Sample size calculation is an integral part in the inferential statistics, in which 
we are estimating a **population parameter** (our interest) from a **sample statistic** 
(data available to us).
We carry out experiment on a finite numbered sample and calculate a summary measure 
out of it (sample statistic, lets say, sample mean) with an intent to estimate the 
unknown population parameter (population mean).

The basic premise of statistical estimation is that as the sample size increases, the sample 
statistic will be reflecting the population parameter more accurately (its variation 
will be less around the population parameter). So, sample size may be titrated to achieve 
the required accuracy in estimation.

The software, PS, whose link is given below will help us in standardising calculation of sample size 
for our research works.

## The software

The intended software is **PS - Power and Sample Size Program**, a freely downloadable 
and distributable software for Windows user made by the Biostatistics Department  of 
Vanderbilt University.  The software is downloadable from <http://biostat.mc.vanderbilt.edu/wiki/Main/PowerSampleSize> 
as `pssetup3.exe` file.

## The study designs for which the software can be used

1. **Case control studies:**  Corrected and uncorrected  chi squared contingency table tests and Fisher's exact test. 
The alternative hypotheses may be specified in terms of **odds ratio** or **exposure prevalence rates**.

2. **Matched case control studies:** McNemar test. The alternative hypothesis may be specified in terms of 
**odds ratio**.

3. **Multiple 2x2 tables:** Mantel Haenszel test. Assume that each 2x2 table consists of cases and controls 
selected from a different stratum that is defined by one or more confounding variables. The odds ratio 
of disease in the the exposed and subjects compared to unexposed subjects is assumed to be equal within 
all the strata. The alternative hypothesis is specified in terms of this **odds ratio**.

4. **Cohort studies with Dichotomous Outcomes:** Independent contingency table tests, McNemar's test. 
The alternative hypotheses may be specified in terms of **relative risks** or **outcome probabilities**.

5. **Linear regression (1 treatment):** Testing the slope of a simple linear regression line. The investigator 
wishes to **detect a regression slope a given magnitude**. The values of the independent variable 
may either be specified by the investigator or determined observationally 
when the study is performed. In the later case, the investigator must estimate the standard deviation of 
the independent variable(s).

6. **Linear regression (2 treatments):** Comparing the slopes and intercepts of two linear regression 
lines. The investigator wishes to determine whether the **difference in slopes and intercepts** is of a 
given magnitude. The values of the independent variables may either be specified by the investigator or determined observationally 
when the study is performed. In the later case, the investigator must estimate the standard deviations of 
the independent variables.

7. **Survival studies:** Evaluating independent cohorts using log rank test. The ratio of number of controls 
per experimental subject may be specified by the investigator. The alternative hypothesis may be specified 
in terms of **hazard ratio** of control subjects relative to experimental subject or **median survival times** 
for the control and experimental subjects.

8. **Continuous response measures in two groups:** Paired and independent t tests. The ratio of number 
of control subjects per experimental subject may be specified by the user.

## How to cite the software

Dupont WD, Plummer WD: 'Power and Sample Size Calculations: A Review and Computer Program', Controlled Clinical 
Trials 1990; 11: 116-28

or

Dupont WD, Plummer WD: 'Power and Sample Size Calculations for studies involving Linear Regression:', Controlled Clinical 
Trials 1998; 19: 589-601


