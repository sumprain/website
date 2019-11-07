+++
title = "Confounding and collinearity"
date = 2018-01-03T19:48:15+05:30
draft = false

# Authors. Comma separated list, e.g. `["Bob Smith", "David Jones"]`.
authors = ["Suman Kumar"]

# Publication type.
# Legend:
# 0 = Uncategorized
# 1 = Conference proceedings
# 2 = Journal
# 3 = Work in progress
# 4 = Technical report
# 5 = Book
# 6 = Book chapter
publication_types = ["4"]

# Publication name and optional abbreviated version.
publication = "Confounding and collinearity"
publication_short = ""

# Abstract and optional shortened version.
abstract = "In this article, I will discuss about the bias introduced in estimation of coefficient of a given explanatory variable due to the presence of confounding factors.  After that, I will try to demonstrate about the effect of variable collinearity on estimation of coefficient."

abstract_short = "In this article, I will discuss about the bias introduced in estimation of coefficient of a given explanatory variable due to the presence of confounding factors.  After that, I will try to demonstrate about the effect of variable collinearity on estimation of coefficient."

# Featured image thumbnail (optional)
image_preview = ""

# Is this a selected publication? (true/false)
selected = true

# Projects (optional).
#   Associate this publication with one or more of your projects.
#   Simply enter the filename (excluding '.md') of your project file in `content/project/`.
#   E.g. `projects = ["deep-learning"]` references `content/project/deep-learning.md`.
projects = []

# Tags (optional).
#   Set `tags = []` for no tags, or use the form `tags = ["A Tag", "Another Tag"]` for one or more tags.
tags = ["R", "statistics", "confounding", "collinearity"]

# Links (optional).
url_pdf = ""
url_preprint = ""
url_code = ""
url_dataset = ""
url_project = ""
url_slides = ""
url_video = ""
url_poster = ""
url_source = ""

# Custom links (optional).
#   Uncomment line below to enable. For multiple links, use the form `[{...}, {...}, {...}]`.
# url_custom = [{name = "Custom Link", url = "http://example.org"}]

# Does this page contain LaTeX math? (true/false)
math = false

# Does this page require source code highlighting? (true/false)
highlight = true

# Featured image
# Place your image in the `static/img/` folder and reference its filename below, e.g. `image = "example.jpg"`.
[header]
image = ""
caption = ""

+++

## Introduction

In this article, I will discuss about the bias introduced in estimation of coefficient of a given explanatory variable due to the presence of confounding factors.  After that, I will try to demonstrate about the effect of variable collinearity on estimation of coefficient.

## Underlying structure

Let us assume that we know about the underlying relationship between the explanatory variable (E), outcome variable (O), confounder (C), variable correlated with E (VE) and variable correlated with O (VO).

Let us assume that the real underlying relation between various variables are as given below:

```
O = 2.E + 2.C + 2.VO + 0.VE + error

E, C and VE are correlated with each other with correlation coef of 0.7 each.

All the random variables are normally distributed, with mean 0 and sd 5.  error is distributed as standard normal distribution.
```
The directed acyclic graph depicts the relationship between the variables.

![Underlying structure](/img/conf_orig.png)

## Samples for simulation

```r
library(MASS)
library(dplyr)
```

```r
set.seed(100)
sam <- mvrnorm(n = 10000, mu = c(0, 0, 0), Sigma = matrix(c(5, 4, 4, 4, 5, 4, 
    4, 4, 5), 3, 3, byrow = T))
df <- data.frame(E = sam[, 1], C = sam[, 2], VE = sam[, 3])
df <- mutate(df, VO = rnorm(10000, 0, 5), O = 2 * E + 2 * C + 2 * VO + rnorm(10000))
```

## Scenarios

### <a name="scen1"></a>Scenario 1 (Only E as covariate, typical univariate analysis)

```r
mod.E <- lm(O ~ E, data = df)
summary(mod.E)
```

```
## 
## Call:
## lm(formula = O ~ E, data = df)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -40.47  -7.05  -0.05   7.01  45.04 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   0.0069     0.1034    0.07     0.95    
## E             3.6420     0.0465   78.40   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 10.3 on 9998 degrees of freedom
## Multiple R-squared:  0.381,	Adjusted R-squared:  0.381 
## F-statistic: 6.15e+03 on 1 and 9998 DF,  p-value: <2e-16
```

```r
confint(mod.E)
```

```
##               2.5 % 97.5 %
## (Intercept) -0.1959 0.2097
## E            3.5509 3.7330
```

The model overestimates coefficient of E, its value is the sum of coef of E and C.  This type of bias is known as **bias due to confounding** and occurs when confounders are not included as co-variates.

### Scenario 2 (E and C as covariates)

This model is incorrect as it does not include the other covariate, VO.

```r
mod.EC <- lm(O ~ E + C, data = df)
summary(mod.EC)
```

```
## 
## Call:
## lm(formula = O ~ E + C, data = df)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -35.99  -6.75  -0.11   6.71  43.32 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  -0.0194     0.0998   -0.19     0.85    
## E             2.0293     0.0741   27.40   <2e-16 ***
## C             2.0231     0.0740   27.35   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 9.98 on 9997 degrees of freedom
## Multiple R-squared:  0.424,	Adjusted R-squared:  0.424 
## F-statistic: 3.68e+03 on 2 and 9997 DF,  p-value: <2e-16
```

```r
confint(mod.EC)
```

```
##              2.5 % 97.5 %
## (Intercept) -0.215 0.1762
## E            1.884 2.1745
## C            1.878 2.1681
```


Due to addition of C, the estimate of E is unbiased. Exclusion of VO has not influenced the coefficient of E.

### Scenario 3 (E and VO as covariates)

E and VO are not correlated. C has been excluded from the model.


```r
mod.EVO <- lm(O ~ E + VO, data = df)
summary(mod.EVO)
```

```
## 
## Call:
## lm(formula = O ~ E + VO, data = df)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -9.638 -1.938  0.005  1.954 12.780 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.03893    0.02884    1.35     0.18    
## E            3.59175    0.01295  277.27   <2e-16 ***
## VO           2.00167    0.00581  344.39   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 2.88 on 9997 degrees of freedom
## Multiple R-squared:  0.952,	Adjusted R-squared:  0.952 
## F-statistic: 9.88e+04 on 2 and 9997 DF,  p-value: <2e-16
```

```r
confint(mod.EVO)
```

```
##                2.5 %  97.5 %
## (Intercept) -0.01761 0.09546
## E            3.56636 3.61714
## VO           1.99028 2.01307
```


There is again overestimation of coefficient of E (as in [scenario 1](#scen1)), due to exclusion of C. Estimate of VO is unbiased as it has no correlation with either E, VO or VE.

### <a name="scen4"></a>Scenario 4 (E and VE as covariates)

E and VE are correlated, but VE is not associated with O.


```r
mod.EVE <- lm(O ~ E + VE, data = df)
summary(mod.EVE)
```

```
## 
## Call:
## lm(formula = O ~ E + VE, data = df)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -41.55  -7.00  -0.02   6.93  43.46 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   0.0106     0.1026     0.1     0.92    
## E             2.8584     0.0761    37.6   <2e-16 ***
## VE            0.9857     0.0761    12.9   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 10.3 on 9997 degrees of freedom
## Multiple R-squared:  0.391,	Adjusted R-squared:  0.391 
## F-statistic: 3.21e+03 on 2 and 9997 DF,  p-value: <2e-16
```

```r
confint(mod.EVE)
```

```
##               2.5 % 97.5 %
## (Intercept) -0.1905 0.2117
## E            2.7093 3.0075
## VE           0.8364 1.1350
```

Estimates of both E and VE are biased. But coefficient of E is less biased than [scenario 1](#scen1). Reason for biasness of coefficient of VE is unknown to me.

### <a name="scen5"></a>Scenario 5 (E, VE, VO, C as covariates) (Complete model)


```r
mod.comp <- lm(O ~ E + C + VO + VE, data = df)
summary(mod.comp)
```

```
## 
## Call:
## lm(formula = O ~ E + C + VO + VE, data = df)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -3.434 -0.677 -0.014  0.678  3.616 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.01297    0.01000    1.30     0.19    
## E            1.98794    0.00826  240.62   <2e-16 ***
## C            1.99993    0.00830  240.81   <2e-16 ***
## VO           2.00033    0.00202  992.51   <2e-16 ***
## VE           0.01214    0.00831    1.46     0.14    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1 on 9995 degrees of freedom
## Multiple R-squared:  0.994,	Adjusted R-squared:  0.994 
## F-statistic: 4.29e+05 on 4 and 9995 DF,  p-value: <2e-16
```

```r
confint(mod.comp)
```

```
##                 2.5 %  97.5 %
## (Intercept) -0.006635 0.03258
## E            1.971745 2.00413
## C            1.983651 2.01621
## VO           1.996381 2.00428
## VE          -0.004159 0.02844
```

Unbiased coefficients of all the covariates.

### Scenario 6 (E, VE, C as covariates)


```r
mod.EVEC <- lm(O ~ E + VE + C, data = df)
summary(mod.EVEC)
```

```
## 
## Call:
## lm(formula = O ~ E + VE + C, data = df)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -36.18  -6.74  -0.08   6.71  43.21 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  -0.0185     0.0998   -0.19     0.85    
## E             1.9892     0.0824   24.13   <2e-16 ***
## VE            0.0920     0.0830    1.11     0.27    
## C             1.9817     0.0829   23.92   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 9.98 on 9996 degrees of freedom
## Multiple R-squared:  0.424,	Adjusted R-squared:  0.424 
## F-statistic: 2.45e+03 on 3 and 9996 DF,  p-value: <2e-16
```

```r
confint(mod.EVEC)
```

```
##                2.5 % 97.5 %
## (Intercept) -0.21411 0.1771
## E            1.82758 2.1507
## VE          -0.07061 0.2546
## C            1.81930 2.1441
```

Addition of C has removed the biasedness from the coefficients of E and VE (see [Scenario 4](#scen4)).

### Scenario 7 (E, VO and C as covariates)

```r
mod.EVOC <- lm(O ~ E + VO + C, data = df)
summary(mod.EVOC)
```

```
## 
## Call:
## lm(formula = O ~ E + VO + C, data = df)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -3.427 -0.676 -0.014  0.680  3.641 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.01285    0.01000    1.29      0.2    
## E            1.99323    0.00742  268.50   <2e-16 ***
## VO           2.00036    0.00202  992.52   <2e-16 ***
## C            2.00539    0.00741  270.45   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1 on 9996 degrees of freedom
## Multiple R-squared:  0.994,	Adjusted R-squared:  0.994 
## F-statistic: 5.72e+05 on 3 and 9996 DF,  p-value: <2e-16
```

```r
confint(mod.EVOC)
```

```
##                 2.5 %  97.5 %
## (Intercept) -0.006751 0.03246
## E            1.978682 2.00779
## VO           1.996410 2.00431
## C            1.990857 2.01993
```


Performance of this model is almost the same as [Scenario 5](#scen5)

## Performance of different scenarios

Models | Complete (E, C, VO, VE) | E, C, VE | E, C, VO | E, C | E, VE | E, VO | E only
-------|-------------------------|----------|----------|------|-------|-------|--------
E (2) | 1.9879, 0.0083 | 1.9892, 0.0824 | 1.9932, 0.0074 | 2.0293, 0.0741 | 2.8584, 0.0761 | 3.5918, 0.013 | 3.642, 0.0465 
C (2) | 1.9999, 0.0083 | 1.9817, 0.0829 | 2.0054, 0.0074 | 2.0231, 0.074 | - | - | - 
VE (0) | 0.0121, 0.0083 | 0.092, 0.083 | - | - | 0.9857, 0.0761 | - | -
VO (2) | 2.0003, 0.002 | - | 2.0004, 0.002 | - | - | 2.0017, 0.0058 | -

E (2): means covariate E and its actual value.

--,--: coefficient, standard error.

## Conclusions

Models with C, will have unbiased estimate of E. Models with correlated variables (E, C, VE) has larger SE due to near collinearity.
