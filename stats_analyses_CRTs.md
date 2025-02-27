---
title: "statistical analyses cluster RCTs"
author: "A.Amstutz"
date: "2023-10-18"
output:
  html_document:
    keep_md: yes
    toc: yes
    toc_float: yes
    code_folding: hide
  pdf_document:
    toc: yes
---

# Analyses of cluster randomized trials (CRTs) including stepped-wedge cluster randomized trials (SW-CRTs)
This is an (incomplete) excerpt from an Inserm France workshop, prepared and taught by the wonderful team in Tours and Bordeaux (Bruno Giraudeau, Laurent Billot, Agnès Caille and several of their PhD students).

## There are two general ways to analyse CRTs
1. On cluster level ("cluster level = unit of analysis")
This is done simply by comparing cluster level data (or aggregating individual level data on cluster level) across arms using simple tests (t-test, Wilcoxon, etc.). Recommended if total number of cluster small, i.e., below ~ 10.
2. On individual level ("individual level = unit of analysis")
This is more powerful and more flexible. And recommended when there are enough clusters.
There are two approaches:
* Using cluster-specific models, i.e., mixed-models (GLMM): Interpretation: effect if an individual moves from a control CLUSTER to intervention CLUSTER. -> conditional effect
* Using population-averaged models, i.e., Generalized estimating equations (GEE): Interpretation: effect if an individual from the target population moves from control to intervention. -> marginal effect
There is an entire literature on when/how to use both of these models, and their benefit/challenges, most use mixed-models. "GEE expert": Liz Turner (https://scholars.duke.edu/person/liz.turner/publications)
3. The same is true for SW-CRT, but more complex. In addition you need to account for time and the correlation structure is more complex. Check out: https://steppedwedgehog.blog/what-is-a-stepped-wedge-trial/

## ICC reporting
1. It is good practice to report the ICC (and its 95%CI) in the results publication of a CRT, at least for the primary outcome, better, for all outcomes. For other trialists to use it, see e.g., ICC database: https://monash-biostat.shinyapps.io/CLOUDbank/
2. It is good practice to report the ICC by arm.
3. There are several ways how to calculate the ICC, most straight-forward way is using one-way ANOVA by group

## Parallel CRT with baseline period
1. "A common enhancement of a simple parallel CRT is to add an assessment of participants’ outcomes in a baseline period (before randomisation). Even if different participants are assessed at baseline and follow-up [i.e. cross-sectional sampling], the fact that they are sampled from the same cluster allows some control for cluster differences." -> https://www.bmj.com/content/360/bmj.k1121
2. This is illustratively shown in the sample size calculator: https://clusterrcts.shinyapps.io/rshinyapp/ (switch between "Parallel" and "Parallel with baseline measure") -> can yield a substantial increase in power! See last chapter below.

# Load packages

```r
library(tidyverse)
library(here)
library(readr)
library(sjPlot) # for tab_model()

library(lmerTest) # GLMM for CRTs with cont outcome: cluster-specific model (conditional)
library(geepack) # GEE for CRTs: population-averaged model (marginal) incl. sandwich estimator and exchangeable correlation structure
library(ICC) # one-way ANOVA for the calculation of the ICC, using mean squares within and between clusters
library(swCRTdesign) # stepped-wedge design plot

library(pwr)
```

# Parallel CRT
We want to evaluate the impact of a 6-hours fasting period prior to extubation in mechanically ventilated intensive care patients. To this end, the AMBROISIE study has been set up. The "ambroisie.csv" database contains some of the variables collected as part of this trial. Patients are included in the study at the time of the medical decision to extubate. The primary endpoint was extubation failure (reintubation or death within 7 days of extubation). Caloric intake on the day before extubation was also recorded.
Sampling frame: cohort sampling (the same people recruited and followed up, but only assessed once, at the end)

## Variable description
1. CENTER : Center
2. PATIENT : Patient number in corresponding center
3. GROUP : Center randomization group
4. BMI : BMI
5. CALBEFORE : Caloric intake the day before extubation (kcal)
6. INTUBATIONJ7 : Reintubation before D7
7. DEATHJ7 : Death on D7
8. UNIVERSITY : University hospitals or not (stratification variable at the cluster level)

## Load data


## Continuous outcome

```r
# reformat
df$GROUP <- as.factor(df$GROUP)
df$CENTER <- as.factor(df$CENTER)

## Use outcome CALBEFORE: An intermediate outcome, on the causal pathway between randomisation and target outcome
# GLMM 
calintake.glmm <- lmer(CALBEFORE ~ (1|CENTER) + GROUP, data = df)
tab_model(calintake.glmm)
```

<table style="border-collapse:collapse; border:none;">
<tr>
<th style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm;  text-align:left; ">&nbsp;</th>
<th colspan="3" style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm; ">CALBEFORE</th>
</tr>
<tr>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  text-align:left; ">Predictors</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">Estimates</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">CI</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">p</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">(Intercept)</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1771.65</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1666.16&nbsp;&ndash;&nbsp;1877.14</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">GROUP2 × Fasting</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">&#45;284.03</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">&#45;436.01&nbsp;&ndash;&nbsp;-132.06</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>
<tr>
<td colspan="4" style="font-weight:bold; text-align:left; padding-top:.8em;">Random Effects</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&sigma;<sup>2</sup></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">341560.36</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&tau;<sub>00</sub> <sub>CENTER</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">25255.34</td>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">ICC</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.07</td>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">N <sub>CENTER</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">22</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm; border-top:1px solid;">Observations</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left; border-top:1px solid;" colspan="3">1130</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">Marginal R<sup>2</sup> / Conditional R<sup>2</sup></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.052 / 0.117</td>
</tr>

</table>

```r
# GEE
calintake.gee <- geeglm(CALBEFORE ~ GROUP, id = CENTER, data = df, corstr = "exchangeable")
tab_model(calintake.gee) # same as GLMM
```

<table style="border-collapse:collapse; border:none;">
<tr>
<th style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm;  text-align:left; ">&nbsp;</th>
<th colspan="3" style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm; ">CALBEFORE</th>
</tr>
<tr>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  text-align:left; ">Predictors</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">Estimates</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">CI</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">p</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">(Intercept)</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1771.52</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1670.23&nbsp;&ndash;&nbsp;1872.82</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">GROUP2 × Fasting</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">&#45;284.14</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">&#45;428.88&nbsp;&ndash;&nbsp;-139.40</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">N <sub>CENTER</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">22</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm; border-top:1px solid;">Observations</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left; border-top:1px solid;" colspan="3">1130</td>
</tr>

</table>

```r
# GLMM, adjusted for BMI
calintake.glmm.bmi <- lmer(CALBEFORE ~ (1|CENTER) + GROUP + BMI, data = df)
# GEE, adjusted for BMI
calintake.gee.bmi <- geeglm(CALBEFORE ~ GROUP + BMI, id = CENTER, data = df, corstr = "exchangeable")

# ICC directly from GLMM model
ICC_BMI_unadj <- 158.9^2/(158.9^2+584.4^2) # but this is misleading, since conditioned on intervention!
ICC_BMI_adj <- 157.7^2/(157.7^2+579.8^2) # but this is misleading, since conditioned on intervention!
# => ICC directly from GLMM without conditioning on intervention
calintake.glmm.uncond <- lmer(CALBEFORE ~ (1|CENTER) + 1,
                   data = df)
ICC_BMI_unadj_uncond <- 210.4^2/(210.4^2+584.5^2) # but this does not provide 95%CI, and is not by group
# However, it is the same ICC as from one-way ANOVA overall
ICCest(x = CENTER, y = CALBEFORE, data = df, alpha = 0.05, CI.type = "THD") # ~ same as above, ICC ~ 0.11
```

```
## $ICC
## [1] 0.1164643
## 
## $LowerCI
## [1] 0.06486855
## 
## $UpperCI
## [1] 0.2252985
## 
## $N
## [1] 22
## 
## $k
## [1] 51.06574
## 
## $varw
## [1] 341767.7
## 
## $vara
## [1] 45050.49
```

```r
# One-way ANOVA ICC by group -> Gold standard to report ICC by group
ICCest(x = CENTER, y = CALBEFORE, data = df[df$GROUP == "1: Maintaining caloric intake",], alpha = 0.05, CI.type = "THD") # ICC control group
```

```
## $ICC
## [1] 0.08037775
## 
## $LowerCI
## [1] 0.03190476
## 
## $UpperCI
## [1] 0.2353509
## 
## $N
## [1] 11
## 
## $k
## [1] 55.70243
## 
## $varw
## [1] 323114.5
## 
## $vara
## [1] 28241.18
```

```r
ICCest(x = CENTER, y = CALBEFORE, data = df[df$GROUP == "2: Fasting",], alpha = 0.05, CI.type = "THD") # ICC intervention group
```

```
## $ICC
## [1] 0.05973571
## 
## $LowerCI
## [1] 0.01896594
## 
## $UpperCI
## [1] 0.1950354
## 
## $N
## [1] 11
## 
## $k
## [1] 45.82144
## 
## $varw
## [1] 364285.3
## 
## $vara
## [1] 23143.33
```

```r
# ICC from GEE model
calintake.gee.uncond <- geeglm(CALBEFORE ~ 1, id = CENTER, data = df, corstr = "exchangeable")
calintake.gee.uncond # See Estimated Correlation Parameters (alpha), the same: ICC ~ 0.11 // 95%CI ?
```

```
## 
## Call:
## geeglm(formula = CALBEFORE ~ 1, data = df, id = CENTER, corstr = "exchangeable")
## 
## Coefficients:
## (Intercept) 
##     1634.43 
## 
## Degrees of Freedom: 1130 Total (i.e. Null);  1129 Residual
## 
## Scale Link:                   identity
## Estimated Scale Parameters:  [1] 384218.7
## 
## Correlation:  Structure = exchangeable    Link = identity 
## Estimated Correlation Parameters:
##     alpha 
## 0.1076625 
## 
## Number of clusters:   22   Maximum cluster size: 70
```

```r
calintake.gee.uncond.cont <- geeglm(CALBEFORE ~ 1, id = CENTER, data = df[df$GROUP == "1: Maintaining caloric intake",], corstr = "exchangeable") # ICC control group
calintake.gee.uncond.int <- geeglm(CALBEFORE ~ 1, id = CENTER, data = df[df$GROUP == "2: Fasting",], corstr = "exchangeable") # ICC intervention group
```

## Binary outcome

```r
## binary outcome: Reintubation or death (the target trial outcome)
df <- df %>% # Create the variable outcome, which is equal to 1 for failure and 0 for success. => "Failure rate"
  mutate(outcome = case_when(INTUBATIONJ7 == 0 & DEATHJ7 == 0 ~ 0,
                             INTUBATIONJ7 == 1 | DEATHJ7 == 1 ~ 1))

# GLMM, the effect of the intervention on the failure rate
outcome.glmm <- glmer(outcome ~ (1|CENTER) + GROUP, data = df, family = "binomial")
tab_model(outcome.glmm)
```

<table style="border-collapse:collapse; border:none;">
<tr>
<th style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm;  text-align:left; ">&nbsp;</th>
<th colspan="3" style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm; ">outcome</th>
</tr>
<tr>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  text-align:left; ">Predictors</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">Odds Ratios</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">CI</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">p</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">(Intercept)</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.21</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.16&nbsp;&ndash;&nbsp;0.26</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">GROUP2 × Fasting</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.03</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.74&nbsp;&ndash;&nbsp;1.43</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.877</td>
</tr>
<tr>
<td colspan="4" style="font-weight:bold; text-align:left; padding-top:.8em;">Random Effects</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&sigma;<sup>2</sup></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">3.29</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&tau;<sub>00</sub> <sub>CENTER</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.02</td>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">ICC</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.01</td>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">N <sub>CENTER</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">22</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm; border-top:1px solid;">Observations</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left; border-top:1px solid;" colspan="3">1130</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">Marginal R<sup>2</sup> / Conditional R<sup>2</sup></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.000 / 0.006</td>
</tr>

</table>

```r
# GEE, the effect of the intervention on the failure rate
outcome.gee <- geeglm(outcome ~ GROUP, id = CENTER, data = df, corstr = "exchangeable", family = "binomial")
tab_model(outcome.gee) # same result as GLMM
```

<table style="border-collapse:collapse; border:none;">
<tr>
<th style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm;  text-align:left; ">&nbsp;</th>
<th colspan="3" style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm; ">outcome</th>
</tr>
<tr>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  text-align:left; ">Predictors</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">Odds Ratios</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">CI</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">p</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">(Intercept)</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.21</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.17&nbsp;&ndash;&nbsp;0.26</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">GROUP2 × Fasting</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.03</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.74&nbsp;&ndash;&nbsp;1.43</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.877</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">N <sub>CENTER</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">22</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm; border-top:1px solid;">Observations</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left; border-top:1px solid;" colspan="3">1130</td>
</tr>

</table>

```r
# GLMM, adjusted for BMI
outcome.glmm.bmi <- glmer(outcome ~ (1|CENTER) + GROUP + BMI, data = df, family = "binomial")
# GEE, adjusted for BMI
outcome.gee.bmi <- geeglm(outcome ~ GROUP + BMI, id = CENTER, data = df, corstr = "exchangeable", family = "binomial")

# Gold standard: Report one-way ANOVA ICC by group // GLMM
ICCest(x = CENTER, y = outcome, data = df[df$GROUP == "1: Maintaining caloric intake",], alpha = 0.05, CI.type = "THD")
```

```
## $ICC
## [1] 0.002250531
## 
## $LowerCI
## [1] -0.008255373
## 
## $UpperCI
## [1] 0.04262951
## 
## $N
## [1] 11
## 
## $k
## [1] 55.70243
## 
## $varw
## [1] 0.142225
## 
## $vara
## [1] 0.0003208038
```

```r
ICCest(x = CENTER, y = outcome, data = df[df$GROUP == "2: Fasting",], alpha = 0.05, CI.type = "THD")
```

```
## $ICC
## [1] 0.008011138
## 
## $LowerCI
## [1] -0.007462032
## 
## $UpperCI
## [1] 0.06602751
## 
## $N
## [1] 11
## 
## $k
## [1] 45.82144
## 
## $varw
## [1] 0.1439024
## 
## $vara
## [1] 0.001162132
```

```r
## overall ICC (but conditioned on the intervention!)
ICCest(x = CENTER, y = outcome, data = df, alpha = 0.05, CI.type = "THD")
```

```
## $ICC
## [1] 0.003726263
## 
## $LowerCI
## [1] -0.00591458
## 
## $UpperCI
## [1] 0.02745922
## 
## $N
## [1] 22
## 
## $k
## [1] 51.06574
## 
## $varw
## [1] 0.142985
## 
## $vara
## [1] 0.0005347924
```

```r
# Gold standard: Report one-way ANOVA ICC by group // GEE
calintake.gee.uncond.cont <- geeglm(outcome ~ 1, id = CENTER, data = df[df$GROUP == "1: Maintaining caloric intake",], corstr = "exchangeable", family = "binomial")
calintake.gee.uncond.int <- geeglm(outcome ~ 1, id = CENTER, data = df[df$GROUP == "2: Fasting",], corstr = "exchangeable", family = "binomial")
```

# SW-CRTs
Trial publication: https://pubmed.ncbi.nlm.nih.gov/30913216/

## Variable description
1. phc_code: Eighteen primary health centre (1 to 18), CLUSTER
2. PHASE: Time Period (1=6months, 2=12months, 3=18months, 4=24months) // "VERTICAL"
3. block: 3 sequences, each includes 6 PHCs, unit of randomization // "HORIZONTAL"
4. TRT:	Treatment allocation
5. primary_event:	primary, binary, outcome (SBP<140mmHg)
6. EQUK_change:	EQUK change from baseline to endline (a secondary, cont, outcome)

## Load data

```r
df <- read_delim("SMART.csv", delim = ";", 
    escape_double = FALSE, trim_ws = TRUE)
```

## Binary outcome (primary outcome)

```r
# reformat
df$primary_event_f <- as.factor(df$primary_event)
df$TRT_f <- as.factor(df$TRT)
df$PHASE_f <- as.factor(df$PHASE)

df$phc_code <- as.factor(df$phc_code) # always a factor
df$phc_code_modif <- as.factor(df$phc_code_modif) # only used for the SWplot

df <- df %>%
  mutate(primary_event_n = case_when(primary_event == "No" ~ 0,
                             primary_event == "Yes" ~ 1))
df <- df %>%
  mutate(TRT_n = case_when(TRT == "Control" ~ 0,
                             TRT == "Intervention" ~ 1))

# SW plot
swPlot(EQUK_change, TRT_n, PHASE_f, phc_code_modif, df, by.wave=FALSE,
       combined.plot=FALSE, 
       choose.tx.pos="bottomright",
       choose.legend.pos="bottom")

# table(df$PHASE,df$TRT)
# table(df$block,df$TRT) # block = sequence = randomised. Important: Different to a parallel CRT (or individual RCT) the randomized group variable is not used in the model

# GLMM
outcome.glmm <- glmer(primary_event_f ~ (1|phc_code) + TRT_f + PHASE_f, data = df, family = "binomial")
tab_model(outcome.glmm)
```

<table style="border-collapse:collapse; border:none;">
<tr>
<th style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm;  text-align:left; ">&nbsp;</th>
<th colspan="3" style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm; ">primary event f</th>
</tr>
<tr>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  text-align:left; ">Predictors</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">Odds Ratios</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">CI</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">p</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">(Intercept)</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.51</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.46&nbsp;&ndash;&nbsp;0.57</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">TRT f [Intervention]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.00</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.86&nbsp;&ndash;&nbsp;1.17</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.955</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">PHASE f [2]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.86</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.63&nbsp;&ndash;&nbsp;2.12</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">PHASE f [3]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.08</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.92&nbsp;&ndash;&nbsp;1.27</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.355</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">PHASE f [4]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.47</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.21&nbsp;&ndash;&nbsp;1.79</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>
<tr>
<td colspan="4" style="font-weight:bold; text-align:left; padding-top:.8em;">Random Effects</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&sigma;<sup>2</sup></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">3.29</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&tau;<sub>00</sub> <sub>phc_code</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.01</td>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">ICC</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.00</td>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">N <sub>phc_code</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">18</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm; border-top:1px solid;">Observations</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left; border-top:1px solid;" colspan="3">8642</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">Marginal R<sup>2</sup> / Conditional R<sup>2</sup></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.018 / 0.023</td>
</tr>

</table>

```r
# GEE - what the authors used in the publication // takes ~ 1min to converge
outcome.gee <- geeglm(primary_event_n ~ TRT_f + PHASE_f, id = phc_code, data = df, corstr = "exchangeable", family = "binomial")
tab_model(outcome.gee) # same as in publication
```

<table style="border-collapse:collapse; border:none;">
<tr>
<th style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm;  text-align:left; ">&nbsp;</th>
<th colspan="3" style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm; ">primary event n</th>
</tr>
<tr>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  text-align:left; ">Predictors</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">Odds Ratios</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">CI</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">p</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">(Intercept)</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.51</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.45&nbsp;&ndash;&nbsp;0.58</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">TRT f [Intervention]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.01</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.80&nbsp;&ndash;&nbsp;1.26</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.961</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">PHASE f [2]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.85</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.54&nbsp;&ndash;&nbsp;2.23</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">PHASE f [3]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.08</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.87&nbsp;&ndash;&nbsp;1.33</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.494</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">PHASE f [4]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.47</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.15&nbsp;&ndash;&nbsp;1.87</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>0.002</strong></td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">N <sub>phc_code</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">18</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm; border-top:1px solid;">Observations</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left; border-top:1px solid;" colspan="3">8642</td>
</tr>

</table>

```r
# Note: think about decaying correlation structure
```

![](stats_analyses_CRTs_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

# Parallel CRT with baseline period
Based on: https://www.bmj.com/content/360/bmj.k1121.long
Using data from PEBRA trial.
1. Trial publication: https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1004150
2. Trial code repo: https://github.com/alainamstutz/pebra 

## There are various ways to do it:
1. Analysis of covariance (ANCOVA): Aggregate outcomes at baseline, and adjusts each individual participant at follow-up for the baseline cluster mean
2. Constrained baseline analysis: Treat outcomes collected at baseline and follow-up as longitudinal, and to use a repeated measures analysis to estimate the effect of the intervention being switched on in one of the randomised groups on the second of these occasions, see design matrix in https://clusterrcts.shinyapps.io/rshinyapp/. Unlike a difference of differences analysis, it assumes that there is no systematic difference between the groups at baseline.

## Load data

```r
df <- readRDS("df_pebra.RData")
```

## Primary model, without baseline period

```r
# ITT model on primary endpoint (viral load), see publication
vs <- glmer(endpoint_reached ~ ARM + (1|USER) + DISTRICT + GENDER, data = df,
              family = "binomial")
```

```
## boundary (singular) fit: see help('isSingular')
```

```r
tab_model(vs)
```

<table style="border-collapse:collapse; border:none;">
<tr>
<th style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm;  text-align:left; ">&nbsp;</th>
<th colspan="3" style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm; ">endpoint reached</th>
</tr>
<tr>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  text-align:left; ">Predictors</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">Odds Ratios</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">CI</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">p</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">(Intercept)</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.93</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.26&nbsp;&ndash;&nbsp;2.97</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>0.003</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">ARM [interv.]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.27</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.79&nbsp;&ndash;&nbsp;2.03</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.327</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">DISTRICT [Leribe]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.71</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.34&nbsp;&ndash;&nbsp;1.49</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.370</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">DISTRICT [MKG]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.59</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.36&nbsp;&ndash;&nbsp;0.98</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>0.040</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">GENDER [male]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.11</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.66&nbsp;&ndash;&nbsp;1.88</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.688</td>
</tr>
<tr>
<td colspan="4" style="font-weight:bold; text-align:left; padding-top:.8em;">Random Effects</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&sigma;<sup>2</sup></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">3.29</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&tau;<sub>00</sub> <sub>USER</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.00</td>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">N <sub>USER</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">20</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm; border-top:1px solid;">Observations</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left; border-top:1px solid;" colspan="3">307</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">Marginal R<sup>2</sup> / Conditional R<sup>2</sup></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.023 / NA</td>
</tr>

</table>

## With baseline period: Analysis of covariance

```r
# Calculate the mean cluster value of the baseline viral load
df$VL_RESULT_baseline <- as.numeric(df$VL_RESULT_baseline)
df <- df %>% # there are several baseline VL variables, to force <20 into 0 is not ideal, but best we can do.
  mutate(baseline_Vl_num = case_when(baseline_Vl_cat == "<20" ~ 0,
                                     baseline_Vl_cat == ">999" ~ VL_RESULT_baseline,
                                     baseline_Vl_cat == "20-999" ~ VL_RESULT_baseline))
cluster_mean <- df %>%
  group_by(USER) %>%
  drop_na(baseline_Vl_num) %>% 
  summarize(baseline_Vl_meanUSER = mean(baseline_Vl_num))
df <- left_join(df, cluster_mean[, c("baseline_Vl_meanUSER", "USER")], by = join_by(USER == USER))

# individual-level ANCOVA with cluster-level adjustment
vs.ancova <- glmer(endpoint_reached ~ DISTRICT + ARM + GENDER + (1|USER) + baseline_Vl_meanUSER, 
              data = df, family = "binomial")
# tab_model(vs.ancova)
# not ideal, due to cluster level aggregation of viral load categories
```

## With baseline period: Constrained baseline analysis
Two possible ways: The first approach assumes that the correlation between two people from the same cluster is the same whether they are sampled in the same period or a different period. The second approach allows the correlation to be weaker between different periods. The method is extremely flexible, is available in cohort or repeated cross section forms, and allows an analysis based on individual level data, with no aggregation needed either at baseline or at follow-up.

```r
# First, reshape dataset to mirror the design
# Duplicate the dataset
df_dup <- rep(list(df), times = 2)
df_dup <- do.call(rbind, df_dup)
# Create a new variable "time" and assign values 0 and 1 to each of the created clones, corresponding to 0=baseline and 1=follow-up
df_dup$time <- rep(0:1, each = nrow(df_dup) / 2)
# Add the baseline VL to the baseline clone, using the same definition as for the outcome
df_dup <- df_dup %>% 
  mutate(baseline_endpoint_reached = case_when(baseline_Vl_cat == "<20" ~ 1,
                                               baseline_Vl_cat == "20-999" | baseline_Vl_cat == ">999" ~ 0))
df_dup <- df_dup %>% 
  mutate(endpoint_reached = case_when(time == 0 ~ baseline_endpoint_reached,
                           TRUE ~ endpoint_reached))
# Create the treatment variable by period (=time)
df_dup <- df_dup %>% 
  mutate(treat = case_when(ARM == "interv." & time == 1 ~ 1,
                           TRUE ~ 0))
# df_dup %>%
#   select(time, USER, IND_ID, ARM, treat, endpoint_reached, baseline_endpoint_reached) %>%
#   View()

# # just for simplicity: take out NA, assign not reached endpoint
# df_dup <- df_dup %>% 
#   mutate(endpoint_reached = case_when(is.na(endpoint_reached) ~ 0,
#                            TRUE ~ endpoint_reached))

# Approach 1: constrained baseline analysis – inflexible correlation structure, assuming a random effect of cluster and a random effect of individual nested within cluster, but no random effect of time nested within cluster (this fits a model where the cluster autocorrelation is assumed to be 1).
vs.constrained.inflex <- glmer(endpoint_reached ~ time + treat + (1|USER) + DISTRICT + GENDER, data = df_dup, family = "binomial")

# Approach 2: constrained baseline analysis – flexible correlation structure, assume instead a random effect of cluster, a random effect of individual nested within cluster, and a random effect of time nested within cluster (this allows the cluster autocorrelation to be less than 1)
vs.constrained.flex <- glmer(endpoint_reached ~ time + treat + (1|USER) + (1 + time | USER) + DISTRICT + GENDER, data = df_dup, family = binomial, control = glmerControl(optimizer = "bobyqa"))

# compare the three models, a) primary, b) constrained baseline inflexible, c) constrained baseline flexible
tab_model(vs) # primary model
```

<table style="border-collapse:collapse; border:none;">
<tr>
<th style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm;  text-align:left; ">&nbsp;</th>
<th colspan="3" style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm; ">endpoint reached</th>
</tr>
<tr>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  text-align:left; ">Predictors</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">Odds Ratios</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">CI</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">p</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">(Intercept)</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.93</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.26&nbsp;&ndash;&nbsp;2.97</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>0.003</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">ARM [interv.]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.27</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.79&nbsp;&ndash;&nbsp;2.03</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.327</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">DISTRICT [Leribe]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.71</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.34&nbsp;&ndash;&nbsp;1.49</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.370</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">DISTRICT [MKG]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.59</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.36&nbsp;&ndash;&nbsp;0.98</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>0.040</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">GENDER [male]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.11</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.66&nbsp;&ndash;&nbsp;1.88</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.688</td>
</tr>
<tr>
<td colspan="4" style="font-weight:bold; text-align:left; padding-top:.8em;">Random Effects</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&sigma;<sup>2</sup></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">3.29</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&tau;<sub>00</sub> <sub>USER</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.00</td>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">N <sub>USER</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">20</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm; border-top:1px solid;">Observations</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left; border-top:1px solid;" colspan="3">307</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">Marginal R<sup>2</sup> / Conditional R<sup>2</sup></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.023 / NA</td>
</tr>

</table>

```r
tab_model(vs.constrained.inflex) # As with a difference of differences analysis, the treatment effect is the regression coefficient for treat
```

<table style="border-collapse:collapse; border:none;">
<tr>
<th style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm;  text-align:left; ">&nbsp;</th>
<th colspan="3" style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm; ">endpoint reached</th>
</tr>
<tr>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  text-align:left; ">Predictors</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">Odds Ratios</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">CI</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">p</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">(Intercept)</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">2.20</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.58&nbsp;&ndash;&nbsp;3.06</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">time</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.96</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.64&nbsp;&ndash;&nbsp;1.45</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.861</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">treat</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.32</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.82&nbsp;&ndash;&nbsp;2.11</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.256</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">DISTRICT [Leribe]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.04</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.60&nbsp;&ndash;&nbsp;1.82</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.884</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">DISTRICT [MKG]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.52</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.36&nbsp;&ndash;&nbsp;0.75</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">GENDER [male]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.78</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.54&nbsp;&ndash;&nbsp;1.14</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.197</td>
</tr>
<tr>
<td colspan="4" style="font-weight:bold; text-align:left; padding-top:.8em;">Random Effects</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&sigma;<sup>2</sup></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">3.29</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&tau;<sub>00</sub> <sub>USER</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.00</td>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">N <sub>USER</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">20</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm; border-top:1px solid;">Observations</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left; border-top:1px solid;" colspan="3">580</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">Marginal R<sup>2</sup> / Conditional R<sup>2</sup></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.038 / NA</td>
</tr>

</table>

```r
tab_model(vs.constrained.flex) # As with a difference of differences analysis, the treatment effect is the regression coefficient for treat
```

<table style="border-collapse:collapse; border:none;">
<tr>
<th style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm;  text-align:left; ">&nbsp;</th>
<th colspan="3" style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm; ">endpoint reached</th>
</tr>
<tr>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  text-align:left; ">Predictors</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">Odds Ratios</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">CI</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">p</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">(Intercept)</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">2.20</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.58&nbsp;&ndash;&nbsp;3.06</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">time</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.96</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.64&nbsp;&ndash;&nbsp;1.45</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.861</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">treat</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.32</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.82&nbsp;&ndash;&nbsp;2.11</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.256</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">DISTRICT [Leribe]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">1.04</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.60&nbsp;&ndash;&nbsp;1.82</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.884</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">DISTRICT [MKG]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.52</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.36&nbsp;&ndash;&nbsp;0.75</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  "><strong>&lt;0.001</strong></td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">GENDER [male]</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.78</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.54&nbsp;&ndash;&nbsp;1.14</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">0.197</td>
</tr>
<tr>
<td colspan="4" style="font-weight:bold; text-align:left; padding-top:.8em;">Random Effects</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&sigma;<sup>2</sup></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">3.29</td>
</tr>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&tau;<sub>00</sub> <sub>USER</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.00</td>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&tau;<sub>00</sub> <sub>USER.1</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.00</td>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&tau;<sub>11</sub> <sub>USER.1.time</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.00</td>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">&rho;<sub>01</sub> <sub>USER.1</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">&nbsp;</td>

<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">N <sub>USER</sub></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">20</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm; border-top:1px solid;">Observations</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left; border-top:1px solid;" colspan="3">580</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">Marginal R<sup>2</sup> / Conditional R<sup>2</sup></td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="3">0.038 / NA</td>
</tr>

</table>

```r
summary(vs.constrained.flex)
```

```
## Generalized linear mixed model fit by maximum likelihood (Laplace
##   Approximation) [glmerMod]
##  Family: binomial  ( logit )
## Formula: endpoint_reached ~ time + treat + (1 | USER) + (1 + time | USER) +  
##     DISTRICT + GENDER
##    Data: df_dup
## Control: glmerControl(optimizer = "bobyqa")
## 
##      AIC      BIC   logLik deviance df.resid 
##    772.4    816.0   -376.2    752.4      570 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.7047 -1.0683  0.6730  0.7624  1.0774 
## 
## Random effects:
##  Groups Name        Variance  Std.Dev.  Corr
##  USER   (Intercept) 0.000e+00 0.000e+00     
##  USER.1 (Intercept) 0.000e+00 0.000e+00     
##         time        1.389e-15 3.727e-08  NaN
## Number of obs: 580, groups:  USER, 20
## 
## Fixed effects:
##                Estimate Std. Error z value Pr(>|z|)    
## (Intercept)     0.78726    0.16951   4.644 3.41e-06 ***
## time           -0.03641    0.20830  -0.175 0.861232    
## treat           0.27466    0.24158   1.137 0.255567    
## DISTRICTLeribe  0.04133    0.28386   0.146 0.884241    
## DISTRICTMKG    -0.65514    0.18522  -3.537 0.000405 ***
## GENDERmale     -0.24482    0.18983  -1.290 0.197164    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) time   treat  DISTRICTL DISTRICTM
## time        -0.474                                  
## treat        0.038 -0.550                           
## DISTRICTLrb -0.348  0.004  0.000                    
## DISTRICTMKG -0.547  0.015 -0.021  0.326             
## GENDERmale  -0.336  0.051 -0.081 -0.030    -0.002   
## optimizer (bobyqa) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

# Sample size calculation for parallel CRT with binary outcome (AND baseline period AND closed cohort design)
Following the guidance of: 
* Hooper Richard et al. Sample size calculation for stepped wedge and other longitudinal cluster randomised trials. Statistics in Medicine. 2016: https://onlinelibrary.wiley.com/doi/10.1002/sim.7028 and 
* Leyrat Clémence et al. Practical considerations for sample size calculation for cluster randomized trials. Journal of Epidemiology and Population Health. 2024: https://www.sciencedirect.com/science/article/pii/S2950433324000090

On the example of the Hair SALON hybrid implementation-effectiveness parallel CRT in Lesotho.

PICO:
* P: Adolescent girls and young women aged 15-30
* I: HIV/SRH package with at least oral PrEP (plus HIVST, oral contraception, menstrual health products, linkage services, etc.)
* C: ”SOC” at hair salon (info/flyer and referral)
* O: Feasibility indicators, but explore clinical effectiveness on PrEP uptake at 6 months, for which we power our sample size calculation

Design features:
* Cluster randomized (hair salons)
* 1:1 randomization, stratified by district
* Binary outcome (self-reported uptake within past week at 6 months)
* Closed cohort
* With baseline period. To increase efficiency, but also avoid recruitment bias, i.e. enrollment before allocation.

Requirements:
* power: 90%
* alpha: 5%
* as few clusters as possible, as few participants as possible

Assumptions:
* Baseline uptake: 20% (see survey data)
* Increase in uptake through intervention: 10%
* ICC: 0.15 (this is rather conservative, according to guidance from CDC guidance re ICC for PrEP uptake as outcome [https://www.cdc.gov/hiv/pdf/research/interventionresearch/compendium/prep/PrEP_Chapter_EBI_Criteria.pdf] and see data from SEARCH [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9169331/], and realistic for a behavioral intervention)
* CAC: 1 (since closed cohort, i.e. same individuals)
* IAC: 0.7 (this is rather conservative; often 0.8 taken; to be explored from longitudinal data in the pilot)


```r
# Define the num. parameters
alpha <- 0.05 # alpha level
power <- 0.90 # power
p_cont <- 0.20 # Proportion of outcome in the control group
p_int <- 0.30 # Proportion of outcome in the intervention group

cluster_size <- 20  # Average cluster size, i.e. participants per cluster

icc <- 0.15  # Intra-cluster correlation coefficient; "how similar outcomes are within the same cluster"
cac <- 1 # cluster auto-correlation; "how similar outcomes are within the same cluster over time"; a cluster autocorrelation less than 1 reflects a situation where individuals sampled from the same cluster at different times have less correlated outcomes than individuals sampled from the same cluster at the same time; however, if closed cohort then CAC = 1, because same individuals.
iac <- 0.7 # individual auto-correlation; correlation between repeated measurements taken from the same individual over time; r=1: Perfect positive autocorrelation — successive measurements are identical; r=0: No autocorrelation — successive measurements are independent. (of course, during usual care, without the intervention)

deff_c <- 1 + (cluster_size - 1) * icc # design effect due to cluster randomization
rcc <- ((cluster_size*icc*cac) + ((1-icc)*iac)) / (1 + ((cluster_size-1)*icc)) # the correlation between Yt1kl and Yt2kl for any k, l, t1 ≠ t2 under model (see formula 9 in Hooper et al): "The only other thing that matters to the variance of the treatment effect estimator is the correlation, r"; the correlation between outcomes at different times within the same cluster (r) is the crucial factor that affects the variance of the treatment effect estimator.
deff_r <- (1-rcc^2) # design effect due to repeated assessment


# First, calculate the sample size for an individual RCT, using pwr, two-sided (effect could go either way)
ss_ind_1arm <- pwr.2p.test(h = ES.h(p_int, p_cont), 
                           sig.level = alpha, 
                           power = power)
ss_ind <- ss_ind_1arm$n * 2
# cat("Total sample size individual RCT, pwr:", ss_ind)
# cat("Total sample size individual RCT, pwr:", round(ss_ind, 0))

## As a comparison, use "manual" formula instead of pwr package
Z_alpha_half <- qnorm(1 - alpha / 2) # translate into Z-distribution -> equals 0.975 (95% CI) 
Z_beta <- qnorm(power)
ss_ind_1arm_man <- ((Z_alpha_half + Z_beta)^2 * (p_int * (1 - p_int) + p_cont * (1 - p_cont))) / (p_int - p_cont)^2
ss_ind_man <- ss_ind_1arm_man * 2
# cat("Total sample size individual RCT, manual:", ss_ind_man) # total sample size
# => same result, use pwr result going forward


# Second, inflate for clustering, following a standard parallel 1:1 CRT without repeated measures/baseline period
ss_crt_standard <- ss_ind * deff_c
cat("Total sample size CRT, standard:", ss_crt_standard) # total sample size for a standard CRT
```

```
## Total sample size CRT, standard: 3006.767
```

```r
n_clus_crt_standard <- ss_crt_standard / cluster_size
cat("Total N clusters CRT, standard:", n_clus_crt_standard) # total number of clusters for a standard CRT (divide by arm or sequence/steps if SWCRT)
```

```
## Total N clusters CRT, standard: 150.3383
```

```r
# Third, add the baseline period and assume a closed cohort, i.e., 1 repeated measure among the same individuals
ss_crt_b_period <- deff_r*deff_c*ss_ind 
cat("Total sample size CRT, baseline period:", ss_crt_b_period) # total sample size for a CRT with baseline period design
```

```
## Total sample size CRT, baseline period: 385.1086
```

```r
n_clus_crt_b_period <- ss_crt_b_period / cluster_size
cat("Total N clusters CRT, baseline period:", n_clus_crt_b_period) # total number of clusters for a CRT with baseline period design (divide by arm or sequence/steps if SWCRT)
```

```
## Total N clusters CRT, baseline period: 19.25543
```
The exact same results were obtained when using the online sample size calculator for CRTs, developed by Hemming et al ! ->  https://clusterrcts.shinyapps.io/rshinyapp/ (assuming an exchangeable correlation structure, i.e. CAC = 1)


# Simulation of datasets and sample sizes for a parallel CRT with binary outcome AND baseline period AND closed cohort
Following the code provided in the supplement of (adapted for binary outcome and R): 
* Hooper Richard et al. Sample size calculation for stepped wedge and other longitudinal cluster randomised trials. Statistics in Medicine. 2016: https://onlinelibrary.wiley.com/doi/10.1002/sim.7028

```r
s_cohortstep_binary <- function(p_cont, p_int, icc, cac, iac, ncluspergrp, cluster_size, nstep) {
  
  # Set seed for reproducibility
  set.seed(102030)
  
  # Calculate log odds for control and intervention
  odds_control <- p_cont / (1 - p_cont)
  odds_intervention <- p_int / (1 - p_int)
  
  # Calculate treatcoeff as the log of the odds ratio
  treatcoeff <- log(odds_intervention / odds_control)
  
  # Number of clusters
  nclus <- nstep * ncluspergrp
  
  # Generate cluster-level random effects
  idclus <- rep(1:nclus, each = cluster_size * (nstep + 1))  # The vector contains integers from 1 to nclus; each integer represents a unique cluster ID.
  group <- rep(rep(1:nstep, each = cluster_size * (nstep + 1)), ncluspergrp) # Assign each cluster to a specific group based on the step in the stepped-wedge design. nstep is the number of steps or time periods. ncluspergrp is the number of clusters per group. 
  rand_clus <- rnorm(nclus, 0, sqrt(icc * cac)) # generates random values drawn from a normal distribution with a mean of 0 and a standard deviation of sqrt(icc * cac). The product icc * cac provides a variance measure that represents the contribution of cluster-level random effects to the overall variability of the outcome (using the square root to get standard deviations), i.e. the unobserved differences between clusters that might influence the outcome.
  
  # Generate time-level random effects
  rand_time <- matrix(rnorm(nclus * (nstep + 1), 0, sqrt(icc * (1 - cac))), nrow = nclus) # generates random effects that vary across both clusters and time periods. nstep + 1 is the total number of time periods (steps) in the study. nclus * (nstep + 1) represents the total number of cluster-time combinations, i.e., the total number of random effects to generate. icc * (1 - cac): portion of the within-cluster variability that is not explained by the cluster autocorrelation, i.e., the variance of the time-level random effects within clusters. => variability not only between clusters but also over time within each cluster! The outcome for a cluster may vary across different time points due to factors not accounted for by the fixed effects (like the intervention or the passage of time), influenced by unmeasured factors.
  
  # Expand to individual level
  id <- rep(1:(nclus * cluster_size), each = ncluspergrp) # nclus * cluster_size = total number of individuals in the entire study.
  rand_char <- rnorm(length(id), 0, sqrt(iac * (1 - icc))) # generates individual-level random effects for each participant in the study. The product iac * (1 - icc) gives the variance of the individual-level random effects, capturing the variability at the individual level after accounting for the intracluster correlation.
  
  # Create time variable and treatment indicator
  time <- rep(1:(nstep + 1), times = cluster_size * ncluspergrp)
  treat <- as.numeric(time >= group) # It assigns a value of 1 if the individual is in the treatment group (i.e., their time step is equal to or greater than their cluster's crossover step) and 0 if they are still in the control group.
  
  # Generate individual-level error term
  rand_err <- rnorm(length(id), 0, sqrt((1 - iac) * (1 - icc))) # The rand_err variable represents the random error at the individual level. This error accounts for variability in the outcome that is not explained by the cluster-level effects, time-level effects, or individual-level characteristics. The unstructured random noise at the individual level. While rand_char represents structured random variability associated with individual-level characteristics (e.g. baseline severity).
  
  # Generate the linear predictor (log-odds)
  # The linear predictor (the log-odds) is the input to this logistic function, determining the probability of success (e.g., PrEP uptake) for each individual based on the treatment and other factors. The linear predictor is the sum of all these terms and represents the combined effect of the treatment, cluster, time, individual characteristics, and random noise on the log-odds of the binary outcome. In logistic regression, the probability p of the binary outcome occurring is computed using the logistic function (see formula below, next line).
  linear_predictor <- treatcoeff * treat + 
                   rep(rnorm(nclus, 0, sqrt(icc * cac)), each = cluster_size * (nstep + 1)) + 
                   as.vector(rnorm(nclus * (nstep + 1), 0, sqrt(icc * (1 - cac)))) + 
                   rnorm(cluster_size * nclus * (nstep + 1), 0, sqrt(iac * (1 - icc))) + 
                   rnorm(cluster_size * nclus * (nstep + 1), 0, sqrt((1 - iac) * (1 - icc)))
  
  # Convert linear predictor to probability using logistic function
  prob <- 1 / (1 + exp(-linear_predictor))
  
  # Generate binary outcome based on the probability
  y <- rbinom(length(prob), size = 1, prob = prob) # For each individual, rbinom() performs a Bernoulli trial (a random experiment with two possible outcomes) with the specified probability of success (from the prob vector), e.g.the rbinom() function will generate a 1 (indicating success) with a 20% probability and a 0 (indicating failure) with a 80% probability.
  
  # Create a data frame
  data <- data.frame(id = 1:length(y), idclus, time, treat, y)
  
  # Fit the GLMM model (logistic regression with random effects)
  library(lme4)
  model <- glmer(y ~ factor(time) + treat + (1 | idclus) + (1 + time | idclus), 
                 data = data, family = binomial, control = glmerControl(optimizer = "bobyqa"))
  
  # Calculate p-value for the treatment effect
  treat_effect <- summary(model)$coefficients["treat", "Estimate"]
  treat_se <- summary(model)$coefficients["treat", "Std. Error"]
  p_value <- 2 * pnorm(-abs(treat_effect / treat_se))
  
  return(list(model = model, p_value = p_value))
}


### 
simulate_power_binary <- function(p_cont, p_int, icc, cac, iac, cluster_size, nstep, p, start, inc) {
  power <- 0
  ncluspergrp <- start
  
  while (power < p) { # The loop continues until the simulated power meets or exceeds the desired power level p.
    results <- replicate(1000, {
      sim_result <- s_cohortstep_binary(p_cont, p_int, icc, cac, iac, ncluspergrp, cluster_size, nstep)
      sim_result$p_value < 0.05 # if below 0.05, then TRUE
    })
    
    power <- mean(results) #proportion of simulations in which the treatment effect was significant (i.e., the proportion of TRUE values in results)
    
    if (power >= p) { # if power reached, the loop breaks, and the function stops increasing the number of clusters per group.
      break
    }
    
    ncluspergrp <- ncluspergrp + inc
  }
  
  return(ncluspergrp)
}



### Example usage
set.seed(102030)

p_value_simulation <- s_cohortstep_binary(
  p_cont = 0.20,                  # Proportion in control group
  p_int = 0.35,                   # Proportion in intervention group
  icc = 0.15,                     # Intracluster correlation coefficient
  cac = 1,                        # Cluster autocorrelation
  iac = 0.8,                      # Individual autocorrelation
  ncluspergrp = 3,                # clusters per group
  cluster_size = 16,              # Number of individuals per cluster
  nstep = 2                       # Number of steps in the stepped-wedge design
)

print(p_value_simulation)
```

```
## $model
## Generalized linear mixed model fit by maximum likelihood (Laplace
##   Approximation) [glmerMod]
##  Family: binomial  ( logit )
## Formula: y ~ factor(time) + treat + (1 | idclus) + (1 + time | idclus)
##    Data: data
##       AIC       BIC    logLik  deviance  df.resid 
##  385.2973  414.6010 -184.6487  369.2973       280 
## Random effects:
##  Groups   Name        Std.Dev. Corr 
##  idclus   (Intercept) 0.0000        
##  idclus.1 (Intercept) 0.7153        
##           time        0.1768   -1.00
## Number of obs: 288, groups:  idclus, 6
## Fixed Effects:
##   (Intercept)  factor(time)2  factor(time)3          treat  
##       0.25085        0.08694        0.16666        0.32891  
## optimizer (bobyqa) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 
## 
## $p_value
## [1] 0.5623465
```

```r
# nclus_required_binary <- simulate_power_binary(
#   p_cont = 0.20,                  # Proportion in control group
#   p_int = 0.35,                   # Proportion in intervention group
#   icc = 0.15,                     # Intracluster correlation coefficient
#   cac = 1,                        # Cluster autocorrelation
#   iac = 0.8,                      # Individual autocorrelation
#   cluster_size = 16,              # Number of individuals per cluster
#   nstep = 2,                      # Number of steps in the stepped-wedge design
#   p = 0.1,                        # Desired power level
#   start = 2,                      # Initial number of clusters per group
#   inc = 1                         # Increment for increasing clusters per group
# )
```

# Sample size calculation for parallel CRT with binary outcome (AND baseline period AND closed cohort design) AND multiarm.

Example: EHR nudges for statin prescription CRT

PICO:
* P: People living with HIV with low-moderate-high CVD risk score and not on a statin
* I: Statin EHR nudges on provider and participant level
* C: SOC: No Statin EHR nudges
* O: Statin prescription rate

Design features:
* Cluster randomized (SHCS sites)
* 1:1:1 randomization (stratified by mean statin prescription rate)
* Binary outcome (statin prescription rate after implementation)
* Closed cohort
* With baseline period. To increase efficiency.
* 3 arms: 1) SOC, 2) provider nudges, 3) provider+participant nudges

Requirements:
* power: 90%
* To avoid adjustment for multiple testing, share alpha for the two comparisons
* as few clusters as possible, as few participants as possible

Assumptions:
* Baseline rate: 25% (see SHCS data)
* Increase through intervention: 10%
* ICC: 0.15 (conservative; behavioral)
* CAC: 1 (since closed cohort, i.e. same individuals)
* IAC: 0.7 (behavioral)


```r
# Define parameters
alpha <- 0.025  # Adjusted alpha for multiple comparisons (split from 0.05)
power <- 0.90  # 90% power
p_cont <- 0.25 # Proportion of outcome in the control group
p_int <- 0.35 # Proportion of outcome in the intervention group

cluster_size <- 33  # Average cluster size
icc <- 0.15  # Intra-cluster correlation coefficient
cac <- 1  # Cluster autocorrelation
iac <- 0.7  # Individual autocorrelation (correlation between measurements within individuals)

# First, calculate the sample size for an individual RCT, using pwr, two-sided (effect could go either way)
ss_ind_1arm <- pwr.2p.test(h = ES.h(p_int, p_cont), 
                           sig.level = alpha, 
                           power = power)
ss_ind <- ss_ind_1arm$n * 3
# cat("Total sample size individual RCT, pwr:", ss_ind)

# Second, inflate for clustering, following a standard parallel 1:1 CRT without repeated measures/baseline period
ss_crt_standard <- ss_ind * deff_c
cat("Total sample size CRT, standard:", ss_crt_standard) # total sample size for a standard CRT for 2 comparisons (shared control)
```

```
## Total sample size CRT, standard: 5982.878
```

```r
n_clus_crt_standard <- ss_crt_standard / cluster_size
cat("Total N clusters CRT, standard:", n_clus_crt_standard) # total number of clusters for a standard CRT for 2 comparisons (shared control)
```

```
## Total N clusters CRT, standard: 181.2993
```

```r
# Third, add the baseline period and assume a closed cohort, i.e., 1 repeated measure among the same individuals
ss_crt_b_period <- deff_r*deff_c*ss_ind 
cat("Total sample size CRT, baseline period:", ss_crt_b_period) # total sample size for a CRT with baseline period design, for 2 comparisons (shared control)
```

```
## Total sample size CRT, baseline period: 766.2908
```

```r
n_clus_crt_b_period <- ss_crt_b_period / cluster_size
cat("Total N clusters CRT, baseline period:", n_clus_crt_b_period) # total number of clusters for a CRT with baseline period design, for 2 comparisons (shared control)
```

```
## Total N clusters CRT, baseline period: 23.22093
```

```r
### since I divided my alpha, I can make two comparisons (Intervention 1 vs Usual care & Intervention 2 vs Usual care)
```


# Sample size calculation for parallel CRT with cont outcome (AND baseline period AND closed cohort design) AND multiarm.

Example: EHR nudges for statin prescription CRT

PICO:
* P: People living with HIV with low-moderate-high CVD risk score and not on a statin
* I: Statin EHR nudges on provider and participant level
* C: SOC: No Statin EHR nudges
* O: nonHDL reduction

Design features:
* Cluster randomized (SHCS sites)
* 1:1:1 randomization (stratified by mean statin prescription rate / mean nonHDL)
* Continuous outcome (15mg/dL nonHDL reduction at 6 months, see guidance)
* Closed cohort
* With baseline period. To increase efficiency.
* 3 arms: 1) SOC, 2) provider nudges, 3) provider+participant nudges

Requirements:
* power: 90%
* To avoid adjustment for multiple testing, share alpha for the two comparisons
* as few clusters as possible, as few participants as possible

Assumptions:
* SD nonHDL: 44mg/dL (see guidance)
* Clinically meaningful difference nonHDL: 15mg/dL
* ICC: 0.15 (conservative)
* CAC: 1 (since closed cohort, i.e. same individuals)
* IAC: 0.6 (lower than above since less behavioural)


```r
# Define parameters
alpha <- 0.025  # Adjusted alpha for multiple comparisons (split from 0.05)
power <- 0.90  # 90% power
mean_diff <- 15  # Detecting at least a 15 mg/dL reduction in nonHDL cholesterol
sd <- 44  # Standard deviation of the outcome (non-HDL cholesterol)

cluster_size <- 13  # Average cluster size
icc <- 0.15  # Intra-cluster correlation coefficient
cac <- 1  # Cluster autocorrelation
iac <- 0.6  # Individual autocorrelation (correlation between measurements within individuals)

# Calculate the design effect due to clustering
deff_c <- 1 + (cluster_size - 1) * icc  # Design effect due to clustering

# Calculate the correlation between outcomes at different time points within the same cluster
rcc <- ((cluster_size * icc * cac) + ((1 - icc) * iac)) / (1 + ((cluster_size - 1) * icc))

# Calculate the design effect due to repeated measures
deff_r <- (1 - rcc^2)  # Design effect due to repeated measures

# First, calculate the sample size for an individual RCT with continuous outcome
# Alpha is set to 0.025 to account for multiple comparisons
ss_ind_1arm <- pwr.t.test(d = mean_diff / sd, 
                          sig.level = alpha, 
                          power = power, 
                          type = "two.sample", 
                          alternative = "two.sided")
ss_ind <- ss_ind_1arm$n * 2  # Total sample size for 2 comparisons (shared control)
# cat("Total sample size individual RCT (per arm):", ss_ind, "\n")

# Second, inflate for clustering, following a standard parallel 1:1 CRT without repeated measures/baseline period
ss_crt_standard <- ss_ind * deff_c
cat("Total sample size CRT, standard:", ss_crt_standard) # total sample size for a standard CRT, for 2 comparisons (shared control)
```

```
## Total sample size CRT, standard: 1203.13
```

```r
n_clus_crt_standard <- ss_crt_standard / cluster_size
cat("Total N clusters CRT, standard:", n_clus_crt_standard) # total number of clusters for a standard CRT, for 2 comparisons
```

```
## Total N clusters CRT, standard: 92.54844
```

```r
# Third, add the baseline period and assume a closed cohort, i.e., 1 repeated measure among the same individuals
ss_crt_b_period <- deff_r*deff_c*ss_ind 
cat("Total sample size CRT, baseline period:", ss_crt_b_period) # total sample size for a CRT with baseline period design, for 2 comparisons (shared control)
```

```
## Total sample size CRT, baseline period: 274.4486
```

```r
n_clus_crt_b_period <- ss_crt_b_period / cluster_size
cat("Total N clusters CRT, baseline period:", n_clus_crt_b_period) # total number of clusters for a CRT with baseline period design, for 2 comparisons (shared control)
```

```
## Total N clusters CRT, baseline period: 21.11143
```

```r
### since I divided my alpha, I can make two comparisons (Intervention 1 vs Usual care & Intervention 2 vs Usual care)
```

# Sample Size calculation for a pilot and feasibilty study, according to "Guidelines for Designing and Evaluating Feasibility Pilot Studies"
(https://pmc.ncbi.nlm.nih.gov/articles/PMC8849521/)
1) The aim of the sample size calculation is to calculate the CI around the single proportion of PrEP uptake within 6 months in our pilot, to inform the larger cluster randomized trial
2) We take into account the clustering of the data, i.e., clusters = stylists with several clients
3) Target population: Young women attending a hair salon in Lesotho, sexually active, HIV-
4) Estimated demand of PrEP: Based on the data from our citizen scientist survey -> 22.1% were interested in taking up PrEP when offered
5) 95% confidence: The probability that our sample results will be within the margin of error of the true population estimate.
6) Margin of error (precision): 10%. The maximum acceptable difference between our sample estimate and the true population parameter. If we find that 22.1% of our participating clients take up PrEP, we can assume that 12.1-32.1% of our target population will do so.
7) Estimate the intracluster correlation coefficient (ICC): The ICC measures the similarity or correlation between responses within the same cluster. A higher ICC (closer to 1) implies more similarity (i.e. more correlation). Based on the ICC, the design effect (or inflation factor), is calculated. We calculated the ICC from our citizen scientist survey data -> 0.06543012
8) Due to operational/budget and piloting reasons: as few clusters as possible


```r
# Define parameters
confidence_level <- 0.95 # confidence level
z <- qnorm(1 - (1 - confidence_level) / 2) # z statistic
d <- 0.10 # margin of error
cluster_size <- 15 # Average cluster size (fixed in our case)
icc <- 0.06543012 # Intracluster correlation coefficient -> from the survey
prop <- 0.221 # Estimated proportion (= 22.1% demand) -> from the survey
deff <- 1 + (cluster_size - 1) * icc # Design Effect
sample_size <- ceiling((z^2 * prop * (1 - prop) * deff) / d^2)
# Number of individual participants
cat("Required sample size:", sample_size)
```

```
## Required sample size: 127
```

```r
# Number of clusters
n_clusters <- sample_size / cluster_size
cat("Required clusters:", round(n_clusters, 0))
```

```
## Required clusters: 8
```

# R21/R33: Sample Size calculation for R21/R33 Vending machines proposal
## according to "Guidelines for Designing and Evaluating Feasibility Pilot Studies"
(https://pmc.ncbi.nlm.nih.gov/articles/PMC8849521/)
1) The aim of the sample size calculation is to calculate the CI around the single proportion of participants at high HIV acquisition risk accessing HIVST at a vending machine and link to the clinic within 3 months
2) We take into account the clustering of the data, i.e., clusters = vending machines with several participants
3) Target population: Participants accessing one of the vending machines at high HIV acquisition risk (epi 3+ risk score)
4) Based on previous data, we know there are ca. 20'000 campus students, 12'000 sexually active and 8000 of these with a high epi risk score. We expect to at least reach 5% of them to access one of our 6 machines (ca. 70 per machine) and among these reaching ca. 35%.
5) confidence interval set to 95%: The probability that our sample results will be within the margin of error of the true population estimate.
6) Margin of error (precision) set to 10%. The maximum acceptable difference between our sample estimate and the true population parameter.
7) Estimate the intracluster correlation coefficient (ICC): The ICC measures the similarity or correlation between responses within the same cluster. A higher ICC (closer to 1) implies more similarity (i.e. more correlation). Based on the ICC, the design effect (or inflation factor), is calculated. There is no reason to assume high correlation.
8) Due to operational/budget and piloting reasons: as few clusters as possible


```r
# Define parameters
confidence_level <- 0.95 # confidence level
z <- qnorm(1 - (1 - confidence_level) / 2) # z statistic
d <- 0.10 # margin of error
cluster_size <- 70 # Average cluster size
icc <- 0.05 # Intracluster correlation coefficient
prop <- 0.35 # Estimated proportion (= 35% reach)
deff <- 1 + (cluster_size - 1) * icc # Design Effect
sample_size <- ceiling((z^2 * prop * (1 - prop) * deff) / d^2)
# Number of individual participants
cat("Required sample size:", sample_size)
```

```
## Required sample size: 389
```

```r
# Number of clusters
n_clusters <- sample_size / cluster_size
cat("Required clusters:", round(n_clusters, 0))
```

```
## Required clusters: 6
```
