# Pancreatic-cancer-subtypes-survival-linear-predictor
Show that there is no significant relationship between cancer subtypes and our linear predictor of survival



## DATA
- 4 subtypes : 1squamous 2immunogenic 3progenitor 4ADEX
- Retrieve by querying the TCGA GDC database (type = paad)
- Merge with the linear predictor of the models of interest (8)

### Objective: to verify that the cancer subtypes are independent of the linear predictor in continuous form (Anova) and dichotomized on Maxstat (Fisher)  

## RESULTS
The ANOVA test makes the following assumptions about the data:

- Independence of observations. Each subject must belong to only one group. There is no relationship between the observations in each group. No repeated measures are allowed for the same participants.

- No significant outliers in any cell of the design

- Normality. The data in each cell of the design should be approximately normally distributed.

- Homogeneity of variances. The variance of the response variable should be equal in each cell of the design.

#### Groups vs cancer_subtypes
![Image of aciduino on protoboard](https://github.com/dinaOuahbi/pancreatic-cancer-subtypes-survival-linear-predictor/blob/main/groups_subtypes.png)

#### Linear predictor vs cancer subtypes _ some statistics 
![Image of aciduino on protoboard](https://github.com/dinaOuahbi/pancreatic-cancer-subtypes-survival-linear-predictor/blob/main/linearPred_subtypes.png)

#### Linear predictor vs cancer subtypes _ Normality
![Image of aciduino on protoboard](https://github.com/dinaOuahbi/pancreatic-cancer-subtypes-survival-linear-predictor/blob/main/norm_hypo.png)
- p values (Shapiro test)
  - 1 : 0.285
  - 2 : 0.0000611
  - 3 : 0.125
  - 4 : 0.00511

#### Linear predictor vs cancer subtypes _ Homogeneity of variances
![Image of aciduino on protoboard](https://github.com/dinaOuahbi/pancreatic-cancer-subtypes-survival-linear-predictor/blob/main/homoG_hypo.png)
- Levene Test : 0.0450

#### Linear predictor vs cancer subtypes _ welch ANOVA
![Image of aciduino on protoboard](https://github.com/dinaOuahbi/pancreatic-cancer-subtypes-survival-linear-predictor/blob/main/welsh_anova.png)





##### thanks for following (OD)
