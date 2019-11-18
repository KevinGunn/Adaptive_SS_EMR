# Adaptive_SS_EMR

The files “SSQ_SNP_code_example.R” and "SSQSNP_allsims_500.R" contain code to replicate the methods of 
"Adaptive Semi-Supervised Inference for Optimal Treatment Decisions with Electronic Medical Record Data" 
by Gunn, Lu, and Song.  An example is provided of the component-wise bias, empirical SE, asymptotic SE, 
coverage probability for the 95% CI, and component-wise RE calculated on simulated data. 

The file contains R functions to perform SSQ-SNP as described in Section 4.2, and estimate the asymptotic SE 
for the regression coefficients using the variance estimator explained in Section 4.3.

Additionally, the R packages MASS and RandomForest may need to be installed to run code for the simulated data. 

A description of the R functions used can be found below.

#### KS_2D_Gauss: 
Performs kernel smoothing on Y with 2 covariates and outputs a vector for imputed data.  

#### cv.h:
Obtains the bandwidth, h, used in the kernel smoothing function, KS_2D_Gauss.

#### hlscv.ks:
Estimates theta_1 and theta_0 and then predicts missing response with Q(x,a,theta).

#### double_cv.ks:
Performs double cross validation which is used to get an estimator for the variance of the parameters, beta.

#### IF_se:


#### IF_OLS_se
