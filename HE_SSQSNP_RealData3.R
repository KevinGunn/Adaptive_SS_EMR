#####################################################
##         Real Data Code on Hypotension Data Set  ##
##                  Kevin Gunn                     ##
##                  5/27/2017                      ##
#####################################################

library("car")
library('ggplot2')
library('MASS')

library("readr")

#n = 300

# Other Seeds.
#### 42173, #### 39433, #### 47516, #### 96279, #### 11946, #### 48974, #### 30738


set.seed(42173)

# Load Dataset
HE_dataset_full_in <- read_csv("EMR-research/SQL/HE_MY_CODE/HYPO_FINAL_DS.csv")
#View(HE_dataset_full_in)

#########################################################################################

#KS_2D_Gauss performs kernel smoothing on Yt and outputs a vector for imputed data.  
KS_2D_Gauss = function(Yt.v, X_impute , X_label, const=c(2,2)){ 
  
  gauss <- function(x) 1/sqrt(2*pi) * exp(-(x^2)/2)
  #Calculate bandwidth first as it is used in every iteration.
  n = length(Yt.v)
  m = dim(X_impute)[1]
  
  h1 = const[1]*n^(-1/5)
  h2 = const[2]*n^(-1/5)
  
  kde1.g = matrix(0 , ncol = m , nrow = n)
  kde2.g = matrix(0 , ncol = m , nrow = n)
  
  C_np_impute = rep(0, m)
  
  for(j in 1:m){
    kde1.g[,j] = (1/h1)*gauss((X_label[,c("X1")] - X_impute[j,c("X1")])/h1)
    kde2.g[,j] = (1/h2)*gauss((X_label[,c("X2")] - X_impute[j,c("X2")])/h2)
    C_np_impute[j] = sum(kde1.g[,j]*kde2.g[,j]*Yt.v) / sum(kde1.g[,j]*kde2.g[,j])
  }
  #C_np_impute[is.nan(C_np_impute)] = 0
  return(C_np_impute)
}

#CV1 function - use directly with KS_2D_Gauss function.
cv.h <- function(Yt.in , x.in , seq.c){
  
  n = length(Yt.in)
  num.folds=5
  bandwidths <- expand.grid(seq.c,seq.c)
  bw_seq_length <- dim(bandwidths)[1]
  fold_MSEs <- matrix(0,nrow=num.folds,
                      ncol=bw_seq_length)
  
  colnames(fold_MSEs) <- 1:bw_seq_length
  case.folds <- rep(1:num.folds,length.out=n)
  case.folds <- sample(case.folds)
  
  for (fold in 1:num.folds) {
    test.rows = which(case.folds==fold)
    x.test = x.in[test.rows,]
    y.test = Yt.in[test.rows]
    x.train = x.in[-test.rows,]
    y.train = Yt.in[-test.rows]
    for (bw in 1:bw_seq_length) {
      mx.hat <- KS_2D_Gauss( y.train , x.test[,c("X1","X2")] , 
                             x.train[,c("X1","X2")] , const = c(bandwidths[bw,1] , bandwidths[bw,2]) ) 
      fold_MSEs[fold,paste(bw)] <- sum((y.test - mx.hat)^2)
    }
  }
  
  CV_MSEs = colMeans(fold_MSEs)
  
  best.bw = bandwidths[which.min(CV_MSEs),]
  return(as.numeric(best.bw))
}

hlscv.ks <- function(Yt.in , x.in , x.impute , const_in, prop_score){
  
  ####
  ## Get mu_k.hat for use in the IF functions
  ####
  n = length(Yt.in)
  num.folds=5
  m = dim(x.impute)[1]
  
  fold_dfs = data.frame()
  imp_mat = matrix(0, ncol=num.folds , nrow = m)
  
  case.folds <- rep(1:num.folds,length.out=n)
  case.folds <- sample(case.folds)
  
  #This for loop gets mx_k.hat. needs to be readjusted with coefficients for k !=k'.
  for (fold in 1:num.folds) {
    test.rows = which(case.folds==fold,arr.ind=TRUE)
    x.test = x.in[test.rows,]
    y.test = Yt.in[test.rows]
    x.train = x.in[-test.rows,]
    y.train = Yt.in[-test.rows]
    
    #kernel smoothing
    mx.hat <- KS_2D_Gauss(y.train , x.test[,c("X1","X2")] , 
                          x.train[,c("X1","X2")] , const = const_in)
    
    mx.hat.imp <- KS_2D_Gauss(y.train , x.impute[,c("X1","X2")] , 
                              x.train[,c("X1","X2")] , const = const_in)
    
    #Add to data set to perform regression on residuals
    offLM = cbind(y.test , mx.hat , x.test , 
                  rep(fold , length(y.test)) )
    #print(head(x.test))
    #print(head(offLM))
    colnames(offLM) <- c("yt.test" , "mx_k.hat","int","X1","X2","fold")
    
    fold_dfs = rbind(fold_dfs , offLM)
    imp_mat[,fold] = mx.hat.imp 
    
  }
  
  #Need to change this part so it fits coefficients for the different folds.
  return_df = fold_dfs[order(as.numeric(rownames(fold_dfs))),]
  offlm.model = lm(return_df$yt.test ~ X1 + X2 , data = return_df , offset = mx_k.hat,
                   weights = prop_score^-1)
  
  beta.off = offlm.model$coefficients
  
  #mu_hat = return_df$mx_k.hat + beta.off%*%t(return_df[c("int","X1","X2")])
  mu_hat = as.vector(rowMeans(imp_mat) + beta.off%*%t(x.impute))
  mu_all = list(mu_hat,beta.off)
  return(mu_all)
}

double_cv.ks1 <- function(Yt.in , x.in , const_in, prop_score){
  
  ####
  ## Get mu_k.hat for use in the IF functions
  ####
  n = length(Yt.in)
  num.folds=5
  
  fold_dfs = data.frame()
  
  case.folds <- rep(1:num.folds,length.out=n)
  case.folds <- sample(case.folds)
  
  #This for loop gets mx_k.hat. needs to be readjusted with coefficients for k !=k'.
  for (fold in 1:num.folds) {
    test.rows = which(case.folds==fold,arr.ind=TRUE)
    x.test = x.in[test.rows,]
    y.test = Yt.in[test.rows]
    x.train = x.in[-test.rows,]
    y.train = Yt.in[-test.rows]
    prop_score.train = prop_score[-test.rows]
    prop_score.test = prop_score[test.rows]
    
    #kernel smoothing
    mx.hat <- KS_2D_Gauss(y.train , x.test[,c("X1","X2")] , 
                          x.train[,c("X1","X2")] , const = const_in)
    
    #Add to data set to perform regression on residuals
    offLM = cbind(y.test , mx.hat , x.test , 
                  rep(fold , length(y.test)) )
    
    colnames(offLM) <- c("yt.test" , "mx_k.hat","int","X1","X2","fold")
    
    fold_dfs = rbind(fold_dfs , offLM )
  }
  
  return_df = fold_dfs[order(as.numeric(rownames(fold_dfs))),]
  #print(return_df)
  mu_bc = c()
  for (k in 1:num.folds){
    # folds when k != k'
    fold_df = subset( return_df , fold==k, select = -fold )
    #print(fold_df)
    fold_lm.model = lm(fold_df$yt.test ~ X1 + X2 , data = fold_df , offset = mx_k.hat )
    beta.off = fold_lm.model$coefficients
    
    # k = k'
    kfold_df = subset( return_df , fold==k, select = -fold )
    mu_hat =  kfold_df$mx_k.hat + beta.off%*%t( kfold_df[,c("int","X1","X2")] )
    mu_bc = c(mu_bc,mu_hat)
  }
  mu_bc = mu_bc[order(as.numeric(rownames(fold_dfs)))]
  return(mu_bc)
}
double_cv.ks0 <- function(Yt.in , x.in , const_in, prop_score){
  
  ####
  ## Get mu_k.hat for use in the IF functions
  ####
  n = length(Yt.in)
  num.folds=5
  
  fold_dfs = data.frame()
  
  case.folds <- rep(1:num.folds,length.out=n)
  case.folds <- sample(case.folds)
  
  #This for loop gets mx_k.hat. needs to be readjusted with coefficients for k !=k'.
  for (fold in 1:num.folds) {
    test.rows = which(case.folds==fold,arr.ind=TRUE)
    x.test = x.in[test.rows,]
    y.test = Yt.in[test.rows]
    x.train = x.in[-test.rows,]
    y.train = Yt.in[-test.rows]
    prop_score.train = prop_score[-test.rows]
    prop_score.test = prop_score[test.rows]
    
    #kernel smoothing
    mx.hat <- KS_2D_Gauss(y.train , x.test[,c("X1","X2")] , 
                          x.train[,c("X1","X2")] , const = const_in)
    
    #Add to data set to perform regression on residuals
    offLM = cbind(y.test , mx.hat , x.test , 
                  rep(fold , length(y.test)) )
    
    colnames(offLM) <- c("yt.test" , "mx_k.hat","int","X1","X2","fold")
    
    fold_dfs = rbind(fold_dfs , offLM )
  }
  
  
  return_df = fold_dfs[order(as.numeric(rownames(fold_dfs))),]
  #print(return_df)
  mu_bc = c()
  for (k in 1:num.folds){
    # folds when k != k'
    fold_df = subset( return_df , fold!=k, select = -fold )
    #print(fold_df)
    fold_lm.model = lm(fold_df$yt.test ~ X1 + X2 , data = fold_df , offset = mx_k.hat )
    beta.off = fold_lm.model$coefficients
    
    # k = k'
    kfold_df = subset( return_df , fold==k, select = -fold )
    mu_hat =  kfold_df$mx_k.hat + beta.off%*%t( kfold_df[,c("int","X1","X2")] )
    mu_bc = c(mu_bc,mu_hat)
  }
  mu_bc = mu_bc[order(as.numeric(rownames(fold_dfs)))]
  return(mu_bc)
}

# Influence Function
IF_sd <- function(XtX.inv,X,Y,A,mv,prop ){
  
  #prop1 vector of propensity scores for patients assigned to trt 1
  #prop0 vector of propensity scores for patients assigned to trt 0
  
  p = dim(XtX.inv)[2]
  varm1 = matrix(0,ncol=p,nrow=p)
  varm0 = matrix(0,ncol=p,nrow=p)
  n = length(Y)
  
  for(i in 1:n){
    
    if(A[i]==1){
      # Variance  of IF
      #print((Y[i] - mv[i]))
      IF1 = X[i,]*(Y[i] - mv[i])*prop[i]^-1
      Var_IF1 = ( IF1%*%t(IF1) ) 
      varm1 = varm1 + Var_IF1
    } else{
      #print((Y[i] - mv[i]))
      IF0 = -X[i,]*(Y[i] - mv[i])*(1-prop[i])^-1 
      Var_IF0 = ( IF0%*%t(IF0) ) 
      varm0 = varm0 + Var_IF0
    } 
    
  }
  #print(varm1+varm0)
  V1 = XtX.inv%*%varm1%*%t(XtX.inv) 
  V0 = XtX.inv%*%varm0%*%t(XtX.inv) 
  SE = sqrt(diag((V1+V0)/n^2 ) )
  #print(V1);print(V0)
  return(SE)
}

IF_ols_sd <- function(XtX.inv,X,Yt,mv ){
  
  p = dim(XtX.inv)[2]
  varm = matrix(0,ncol=p,nrow=p)
  n = length(Yt)
  
  for(i in 1:n){
    
    # Variance  of IF
    IF = X[i,]*(Yt[i] - mv[i])
    Var_IF = ( IF%*%t(IF) ) 
    varm = varm + Var_IF
    
  }
  #print(varm)
  V = XtX.inv%*%varm%*%t(XtX.inv) 
  SE = sqrt(diag((V)/n^2 ) )
  return(SE)
}

full_IF_ols_sd <- function(XtX.inv,X,Yt,mv ){
  
  p = dim(XtX.inv)[2]
  varm = matrix(0,ncol=p,nrow=p)
  n = length(Yt)
  
  for(i in 1:n){
    
    # Variance  of IF
    IF = as.matrix( X[i,]*(Yt[i] - mv[i]) )
    Var_IF =  t(IF) %*% IF 
    varm = varm + Var_IF
    
  }
  #print(varm)
  V = XtX.inv%*%varm%*%t(XtX.inv) 
  SE = sqrt(diag((V)/n^2 ) )
  return(SE)
}

###########################################################################################

# Initial Clean up.
head(HE_dataset_full_in)

# Make elixhauser scores numeric.
HE_dataset_full_in[,28:30] <- lapply(HE_dataset_full_in[,28:30], as.numeric)

# Remove patients that recieved both.
#HE_dataset_full_in <- subset(HE_dataset_full_in,!( ( iv_fluid_ind ==1) & ( vp_ind ==1 ) ) )
#HE_dataset_full_in <- HE_dataset_full_in[which(!is.na(HE_dataset_full_in$mean_resp)),] 

# Create missing indicator for each patient.
missing_ind_resp <- ifelse( is.na(HE_dataset_full_in$rise_creat) , 1 , 0)

missing_ind_trt <- ifelse(( ( HE_dataset_full_in$iv_fluid_ind ==0) & ( HE_dataset_full_in$vp_ind ==0 ) ), 1 , 0)

# Include patients who recieved both as patients that are not in the IV or VP cohort.
missing_ind_trt2 <- ifelse(( ( HE_dataset_full_in$iv_fluid_ind ==1) & ( HE_dataset_full_in$vp_ind ==1 ) ), 1 , 0)

missing_ind_count <- missing_ind_resp + missing_ind_trt + missing_ind_trt2

missing_ind <- ifelse(missing_ind_count>0, 1, 0) # 2118 with missing data.

###########################
## Dataset with missing indicator ##
HE_Cohort_Fluid_VP <- cbind(HE_dataset_full_in, missing_ind)

############################################################################################################
# Normalize covariates.

norm_tuv <- scale(HE_Cohort_Fluid_VP$tot_urine_vol)

norm_pre_cr <-  scale(HE_Cohort_Fluid_VP$pre_creat)

norm_mean_map <-  scale(HE_Cohort_Fluid_VP$mean_map)

norm_mean_spo2 <-  scale(HE_Cohort_Fluid_VP$mean_spo2)

norm_mean_resp <-  scale(HE_Cohort_Fluid_VP$mean_resp)

norm_mean_heart_rate <-  scale(HE_Cohort_Fluid_VP$mean_heart_rate)

norm_age <-  scale(HE_Cohort_Fluid_VP$age_numeric)

log_elix <- log(-min(HE_Cohort_Fluid_VP$elixhauser_vanwalraven) + 
                  1 + HE_Cohort_Fluid_VP$elixhauser_vanwalraven )

log_saps <- log(HE_Cohort_Fluid_VP$saps)

###########################################################################################################

## Dataset with scaled covariates.
HE_Cohort_Fluid_VP_scale <- cbind(HE_Cohort_Fluid_VP, norm_pre_cr, norm_tuv, norm_mean_map, norm_mean_resp,
                                  norm_mean_spo2, norm_mean_heart_rate, norm_age, log_elix, log_saps )

# service indicator -> did patient recieve surgical or other services?
service_ind <- ifelse(HE_Cohort_Fluid_VP_scale$curr_serv 
                      %in% c("CSURG","SURG","VSURG","TSURG","PSURG"),1,0)

# Gender indicator 
gender_ind <- ifelse(HE_Cohort_Fluid_VP_scale$gender=='M',1,0)

HE_Cohort_Fluid_VP_all <- cbind(HE_Cohort_Fluid_VP_scale,service_ind, gender_ind)

# Seperate two treatment groups.
IV_Cohort <- subset( HE_Cohort_Fluid_VP_all, iv_fluid_ind==1 &  vp_ind==0 & missing_ind == 0 )
VP_Cohort <- subset( HE_Cohort_Fluid_VP_all, vp_ind==1 & iv_fluid_ind==0 & missing_ind == 0)


trt_set1 <- rbind(IV_Cohort, VP_Cohort)

#trt_ind = 1 if Vasoactive therapy fluid, trt_ind = 0 if IV_Fluid.
trt_ind <- c( rep(0,dim(IV_Cohort)[1]) , rep(1,dim(VP_Cohort)[1]) ) 

trt_set <- cbind(trt_set1, trt_ind)

###########################################################################################################
###########################################################################################################
# Missing Data set

### NOTE: HAVE TO DROP PATIENTS W/O BASELINE CREATININE MEASUREMENTS.
unlabeled_set <- subset(HE_Cohort_Fluid_VP_all, missing_ind == 1 )

############################################################################################################

# Fully supervised approach with entire treated cohort.

# Use all data to fit a propensity score.
trt_propensity_all <- glm(trt_ind ~  norm_pre_cr + norm_tuv + norm_mean_map + norm_age + 
                            norm_mean_heart_rate + log_elix + log_saps + gender_ind + service_ind
                          ,data=trt_set,family = binomial)

prop_score_all <-  trt_propensity_all$fitted.values


# Decrease in serum creatinine = good kidney function.
Yt <- (trt_set$rise_creat*(trt_set$trt_ind - prop_score_all)) / 
  (prop_score_all*(1-prop_score_all)) 

Lin_ds <- as.data.frame(cbind(Yt, trt_set))

Lin_Reg <- lm(-Yt ~ norm_pre_cr + norm_tuv + norm_mean_map + norm_age + 
                norm_mean_heart_rate + 
                log_elix + log_saps + gender_ind + service_ind ,data=Lin_ds )


beta.Yt <- Lin_Reg$coefficients
beta.Yt.full <- beta.Yt

int.trt <- rep(1,dim(Lin_ds)[1])

X.trt <- cbind(int.trt, Lin_ds$norm_pre_cr, Lin_ds$norm_tuv, Lin_ds$norm_mean_map, Lin_ds$norm_age,
               Lin_ds$norm_mean_heart_rate, Lin_ds$log_elix, Lin_ds$log_saps,
               Lin_ds$gender_ind, Lin_ds$service_ind)

X.trt <- data.frame(apply(X.trt, 2, as.numeric))

trt_xtx = t(as.matrix(X.trt)) %*% as.matrix(X.trt) / dim(X.trt)[1]
Lambda_trt_inv = solve(trt_xtx)

# IF estimation
IF_sd_OLS <- full_IF_ols_sd(Lambda_trt_inv, X.trt, -Yt, Lin_Reg$fitted.values)
print(IF_sd_OLS)

# p-value calculation
#pvalue2sided=2*pnorm(-abs(z))

z = Lin_Reg$coefficients / IF_sd_OLS
pvalue_tr = 2*pnorm(-abs(Lin_Reg$coefficients / IF_sd_OLS))

# Correlation of covariance matrix.
cor(X.trt[,-1])
##################################################################################################################

Lin_Reg2 <- lm(-Yt ~ norm_mean_map + norm_mean_heart_rate ,data=Lin_ds )


beta.Yt.2 <- Lin_Reg2$coefficients

X.trt2 <- cbind(int.trt, Lin_ds$norm_mean_map, Lin_ds$norm_mean_heart_rate)

X.trt2 <- data.frame(apply(X.trt2, 2, as.numeric))

trt_xtx2 = t(as.matrix(X.trt2)) %*% as.matrix(X.trt2) / dim(X.trt2)[1]
Lambda_trt_inv2 = solve(trt_xtx2)

# IF estimation
IF_sd_OLS2 <- full_IF_ols_sd(Lambda_trt_inv2, X.trt2, -Yt, Lin_Reg2$fitted.values)
print(IF_sd_OLS2)

# p-value calculation
#pvalue2sided=2*pnorm(-abs(z))

z = Lin_Reg2$coefficients / IF_sd_OLS2
pvalue_tr2 = 2*pnorm(-abs(Lin_Reg2$coefficients / IF_sd_OLS2))
pvalue_tr2

# Correlation of covariance matrix.
cor(X.trt2[,-1])


##################################################################################################################

# Now perform Semi-Supervised method.

##################################################################################################################
# subset 300 from labeled set.
n=300
#n=350  

# VP trt.
VP_count = length(which(trt_set$trt_ind==1))
# IV fluid trt
IV_count = length(which(trt_set$trt_ind==0))
#proportions
IV_prop = IV_count / (IV_count + VP_count)
VP_prop = VP_count / (VP_count + IV_count)

# subsample
IV_samp_num <- sample(dim(IV_Cohort)[1], size = floor(n*IV_prop))
VP_samp_num <- sample(dim(VP_Cohort)[1], size = ceiling(n*VP_prop))

# Simulated complete data.
IV_subset <- IV_Cohort[IV_samp_num, ]
VP_subset <- VP_Cohort[VP_samp_num, ]

# left out data.
IV_out <- IV_Cohort[-IV_samp_num, ]
VP_out <- VP_Cohort[-VP_samp_num, ]
hold_set <- cbind( rep(1,dim(IV_out)[1]+dim(VP_out)[1]), c( rep(0,dim(IV_out)[1]) , rep(1,dim(VP_out)[1]) )  , 
                   rbind(IV_out,VP_out) )

# Updated unlabeled set.
unlabeled_ds <- rbind(unlabeled_set, IV_out, VP_out )
int_unlabeled <- rep( 1, dim(unlabeled_ds)[1] )
unlabeled_ds_int <- cbind( unlabeled_ds, int = int_unlabeled )

#unlabeled set
unlabeled_SS <- as.matrix(unlabeled_ds_int[,c(47,40,43)])

colnames(unlabeled_SS) <- c("int","X1","X2")
rownames(unlabeled_SS) <- c()

## Supervised Regression with complete cases.
trt0_num = length(rep(0,dim(IV_subset)[1]))
trt_ind <- c( rep(0,dim(IV_subset)[1]) , rep(1,dim(VP_subset)[1]) ) 
trt_set1 <- rbind(IV_subset,VP_subset)
trt_set_in <- cbind(trt_set1, trt_ind)

# Use all data to fit a propensity score.
trt_propensity <- glm(trt_ind ~  norm_pre_cr + norm_tuv + norm_mean_map + norm_age + norm_mean_heart_rate + 
                        log_elix + log_saps + gender_ind + service_ind
                      ,data=trt_set_in,family = binomial)
prop_score <-  trt_propensity$fitted.values

## -------- ## --------- ## -------- ## -------- ## ---------- ##
# All variables on subset of 300.
Yt <- (trt_set_in$rise_creat*(trt_set_in$trt_ind - prop_score)) / (prop_score*(1-prop_score)) 

Lin_Reg <- lm(-Yt ~ norm_pre_cr + norm_tuv + norm_mean_map + norm_age + norm_mean_heart_rate + 
                log_elix + log_saps + gender_ind + service_ind ,data=trt_set_in )
summary(Lin_Reg)

## ------------ ## -------------- ## -------------- ## ---------- ##

# Treatment cohort subsets.
IV_subset_int <- cbind(IV_subset, rep(1,dim(IV_subset)[1]))
VP_subset_int <- cbind(VP_subset, rep(1,dim(VP_subset)[1]))

# Make subsets of dataset to work in semisupervised model.
IV_SS <- as.matrix(cbind(IV_subset_int[,c(23,47,40,43)],prop_score_in = prop_score[1:trt0_num] ) )
VP_SS <- as.matrix(cbind( VP_subset_int[,c(23,47,40,43)],prop_score_in = prop_score[-(1:trt0_num)] ) )
colnames(IV_SS) <- colnames(VP_SS) <- c("chg_crt","int","X1","X2","prop_score")
rownames(IV_SS) <- rownames(VP_SS) <- c()

trt_set_in_2covs <- as.data.frame(cbind(rbind(IV_SS,VP_SS), trt_ind))

hc1 = cv.h(VP_SS[,c("chg_crt")], x.in = VP_SS[,c("X1","X2")],seq.c = seq(0.1,100,2) )
hc0 = cv.h(IV_SS[,c("chg_crt")], x.in= IV_SS[,c("X1","X2")],seq.c = seq(0.1,100,2) )
print(hc1);print(hc0)

# Same approach as simulations
mx.hat_ls1 = hlscv.ks(VP_SS[,c("chg_crt")],x.in= VP_SS[,c("int","X1","X2")],x.impute=unlabeled_SS[,c("int","X1","X2")]
                      ,const_in=hc1, prop_score = VP_SS[,c("prop_score")])[[1]]

mx.hat_ls0 = hlscv.ks(IV_SS[,c("chg_crt")] , x.in= IV_SS[,c("int","X1","X2")],x.impute=unlabeled_SS[,c("int","X1","X2")]
                      , const_in = hc0, prop_score = 1 - IV_SS[,c("prop_score")])[[1]]

CX_ls = mx.hat_ls1 - mx.hat_ls0

X.imp = as.data.frame(cbind(CX_ls, unlabeled_SS))

colnames(X.imp) <- c("CX_ls","int","norm_mean_map","norm_mean_heart_rate")

Reg.CX_ls = lm(-CX_ls ~ norm_mean_map + norm_mean_heart_rate , data = X.imp)
beta.CX_ls = Reg.CX_ls$coefficients
summary(Reg.CX_ls)
leveragePlots(Reg.CX_ls)
outlierTest(Reg.CX_ls)


# Decrease in serum creatinine = good kidney function.
Yt <- (trt_set_in_2covs$chg_crt*(trt_set_in_2covs$trt_ind - trt_set_in_2covs$prop_score)) / 
  (trt_set_in_2covs$prop_score*(1-trt_set_in_2covs$prop_score)) 

Lin_ds_subset <- as.data.frame(cbind(Yt, trt_set_in_2covs))

Lin_Reg <- lm(-Yt ~ X1 + X2 ,data=Lin_ds_subset )
summary(Lin_Reg)

beta.Yt = Lin_Reg$coefficients

###################################################################################################################
# Variance Estimation

Y_A1 = Lin_ds_subset[which(Lin_ds_subset$trt_ind == 1),2]
X.L1 = as.matrix(Lin_ds_subset[which(Lin_ds_subset$trt_ind == 1),c(3,4,5)]) 
colnames(X.L1) <- c("int","X1","X2")

Y_A0 = Lin_ds_subset[which(Lin_ds_subset$trt_ind == 0),2]
X.L0 = as.matrix(Lin_ds_subset[which(Lin_ds_subset$trt_ind == 0),c(3,4,5)])
colnames(X.L0) <- c("int","X1","X2")

prop1 = prop_score[trt_ind==1]; prop0 = prop_score[trt_ind==0]

# Inverse of XtX
X.L = rbind(X.L0, X.L1)
Y.if = c(Y_A0, Y_A1)
prop = c(prop0,prop1)
X.m = as.matrix(X.imp[,-1])

XL_inv = t(as.matrix(X.L)) %*% as.matrix(X.L) / dim(X.L)[1]
Lambda_n_inv = solve(XL_inv)


Xm_inv = t(as.matrix(X.m)) %*% as.matrix(X.m) / dim(X.m)[1]
Lambda_N_inv = solve(Xm_inv)

X.all = rbind(X.m, X.L )
full_xtx = t(as.matrix(X.all)) %*% as.matrix(X.all) / dim(X.all)[1]
Lambda_all_inv = solve(full_xtx)

# IF estimation
# double estimates for IF.

mv1 = double_cv.ks1(-Y_A1, X.L1, hc1, prop1)
mv0 = double_cv.ks0(-Y_A0, X.L0, hc0, 1-prop0)

IF_sd_CXls = IF_sd(Lambda_all_inv, X.L, -Y.if, trt_ind, c(mv0, mv1) , prop )
print(IF_sd_CXls)

IF_sd_OLS = IF_ols_sd(Lambda_all_inv, X.L, -Yt, Lin_Reg$fitted.values)
print(IF_sd_OLS)

# p-value calculation
#pvalue2sided=2*pnorm(-abs(z))

beta.CX_ls / IF_sd_CXls
pvalue_ss = 2*pnorm(-abs(beta.CX_ls / IF_sd_CXls))

Lin_Reg$coefficients / IF_sd_OLS
pvalue_tr = 2*pnorm(-abs(Lin_Reg$coefficients / IF_sd_OLS))

# RE
IF_sd_OLS^2 / IF_sd_CXls^2

#CI
beta.CX_ls + 1.96*IF_sd_CXls
beta.CX_ls - 1.96*IF_sd_CXls

beta.Yt + 1.96*IF_sd_OLS
beta.Yt - 1.96*IF_sd_OLS

# output of p-values
pvalue_ss

pvalue_tr
###################################################################################################################

# Treatment recommendation
OLS_trt <- ifelse(beta.Yt%*%t(X.all[,c("int","norm_mean_map","norm_mean_heart_rate")]) > 0 ,1,0)
SS_trt <- ifelse(beta.CX_ls%*%t(X.all[,c("int","norm_mean_map","norm_mean_heart_rate")]) > 0 ,1,0)

OLS_trt_full <- ifelse(beta.Yt.full %*% t(X.trt) > 0 ,1,0)
sum(OLS_trt_full)
SS_trt_onlytrt <- ifelse(beta.CX_ls %*% t(X.trt[,c(1,4,6)]) > 0 ,1,0)

# Distribution Tests.
ks.test(X.L[,2],X.imp[,3])
ks.test(X.L[,3],X.imp[,4])

wilcox.test(trt_set_in$mean_map,unlabeled_ds_int$mean_map)
wilcox.test(trt_set_in$mean_heart_rate,unlabeled_ds_int$mean_heart_rate)

ks.test(trt_set_in$mean_map,unlabeled_ds_int$mean_map)
ks.test(trt_set_in$mean_heart_rate,unlabeled_ds_int$mean_heart_rate)

mean(trt_set_in$mean_map);sd(trt_set_in$mean_map)
mean(unlabeled_ds_int$mean_map);sd(unlabeled_ds_int$mean_map)

mean(trt_set_in$mean_heart_rate);sd(trt_set_in$mean_heart_rate)
mean(unlabeled_ds_int$mean_heart_rate);sd(unlabeled_ds_int$mean_heart_rate)

##################################################################################################

sum((OLS_trt)*(SS_trt))

sum((1-OLS_trt)*(1-SS_trt))

sum((1-OLS_trt)*SS_trt)

sum(OLS_trt*(1-SS_trt))

table_22 = cbind(c(sum((OLS_trt)*(SS_trt)),sum(OLS_trt*(1-SS_trt))), 
                 c(sum((1-OLS_trt)*(1-SS_trt)), sum((1-OLS_trt)*SS_trt)))

##################################################################################################
# Compare with OLS_FULL
sum((OLS_trt_full)*(SS_trt_onlytrt))

sum((1-OLS_trt_full)*(1-SS_trt_onlytrt))

sum((1-OLS_trt_full)*SS_trt_onlytrt)

sum(OLS_trt_full*(1-SS_trt_onlytrt))

######################################################################################################
#qplot(prop_score, binwidth = 0.03, geom = "histogram", fill=I("red"),xlab="Propensity Score")


wilcox.test(trt_set_in$mean_heart_rate,unlabeled_ds_int$mean_heart_rate)

beta.Yt.2
IF_sd_OLS2
pvalue_tr2


IF_sd_OLS^2 / IF_sd_CXls^2

pvalue_ss

print(dim(trt_set_in));print(dim(unlabeled_ds_int))
