#####################################################
##          Missing Data OTR Simulation            ##
##                  Kevin Gunn                     ##
##                  12/15/2018                     ##
#####################################################

library('MASS')
library("randomForest")

set.seed(12345)

#########################################################################################
#########################################################################################
#KS_2D_Gauss performs kernel smoothing on Y and outputs a vector for imputed data.  
KS_2D_Gauss = function(Yt.v, X_impute , X_label, const=2){ 
  
  gauss <- function(x) 1/sqrt(2*pi) * exp(-(x^2)/2)
  #Calculate bandwidth first as it is used in every iteration.
  n = length(Yt.v)
  m = dim(X_impute)[1]
  
  h1 = const*n^(-1/5)
  h2 = const*n^(-1/5)
  
  kde1.g = matrix(0 , ncol = m , nrow = n)
  kde2.g = matrix(0 , ncol = m , nrow = n)
  
  C_np_impute = rep(0, m)
  
  for(j in 1:m){
    kde1.g[,j] = (1/h1)*gauss((X_label$X1 - X_impute$X1[j])/h1)
    kde2.g[,j] = (1/h2)*gauss((X_label$X2 - X_impute$X2[j])/h2)
    C_np_impute[j] = sum(kde1.g[,j]*kde2.g[,j]*Yt.v) / sum(kde1.g[,j]*kde2.g[,j])
  }
  C_np_impute[is.nan(C_np_impute)] = 0
  return(C_np_impute)
}

#CV.h function - use directly with KS_2D_Gauss function.
cv.h <- function(Yt.in , x.in , seq.c){
  
  n = length(Yt.in)
  num.folds=10
  bandwidths <- seq.c
  fold_MSEs <- matrix(0,nrow=num.folds,
                      ncol=length(bandwidths))
  
  colnames(fold_MSEs) = bandwidths
  case.folds <- rep(1:num.folds,length.out=n)
  case.folds <- sample(case.folds)
  
  for (fold in 1:num.folds) {
    test.rows = which(case.folds==fold)
    x.test = x.in[test.rows,-1]
    y.test = Yt.in[test.rows]
    x.train = x.in[-test.rows,-1]
    y.train = Yt.in[-test.rows]
    for (bw in bandwidths) {
      mx.hat <- KS_2D_Gauss(y.train , x.test[c("X1","X2")] , 
                            x.train[c("X1","X2")] , const = bw)
      fold_MSEs[fold,paste(bw)] <- sum((y.test - mx.hat)^2)
    }
  }
  
  CV_MSEs = colMeans(fold_MSEs)
  best.bw = bandwidths[which.min(CV_MSEs)]
  return(best.bw)
}

# hlscv.ks estimates theta_1 and theta_0 and then predicts missing response with Q(x,a,theta). 
# This function implements SSQ-SNP.
hlscv.ks <- function(Yt.in , x.in , x.impute , const_in, prop_score){
  
  ####
  ## Obtain Q(X,a,theta)
  ####
  n = length(Yt.in)
  num.folds=5
  m = dim(x.impute)[1]
  
  fold_dfs = data.frame()
  imp_mat = matrix(0, ncol=num.folds , nrow = m)
  
  case.folds <- rep(1:num.folds,length.out=n)
  case.folds <- sample(case.folds)
  
  #This for loop gets mx_k.hat. 
  for (fold in 1:num.folds) {
    test.rows = which(case.folds==fold,arr.ind=TRUE)
    x.test = x.in[test.rows,]
    y.test = Yt.in[test.rows]
    x.train = x.in[-test.rows,]
    y.train = Yt.in[-test.rows]
    
    #kernel smoothing
    mx.hat <- KS_2D_Gauss(y.train , x.test[c("X1","X2")] , 
                          x.train[c("X1","X2")] , const = const_in)
    
    mx.hat.imp <- KS_2D_Gauss(y.train , x.impute[c("X1","X2")] , 
                              x.train[c("X1","X2")] , const = const_in)
    
    #Add to data set to perform regression on residuals
    offLM = cbind(y.test , mx.hat , x.test , 
                  rep(fold , length(y.test)) )
   
    colnames(offLM) <- c("yt.test" , "mx_k.hat","int","X1","X2","fold")
    
    fold_dfs = rbind(fold_dfs , offLM)
    imp_mat[,fold] = mx.hat.imp 
    
  }
  
 
  return_df = fold_dfs[order(as.numeric(rownames(fold_dfs))),]
  offlm.model = lm(return_df$yt.test ~ X1 + X2 , data = return_df , offset = mx_k.hat,
                   weights = prop_score^-1)
  
  beta.off = offlm.model$coefficients
  
  mu_hat = as.vector(rowMeans(imp_mat) + beta.off%*%t(x.impute))
  mu_all = list(mu_hat,beta.off)
  return(mu_all)
}

# Functions to obtain double-CV estimators for inference.
double_cv.ks1 <- function(Yt.in , x.in , const_in, prop_score){
  
  ####
  ## Get Q_k(X,a,1) for use in the IF functions.
  ####
  n = length(Yt.in)
  num.folds=5
  
  fold_dfs = data.frame()
 
  case.folds <- rep(1:num.folds,length.out=n)
  case.folds <- sample(case.folds)
  
  #This for loop gets mx_k.hat.
  for (fold in 1:num.folds) {
    test.rows = which(case.folds==fold,arr.ind=TRUE)
    x.test = x.in[test.rows,]
    y.test = Yt.in[test.rows]
    x.train = x.in[-test.rows,]
    y.train = Yt.in[-test.rows]
    prop_score.train = prop_score[-test.rows]
    prop_score.test = prop_score[test.rows]
    
    #kernel smoothing
    mx.hat <- KS_2D_Gauss(y.train , x.test[c("X1","X2")] , 
                          x.train[c("X1","X2")] , const = const_in)
    
    #Add to data set to perform regression on residuals
    offLM = cbind(y.test , mx.hat , x.test , 
                  rep(fold , length(y.test)), prop_score.test )
    
    colnames(offLM) <- c("yt.test" , "mx_k.hat","int","X1","X2","fold","prop_k")
    
    fold_dfs = rbind(fold_dfs , offLM )
  }
  
  return_df = fold_dfs[order(as.numeric(rownames(fold_dfs))),]
  
  mu_bc = c()
  for (k in 1:num.folds){
    # folds when k != k'
    fold_df = subset( return_df , fold!=k, select = -fold )
    
    fold_lm.model = lm(fold_df$yt.test ~ X1 + X2 , data = fold_df , offset = mx_k.hat,
                       weights = prop_k^-1)
    beta.off = fold_lm.model$coefficients
    
    # k = k'
    kfold_df = subset( return_df , fold==k, select = -fold )
    mu_hat =  kfold_df$mx_k.hat + beta.off%*%t( kfold_df[c("int","X1","X2")] )
    mu_bc = c(mu_bc,mu_hat)
  }
  mu_bc = mu_bc[order(as.numeric(rownames(fold_dfs)))]
  return(mu_bc)
}
double_cv.ks0 <- function(Yt.in , x.in , const_in, prop_score){
  
  ####
  ## Get Q_k(X,a,0) for use in the IF functions
  ####
  n = length(Yt.in)
  num.folds=5
  
  fold_dfs = data.frame()
  
  case.folds <- rep(1:num.folds,length.out=n)
  case.folds <- sample(case.folds)
  
  #This for loop gets mx_k.hat.
  for (fold in 1:num.folds) {
    test.rows = which(case.folds==fold,arr.ind=TRUE)
    x.test = x.in[test.rows,]
    y.test = Yt.in[test.rows]
    x.train = x.in[-test.rows,]
    y.train = Yt.in[-test.rows]
    prop_score.train = prop_score[-test.rows]
    prop_score.test = prop_score[test.rows]
    
    #kernel smoothing
    mx.hat <- KS_2D_Gauss(y.train , x.test[c("X1","X2")] , 
                          x.train[c("X1","X2")] , const = const_in)
    
    #Add to data set to perform regression on residuals
    offLM = cbind(y.test , mx.hat , x.test , 
                  rep(fold , length(y.test)), prop_score.test )
    
    colnames(offLM) <- c("yt.test" , "mx_k.hat","int","X1","X2","fold","prop_k")
    
    fold_dfs = rbind(fold_dfs , offLM )
  }
  
  return_df = fold_dfs[order(as.numeric(rownames(fold_dfs))),]
  
  mu_bc = c()
  for (k in 1:num.folds){
    # folds when k != k'
    fold_df = subset( return_df , fold!=k, select = -fold )
   
    fold_lm.model = lm(fold_df$yt.test ~ X1 + X2 , data = fold_df , offset = mx_k.hat,
                       weights = prop_k^-1)
    beta.off = fold_lm.model$coefficients
    
    # k = k'
    kfold_df = subset( return_df , fold==k, select = -fold )
    mu_hat =  kfold_df$mx_k.hat + beta.off%*%t( kfold_df[c("int","X1","X2")] )
    mu_bc = c(mu_bc,mu_hat)
  }
  mu_bc = mu_bc[order(as.numeric(rownames(fold_dfs)))]
  return(mu_bc)
}
###########################################################################################

# Influence Functions
# IF for SSQ-SNP.
IF_se <- function(XtX.inv,X,Y,A,mv,prop ){
  
  #prop1 vector of propensity scores for patients assigned to trt 1
  #prop0 vector of propensity scores for patients assigned to trt 0
  
  p = dim(XtX.inv)[2]
  varm1 = matrix(0,ncol=p,nrow=p)
  varm0 = matrix(0,ncol=p,nrow=p)
  n = length(Y)
  
  for(i in 1:n){
    
    if(A[i]==1){
      # IF function
      IF1 = X[i,]*(Y[i]- mv[i])*prop[i]^-1
      Var_IF1 = ( IF1%*%t(IF1) ) 
      varm1 = varm1 + Var_IF1
    } else{
      
      IF0 = -X[i,]*(Y[i] - mv[i])*(1-prop[i])^-1 
      Var_IF0 = ( IF0%*%t(IF0) ) 
      varm0 = varm0 + Var_IF0
    } 
    
  }
  
  V1 = XtX.inv%*%varm1%*%t(XtX.inv) 
  V0 = XtX.inv%*%varm0%*%t(XtX.inv) 
  SE = sqrt(diag((V1+V0)/n^2 ) )
  
  return(SE)
}

# IF for TR-OLS.
IF_OLS_se <- function(XtX.inv,X,Yt,bx.in ){
  
  #prop1 vector of propensity scores for patients assigned to trt 1
  #prop0 vector of propensity scores for patients assigned to trt 0
  
  p = dim(XtX.inv)[2]
  varm = matrix(0,ncol=p,nrow=p)
  n = length(Yt)
  
  for(i in 1:n){
    
    # Variance  of IF
    IF = X[i,]*(Yt[i] - bx.in[i])
    Var_IF = ( IF%*%t(IF) ) 
    varm = varm + Var_IF
    
  }
  
  V = XtX.inv%*%varm%*%t(XtX.inv) 
  SE = sqrt(diag(V/n^2 ) )

  return(SE)
}

#########################################################################################

# Monte Carlo Simulation to get true values of regression coefficients
MC_coeff <-function(model = "L2L", n.samp = 300000,p,g=NULL,theta,lambda,beta,alpha){
  
  
  X = mvrnorm(n=n.samp , mu = rep(0,p), Sigma = diag(p))
  
  colnames(X) <- paste("X",1:p,sep="")
  
  logit_coefs <- rep(1,ncol(X))
  true.prop = plogis( 0.5*X[,1] - 0.5*X[,2] )
  
  # Treatment Allocation.
  A = sapply(true.prop, function(x){ rbinom(1,1,x) })
  
  #Generate Y based on the given true model.
  if(model == "CL"){
    Y = as.vector((X%*%theta)^3 + A*lambda*(X%*%beta) + rnorm(n.samp,0,1))
    mu_X = (X%*%theta)^3; CX = X%*%beta
  }else if(model == "L2L"){
    Y = as.vector( (X%*%alpha)*(1+X%*%theta) + A*lambda*(X%*%beta) + rnorm(n.samp,0,1))
    mu_X = (X%*%alpha)*(1+X%*%theta); CX = X%*%beta
  }else if(model == "L2NL3"){
    Y = as.vector(X%*%alpha*(1+X%*%theta) + A*lambda*((X%*%g)^3) + rnorm(n.samp,0,1))
    mu_X = (X%*%alpha)*(1+X%*%theta); CX = (X%*%g)^3
  }else if(model == "CNL3"){
    Y = as.vector((X%*%theta)^3 + A*lambda*((X%*%g)^3) + rnorm(n.samp,0,1))
    mu_X = (X%*%theta)^3; CX = (X%*%g)^3
  }else if(model == "L2S"){
    Y = as.vector(X%*%alpha*(1+X%*%theta) + A*lambda*(sin(X%*%beta) ) + rnorm(n.samp,0,1))
    mu_X = (X%*%alpha)*(1+X%*%theta); CX = sin(X%*%beta)
  }else if(model == "CS"){
    Y = as.vector((X%*%theta)^3 + A*lambda*( sin(X%*%beta) ) + rnorm(n.samp,0,1))
    mu_X = (X%*%theta)^3; CX = sin(X%*%beta)
  }
  
  
  Yt.trueP = Y*(A - true.prop) / (true.prop*(1 - true.prop))
  Yt = Yt.trueP
  
  
  formula.linReg = as.formula(paste("CX~", paste(colnames(X), collapse="+")))
  regM = lm(formula.linReg, data=as.data.frame(X))
  
  
  V_0 = mu_X + ifelse(CX>0,1,0)*CX
  V_b0 = mu_X + ifelse(cbind(rep(1,dim(X)[1]),X)%*%regM$coefficients>0,1,0)*CX
  
  ret_list = list(coefs = regM$coefficients , ux = mu_X , X = X, CX=CX,V_0=V_0,V_b0=V_b0 )
  return(ret_list)
}

#########################################################################################

#Simulation function

#########################################################################################


# example input: n = 500;N=5000;p=2;mu_0 = 1;g=1;g2=1;lambda=1;model = "CL";beta = c(1,1);num.sims=50

Sim_OTR <- function(n , N , p=2, model = "SL" ,lambda=1,beta=NULL,g=NULL,num.sims=500){
  
  theta = rep(c(0.5,0.5),floor(p/2))
  alpha = rep(c(0.75,0.75),floor(p/2))
  
  MSE.ols = rep(0,num.sims)
  MSE.np1 = rep(0,num.sims)
  MSE.np2 = rep(0,num.sims)
  
  MSE.ols_b0 = rep(0,num.sims)
  MSE.ols_b1 = rep(0,num.sims)
  MSE.ols_b2 = rep(0,num.sims)
  
  MSE.CXls_b0 = rep(0,num.sims)
  MSE.CXls_b1 = rep(0,num.sims)
  MSE.CXls_b2 = rep(0,num.sims)
  
  #Store beta coefficients.
  Beta_mat_OLS = matrix(0 , ncol=p+1 , nrow = num.sims)
  
  Beta_mat1 = matrix(0 , ncol=p+1 , nrow = num.sims)
  imput_est_mat1.CX = matrix(0 , ncol = 4 , nrow= num.sims)
  
  Beta_mat2 = matrix(0 , ncol=p+1 , nrow = num.sims)
  imput_est_mat2.CX = matrix(0 , ncol = 4 , nrow= num.sims)
  
  #Store ASE 
  IF.CXls_mat = matrix(0 , ncol=p+1 , nrow = num.sims)
  IF.RF_mat = matrix(0 , ncol=p+1 , nrow = num.sims)
  IF.OLS_mat = matrix(0 , ncol=p+1 , nrow = num.sims)
  #CP
  CP_CX.mat = matrix(0 , ncol=p+1 , nrow = num.sims)
  CP_CX.mat2 = matrix(0 , ncol=p+1 , nrow = num.sims)
  CP_CX.mat3 = matrix(0 , ncol=p+1 , nrow = num.sims)
  
  #pcd
  pcd_sim.RF = rep(0,num.sims)
  pcd_sim.CXLS = rep(0,num.sims)
  pcd_sim.OLS = rep(0,num.sims)
  # Value functions
  VF.RF = rep(0,num.sims)
  VF.CXLS = rep(0,num.sims)
  VF.OLS = rep(0,num.sims)
  VF.true = rep(0,num.sims)
  
  #Get beta.t
  MC_list = MC_coeff(model = model, n.samp = 500000,p=p,
                     g=g,theta=theta,lambda=lambda,beta=beta,alpha = alpha)
  beta.t = MC_list$coefs
  print(beta.t)
  
  for(sim in 1:num.sims){
    print(sim)
    X = mvrnorm(n=n+N , mu = rep(0,p), Sigma = diag(1,p))
    X.m = apply(X,1,function(x){ifelse(x > 5,5,x)})
    X = apply(X.m,1,function(x){ifelse(x < -5,-5,x)})
    
    colnames(X) <- paste("X",1:p,sep="")
    
    #Impute incomplete cases with Kernel Smoothing and Random Forests.
    X.imp = as.data.frame(X[(n+1):(N+n),]); X.label = as.data.frame(X[1:n,])
    
    logit_coefs <- rep(1,ncol(X))
    true.prop = plogis( 0.5*X[,1] - 0.5*X[,2] )
    # Treatment Allocation.
    A = sapply(true.prop, function(x){ rbinom(1,1,x) })
    
    #Generate Y based on the given true model.
    if(model == "CL"){
      Y = as.vector((X%*%theta)^3 + A*lambda*(X%*%beta) + rnorm(n+N,0,1))
      CX = X%*%beta
    }else if(model == "L2L"){
      Y = as.vector( X%*%alpha*(1+X%*%theta) + A*lambda*(X%*%beta) + rnorm(n+N,0,1))
      CX = X%*%beta
    }else if(model == "L2NL3"){
      Y = as.vector(X%*%alpha*(1+X%*%theta) + A*lambda*((X%*%g)^3) + rnorm(n+N,0,1))
      CX = (X%*%g)^3
    }else if(model == "CNL3"){
      Y = as.vector((X%*%theta)^3 + A*lambda*((X%*%g)^3) + rnorm(n+N,0,1))
      CX = (X%*%g)^3
    }else if(model == "L2S"){
      Y = as.vector(X%*%alpha*(1+X%*%theta) + A*lambda*(sin(X%*%beta) ) + rnorm(n+N,0,1))
      CX = sin(X%*%beta)
    }else if(model == "CS"){
      Y = as.vector((X%*%theta)^3 + A*lambda*( sin(X%*%beta) ) + rnorm(n+N,0,1))
      CX = sin(X%*%beta)
    }
    
    #need these to get propensity score.
    X.known = X[1:n,]
    A.known = A[1:n]
    Y.known = Y[1:n]
    
    prop = true.prop[1:n]
    
    #values for response variable in transformed outcome linear regression.
    Yt = Y.known*(A.known - prop) / (prop*(1 - prop))
    
    # Estimate Y_A1 and Y_A0 separately. 
    Y_A1 = Y.known[which(A.known==1)];Y_A0 = Y.known[which(A.known==0)];prop0 = prop[which(A.known==0)]
    X.labelA1 = X.label[A.known==1,];X.labelA0 = X.label[A.known==0,];prop1 = prop[which(A.known==1)]
    
    ###################################################################################################
    # Random Forest approach
    
    #Predict missing Y.
    formula.YA1 = as.formula(paste("Y_A1~", paste(colnames(X.label), collapse="+")))
    formula.YA0 = as.formula(paste("Y_A0~", paste(colnames(X.label), collapse="+")))
    
    #Random Forest Model
    RF.model_A1 = randomForest(formula.YA1 , data = X.labelA1)
    RF.model_A0 = randomForest(formula.YA0 , data = X.labelA0)
    mx.hat_A1 = predict(RF.model_A1,newdata = X.imp)
    mx.hat_A0 = predict(RF.model_A0,newdata = X.imp)
    CX_RF = mx.hat_A1 - mx.hat_A0
    
    ##################################################################################################
    
    #Kernel Smoothing approaches.
    X.l_1 = cbind(rep(1,dim(X.labelA1)[1]),X.labelA1);X.l_0 = cbind(rep(1,dim(X.labelA0)[1]),X.labelA0)
    X.i_KS = cbind(rep(1,N),X.imp)
    X.l_KS = cbind(rep(1,n),X.label)
    
    #Get h by grid search
    if(sim==1){
      hc1 = cv.h(Y_A1, X.l_1,seq.c = seq(0.1,5,0.2)) #0.75
      hc0 = cv.h(Y_A0, X.l_0,seq.c = seq(0.1,5,0.2))
      print(hc1);print(hc0)
    }
    
    #Predictions with bias reduction by least squares
    mx.hat_ls1 = hlscv.ks(Y_A1,x.in= X.l_1,x.impute=X.i_KS,const_in=hc1, prop_score = prop1 )[[1]]
    mx.hat_ls0 = hlscv.ks(Y_A0 , X.l_0, X.i_KS, const_in = hc0, prop_score = 1 - prop0 )[[1]]
    CX_ls = mx.hat_ls1 - mx.hat_ls0
    
    ####################################################################################################
    #Regression predicitions.
    formula.Yt = as.formula(paste("Yt~", paste(colnames(X.label), collapse="+")))
    Reg.OLS = lm(formula.Yt , data = X.label)
    beta.OLS = Reg.OLS$coefficients
    
    formula.CX_RF  = as.formula(paste("CX_RF ~", paste(colnames(X.label), collapse="+")))
    Reg.CX_RF = lm(formula.CX_RF  , data = X.imp)
    beta.CX_RF = Reg.CX_RF$coefficients
    
    formula.CX_ls  = as.formula(paste("CX_ls ~", paste(colnames(X.label), collapse="+")))
    Reg.CX_ls = lm(formula.CX_ls  , data = X.imp)
    beta.CX_ls = Reg.CX_ls$coefficients
    
    # Relative Efficiency
    MSE.ols[sim] = sum((beta.OLS - beta.t)^2)
    MSE.np1[sim] = sum((beta.CX_RF - beta.t)^2)
    MSE.np2[sim] = sum((beta.CX_ls - beta.t)^2)
    
    # Individual regression coefficients RE.
    MSE.ols_b0[sim] = sum((beta.OLS[1] - beta.t[1])^2)
    MSE.ols_b1[sim] = sum((beta.OLS[2] - beta.t[2])^2)
    MSE.ols_b2[sim] = sum((beta.OLS[3] - beta.t[3])^2)
    
    MSE.CXls_b0[sim] = sum((beta.CX_ls[1] - beta.t[1])^2)
    MSE.CXls_b1[sim] = sum((beta.CX_ls[2] - beta.t[2])^2)
    MSE.CXls_b2[sim] = sum((beta.CX_ls[3] - beta.t[3])^2)
    
    # Record regression coefficients
    Beta_mat_OLS[sim,] = beta.OLS
    Beta_mat1[sim,] = beta.CX_RF
    Beta_mat2[sim,] = beta.CX_ls
     
    #Record imputation value
    imput_est_mat1.CX[sim,] = c(Q1_bias=quantile(CX_RF)[2]-quantile(CX[-(1:n)])[2] ,mean_bias=mean(CX_RF) - mean(CX[-(1:n)]) , 
                                Q3_bias=quantile(CX_RF)[4]-quantile(CX[-(1:n)])[4],MSE=mean((CX_RF - CX[-(1:n)])^2))
    imput_est_mat2.CX[sim,] = c(Q1_bias=quantile(CX_ls)[2]-quantile(CX[-(1:n)])[2] ,mean_bias=mean(CX_ls) - mean(CX[-(1:n)]) , 
                                Q3_bias=quantile(CX_ls)[4]-quantile(CX[-(1:n)])[4],MSE=mean((CX_ls - CX[-(1:n)])^2))
    
    ##########################################################################
    # Variance Estimation
    # First get double cv mu_hat estimates
    dcv.mx1 = double_cv.ks1(Y_A1,x.in= X.l_1,const_in=hc1, prop_score = prop1 )
    dcv.mx0 = double_cv.ks0(Y_A0 , X.l_0, const_in = hc0, prop_score = 1 - prop0 )
    
    n1 = length(Y_A1)
    n0 = length(Y_A0)
    
    full_X.imp = as.matrix(cbind(rep(1,N),X.imp))
    names(X.l_1) <- names(X.l_0) 
    X.l_1mat = as.matrix(X.l_1);X.l_0mat = as.matrix(X.l_0)
    
    X.KS = as.matrix(X.i_KS)
    XtX_sum= matrix(0, ncol=p+1, nrow = p+1)
    for(a in 1:dim(X.KS)[1]){
      XtX = X.KS[a,] %*% t(X.KS[a,])
      XtX_sum = XtX_sum + XtX
    }
    GammaInv = solve(XtX_sum / dim(X.i_KS)[1])
    
    #print(GammaInv)
    
    full_XL = as.matrix(cbind(rep(1,n),X.label))
    XtXn_sum = matrix(0, ncol=p+1, nrow = p+1)
    for(a in 1:dim(full_XL)[1]){
      XtXn = full_XL[a,] %*% t(full_XL[a,])
      XtXn_sum = XtXn_sum + XtXn
    }
    GammaInv_n = solve(XtXn_sum / dim(full_XL)[1])
    
    X.in = rbind( as.matrix(X.l_1) , as.matrix(X.l_0) );Y.in = c(Y_A1,Y_A0);prop.in = c(prop1,prop0)
    mv.input = c(dcv.mx1,dcv.mx0); A.in = c( rep(1,length(Y_A1)) , rep(0,length(Y_A0)) ) 
    RF.in = c(RF.model_A1$predicted,RF.model_A0$predicted)
   
    ####################################################################################################
    
    print(beta.CX_ls)
    IF.CXls_sd = IF_se(GammaInv, X.in, Y.in, A.in, mv.input,prop.in )
    print(IF.CXls_sd)
    
    IF.CXls_mat[sim,] = IF.CXls_sd
    
    IF.RF_sd = IF_se(GammaInv, X.in, Y.in, A.in, RF.in, prop.in )
    IF.RF_mat[sim,] = IF.RF_sd
    print(beta.CX_RF)
    print(IF.RF_sd)
    
    X.L = as.matrix(cbind(rep(1,length(Yt)),X.label))
    IF.OLS_sd = IF_OLS_se(GammaInv_n, X.L, Yt=Yt, bx.in = Reg.OLS$fitted.values )
    IF.OLS_mat[sim,] = IF.OLS_sd
    
    ##############################################################
    #Coverage probability
    #OLS
    z = qnorm(0.975)
    #CX_ls RAL estimator
    lower.cx = beta.CX_ls - z*IF.CXls_sd
    upper.cx = beta.CX_ls + z*IF.CXls_sd
    CP_CX.mat[sim,] = 1*((lower.cx <= beta.t) & (beta.t <= upper.cx))
    print(1*((lower.cx <= beta.t) & (beta.t <= upper.cx)))
    
    lower.cx_RF = beta.CX_RF - z*IF.RF_sd
    upper.cx_RF = beta.CX_RF + z*IF.RF_sd
    CP_CX.mat2[sim,] = 1*((lower.cx_RF <= beta.t) & (beta.t <= upper.cx_RF))
    print( 1*((lower.cx_RF <= beta.t) & (beta.t <= upper.cx_RF)))
    
    lower.ols = beta.OLS - z*IF.OLS_sd
    upper.ols = beta.OLS + z*IF.OLS_sd
    CP_CX.mat3[sim,] = 1*((lower.ols <= beta.t) & (beta.t <= upper.ols))
    print( 1*((lower.ols <= beta.t) & (beta.t <= upper.ols)))
    
    ###############################################################
    #pcd for the current simulation.
    int_X = cbind(rep(1,dim(X)[1]),X)
    BX_RF = beta.CX_RF%*%t(int_X)
    BX_CXLS = beta.CX_ls%*%t(int_X)
    BX_OLS = beta.OLS%*%t(int_X)
    BX_true = beta.t%*%t(int_X)
    
    Ind_RF = sapply(BX_RF,function(x)ifelse((x>0),1,0))
    Ind_CXLS = sapply(BX_CXLS,function(x)ifelse((x>0),1,0))
    Ind_OLS = sapply(BX_OLS,function(x)ifelse((x>0),1,0))
    Ind_true = sapply(BX_true,function(x)ifelse((x>0),1,0))
    
    pcd_sim.RF[sim] = 1 - mean(abs(Ind_RF - Ind_true))
    pcd_sim.CXLS[sim] = 1 - mean(abs(Ind_CXLS - Ind_true))
    pcd_sim.OLS[sim] = 1 - mean(abs(Ind_OLS - Ind_true))
    
    #VF 
    X_MC = MC_list$X
    int_MCX = cbind(rep(1,dim(X_MC)[1]),X_MC)
    MC_RF = beta.CX_RF%*%t(int_MCX)
    MC_CXLS = beta.CX_ls%*%t(int_MCX)
    MC_OLS = beta.OLS%*%t(int_MCX)
    
    
    Ind_MC_RF = sapply(MC_RF,function(x)ifelse((x>0),1,0))
    Ind_MC_CXLS = sapply(MC_CXLS,function(x)ifelse((x>0),1,0))
    Ind_MC_OLS = sapply(MC_OLS,function(x)ifelse((x>0),1,0))
    
    VF.RF[sim] = mean(as.vector(MC_list$ux + Ind_MC_RF*lambda*MC_list$CX))
    VF.CXLS[sim] = mean(as.vector(MC_list$ux + Ind_MC_CXLS*lambda*MC_list$CX))
    VF.OLS[sim] = mean(as.vector(MC_list$ux + Ind_MC_OLS*lambda*MC_list$CX))
  }
  
  #pcd
  pcd.RF.mean = mean(pcd_sim.RF);pcd.RF.sd = sd(pcd_sim.RF)
  pcd.CXLS.mean = mean(pcd_sim.CXLS);pcd.CXLS.sd = sd(pcd_sim.CXLS)
  pcd.OLS.mean = mean(pcd_sim.OLS);pcd.OLS.sd = sd(pcd_sim.OLS)
  #VF
  VF.RF.mean = mean(VF.RF);VF.RF.sd = sd(VF.RF)
  VF.CXLS.mean = mean(VF.CXLS);VF.CXLS.sd = sd(VF.CXLS)
  VF.OLS.mean = mean(VF.OLS);VF.OLS.sd = sd(VF.OLS )
  
  #mean
  betaOLS.mean = apply(Beta_mat_OLS,2,mean)
  beta.CX_RF.mean = apply(Beta_mat1,2,mean)
  beta.CX_ls.mean = apply(Beta_mat2,2,mean)
  
  #Empirical sd
  betaOLS.sd = apply(Beta_mat_OLS,2,sd)
  beta.CX_RF.sd = apply(Beta_mat1,2,sd)
  beta.CX_ls.sd = apply(Beta_mat2,2,sd)
  
  #Average IF/boot sd
  IF.CXls_avg = apply(IF.CXls_mat,2,mean)
  
  #Average IF/boot sd
  IF.RF_avg = apply(IF.RF_mat,2,mean)
  
  IF.OLS_avg = apply(IF.OLS_mat,2,mean)
  
  #CP
  CP_CX = apply(CP_CX.mat,2,mean)
  CP_CX2 = apply(CP_CX.mat2,2,mean)
  CP_CX3 = apply(CP_CX.mat3,2,mean)
  
  #Relative Efficiency
  RE.RF = mean(MSE.ols)/mean(MSE.np1)
  
  RE.CX_ls = mean(MSE.ols) / mean(MSE.np2)
  
  # RE of individual coefficients.
  RE.CX_ls_b0 = mean(MSE.ols_b0) / mean(MSE.CXls_b0)
  RE.CX_ls_b1 = mean(MSE.ols_b1) / mean(MSE.CXls_b1)
  RE.CX_ls_b2 = mean(MSE.ols_b2) / mean(MSE.CXls_b2)
  
  RE_coefs = c(RE.CX_ls_b0, RE.CX_ls_b1, RE.CX_ls_b2)
  ##########################################################################################
  print(hist(Beta_mat2[,1]))
  print(hist(Beta_mat2[,2]))
  print(hist(Beta_mat2[,3]))
  
  RF_2_imput = as.matrix(apply(imput_est_mat1.CX,2,mean),ncol=4)
  rownames(RF_2_imput)<-c("Q1_bias","mean_bias","Q3_bias","MSE")
  CXls_imput = as.matrix(apply(imput_est_mat2.CX,2,mean),ncol=4)
  rownames(CXls_imput)<-c("Q1_bias","mean_bias","Q3_bias","MSE")
  
  RE = list( RE_RF = RE.RF,
             RE_CXls = RE.CX_ls,
             RF_2_imput = RF_2_imput,
             CXls_imput = CXls_imput,
             bias_sd_OLS = cbind(beta.t,bias=betaOLS.mean-beta.t,Esd=betaOLS.sd,
                                 Asd = IF.OLS_avg,CP = CP_CX3),
             bias_sd_RF = cbind(beta.t,bias=beta.CX_RF.mean-beta.t ,sd=beta.CX_RF.sd),
             bias_sd_CXls = cbind(beta.t,bias=beta.CX_ls.mean-beta.t ,Esd=beta.CX_ls.sd,Asd_IF = IF.CXls_avg,CP=CP_CX,
                                  RE_beta = RE_coefs),
             pcd_OLS = c(pcd.OLS.mean,pcd.OLS.sd),
             VF_OLS = c(VF.OLS.mean,VF.OLS.sd),
             pcd_RF = c(pcd.RF.mean, pcd.RF.sd),
             VF_RF = c(VF.RF.mean,VF.RF.sd),
             pcd_CXls = c(pcd.CXLS.mean, pcd.CXLS.sd),
             VF_Cxls = c(VF.CXLS.mean,VF.CXLS.sd),
             VF_0 = c(mean_V0 = mean(MC_list$V_0),sd_V0 = sd(MC_list$V_0)),
             VF_B0 = c(mean_VB0 = mean(MC_list$V_b0),sd_VB0 = sd(MC_list$V_b0)),
             PCD_vals_RF = pcd_sim.RF,
             PCD_vals_CXLS = pcd_sim.CXLS,
             PCD_vals_oLS = pcd_sim.OLS,
             V_RF = VF.RF,
             V_CXLS = VF.CXLS,
             V_OLS = VF.OLS,
             Beta_CX = Beta_mat2
             
  )
  
  return(RE)
}

#p=10
beta.in = rep(1,2);g.in = c(0.3,0.6)
###################################################################################

Sim_CL = Sim_OTR(n=500,N=5000, p=2, model = "CL",beta = beta.in,num.sims=500)

Sim_L2L = Sim_OTR(n=500,N=5000, p=2, model = "L2L",beta = beta.in,num.sims = 500)

Sim_L2NL3 = Sim_OTR(n=500,N=5000, p=2, model = "L2NL3",beta=beta.in,g = g.in,num.sims=500)

Sim_CNL3 = Sim_OTR(n=500,N=5000, p=2, model = "CNL3",beta=NULL,g = g.in,num.sims= 500)

Sim_L2S = Sim_OTR(n=500,N=5000, p=2, model = "L2S",beta=beta.in, num.sims=500)

Sim_CS = Sim_OTR(n=500,N=5000, p=2, model = "CS",beta=beta.in, num.sims=500)

capture.output(Sim_L2L,file="MCAR_L2L_500.txt")
capture.output(Sim_CL,file="MCAR_CL_500.txt")
capture.output(Sim_L2NL3,file="MCAR_L2NL3_500.txt")
capture.output(Sim_CNL3,file="MCAR_CNL3_500.txt")
capture.output(Sim_L2S,file="MCAR_L2S_500.txt")
capture.output(Sim_CS,file="MCAR_CS_500.txt")


