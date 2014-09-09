#' Multinomial Logistic Regression with binary and real-valued features.
#' 
#' @param 
#' @param y is factor vector with classes and the same length of x.
#' @param method is the hill-climbing method used to optimize the function.
#' @return A model to predict.
#' @examples


maximumentropy <- function(x, y=NULL, data=NULL, iteration=NULL, verbose=TRUE, method=c("L-BFGS-B", "GIS", "CG", "BFGS"), addslack=FALSE, normalize=FALSE){
  
  # Using formula
  if(class(x) == "formula" && !is.null(data)){
    data <- as.data.frame(data)
    df_data <- model.frame(x, data = data)
    y <- df_data[[1]]
    x <- df_data[2:ncol(df_data)]
  }
  # if addslack feature 
  if(addslack){
    x$slack <- rep(1,nrow(x))
  }
  # if normalize
  if(normalize){
    x <- normalize_x(x);
  }
  y <- as.factor(y)
  x <- as.matrix(x);
  
  
  if(verbose){cat("Treat data ... ","\n")}
  
  features_names <- colnames(x);
  
  nb_observations <- nrow(x); #N
  nb_classes <- length(levels(y)); #C
  nb_features <- ncol(x); #M
  featureSize <- nb_classes*nb_features;
  
  if(verbose){
  cat("Number of classes",nb_classes,"\n")
  cat("Number of features",nb_features,"\n")
  cat("Number of observations",nb_observations,"\n")
  cat("featureSize",featureSize,"\n")
  }
  
  # Generate Lambda
  pars <- rep(1,featureSize)
  obsCounts <- getObsCounts(x,y,verbose)

  if(verbose){
    cat("Training with",method[1],"...\n")
  }
  if(method == "L-BFGS-B" || method == "BFGS" || method == "CG"){
    model <- 1;
    if(is.null(iteration)){
      model <- optim(pars,objectiveFunction,getGradients,x=x,y=y,obsCounts=obsCounts,method=method[1],control=list(fnscale=-1,trace=verbose,REPORT=1,maxit=10000))
    }
    else{
      model <- optim(pars,objectiveFunction,getGradients,x=x,y=y,obsCounts=obsCounts,method=method[1],control=list(fnscale=-1,trace=verbose,REPORT=1,maxit=iteration))
    }
    model <- list("levels" = levels(y), "lambda"=model$par, "value" = model$value, "features" = features_names, "convergence"=model$convergence, "method" = method[1],"iteration" = model$counts[2])
  }
  else if(method == "GIS"){
    GIS_object <- GIS(x=x,y=y,obsCounts=obsCounts, iteration=iteration, verbose=verbose)
    model <- list("levels" = levels(y), "lambda"=GIS_object$lambda, "value"=GIS_object$value, "features" = features_names, "convergence" = GIS_object$convergence, "method"= method[1], "iteration" = GIS_object$iteration)
  }
  class(model) <- "maximumentropy"
  return(invisible(model))
}


GIS <- function(x=x,y=y,obsCounts=obsCounts,iteration=NULL,verbose=FALSE){
  nb_observations <- nrow(x); #N
  nb_classes <- length(levels(y)); #C
  nb_features <- ncol(x); #M
  featureSize <- nb_classes*nb_features;
  lambda <- matrix(rep(0,featureSize),ncol=nb_features)
  C <- max(margin.table(obsCounts,1))
  #C <- max(obsCounts)
  # the larger the value of C, the smaller the step size. (Malouf, 2002)
  cat("C:",C,"\n")
  
  #--- Iteration ---#
  currIteration <- 1;
  convergence <- 1;
  previousLog <- 0;
  flag <- 1;
  diff_log <- 1;
  while( ((currIteration <= iteration) || is.null(iteration)) && flag !=0 ){
    
    log_likelihood = objectiveFunction(lambda,x,y,obsCounts)
    
    # Step 9 : getExpectedCount
    conditional_probability <- getConditionalProbability(x,y,lambda)
    expCounts <- getExpCounts(x, y, conditional_probability);

    # Step 8 : update Lambda
      
      new_lambda <- matrix(rep(0,featureSize),ncol=nb_features)
      for(i in 1:nb_classes){
        for(j in 1:nb_features){
          if(obsCounts[i,j] == 0){
            new_lambda[i,j] <- 0;
          }else{
            new_lambda[i,j] <- lambda[i,j] + ((log(obsCounts[i,j]) - log(expCounts[i,j])) * (1/C))
#             new_lambda[i,j] <- lambda[i,j] + (log(obsCounts[i,j]/expCounts[i,j]) * (1/C))
          }
        }
      }
    lambda <- new_lambda;
  
    if(currIteration == 1){
      previousLog <- log_likelihood;
    }
    else{
      diff_log <- abs(log_likelihood - previousLog)/abs(log_likelihood);
      if( diff_log < 1e-7 ){
        if(verbose){        
          cat("Convergence!")
        }
        flag <- 0;
        convergence <- 0;
      }
      
      if(previousLog > log_likelihood){
        flag <- 0;
        convergence <- 53;
        cat("Divergence!\n")
      }
      previousLog <- log_likelihood;
      
    }
    if(verbose){
      cat(currIteration,"logLikelihood:",log_likelihood,"difflog",diff_log, "\n")
    }
  currIteration <- currIteration +1;
  }
  #-- endOfIteration --# 
  GIS_object <- list("lambda"=lambda,"value"=log_likelihood, "convergence"=convergence, "iteration" = currIteration)
}


objectiveFunction <- function(pars,x,y,obsCounts){
  
  pars <- unlist(pars)
  loglikelihood <- getConditionalProbabilityObjective(x,y,pars)
  loglikelihood
}

getConditionalProbabilityObjective <- function(x,y,lambda){
  
  nb_features <- ncol(x); #M
  if(class(lambda) != "matrix"){
    lambda <- matrix(lambda,ncol=nb_features)
  }
  # Step 1 : some variables
  nb_classes <- length(levels(y)); #C
  nb_features <- ncol(x); #M
  nb_observations <- nrow(x); #N
  featureSize <- nb_classes * nb_features;
  x <- as.matrix(x)
  
  # Step 4 : Calculte first term
  # Iteration on Training Set
  inner_product_size <- nb_observations * nb_classes;
  inner_product <- matrix(rep(0,inner_product_size),ncol=nb_classes)
  for(currInstance in 1:nb_observations){
    for(currClass in 1:nb_classes){
      inner_product[currInstance,currClass] <- lambda[currClass,] %*% x[currInstance,]
    }
  }
  
    # Step 5 : Calculte second term (LogSumExp trick)
    Z_constante <- numeric(nb_observations)
    for(i in 1:nrow(inner_product)){
      A <- max(inner_product[i,])
      A <- A - 700;
      if(A < 0){
        A <- 0;
      }
      else{
      }
      tmp_total <- 0;
      for(j in 1:nb_classes){
        tmp_total <- tmp_total + exp(inner_product[i,j] - A)
      }
      Z_constante[i] <- log(tmp_total)+A;
    }
  
  # Step 6 : get Model Probability
  conditional_probability_size <- nb_observations * nb_classes;
  conditional_probability = matrix(rep(0,conditional_probability_size),ncol=nb_classes)
  
  log_likelihood <- 0;
  for(i in 1:nrow(conditional_probability)){
    log_likelihood <-  log_likelihood + (inner_product[i,y[i]] - Z_constante[i]);
  }
  log_likelihood
  
}

getConditionalProbability <- function(x,y,lambda){
  
  nb_features <- ncol(x); #M
  if(class(lambda) != "matrix"){
    lambda <- matrix(lambda,ncol=nb_features)
  }
  # Step 1 : some variables
  nb_classes <- length(levels(y)); #C
  nb_features <- ncol(x); #M
  nb_observations <- nrow(x); #N
  featureSize <- nb_classes * nb_features;
  x <- as.matrix(x)
  
  # Step 4 : Calculte first term
  #Iteration on Training Set
  inner_product_size <- nb_observations * nb_classes;
  inner_product <- matrix(rep(0,inner_product_size),ncol=nb_classes)
  for(currInstance in 1:nb_observations){
    for(currClass in 1:nb_classes){
      inner_product[currInstance,currClass] <- lambda[currClass,] %*% x[currInstance,]
    }
  }
  
  # Step 6 : get Model Probability
  conditional_probability_size <- nb_observations * nb_classes;
  conditional_probability = matrix(rep(0,conditional_probability_size),ncol=nb_classes)
  
  for(i in 1:nrow(conditional_probability)){
    # Avoid overflow and underflow.
    # See implementation of Tsuakora (MaxEnt C++)
    pmax <- max(inner_product[i,]);
    sum <- 0;
    offset <- max(c(0.0, pmax - 700));
    for(j in 1:nb_classes){

      pow <- inner_product[i,j] - offset;
      prod <- exp(pow)
      conditional_probability[i,j] <- exp(pow)
      sum <- sum + prod;
#     conditional_probability[i,j] <- exp(inner_product[i,j] - Z_constante[i]);
    }
    conditional_probability[i,] <- conditional_probability[i,] / sum;
  }
  

  vec <- sum(margin.table(conditional_probability,1));
  vec <- round(vec,digits=2)
  if(vec != as.numeric(nb_observations) ){
    cat("Erreur : ","\n")
    print(sprintf("observations %d\n", as.numeric(nb_observations) ))
    print(sprintf("vec %f\n", vec ))
    cat(class(vec))
    cat(class(as.numeric(nb_observations)))
    visualize(conditional_probability)
  }
  conditional_probability
}



getGradients <- function(pars,x=x,y=y,obsCounts=obsCounts){
  
  conditional_probability <- getConditionalProbability(x,y,pars)
  expCounts <- getExpCounts(x, y, conditional_probability)
  gradients <- obsCounts - expCounts;
  gradients <- as.vector(gradients)
  
  gradients
  
}

getExpCounts <- function(x, y, conditional_probability){
  nb_classes <- length(levels(y)); #C
  nb_features <- ncol(x); #M
  featureSize <- nb_classes * nb_features;
  
  # ExpCount Matrix : C * M
  expCounts <- matrix(rep(0,featureSize),ncol=nb_features)
  
  for(i in 1:nb_classes){
    for(j in 1:nb_features){
       expCounts[i,j] <-  conditional_probability[,i] %*% x[,j];
       if(is.na(expCounts[i,j]) || is.infinite(expCounts[i,j])){
         cat("ERROR IN GRADIENT")
       }
    }
  }
  expCounts
}


getObsCounts <- function(x,y,verbose=TRUE){
  
  # TODO : table margins
  if(verbose){
    cat("Count obsCounts ...", "\n")
  }
  nb_classes <- length(levels(y)); #C
  nb_features <- ncol(x); #M
  featureSize <- nb_classes * nb_features;
  
  # obsCount Matrix : C * M
  obsCount <- matrix(rep(0,featureSize),ncol=nb_features)
  
  for(i in 1:nrow(x)){
    current_class <- y[i];
    for(j in 1:ncol(x) ){
      obsCount[current_class,j] <- obsCount[current_class,j] + x[i,j]
    }
  }
  obsCount;
}



summary.maximumentropy <- function(model){
  cat("Class:", class(model))
  cat("\n\n")
  cat("Method:", model$method,"\n")
  cat("Value: ",model$value, "\n")
  cat("Convergence:", model$convergence,"\n")
  cat("Iterations:", model$iteration,"\n")
  cat("\n")
  cat("Lambda:\n")
  featureSize = length(model$features) * length(model$levels)
  lambda <- matrix(model$lambda,ncol=length(model$features))
  features_names <- model$features;
  rownames(lambda) <- model$levels;
  colnames(lambda) <- features_names;
  print(lambda)
  cat("\n")
  cat("Levels:\n");
  print(model$levels);
  cat("\n")
}

predict.maximumentropy <- function (model, x)  {
  lambda <- model$lambda;
  nb_classes <- length(model$levels)
  internal_levels <- model$levels;
  nb_features <- length(lambda) / nb_classes;
  
  if(ncol(x) < nb_features){
    cat("Add Slack feature")
    x$slack <- rep(1,nrow(x)) #slack features
  }
  
  x <- as.matrix(x);
  y <- as.factor(y);
  nb_observations <- nrow(x); #N
  if(class(lambda) != "matrix"){
    lambda <- matrix(lambda,ncol=nb_features)
  }
  conditional_probability <- getConditionalProbability(x,y,lambda)
  
  predict_y <- factor();
  for(i in 1:nrow(conditional_probability)){
    
    argmax_p <- which.max( conditional_probability[i,] )
    predict_y <- c(predict_y,internal_levels[argmax_p])
  }
  predict_y <- as.factor(predict_y)
  return(invisible(predict_y));
}


evaluate <- function(y, predict_y){
  acc <- 0;
  for(i in 1:length(y)){
    if(y[i] == predict_y[i]){
      acc <- acc +1;
    }
  }
  acc <- (acc / length(y)) * 100
  cat('Accuracy',acc)
}


normalize_x <- function(x){
  # Step 1 : if real values, features scaling.
  x <- as.matrix(x)
  cat("Feature Scaling ...","\n")
  for(i in 1:ncol(x)){
    currVector <- x[, i]
    if(max(currVector) - min(currVector) != 0){
      cat("Diff de 0");
      currVector <- (currVector - min(currVector))/( max(currVector) - min(currVector) )
    }
    x[, i] <- currVector;
    cat("feature ",i,"/",ncol(x),"\n")
  }
  x <- as.matrix(x,row.names=FALSE);
  cat("End of Feature Scaling","\n")
  x
}


# for_new <- function(new,y){
#   acc <- 0;
#   new <- new[,2:ncol(new)]
#   for(i in 1:nrow(new)){
#     argmax_p <- which.max( new[i,] )
#     if(argmax_p == as.integer(y[i])){
#       acc <- acc+1;
#     }
#   }
#   acc <- (acc/150)*100;
#   acc <- round(acc,digits=2)
#   cat("Accuracy",acc)
# }


