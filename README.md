# AMI

[AMI](https://github.com/jfaucon/AMI) (Another Maxent Implementation) is an R implementation of multinomial logistic regression, also known as Maximum Entropy classifier. This implementation deals with binary and real-valued features and uses standard R function 'optim()' to maximize the objective function. It is possible to use different iterative methods to compute lambda : BFGS, LM-BFGS-B, Conjugate Gradient (CG) and Generalized Iterative Scaling (GIS). This R package can be modified under the terms of the GNU GPL.


## Usage

* General usage:

        maximumentropy(x, y=NULL, data=NULL, iteration=NULL, 
        method=c("L-BFGS-B", "GIS", "CG", "BFGS"), verbose=TRUE, normalize=FALSE)


* Example with formula object:

        model <- maximumentropy(Species ~ ., data=iris, iteration=100)
        predict(model, new_dataset)
        

* Example with data.frame or matrix:

        x <- subset(iris, select=-Species)
        y <- iris$Species
        model <- maximumentropy(x, y)


* Prediction and evaluation:

        predicted_y <- predict(model, new_dataset)
        evaluate(true_y, predicted_y) 
                
        
* See information about a model:

        summary(model) #summary of model
        # or
        model$value #final value of log-likelihood 
        model$lambda #computed lambda 
        
        model$levels #discrete classes Y
        model$features #names of features x
        
        model$method #name of iterative method
        model$iteration #number of iterations 
        model$convergence #0=converged, 1=max number of iterations reached, ... (see ?optim)
        

 
        
        
        
## Comments

* By default, the MaxEnt is trained with LM-BFGS-B. The GIS is the slower of all iterative methods.

* By default, the maximization of the objective function is done until convergence.

* The y vector must be a factor vector.

* It is possible to normalize by scaling data between 0 and 1 ('normalize' argument).



## Details

* Actual version: 0.4
* Things to do:s
    * Write a documentation
    * Add L1 & L2 regularization
    * Add cut-off for rare events
    * Add feature selection
    * Add automatic parallelization
