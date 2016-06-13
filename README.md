# AMI

[AMI](https://github.com/jfaucon/AMI) (Another Maxent Implementation) is a simple R implementation of multinomial logistic regression, also known as Maximum Entropy classifier. This implementation deals with binary and real-valued features and uses standard R function `optim(.)` to maximize the objective function. It is possible to use different iterative methods to compute lambda : BFGS, Conjugate Gradient (CG) and Generalized Iterative Scaling (GIS). This R package can be modified under the terms of the GNU GPL.

Please, note that this package has a educational objective. For an intensive and scientific use, please see the R package of Timothy P. Jurka (https://github.com/timjurka/maxent), which uses under the hood the MaxEnt C++ library of Tsuruoka (http://www.logos.ic.i.u-tokyo.ac.jp/~tsuruoka/maxent/).


## Usage

* General usage:

```R
maximumentropy(x, y=NULL, data=NULL, iteration=NULL,
method=c("L-BFGS-B", "GIS", "CG", "BFGS"), verbose=TRUE, normalize=FALSE)
```


* Example with formula object:

```R
model <- maximumentropy(Species ~ ., data=iris, iteration=100)
predict(model, new_dataset)
```


* Example with data.frame or matrix:

```R
x <- subset(iris, select=-Species)
y <- iris$Species
model <- maximumentropy(x, y)
```


* Prediction and evaluation:

```R
predicted_y <- predict(model, new_dataset)
evaluate(true_y, predicted_y)
```

* See information about a model:

```R
summary(model) #summary of model
# or
model$value #final value of log-likelihood
model$lambda #computed lambda

model$levels #discrete classes Y
model$features #names of features x

model$method #name of iterative method
model$iteration #number of iterations
model$convergence #0=converged, 1=max number of iterations reached, ... (see ?optim)
```



## Comments

* By default, the MaxEnt is trained with a BFGS. The GIS is the slower of all iterative methods.

* By default, the maximization of the objective function is done until convergence.

* The `y` vector must be a factor vector.

