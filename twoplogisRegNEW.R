#rm(list = ls())

source("E:\\Paper under review\\TPLogistic\\Rcode\\twoplogis.R")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  This is the derivative of logistic
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


XUtplog = function(y, mylocat, myscale, myskew){

zedd <- (y - mylocat) / myscale # Is a n x NOS matrix

# cond1 <- (zedd <= 0)
cond2 <- (zedd > 0)

dl.dlocat <- 2 * tanh(zedd / (4 * myskew)) / (myskew) # cond1
dl.dlocat[cond2] <- (2 * tanh(zedd / (4 * (1-myskew))) / ( 1-myskew))[cond2]
dl.dlocat <- dl.dlocat / (4 * myscale)

dl.dlocat
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  generating Data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = 20

myskew = 0.25; myscale = 2

yy = rtwoplogis(nn, 0, myscale, .5)
x = 1:nn/10 ; X = cbind(1, x)

beta = c(100, 30)
mu = X %*% beta
y = mu + yy


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  TPL.logisticEst
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 logisticREst = function(y, x, myscale, myskew, Tol = 1e-5){

 x0 <- cbind(1, x)
qx0 <- qr(x0)
  ## compute (x'x)^(-1) x'y
betao = solve.qr(qx0, y)
# betao = coef(lm(y~x))
Bo = betao
X = cbind(1, x)
muo = X %*% betao

mun = X %*% (betao +.1)

BWo = (betao)

BWn = (betao +.1)



while(any(abs(BWo - BWn) > Tol)){
BWo = BWn
muo = mun
XU  = XUtplog(y, muo, myscale, myskew)
XtU = t(X) %*% XU

EE = y - mun
WW =  matrix(rep(as.vector(XU / (EE)), ncol(X)), nc = ncol(X)) * X
#print("WW")
#print(WW)
temp15 = (t(y) %*% (XU / (EE)))
BWn = solve(t(WW) %*% (X)) %*% (t(WW) %*% y)


cond1 <- (EE <= 0)
cc = 0.324
sighat = cc * (sum(abs(EE[cond1]/myskew)) +
                                  sum(abs(EE[!cond1]/(1-myskew))))/(nn * 0.9)


Temp20 = myskew*(1 - myskew)
Temp21 = 4 * ((1 - 3 * Temp20) * ((pi^2 / 3) - log(2)^2) + Temp20 * log(2)^2)
v = Temp21 * sighat^2

VV  = v * solve(t(WW) %*% (X)) %*% (t(WW) %*% WW) %*% solve(t(X) %*% (WW))

temp12 = 0.25 * solve(t(X) %*% X) %*% XtU


mun = X %*% BWn
}

df <- nrow(X) - ncol(X)
 list(coefficients = BWn, vcov = VV,
        sigma = sighat,  df = df, p = myskew)}

BWn = logisticREst(y, x,  myscale = myscale, myskew = .75)
BWn


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  TPL.logisticEst
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   TPLogR <- function(x, ...) UseMethod("TPLogR")
#=========================================================================
#         ptpn
#=========================================================================
TPLogR.default <- function(formula, data, weights,  method = "qr",
    model = TRUE, x = FALSE, y = FALSE, qr = TRUE,  contrasts = NULL,
    offset, myskew, ...) {

  x <- as.matrix(x)
  y <- as.numeric(y)


    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", #"subset", "weights", "na.action",
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (method == "model.frame")
        return(mf)

    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")



    offset <- as.vector(model.offset(mf))

    if (is.empty.model(mt)) {
        x <- NULL
        est <- list(coefficients = if (is.matrix(y)) matrix(, 0,
            3) else numeric(), residuals = y, fitted.values = 0 * y,
             rank = 0L )
        if (!is.null(offset)) {
            est$fitted.values <- offset
            est$residuals <- y - offset
        }
    }
    else {
        x <- model.matrix(mt, mf, contrasts)
         est <- logisticREst(y, x, myscale, myskew, Tol = 1e-5)
        }
  ##############################################

  est$fitted.values <- as.vector(x %*% est$coefficients)
  est$residuals <- y - est$fitted.values
  est$p <- p
  #est$call <- match.call()
  class(est) <- "TPLogR"

    est$offset <- offset
    est$contrasts <- attr(x, "contrasts")
    est$xlevels <- .getXlevels(mt, mf)
    est$call <- cl
    est$terms <- mt
    if (model)
        est$model <- mf
    if (ret.x)
        est$x <- x
    if (ret.y)
        est$y <- y
    if (!qr)
        est$qr <- NULL

  est
}

#=========================================================================
#         ptpn
#=========================================================================

print.TPLogR <- function(x, ...) {
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients)
}

#=========================================================================
#         Arash has to change the tests. The test is not correct
#=========================================================================
summary.TPLogR <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  tval <- coef(object) / se
  TAB <- cbind(
    Estimate = coef(object), StdErr = se, t.value = tval,
    if (object$p == 0.5) p.value = 2 * pt(-abs(tval), df = object$df)
     else p.value = 2 * pnorm(-abs(tval))
         )
  res <- list(call = object$call,  coefficients = TAB)
  class(res) <- "summary.TPLogR"
  res
}

#=========================================================================
#         ptpn
#=========================================================================

print.summary.TPLogR <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.value = TRUE, has.Pvalue = TRUE)
}

#=========================================================================
#         ptpn
#=========================================================================
TPLogR.formula <- function(formula, data = list(), ...) {
  mf <- model.frame(formula = formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  est <- TPLogR.default(x, y, myskew, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}

#=========================================================================
#         ptpn
#=========================================================================
predict.TPLogR <- function(object, newdata = NULL, ...) {
  if(is.null(newdata)) y <- fitted(object) else{
   if(!is.null(object$formula)){
  ## model has been fitted using formula interface
  x <- model.matrix(object$formula, newdata)
   }      else{x <- newdata}
  y <- as.vector(x %*% coef(object))
   }
  y
}