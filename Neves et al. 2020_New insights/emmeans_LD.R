emmeans_LD <- function(m0,mf,cov.list,p=0.5) {
  # Find reference grids
  suppressMessages(em0 <- emmeans(m0,mf,at=cov.list[[1]]))
  suppressMessages(em1 <- emmeans(m0,mf,at=cov.list[[2]]))
  # Find offset of intercept
  tran <- em0@misc$tran
  if ((p <= 0)|(p >= 1)) stop("p must be strictly between 0 and 1")
  if (is.null(tran)) stop("model object must include a link function")
  intercept <- switch(tran,
                      probit = qnorm(p),
                      logit = log(p/(1-p)),
                      cloglog = log(-log(1-p)),
                      stop("link function must be probit, logit, or cloglog"))
  # Find indices of non-estimable LD's
  nbasis0 <- em0@nbasis; if (all(is.na(nbasis0))) nbasis0 <- rep(0,length(em0@bhat)) 
  nbasis1 <- em1@nbasis; if (all(is.na(nbasis1))) nbasis1 <- rep(0,length(em1@bhat)) 
  ii <- apply(cbind(em0@linfct %*% nbasis0,em1@linfct %*% nbasis1),1,function(x) {any(x!=0)})
  # Find linear transformations
  A <- em0@linfct
  B <- em1@linfct - em0@linfct
  # Find estimates: Set NA's in hat.theta to zero 
  hat.theta <- em0@bhat; hat.theta[is.na(hat.theta)] <- 0
  alpha  <- c(A%*%hat.theta)
  beta   <- c(B%*%hat.theta)
  LD     <- (intercept-alpha)/beta
  LD[ii] <- NA
  # Find variance matrix
  V <- vcov(m0); V[is.na(V)] <- 0
  H <- cbind(diag(x=-1/beta),diag((alpha-intercept)/(beta^2)))%*%rbind(A,B)
  V <- H%*%V%*%t(H)
  # Make emmGrid object by altering em0
  em <- em0
  em@model.info <- list(call = match.call(), xlev = em0@levels) 
  em@bhat   <- LD
  em@V      <- V[!ii,!ii]
  em@linfct <- diag(1,nrow=length(LD))
  if (sum(ii)==0) {
    em@nbasis <- matrix(NA,1,1)
  } else {
    em@nbasis <- matrix(as.numeric(ii),length(ii),1)
  }
  em@dffun <- function(k,dfargs) {Inf}
  em@dfargs <- list()
  em@misc <- list(estName = "estimate", estType = "prediction", 
                  infer = c(TRUE, FALSE), level = 0.95, adjust = "none", 
                  famSize = nrow(em@linfct), avgd.over = character(0), 
                  pri.vars = names(em@grid), 
                  methDesc = "emmobj")
  # return result
  return(em)
}
