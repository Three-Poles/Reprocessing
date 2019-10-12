## density
dtnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cdtnorm", x, mean, sd, left, right, log))
  if(is.matrix(x)) {
    rval <- matrix(rval, ncol = ncol(x), nrow = nrow(x))
    colnames(rval) <- colnames(x)
    rownames(rval) <- rownames(x)
  }
  return(rval)
}


## distribution function
ptnorm <- function(q, mean = 0, sd = 1, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  input <- data.frame(q = as.numeric(q), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cptnorm", q, mean, sd, left, right, lower.tail, log.p))
  if(is.matrix(q)) {
    rval <- matrix(rval, ncol = ncol(q), nrow = nrow(q))
    colnames(rval) <- colnames(q)
    rownames(rval) <- rownames(q)
  }
  return(rval)
}

## quantiles
qtnorm <- function(p, mean = 0, sd = 1, left = -Inf, right = Inf,
  lower.tail = TRUE, log.p = FALSE) {
  if(log.p) p <- exp(p) 
  lower <- if(lower.tail) left else right
  upper <- if(lower.tail) right else left
  p <- pnorm((lower-mean)/sd, lower.tail = lower.tail) * (1 - p) + 
    p*pnorm((upper - mean)/sd, lower.tail = lower.tail)
  rval <- qnorm(p, lower.tail = lower.tail)*sd + mean
  if(is.matrix(p)) {
    rval <- matrix(rval, ncol = ncol(p), nrow = nrow(p))
    colnames(rval) <- colnames(p)
    rownames(rval) <- rownames(p)
  }
  return(rval)
}

## random numbers
rtnorm <- function(n, mean = 0, sd = 1, left = -Inf, right = Inf) {
  qtnorm(runif(n), mean, sd, left = left, right = right)
}



## scores
stnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf,
  which = c("mu", "sigma")) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  
  for(w in which) {
    if(w == "mu")
      score2 <- with(input, .Call("stnorm_mu", x, mean, sd, left, right))
    if(w == "sigma")
      score2 <- with(input, .Call("stnorm_sigma", x, mean, sd, left, right))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

## Hessian
htnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf,
  which = c("mu", "sigma")) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  for(w in which) {       
    if(w == "mu")         
      hess[[w]] <- with(input, .Call("htnorm_mu", x, mean, sd, left, right))  
    if(w == "sigma")
      hess[[w]] <- with(input, .Call("htnorm_sigma", x, mean, sd, left, right))  
    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- with(input, .Call("htnorm_musigma", x, mean, sd, left, right))  
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}


## Expectation
etnorm <- function(mean = 0, sd = 1, left = -Inf, right = Inf) {
  rmm <- (right-mean)/sd
  lmm <- (left-mean)/sd
  pncens <- pnorm(rmm)-pnorm(lmm)
  rval <- mean + sd*(dnorm(lmm) - dnorm(rmm))/pncens
  rval
}
