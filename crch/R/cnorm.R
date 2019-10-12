## density
dcnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cdcnorm", x, mean, sd, left, right, log))
  if(is.matrix(x)) {
    rval <- matrix(rval, ncol = ncol(x), nrow = nrow(x))
    colnames(rval) <- colnames(x)
    rownames(rval) <- rownames(x)
  }
  return(rval)
}

## distribution function
pcnorm <- function(q, mean = 0, sd = 1, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  input <- data.frame(q = as.numeric(q), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cpcnorm", q, mean, sd, left, right, lower.tail, log.p))
  if(is.matrix(q)) {
    rval <- matrix(rval, ncol = ncol(q), nrow = nrow(q))
    colnames(rval) <- colnames(q)
    rownames(rval) <- rownames(q)
  }
  return(rval)
}

## random numbers
rcnorm <- function(n, mean = 0, sd = 1, left = -Inf, right = Inf) {
  rval <- rnorm(n) * sd + mean
  pmax(pmin(rval, right), left)
}

## quantiles
qcnorm <- function(p, mean = 0, sd = 1, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  rval <- qnorm(p, lower.tail = lower.tail, log.p = log.p) * sd + mean
  rval <- pmax(pmin(rval, right), left)
  if(is.matrix(p)) {
    rval <- matrix(rval, ncol = ncol(p), nrow = nrow(p))
    colnames(rval) <- colnames(p)
    rownames(rval) <- rownames(p)
  }
  return(rval)
}

## scores
scnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf,
  which = c("mu", "sigma")) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  
  for(w in which) {
    if(w == "mu")
      score2 <- with(input, .Call("scnorm_mu", x, mean, sd, left, right))
    if(w == "sigma")
      score2 <- with(input, .Call("scnorm_sigma", x, mean, sd, left, right))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

## Hessian
hcnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, 
  which = c("mu", "sigma")) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  for(w in which) {       
    if(w == "mu")         
      hess[[w]] <- with(input, .Call("hcnorm_mu", x, mean, sd, left, right))  
    if(w == "sigma")
      hess[[w]] <- with(input, .Call("hcnorm_sigma", x, mean, sd, left, right))  
    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- with(input, .Call("hcnorm_musigma", x, mean, sd, left, right))  
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}


## Expectation
ecnorm <- function(mean = 0, sd = 1, left = -Inf, right = Inf) {
  rmm <- (right-mean)/sd
  lmm <- (left-mean)/sd
  pncens <- pnorm(rmm)-pnorm(lmm)
  pncens*etnorm(mean = mean, sd = sd, left = left, right = right) + 
    pnorm(lmm)*left^(is.finite(left)) + 
    pnorm(rmm, lower.tail = FALSE)*right^(is.finite(right))
}
