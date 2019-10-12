## density
dtlogis <- function(x, location = 0, scale = 1, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), location = as.numeric(location), scale = as.numeric(scale), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cdtlogis", x, location, scale, left, right, log))
  if(is.matrix(x)) {
    rval <- matrix(rval, ncol = ncol(x), nrow = nrow(x))
    colnames(rval) <- colnames(x)
    rownames(rval) <- rownames(x)
  }
  return(rval)
}


## distribution function
ptlogis <- function(q, location = 0, scale = 1, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  input <- data.frame(q = as.numeric(q), location = as.numeric(location), scale = as.numeric(scale), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cptlogis", q, location, scale, left, right, lower.tail, log.p))
  if(is.matrix(q)) {
    rval <- matrix(rval, ncol = ncol(q), nrow = nrow(q))
    colnames(rval) <- colnames(q)
    rownames(rval) <- rownames(q)
  }
  return(rval)
}

## quantiles
qtlogis <- function(p, location = 0, scale = 1, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  if(log.p) p <- exp(p) 
  lower <- if(lower.tail) left else right
  upper <- if(lower.tail) right else left
  p <- plogis((lower-location)/scale, lower.tail = lower.tail) * (1 - p) + 
    p*plogis((upper - location)/scale, lower.tail = lower.tail)
  rval <- qlogis(p, lower.tail = lower.tail)*scale + location
  if(is.matrix(p)) {
    rval <- matrix(rval, ncol = ncol(p), nrow = nrow(p))
    colnames(rval) <- colnames(p)
    rownames(rval) <- rownames(p)
  }
  return(rval)
}

## random numbers
rtlogis <- function(n, location = 0, scale = 1, left = -Inf, right = Inf) {
  qtlogis(runif(n), location, scale, left = left, right = right)
}

## scores
stlogis <- function(x, location = 0, scale = 1, left = -Inf, right = Inf, 
  which = c("mu", "sigma")) {
  input <- data.frame(x = as.numeric(x), location = as.numeric(location), scale = as.numeric(scale), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  
  for(w in which) {
    if(w == "mu")
      score2 <- with(input, .Call("stlogis_mu", x, location, scale, left, right))
    if(w == "sigma")
      score2 <- with(input, .Call("stlogis_sigma", x, location, scale, left, right))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

## Hessian
htlogis <- function(x, location = 0, scale = 1, left = -Inf, right = Inf, 
  which = c("mu", "sigma")) {
  input <- data.frame(x = as.numeric(x), location = as.numeric(location), scale = as.numeric(scale), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  for(w in which) {       
    if(w == "mu")         
      hess[[w]] <- with(input, .Call("htlogis_mu", x, location, scale, left, right))  
    if(w == "sigma")
      hess[[w]] <- with(input, .Call("htlogis_sigma", x, location, scale, left, right))  
    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- with(input, .Call("htlogis_musigma", x, location, scale, left, right))  
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}

## Expectation
etlogis <- function(location = 0, scale = 1, left = -Inf, right = Inf) {
  rmm <- (right-location)/scale
  lmm <- (left-location)/scale
  pncens <- plogis(rmm)-plogis(lmm)
  rval <- location + scale*(rmm*plogis(rmm) - log(1+exp(rmm)) - 
    lmm*plogis(lmm) + log(1+exp(lmm))) /pncens
  rval
}
