## density
dclogis <- function(x, location = 0, scale = 1, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), location = as.numeric(location), scale = as.numeric(scale), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cdclogis", x, location, scale, left, right, log))
  if(is.matrix(x)) {
    rval <- matrix(rval, ncol = ncol(x), nrow = nrow(x))
    colnames(rval) <- colnames(x)
    rownames(rval) <- rownames(x)
  }
  return(rval)
}

## distribution function
pclogis <- function(q, location = 0, scale = 1, left = -Inf, right = Inf,
  lower.tail = TRUE, log.p = FALSE) {
  input <- data.frame(q = as.numeric(q), location = as.numeric(location), scale = as.numeric(scale), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cpclogis", q, location, scale, left, right, lower.tail, log.p))
  if(is.matrix(q)) {
    rval <- matrix(rval, ncol = ncol(q), nrow = nrow(q))
    colnames(rval) <- colnames(q)
    rownames(rval) <- rownames(q)
  }
  return(rval)
}

## random numbers
rclogis <- function(n, location = 0, scale = 1, left = -Inf, right = Inf) {
  rval <- rlogis(n) * scale + location
  pmax(pmin(rval, right), left)
}

## quantiles
qclogis <- function(p, location = 0, scale = 1, left = -Inf, right = Inf,
  lower.tail = TRUE, log.p = FALSE) {
  rval <- qlogis(p, lower.tail = lower.tail, log.p = log.p) * scale + location
  rval <- pmax(pmin(rval, right), left)
  if(is.matrix(p)) {
    rval <- matrix(rval, ncol = ncol(p), nrow = nrow(p))
    colnames(rval) <- colnames(p)
    rownames(rval) <- rownames(p)
  }
  return(rval)
}

## scores
sclogis <- function(x, location = 0, scale = 1, left = -Inf, right = Inf,
  which = c("mu", "sigma")) {
  input <- data.frame(x = as.numeric(x), location = as.numeric(location), scale = as.numeric(scale), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  
  for(w in which) {
    if(w == "mu")
      score2 <- with(input, .Call("sclogis_mu", x, location, scale, left, right))
    if(w == "sigma")
      score2 <- with(input, .Call("sclogis_sigma", x, location, scale, left, right))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

# Hessian
hclogis <- function(x, location = 0, scale = 1, left = -Inf, right = Inf, 
  which = c("mu", "sigma")) {
  input <- data.frame(x = as.numeric(x), location = as.numeric(location), scale = as.numeric(scale), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  for(w in which) {       
    if(w == "mu")         
      hess[[w]] <- with(input, .Call("hclogis_mu", x, location, scale, left, right))  
    if(w == "sigma")
      hess[[w]] <- with(input, .Call("hclogis_sigma", x, location, scale, left, right))  
    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- with(input, .Call("hclogis_musigma", x, location, scale, left, right))  
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}

## Expectation
eclogis <- function(location = 0, scale = 1, left = -Inf, right = Inf) {
  rmm <- (right-location)/scale
  lmm <- (left-location)/scale
  pncens <- plogis(rmm)-plogis(lmm)
  pncens*etlogis(location = location, scale = scale, left = left, right = right) + 
    plogis(lmm)*left^(is.finite(left)) + 
    plogis(rmm, lower.tail = FALSE)*right^(is.finite(right))
}
