crch.stabsel <- function( formula, data, ..., nu = 0.1, q, B = 100, thr = .9, maxit = 2000, data_percentage = .5 ) {

   
  ## Checking threshold. Thresholds have to be >0.5 and <1.0
  if ( thr <= 0.5 | thr >=1 )
    stop("Parameter \"thr\" for crch.stabsel has to be > 0.5 and < 1.0!")
  control = crch::crch.boost( maxvar = q, maxit = maxit, nu = nu )

  ## Check if user has set left/right. Left/right can only be numeric
  ## not vectors in the case of stability selection at the moment.
  dotargs <- list( ... )
  for ( n in c("left","right") ) {
    if ( ! n %in% names(dotargs) ) next
    if ( length(dotargs[[n]]) != 1 )
        stop(sprintf("Argument \"%s\" has to be of length 1 for crch.stabsel",n))
  }

  ## Number of data
  N <- nrow(data)

  ## Calling crch
  selected <- list()
  for ( i in 1:B ) {
    mod           <- crch::crch( formula, data = data[sample(1:N,ceiling(N*data_percentage),replace=FALSE),],
                     ..., control = control )
    selected[[i]] <- coef(mod)
  }

  ## Taking selected summary
  tabsel <- table(names(unlist(selected)))
  tabsel <- tabsel[order(tabsel,decreasing=FALSE)]

  ## Comupte per-family-error-rate
  p <- length( attr(mod$terms$scale,   "term.labels") ) +
       length( attr(mod$terms$location,"term.labels") )
  PFER <- (q^2) / ((2*thr - 1)*p)

  ## Generate new formula
  taken <- tabsel[tabsel/B >= thr]
  f <- list( response = as.character(attr(Formula::Formula(mod$formula),"lhs")[[1]]),
             location = names(taken)[grep("^[^\\(]",names(taken))],
             scale    = gsub("^\\(scale\\)_","",names(taken)[grep("^\\(scale\\)_[^\\(]",names(taken))]) )

  if ( length(f$scale)==0 ) {
    f <- sprintf("%s ~ %s",f$response,paste(f$location,collapse=" + "))
  } else {
    f <- sprintf("%s ~ %s | %s",f$response,paste(f$location,collapse=" + "),paste(f$scale,collapse=" + "))
  }

  ## Store censoring, store unique!
  cens <- mod$cens; for ( i in 1:length(cens) ) cens[[i]] <- unique(cens[[i]])

  ## Return
  rval <- list("table"       = tabsel,
               "formula.org" = formula,
               "formula.new" = as.formula(f),
               "family"      = list(dist=mod$dist, cens=cens, truncated=mod$truncated ),
               "parameter"   = list("q" = q, "B" = B, "thr" = thr,
                                    "p" = p, "PFER" = PFER))
  class(rval) <- c("stabsel.crch", "list")
  return(rval)
}


## Small plotting method for S3 object stabsel.crch
## returned from crch.stabsel.
plot.stabsel.crch <- function(x, show = NULL,
   pal = function(n) gray.colors(n, start = 0.9, end = 0.6), main, ...) {

   ## Keep user parameter settings for plotting windows
   hold <- par( no.readonly = TRUE )
   on.exit( par(hold) )

   tabsel <- x$table
   thr    <- x$parameter$thr
   B      <- x$parameter$B
   models <- names(x$formula.new)

   n <- length(tabsel)
   start <- ifelse(is.null(show), 1 , n - show + 1)
   modelID <- factor( as.numeric(grepl("^\\(scale\\)",names(x$table))),
                  levels=0:1, labels=c("location","scale") )

   col <- pal(nlevels(modelID))

   if ( missing(main) ) main <- "Stability Selection Summary"

   par(mar = c(5, 12, 4, 2) + .1)
   bp <- barplot(tabsel[start:n], horiz = TRUE, las = 1,
                 col = col[modelID[start:n]], main=main, ...)
   abline(v = thr*B, col = 1, lty = 3, lwd = 2)
   legend("bottomright", fill = col, legend = levels(modelID), bg = "white")
   title(xlab = "Frequency")

   invisible(bp)
}


## Small print method for S3 object stabsel.crch
## returned from crch.stabsel.
print.stabsel.crch <- function(x, ...) {
    cat("\n  crch Stability Selection Info\n\n")
    cat("  Family information\n")
    cat(sprintf("   Distribution:       %s\n", x$family$dist))
    cat(sprintf("   Truncated:          %s\n", x$family$truncated  ))
    cat(sprintf("   Left-censoring:     %s\n", x$family$cens$left  ))
    cat(sprintf("   Right-censoring:    %s\n", x$family$cens$right ))

    cat("\n  Selection parameters\n")
    cat(sprintf("   Number of parameters (q):      %4d\n", x$parameter$q))
    cat(sprintf("   Number of iterations (B):      %4d\n", x$parameter$B))
    cat(sprintf("   Selection threshold (thr):     %4.2f\n", x$parameter$thr))
    cat(sprintf("   Total number of parameter (p): %4d\n", x$parameter$p))

    cat("\n  Per family error rate\n")
    cat(sprintf("   Given p/q/thr the expected value of falsely\n"))
    cat(sprintf("   chosen parameters is:          %4.2f\n", x$parameter$PFER))

    cat("\n  Selected formula:\n")
    cat("   "); print(x$formula.new)
}

