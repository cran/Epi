poisreg <- function (link = "log")
{
    linktemp <- substitute(link)
    if (!is.character(linktemp)) linktemp <- deparse(linktemp)
    okLinks <- c("log", "identity", "sqrt")
    if (linktemp %in% okLinks)
	stats <- make.link(linktemp)
    else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    } else {
        ## what else shall we allow?  At least objects of class link-glm.
        if(inherits(link, "link-glm")) {
            stats <- link
            if(!is.null(stats$name)) linktemp <- stats$name
        } else {
            stop(gettextf('link "%s" not available for poisreg family; available links are %s',
			  linktemp, paste(sQuote(okLinks), collapse =", ")),
		 domain = NA)
        }
    }
    variance <- function(mu) mu
    validmu <- function(mu) all(is.finite(mu)) && all(mu>0)
    dev.resids <- function(y, mu, wt)
    { ## faster than  2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
	r <- mu*wt
	p <- which(y > 0)
	r[p] <- (wt * (y*log(y/mu) - (y - mu)))[p]
	2*r
    }
    aic <- function(y, n, mu, wt, dev) {
	-2*sum(ifelse(n > 0, (wt/n), 0)*dpois(round(y*n), mu*n, log=TRUE))
    }
    initialize <- expression({
        if (NCOL(y) == 1) {
            n <- rep.int(1, nobs)
            y[weights == 0] <- 0
            if (any(y < 0)) {
                stop("y values must be >= 0")
            }
            m <- weights * y
            if (any(abs(m - round(m)) > 0.001)) {
                warning("non-integer #successes in poisreg glm!")
            }
            mustart <- m + 0.1
        }
        else if (NCOL(y) == 2) {
            if (any(y[,1] < 0)) {
                stop("negative values not allowed for the 'poisreg' family")
            }
            if (any(y[,2] < 0)) {
                stop("negative time not allowed for the 'poisreg' family")
            }
            if(any(y[,1] > 0 & y[,2] == 0)) {
                stop("non-zero counts in zero time in a poisreg glm!")
            }
            if(any(abs(y[,1] - round(y[,1])) > 1e-3)) {
                warning("non-integer counts in a poisreg glm!")
            }

            n <- y[,2]
            y <- ifelse(n == 0, 0, y[, 1]/n)
            
            weights <- weights * n
            mustart <- y + 0.1
        }
        else {
            stop("for the 'poisreg' family, y must be a 2 column matrix where col 1 is no. events and col 2 is time")
        }
    })
    simfun <- function(object, nsim) {
        wts <- object$prior.weights
        ftd <- fitted(object)
        rpois(nsim*length(ftd), ftd*wts)
    }
    structure(list(family = "poisson", ##Fool summary.glm
		   link = linktemp,
		   linkfun = stats$linkfun,
		   linkinv = stats$linkinv,
		   variance = variance,
		   dev.resids = dev.resids,
		   aic = aic,
		   mu.eta = stats$mu.eta,
		   initialize = initialize,
		   validmu = validmu,
		   valideta = stats$valideta,
                   simulate = simfun),
	      class = "family")
}
