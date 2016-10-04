# -------------------------------------------------------------
# extension of ns() from splines, to allow for clamping.

ns.ld <- function(x, df = NULL, knots = NULL, intercept = FALSE,
                  Boundary.knots = range(x), fixsl=c(FALSE,FALSE) )
{
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if(nas <- any(nax))    x <- x[!nax]


  if(!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > Boundary.knots[2L])
  }

  else outside <- FALSE

  if(!missing(df) && missing(knots)) {
    ## df = number(interior knots) + 1 + intercept
    nIknots <- df - 1 - intercept
    if(nIknots < 0) {
      nIknots <- 0
      warning("'df' was too small; have used ", 1 + intercept)
    }

    knots <- if(nIknots > 0) {
      knots <- seq.int(0, 1,
                       length.out = nIknots + 2L)[-c(1L, nIknots + 2L)]
      stats::quantile(x[!outside], knots)
    } ## else  NULL

  } else nIknots <- length(knots)
  Aknots <- sort(c(rep(Boundary.knots, 4L), knots))

  if(any(outside)) {
    basis <- array(0, c(length(x), nIknots + 4L))

    if(any(ol)) {
      k.pivot <- Boundary.knots[1L]
      xl <- cbind(1, x[ol] - k.pivot)
      tt <- splines::spline.des(Aknots, rep(k.pivot, 2L), 4, c(0, 1))$design
      basis[ol,  ] <- xl %*% tt
    }

    if(any(or)) {
      k.pivot <- Boundary.knots[2L]
      xr <- cbind(1, x[or] - k.pivot)
      tt <- splines::spline.des(Aknots, rep(k.pivot, 2L), 4, c(0, 1))$design
      basis[or,  ] <- xr %*% tt
    }


    if(any(inside <- !outside))
      basis[inside,  ] <- splines::spline.des(Aknots, x[inside], 4)$design
  }

  else basis <- splines::spline.des(Aknots, x, 4)$design

# Fix the clamping
if( !is.logical(fixsl) ) warning( "fixsl elements must be of mode logical" )
# Only the 4th parameter affected, should be either 1 or 2 in the two positions
const <- splines::spline.des( Aknots, Boundary.knots, 4, c(2-fixsl[1],2-fixsl[2]) )$design
  
  if(!intercept) {
    const <- const[, -1 , drop = FALSE]
    basis <- basis[, -1 , drop = FALSE]
  }

  qr.const <- qr(t(const))
  basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[,  - (1L:2L), drop = FALSE])
  n.col <- ncol(basis)

  if(nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }

  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = 3, knots = if(is.null(knots)) numeric(0) else knots,
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("cns", "basis", "matrix")
  basis
}

#-------------------------------------------------------------------------------
# An extension of ns that automatically  takes the smallest and largest knots
# as boundary knots without further ado, but also allows centering
# around a reference and detrending by means of projection, as well as
# clamping.

Ns <- function( x, ref = NULL,
                    df = NULL,
                 knots = NULL,
             intercept = FALSE,
        Boundary.knots = NULL,
                 fixsl = c(FALSE,FALSE),
               detrend = FALSE )
{
  ## Check sensibility of arguments
  if( !is.null(ref) ) {
    if( !is.vector(ref) )
      stop( "Argument 'ref' must be a scalar, but it is a ", class(ref), "." )
    if( is.vector(ref) & length(ref)>1 )
      stop( "Argument 'ref' must be a scalar, but has length ", length(ref), "." )
    if( intercept ) {
      warning( "ref= specified, hence intercept=TRUE is ignored")
      intercept <- FALSE
      }
    }
  ## Detrending required?
  if( any(detrend>0) ) { # covers both logical and vector
    if( any(detrend<0) )
      stop( "Some elements of weight are <0, e.g. no",
            (ww <- which(detrend<0))[1:min(5,length(ww))], "." )
    if( !(length(detrend) %in% c(1,length(x))) ) {
      warning( "Weights in inner product diagonal matrix set to 1")
      weight <- rep(1,length(x))
      }
    else weight <- if( is.numeric(detrend) ) detrend else rep(1,length(x))
    detrend <- TRUE
    }
  if( detrend & intercept ) {
    warning( "detrend= specified, hence intercept=TRUE is ignored")
    intercept <- FALSE
    }
  if( detrend & any(!is.na(fixsl)) ) {
    warning( "detrend= specified, hence fixsl argument is ignored")
    fixsl=c(NA,NA)
    }
  ## Here is the specification of the spline basis
  ## df= specified
  if( !is.null(df) )
    MM <- ns.ld( x, df = df, intercept = (intercept & is.null(ref)), fixsl = fixsl )
  else
    ## knots= specified
  {
    if( is.null( Boundary.knots ) )
    {
      if( !is.null( knots ) )
      {
        knots <- sort( unique( knots ) )
        ok <- c(1,length(knots))
        Boundary.knots <- knots[ok]
        knots <- knots[-ok]
      }
    }
    MM <- ns.ld( x, knots = knots,
           Boundary.knots = Boundary.knots,
                intercept = (intercept & is.null(ref)),
                    fixsl = fixsl )
  }
  ## Reference point specified ?
  if( !is.null(ref) )
  {
    MM <- MM - ns.ld( rep(ref,length(x)),
                      knots = attr(MM,"knots"),
             Boundary.knots = attr(MM,"Boundary.knots"),
                      fixsl = fixsl )
  }
    
  ## Detrending required ?
  if( detrend )
  {
    DD <- detrend( MM, x, weight=weight )
    ## NOTE: detrend does not preserve attributes
    for( aa in c("degree","knots","Boundary.knots","intercept","class") )
       attr( DD, aa ) <- attr( MM, aa )
    attr( DD, "detrend" ) <- TRUE
    attr( DD, "proj.wt" ) <- weight
    MM <- DD
  }
  return( MM )
}
