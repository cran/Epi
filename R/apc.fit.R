apc.fit <-
function( data,
             A,
             P,
             D,
             Y,
         ref.c,
         ref.p,
           dist = c("poisson","binomial"),
          model = c("ns","bs","ls","factor"),
        dr.extr = "Y",
           parm = c("ACP","APC","AdCP","AdPC","Ad-P-C","Ad-C-P","AC-P","AP-C"),
           npar = c( A=5, P=5, C=5 ),
          scale = 1,
          alpha = 0.05,
      print.AOV = TRUE )
{
 dist <- match.arg(dist)
model <- match.arg(model)
drtyp <- deparse(substitute(dr.extr))
 parm <- toupper(match.arg(parm))

has.data <- !missing( data )
has.pref <- !missing( ref.p )
has.cref <- !missing( ref.c )

if(has.data)
  {
  if (length(match(c("A", "P", "D", "Y"), names(data))) != 4)
  stop("Data frame ", deparse(substitute(data)),
       " has columns:\n", names(data),
       "\nmust have variables:\n", "A (age), P (period), D (cases) and Y (person-time)")
  data <- data[,c("A","P","D","Y")]
  data <- data[complete.cases(data),]
  A <- data$A
  P <- data$P
  D <- data$D
  Y <- data$Y
  }
else
  {
  nm <- logical(4)
  nm[1] <- missing(A)
  nm[2] <- missing(P)
  nm[3] <- missing(D)
  nm[4] <- missing(Y)
  if (any(nm))
      stop("Variable", if (sum(nm) > 1)
          "s", paste(c(" A", " P", " D", " Y")[nm], collapse = ","),
          " missing from input")
  if( diff(range( lv <-
                  sapply( list(A = A, P = P, D = D, Y = Y),
                          length) ) ) != 0 )
      stop( "\nLengths of variables (", paste(paste(names(lv),
            lv, sep = ":"), collapse = ", "), ") are not the same." )
  }

# Utility to compute where the median amount y is on the x scale
med <-
function(x, y)
{
o <- order(x)
a <- y[o]
names(a) <- x[o]
return( as.numeric(names(a[cumsum(a)/sum(a) > 0.5][1])) )
}
# Set the reference points on the period and cohort scales
p0 <- ifelse( has.pref, ref.p, med(P  , D) )
c0 <- ifelse( has.cref, ref.c, med(P-A, D) )
 
# Number of parameters in the spline modeling 
if( is.list(npar) & length(npar)<3 )
  stop("npar given as a list - should have length 3! \n")
if( !is.list(npar) & length(npar)!=3 )
  {
  npar <- rep(npar, 3)[1:3]
  names(npar) = c("A","P","C")
  cat("NOTE: npar is specified as:\n") ; print( npar )
  }

if( is.null(names(npar)) ) names(npar) <- c("A", "P", "C")
# Labeling of confidence intervals
lu <- paste(formatC( c(alpha/2, 1 - alpha/2) * 100,
                     format = "f", digits = 1), "%", sep = "")

if( is.list(model) ) {
    if (!all(sapply(model, is.function)))
        stop("'model' is a list, but not all elements are functions as they should be.")
    if ((lmod <- length(model)) < 3)
        stop("'model' is a list, with", lmod, "elements, it should have three.")
    if (is.null(names(model)))
        names(model) <- c("A", "P", "C")
    MA <- model[["A"]](A)
    MP <- model[["P"]](P)
    MC <- model[["C"]](P-A)
    Rp <- model[["P"]](p0)
    Rc <- model[["C"]](c0)
}
else {
    if (model == "factor") {
        MA <- model.matrix(~factor(A) - 1)
        MP <- model.matrix(~factor(P) - 1)
        MC <- model.matrix(~factor(P - A) - 1)
        Rp <- MP[abs(P - p0) == min(abs(P - p0)), , drop = FALSE][1, ]
        Rc <- MC[abs(P - A - c0) == min(abs(P - A - c0)), , drop = FALSE][1, ]
    }
    if (model == "ns") {
        # is npar a list
        knl <- is.list( npar )
        # if scalar expand
        if( !knl & length(npar)==1 ) npar <- rep( npar, 3 )
        # if no names, provide them
        if( is.null(names(npar)) ) names(npar) <- c("A","P","C")
        # if names too long or wrong case, rectify
        names( npar ) <- toupper( substr(names(npar),1,1) )
        # if not a list make it one with the correct knots
        if( !knl ){
          nkn <- npar
          eqp <- function(n) (1:n-0.5)/n
          npar <- list( A = quantile( rep(  A,D), probs=eqp(nkn["A"]) ),
                        P = quantile( rep(P  ,D), probs=eqp(nkn["P"]) ),
                        C = quantile( rep(P-A,D), probs=eqp(nkn["C"]) ) )                       
          }
        MA <- Ns(  A, knots = npar[["A"]] )
        MP <- Ns(P  , knots = npar[["P"]] )
        MC <- Ns(P-A, knots = npar[["C"]] )
        Rp <- ns(p0, knots = attr(MP,"knots"),
            Boundary.knots = attr(MP,"Boundary.knots"))
        Rc <- ns(c0, knots = attr(MC,"knots"),
            Boundary.knots = attr(MC,"Boundary.knots"))
        Knots <- list( Age = sort(c(attr(MA,"knots"),
                                    attr(MA,"Boundary.knots"))),
                       Per = sort(c(attr(MP,"knots"),
                                    attr(MP,"Boundary.knots"))),
                       Coh = sort(c(attr(MC,"knots"),
                                    attr(MC,"Boundary.knots"))))
    }
    if (model %in% c("bs", "ls")) {
        deg <- switch(model, ls = 1, bs = 3)
        knl <- is.list(npar)
        if (knl) nk <- sapply(npar, length)
        MA <- if (knl) bs(A, knots = npar[["A"]][-c(1,nk[1])],
                    Boundary.knots = npar[["A"]][ c(1,nk[1])], degree = deg)
              else     bs(A, df = npar[["A"]], degree = deg)
        MP <- if (knl) bs(P, knots = npar[["P"]][-c(1,nk[2])],
                    Boundary.knots = npar[["P"]][ c(1,nk[2])], degree = deg)
              else     bs(P, df = npar[["P"]], degree = deg)
        MC <- if (knl) bs(P - A, knots = npar[["C"]][-c(1, nk[3])],
                        Boundary.knots = npar[["C"]][ c(1, nk[3])], degree = deg)
              else     bs(P - A, df = npar[["C"]], degree = deg)
        Rp <- bs(p0, knots = attr(MP,"knots"),
            Boundary.knots = attr(MP,"Boundary.knots"),
                    degree = attr(MP,"degree"))
        Rc <- bs(c0, knots = attr(MC,"knots"),
            Boundary.knots = attr(MC,"Boundary.knots"),
                    degree = attr(MC,"degree"))
        Knots <- list(Age = sort(c(attr(MA,"knots"),attr(MA,"Boundary.knots"))),
                      Per = sort(c(attr(MP,"knots"),attr(MP,"Boundary.knots"))),
                      Coh = sort(c(attr(MC,"knots"),attr(MC,"Boundary.knots"))))
    }
}
if (tolower(substr(dist, 1, 2)) == "po") {
    m.APC <- glm(D ~ MA + I(P - p0) + MP + MC,
                 offset = log(Y), family = poisson)
    Dist <- "Poisson with log(Y) offset"
}
is.bin <- FALSE
if (is.bin <- tolower(substr(dist, 1, 3)) %in% c("bin")) {
    m.APC <- glm(cbind(D, Y - D) ~ MA + I(P - p0) + MP + MC,
                 family = binomial)
    Dist <- "Binomial regression (logistic) of D/Y"
}
m.AP <- update(m.APC, . ~ . - MC)
m.AC <- update(m.APC, . ~ . - MP)
m.Ad <- update(m.AP , . ~ . - MP)
m.A  <- update(m.Ad , . ~ . - I(P - p0))
m.0  <- update(m.A  , . ~ . - MA)
AOV  <- anova(m.A, m.Ad, m.AC, m.APC, m.AP, m.Ad, test = "Chisq")
colnames(AOV)[1:4] <- c("Mod. df.","Mod. dev.",
                        "Test df.","Test dev.")
AOV <- abs(AOV)
AOV <- cbind(Model = c("Age",
                       "Age-drift",
                       "Age-Cohort",
                       "Age-Period-Cohort",
                       "Age-Period",
                       "Age-drift"),
               AIC = c(AIC(m.A),
                       AIC(m.Ad),
                       AIC(m.AC),
                       AIC(m.APC),
                       AIC(m.AP),
                       AIC(m.Ad)), 
               AOV,
     'Test dev/df' = AOV[,"Test dev."]/AOV[,"Test df."],
        'H0      ' = c("","zero drift ",  
                          "Coh eff|dr.",   
                          "Per eff|Coh",  
                          "Coh eff|Per",  
                          "Per eff|dr."))
A.pt  <- unique(A)
A.pos <- match(A.pt, A)
P.pt  <- unique(P)
P.pos <- match(P.pt, P)
C.pt  <- unique(P - A)
C.pos <- match(C.pt, P - A)
MA <- cbind(1, MA)

# Determine the inner product (diagonal) for projection
if (!mode(dr.extr) %in% c("character", "numeric"))
    stop("\"dr.extr\" must be of mode \"character\" or \"numeric\".\n")

if (is.character(dr.extr))
   { 
   wt <- rep(1, length(D) )
   drtyp <- "1-weights"
   if( toupper(substr(dr.extr, 1, 1)) %in% c("T","D") )
     { wt <- D
       drtyp <- "D-weights" } else
   if( toupper(substr(dr.extr, 1, 1)) %in% c("L","R") )
     { wt <- (Y^2)/D
       drtyp <- "Y^2/D-weights" } else
   if( toupper(substr(dr.extr, 1, 1)) %in% c("Y") )
     { wt <- Y
       drtyp <- "Y-weights" }
   }
if ( is.numeric(dr.extr) )
   {
   if( length(dr.extr)==1 )
     { wt <- D + dr.extr*Y
       drtyp <- paste("D+",dr.extr,"*Y weights",sep="") }
   if( length(dr.extr)==nrow(data) )
     { wt <- dr.extr
       if( any(wt<0) ) stop("dr.extr must be non-negative")
       drtyp <- "extn-weights" }
   }
Rp <- matrix(Rp, nrow = 1)
Rc <- matrix(Rc, nrow = 1)
xP <- Epi::detrend(rbind(Rp, MP), c(p0, P  ), weight = c(0, wt))
xC <- Epi::detrend(rbind(Rc, MC), c(c0, P-A), weight = c(0, wt))
 
MPr <- xP[-1,,drop=FALSE] - has.pref * xP[rep(1, nrow(MP)),,drop=FALSE]
MCr <- xC[-1,,drop=FALSE] - has.cref * xC[rep(1, nrow(MC)),,drop=FALSE]
if (length(grep("-", parm)) == 0) {
    if (parm %in% c("ADPC", "ADCP", "APC", "ACP"))
        m.APC <- update(m.0, . ~ . - 1 + MA + I(P - p0) + MPr + MCr)
    drift <- rbind( ci.exp(m.APC, subset = "I\\(", alpha = alpha),
                    ci.exp(m.Ad , subset = "I\\(", alpha = alpha) )
    rownames(drift) <- c(paste("APC (",drtyp,")",sep=""), "A-d")
    if (parm == "ADCP")
        m.APC <- update(m.0, . ~ . - 1 + MA + I(P - A - c0) + MPr + MCr)
    if (parm == "APC") {
        MPr <- cbind(P - p0, MPr)
        m.APC <- update(m.0, . ~ . - 1 + MA + MPr + MCr)
    }
    if (parm == "ACP") {
        MCr <- cbind(P - A - c0, MCr)
        m.APC <- update(m.0, . ~ . - 1 + MA + MPr + MCr)
    }
    Age <- cbind(Age = A.pt, ci.exp(m.APC, subset = "MA",
        ctr.mat = MA[A.pos,,drop=FALSE], alpha = alpha))[order(A.pt),]
    Per <- cbind(Per = P.pt, ci.exp(m.APC, subset = "MPr",
        ctr.mat = MPr[P.pos,,drop=FALSE], alpha = alpha))[order(P.pt),]
    Coh <- cbind(Coh = C.pt, ci.exp(m.APC, subset = "MCr",
        ctr.mat = MCr[C.pos,,drop=FALSE], alpha = alpha))[order(C.pt),]
    colnames(Age)[-1] <- c("Rate", lu)
    colnames(Per)[-1] <- c("P-RR", lu)
    colnames(Coh)[-1] <- c("C-RR", lu)
    Type <- paste("ML of APC-model", Dist, ": (", parm, "):\n")
    Model <- m.APC
}
else {
    adc <- update(m.0, . ~ . - 1 + MA + I(P - A - c0))
    adp <- update(m.0, . ~ . - 1 + MA + I(P - p0))
    drift <- ci.exp(adc, subset = "I\\(")
    rownames(drift) <- "A-d"
    xP <- cbind(1, P - p0, MPr)
    xC <- cbind(1, P - A - c0, MCr)
    lP <- cbind(P - p0, MPr)
    lC <- cbind(P - A - c0, MCr)
    if (parm == "AD-C-P") {
        rc <- update(m.0, . ~ . - 1 + xC, offset = predict(adc, type = "link"))
        rp <- update(m.0, . ~ . - 1 + xP, offset = predict(adc, type = "link"))
        A.eff <- ci.exp(adc, subset = "MA", ctr.mat = MA[A.pos,], alpha = alpha)
        C.eff <- ci.exp( rc, subset = "xC", ctr.mat = xC[C.pos,], alpha = alpha)
        P.eff <- ci.exp( rp, subset = "xP", ctr.mat = xP[P.pos,], alpha = alpha)
        Model <- list( adc, rc, rp )
    }
    else if (parm == "AD-P-C") {
        rp <- update(m.0, . ~ . - 1 + xP, offset = predict(adp,type = "link"))
        rc <- update(m.0, . ~ . - 1 + xC, offset = predict(rp,type = "link"))
        A.eff <- ci.exp(adp, subset = "MA", ctr.mat = MA[A.pos,], alpha = alpha)
        P.eff <- ci.exp(rp, subset = "xP", ctr.mat = xP[P.pos,], alpha = alpha)
        C.eff <- ci.exp(rc, subset = "xC", ctr.mat = xC[C.pos,], alpha = alpha)
        Model <- list( adp, rp, rc )
    }
    else if (parm == "AC-P") {
        ac <- update(m.0, . ~ . - 1 + MA + lC)
        rp <- update(m.0, . ~ . - 1 + xP, offset = predict(ac,type = "link"))
        A.eff <- ci.exp(ac, subset = "MA", ctr.mat = MA[A.pos,], alpha = alpha)
        C.eff <- ci.exp(ac, subset = "lC", ctr.mat = lC[C.pos,], alpha = alpha)
        P.eff <- ci.exp(rp, subset = "xP", ctr.mat = xP[P.pos,], alpha = alpha)
        Model <- list( ac, rp )
    }
    else if (parm == "AP-C") {
        ap <- update(m.0, . ~ . - 1 + MA + lP)
        rc <- update(m.0, . ~ . - 1 + xC, offset = predict(ap,type = "link"))
        A.eff <- ci.exp(ap, subset = "MA", ctr.mat = MA[A.pos,], alpha = alpha)
        P.eff <- ci.exp(ap, subset = "lP", ctr.mat = lP[P.pos,], alpha = alpha)
        C.eff <- ci.exp(rc, subset = "xC", ctr.mat = xC[C.pos,], alpha = alpha)
        Model <- list( ap, rc )
    }
    Age <- cbind(Age = A.pt, A.eff)[order(A.pt),]
    Per <- cbind(Per = P.pt, P.eff)[order(P.pt),]
    Coh <- cbind(Cph = C.pt, C.eff)[order(C.pt),]
    colnames(Age)[-1] <- c("A.eff", lu)
    colnames(Per)[-1] <- c("P.eff", lu)
    colnames(Coh)[-1] <- c("C.eff", lu)
    Type <- paste("Sequential modelling", Dist, ": (", parm, "):\n")
}
# If the model was binomial we convert to probabilities
o2p <- function(o) o/(1+o)
if( is.bin ) Age[,-1] <- o2p(Age[,-1]) 

res <- list(Type = Type,
           Model = Model,
             Age = Age,
             Per = Per,
             Coh = Coh,
           Drift = drift,
             Ref = c(Per = if ( parm %in%
                                c("APC","ADPC","Ad-P-C","AP-C") ) p0 else NA,
                     Coh = if ( parm %in%
                                c("ACP","ADCP","Ad-C-P","AC-P") ) c0 else NA ),
           Anova = AOV)
# If a spline model is used, add a "Knots" component to the apc-object
if (model %in% c("ns", "bs"))
    res <- c(res, list(Knots = Knots))
res$Age[, -1] <- res$Age[, -1] * scale
if (print.AOV) {
    print(res$Type)
    print(res$Anova)
}
# Print warnings about reference points:
if( !has.pref & parm %in% c("APC","ADPC") )
    cat( "No reference period given; ",
         "reference period for age-effects is chosen as\n",
         "the median date of event: ", p0, ".\n" )
if( !has.cref & parm %in% c("ACP","ADCP") )
    cat( "No reference cohort given; ",
         "reference cohort for age-effects is chosen as\n",
         "the median date of birth for persons  with event: ", c0, ".\n" )
class(res) <- "apc"
invisible(res)
}
