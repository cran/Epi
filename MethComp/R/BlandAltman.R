BlandAltman  <-
function(x, y,
         x.name=NULL,
         y.name=NULL,
         maintit="",
         cex=1,
         pch=16,
         col.points="black",
         col.lines="blue",
         limx=NULL,
         limy=NULL,
         ymax=NULL,
         eqax=FALSE,
         xlab=NULL,
         ylab=NULL,
         print=TRUE,
         conf.int=FALSE,
         reg.line=FALSE,
         digits=2,
         mult=FALSE,
         alpha=0.05,
         ... )
{
# Get names of supplied variables
x.nam <- deparse( substitute( x ) )
y.nam <- deparse( substitute( y ) )

# Check lengths of the supplied variables
if( length(x) != length(y) )
  stop( "\nx and y must have the same length:\n",
        "length(", x.nam, ")=", length(x), " and ",
        "length(", y.nam, ")=", length(y), " !")

# Get the naming of the variables and axes
if( is.null( x.name ) ) x.name <- x.nam
if( is.null( y.name ) ) y.name <- y.nam
if( mult )
  {
  x <- log(x)
  y <- log(y)
  if( is.null( xlab ) ) xlab <- paste( "Geometric mean( ",x.name, " , ", y.name, " )",sep="" )
  if( is.null( ylab ) ) ylab <- paste( x.name, "/", y.name )
  }
else
  {
  if( is.null( xlab ) ) xlab <- paste( "(", x.name, "+", y.name, ") / 2" )
  if( is.null( ylab ) ) ylab <- paste( x.name, "-", y.name )
  }

# The actual calculations
difference <- x-y								               # vector of differences
average    <- (x+y)/2                          # vector of means
n <- sum(!is.na(difference))						       # number of 'observations'
tvalue <- ifelse( missing(alpha), 2, qt(1-alpha/2,n-1)*(n+1)/n )
difference.mean <- mean(difference,na.rm=TRUE) # mean difference
difference.sd   <-   sd(difference,na.rm=TRUE) # SD of differences
al <- tvalue*difference.sd
upper.agreement.limit <- difference.mean+al	   # agreement limits
lower.agreement.limit <- difference.mean-al
difference.se <- difference.sd/sqrt(n)		     # standard error of the mean
al.se <- difference.sd*sqrt(3)/sqrt(n)		     # standard error of the agreement limit
pvalue <- pt(abs(difference.mean/difference.se),n-1,low=FALSE)*2
                                                # p value for H0: mean(diff)=0
difference.mean.ci <- difference.se*tvalue
al.ci <- al.se*tvalue
upper.agreement.limit.ci <- c(upper.agreement.limit-al.ci,
                              upper.agreement.limit+al.ci)
lower.agreement.limit.ci <- c(lower.agreement.limit-al.ci,
                              lower.agreement.limit+al.ci)

# Collect results in a matrix
res <- cbind( c( difference.mean,
                 difference.mean-difference.mean.ci,
                 difference.mean+difference.mean.ci ),
              c( lower.agreement.limit,
                 lower.agreement.limit.ci ),
              c( upper.agreement.limit,
                 upper.agreement.limit.ci ) )
colnames( res ) <- c(ylab,"lower limit","upper limit")
rownames( res ) <- c("Estimate", paste( round(100*   alpha/2 ,1), "%" ),
                                 paste( round(100*(1-alpha/2),1), "%" ) )
# names( dimnames(res) ) <- c( paste( "LoA for ",
#                                     x.name,
#                                     if( mult ) "/" else "-",
#                                     y.name, sep="" ), "" )

# The x and the y limits of the plot (limx, limy)
if( is.null(limx) ) limx <- range( average, na.rm=TRUE )
if( is.null(ymax) ) ymax <- max( c( abs( difference ),
                                      abs( res ),
                                      diff( limx )/2 ), na.rm=TRUE )
# Should the axes be of equal size?
if( eqax )
  {
  maxax <- max( diff(limx), 2*ymax )
  limy <- c(-1,1)*ifelse( maxax == diff(limx), maxax/2, ymax )
  if( maxax != diff(limx) ) limx <- mean(limx) + c(-1,1)*ymax
  }
else if( is.null( limy ) ) limy <- if( is.null(ymax) ) range( difference )
                                   else ymax * c(-1,1)
# If on a log-scale, transform back to display the results
if( mult )
  {
  res <- exp( res )
  average <- exp(average)
  difference <- exp(difference)
  limx <- exp(limx)
  limy <- exp(limy)
  }

# A function that gives the coordinates of the
# point (xf,yf) from ll corner in the current plot.
# if xf or yf are > 1 they are considered percentages
"cnr" <-
function( xf, yf )
{
cn <- par()$usr
xf <- ifelse( xf>1, xf/100, xf )
yf <- ifelse( yf>1, yf/100, yf )
xx <- ( 1 - xf ) * cn[1] + xf * cn[2]
yy <- ( 1 - yf ) * cn[3] + yf * cn[4]
if ( par()$xlog ) xx <- 10^xx
if ( par()$ylog ) yy <- 10^yy
list( x=xx, y=yy )
}

# Plot
plot.default( average, difference, type="n",
              log=if( mult ) "xy" else "",
              xlim=limx, ylim=limy,
              xlab=xlab, ylab=ylab, main=maintit, ... )
if( reg.line )
  {
  if( mult ) m0 <- lm( log10(difference) ~ log10(average) )
  else       m0 <- lm(       difference  ~       average  )
  # It must be log10, because those are the units the
  # plot is referred to when using abline
  alfa  <- coef(m0)[1]
  beta  <- coef(m0)[2]
  sigma <- summary(m0)$sigma
  abline( alfa             , beta, lwd=2, col=col.lines )
  abline( alfa+tvalue*sigma, beta, col=col.lines )
  abline( alfa-tvalue*sigma, beta, col=col.lines )
  if( is.numeric( reg.line ) )
    {
    # Write the regression equations based on regression of Differences on
    # to character objects for printing / plotting
    if( mult )
      {
      dif.eq <-
      paste( x.name,"/", y.name, " = ",
             formatC( 10^alfa, format="f", digits=reg.line ),
             "(", x.name, "*", y.name, ")^",
             formatC( beta/2, format="f", digits=reg.line ),
             " (", paste((1-alpha)*100), "% err.fact.",
             formatC( 10^(sigma*tvalue), format="f", digits=reg.line ),
             ")\n", sep="" )
      y.x.eq <-
      paste( y.name, " = ",
             formatC( 10^(alfa/(1-beta/2)), format="f", digits=reg.line ),
             "(", x.name, ")^",
             formatC( (1+beta/2)/(1-beta/2), format="f", digits=reg.line ),
             " (", paste((1-alpha)*100), "% err.fact.",
             formatC( 10^(sigma*tvalue/(1-beta/2)), format="f", digits=reg.line ),
             ")", sep="" )
      x.y.eq <-
      paste( x.name, " = ",
             formatC( 10^(-alfa/(1+beta/2)), format="f", digits=reg.line ),
             "(", y.name, ")^",
             formatC( (1-beta/2)/(1+beta/2), format="f", digits=reg.line ),
             " (", paste((1-alpha)*100), "% err.fact.",
             formatC( 10^(sigma*tvalue/(1+beta/2)), format="f", digits=reg.line ),
             ")", sep="" )
      }
    else
      {
      dif.eq <-
      paste( x.name,"-", y.name, " = ",
             formatC( alfa, format="f", digits=reg.line ),
             if(beta>0) " + " else " - ",
             formatC( abs(beta), format="f", digits=reg.line ),
             " (", x.name, "+", y.name, ")/2",
             " (", paste((1-alpha)*100), "% p.i.: +/-",
             formatC( sigma*tvalue, format="f", digits=reg.line ),
             ")", sep="" )
      x.y.eq <-
      paste( x.name, " = ",
             formatC( alfa/(1-beta/2), format="f", digits=reg.line ),
             " + ",
             formatC( (1+beta/2)/(1-beta/2), format="f", digits=reg.line ),
             " ", y.name,
             " (", paste((1-alpha)*100), "% p.i.: +/-",
             formatC( sigma*tvalue/(1-beta/2), format="f", digits=reg.line ),
             ")", sep="" )
      y.x.eq <-
      paste( y.name, " = ",
             formatC( -alfa/(1+beta/2), format="f", digits=reg.line ),
             " + ",
             formatC( (1-beta/2)/(1+beta/2), format="f", digits=reg.line ),
             " ", x.name,
             " (", paste((1-alpha)*100), "% p.i.: +/-",
             formatC( sigma*tvalue/(1+beta/2), format="f", digits=reg.line ),
             ")", sep="" )
      }
    text( cnr( 95, 5 ), dif.eq, adj=1 )
    text( cnr( 95,95 ), y.x.eq, adj=1 )
    text( cnr( 95,90 ), x.y.eq, adj=1 )
    }
  }
abline(h=res[1,],col=ifelse(reg.line,"transparent",col.lines),lwd=2)
if( conf.int )
  abline(h=res[2:3,],lty="24",col=ifelse(reg.line,"transparent",col.lines))
points( average, difference,cex=cex,pch=pch,col=col.points )
axis( side=4, at=res[1,],
      col.axis=ifelse(reg.line,"transparent",col.lines),
           col=ifelse(reg.line,"transparent",col.lines),
      labels=formatC( res[1,], format="f", digits=digits ), las=1 )
abline( h=ifelse(mult,1,0) )
box(bty=par("bty"))


# Print results
if( print )
  {
  cat( "\nLimits of agreement",
       if( conf.int ) paste( "with", round(100*(1-alpha),1),
                             "% confidence intervals"), ":\n", sep="" )
  if( conf.int ) print( res ) else print( res[1,] )
  cat("\n")
  if( is.numeric( reg.line ) )
    {
    cat( dif.eq, "\n",
         paste(
     "res.sd =", formatC( summary(m0)$sigma    , format="f", digits=reg.line ),
"   se(beta) =", formatC( summary(m0)$coef[2,2], format="f", digits=reg.line ),
        ", P =", formatC( summary(m0)$coef[2,4], format="f", digits=4 ) ), "\n\n",
        y.x.eq, "\n",
        x.y.eq, "\n" )
    }
  }

# Return list of relevant results
invisible( list( lim.agree=res,
                 p.value=pvalue,
                 if( reg.line ) reg=1:3  ) )
}

