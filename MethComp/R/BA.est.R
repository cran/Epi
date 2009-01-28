BA.est <-
function( data,
        linked = TRUE,
          exch = !linked,
            MI = TRUE # To fit the model with replicates method by item interaction
      # , plot = TRUE # To be used for calling some kind of plotting function,
                      # possibly a variant of plotMeth or generalized version
                      # of BA.plot.
        )
{
# Check that data has item, method and repl
rq.nam <- c("meth","item","repl","y")
if( sum( !is.na( wh <- match( rq.nam, names( data ) ) ) ) < 3 ) stop(
"\nThe supplied dataframe misses columns named ", rq.nam[is.na(wh)], ".\n" )
if( sum( !is.na( wh <- match( rq.nam, names( data ) ) ) ) == 3 ) stop(
"\nThe supplied dataframe misses the column named ", rq.nam[is.na(wh)], ".\n" )

# Package needed for the fitting of the models
require( nlme )

# Sort out a possible discrepancy if exch is given
linked <- !exch

# Make sure that all variables are factors and only has levels actually present
data$meth <- factor( data$meth )
data$item <- factor( data$item )
data$repl <- factor( data$repl )

# More than two methods?
Nm   <- nlevels( data$meth )
Mnam <-  levels( data$meth )

if( MI )
  {
  if( linked )
    {
    if( Nm ==2 )
      m1 <- lme( y ~ meth + item,
                 random=list( item = pdIdent( ~ meth-1 ),
                              repl = ~1 ),
                 weights = varIdent( form = ~1 | meth ),
                 data = data,
               control = lmeControl(returnObject=TRUE) )
    if( Nm > 2 )
      m1 <- lme( y ~ meth + item,
                 random = list( item = pdDiag( ~ meth-1 ),
                              repl = ~1 ),
                 weights = varIdent( form = ~1 | meth ),
                 data = data,
               control = lmeControl(returnObject=TRUE) )
    }
  else
      m1 <- lme( y ~ meth + item,
                 random = list( item = pdDiag( ~ meth-1 ) ),
                 weights = varIdent( form = ~1 | meth ),
                 data = data,
               control = lmeControl(returnObject=TRUE) )
  }
else
  {
  if( linked )
    m1 <- lme( y ~ meth + item,
               random = list( item = pdIdent( ~ repl-1 ) ),
               weights = varIdent( form = ~1 | meth ),
               data = data,
               control = lmeControl(returnObject=TRUE) )

  else
    m1 <- lme( y ~ meth + item,
               random = ~ 1 | one,
               weights = varIdent( form = ~1 | meth ),
               data = cbind(data,one=1),
               control = lmeControl(returnObject=TRUE) )
  }
# Fixed effects estimates
summ <- summary( m1 )$tT
# Extract the estimates of the biases
bias <- c( 0, summ[grep("meth",rownames(summ)),1] )
names( bias ) <- levels( data$meth )

# The two-way random interactions
vc <- VarCorr( m1 )
if(           MI )   tau <- as.numeric( vc[grep("meth",rownames(vc)),2] )
if( linked &  MI ) omega <- as.numeric( vc[grep("Inte",rownames(vc)),2] )
if( linked & !MI ) omega <- as.numeric( vc[grep("repl",rownames(vc)),2][1] )

# The residual variances
sig <- attr(m1$residuals,"std")
sigma <- tapply( sig, names(sig), unique )
vc <- c( if(MI) tau, if(linked) omega, sigma )
names( vc ) <- c( if(MI) paste("MxI",levels(data$meth),sep="."),
                  if(linked) "IxR",
                  paste("res",Mnam,sep=".") )
vcmp <- cbind( vc, vc^2 )
colnames( vcmp ) <- c("SD","Var")
                  
# The limits of agreement
LoA <- matrix( NA, Nm*(Nm+1)/2, 4 )
colnames( LoA ) <- c("Mean","Lower","Upper", "SD")
rownames( LoA ) <- 1:nrow(LoA)
row.no <- 0
for( i in 1:Nm ) for( j in 1:i )
{
  row.no <- row.no + 1
  rownames( LoA )[row.no] <- paste( Mnam[i], "-", Mnam[j], " " )
  LoA[row.no,1] <- bias[i] - bias[j]
  pred.var <- sigma[i]^2 + sigma[j]^2
  if( i!=j & MI ) pred.var <- pred.var + tau[i]^2 + tau[j]^2
  LoA[row.no,4] <- sqrt( pred.var )
  LoA[row.no,2] <- LoA[row.no,1] - 2.00*LoA[row.no,4]
  LoA[row.no,3] <- LoA[row.no,1] + 2.00*LoA[row.no,4]
}
diags <- cumsum(1:Nm)
RC <- cbind(SD=LoA[diags,4],Rep.coef=2.00*sqrt(2)*LoA[diags,4])
rownames( RC ) <- Mnam
list( Bias = bias,
  Var.comp = vcmp,
       LoA = LoA[-diags,],
  Rep.coef = RC )
}