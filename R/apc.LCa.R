###################################################################################
### apc-LCa comparison
apc.LCa <-
function( data,               # cohort reference for the interactions
   keep.models = FALSE, ... )
{
models <- c("APa","ACa","APaC","APCa","APaCa" )
LCa.list <- list()
length( LCa.list ) <- 5
names( LCa.list ) <- models
for( mod in models ){ cat( mod, ":\n" )
   LCa.list[[mod]] <- LCa.fit( data = data,
                              model = mod, ... ) }
APC <- apc.fit( data, npar = list( A=LCa.list[[2]]$knots$a.kn,
                                   P=LCa.list[[1]]$knots$p.kn,
                                   C=LCa.list[[1]]$knots$c.kn ) )
dev <- c( APC$Anova[c(2,5,3,4),2],
          sapply( LCa.list, function(x) x$deviance ) )
df  <- c( APC$Anova[c(2,5,3,4),1],
          sapply( LCa.list, function(x) x$df ) )
names(dev)[1:4] <- names(df)[1:4] <-
gsub( "rift","", gsub("eriod","", gsub("ohort","", gsub("-","",
gsub( "ge", "", rownames(APC$Anova)[c(2,5,3,4)])))))
if( keep.models ) return( list( dev = cbind( dev, df ),
                                apc = APC,
                                LCa = LCa.list ) )
             else return( cbind( dev, df ) )
}

show.apc.LCa <-
function( x, dev.scale=TRUE, top="Ad", ... )
{
if( is.list(x) ) x <- x[[1]]     
TM <- matrix(NA,9,9)
rownames( TM ) <- colnames( TM ) <-
paste( rownames(x), "\n", formatC(x[,1],format="f",digits=1) )
TM[1,2:3] <-
TM[2,c(4,5)] <-
TM[3,c(4,6)] <-
TM[4,c(7,8)] <-
TM[5,7] <-
TM[6,8] <-
TM[c(7,8),9] <- 1
TM
bp <- list( x=c(50,30,70,50,10,90,30,70,50),
            y=c(90,70,70,50,50,50,30,30,10) )
if( !dev.scale )
  boxes.matrix( TM, boxpos=bp, hm=1.5, wm=1.5, ... )
else {
  bp$y <- 5+90*(pmin(x[top,1],x[,1])-x[9,1])/(x[top,1]-x[9,1])
  boxes.matrix( TM, boxpos=bp, hm=1.5, wm=1.5, ... )
  }
}
