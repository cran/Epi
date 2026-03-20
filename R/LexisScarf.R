LexisScarf <-
function(nw = 4,   # fields around
         nl = 20,  # fields along
         nm = 23,  # no. stitches per field
         np = 30,  # no. persons
      scale = 1.5, # maximal FU in units of nw
        clr = c("red", "blue", "limegreen")
         )
{
# everything is measured in units of stitches
wd <- nm * nw
lg <- nm * nl + 1
mm <- matrix(0, wd, lg)
x <- col(mm)
y <- row(mm)
d <- mm

# here is the the grid
d[, 0:nl * nm + 1] <- 1
d[0:(nw - 1) * nm + 1, ] <- 1

# pdf("lexis-scarf.pdf", height = 10, width = 50)
# plot the grid
par(mar = rep(0, 4))
plot(as.vector(x),
     as.vector(y), pch = 16,
     col = gray((c(8, 4) / 10)[as.vector(d) + 1])
     )
# entry points
ew <- floor(runif(np) * wd)
el <- floor(runif(np) * lg)
# life lengths
ll <- floor(runif(np) * wd * scale)
cl <- clr[sample(1:length(clr), np, replace = TRUE)]
# function to generate life lines
LL <- function(xx, # x (length) at start
               yy, # y (width) at start
               ll) # length of FU
      {
      # folding the life lines
      tt <- list(y = ifelse((yy - 0:ll) <= 0,
                       wd + (yy - 0:ll),
                             yy - 0:ll),
                 x = xx + 0:ll)
      # cut the parts outside the scarf
      del <- tt$y < 0  |
             tt$y > wd |
             tt$x < 0  |
             tt$x > lg
      tt$x <- tt$x[!del]
      tt$y <- tt$y[!del]
      tt
      }
# points as life lines for the persons
for(i in 1:np)
    points(LL(el[i],
              ew[i],
              ll[i]),
        col = cl[i], pch = 16)
NULL
# dev.off()
}
