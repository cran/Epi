stackedCIF <- 
    function( x,
          group = 1,
            col = "black",
           fill = "white",
           ylim = c(0,1), 
           xlab = "Time",
           ylab = "Cumulative incidence", ...) 
{
    ne <- ncol(x$pstate) - 1
    ng <- length(x$n)
    if (group %in% 1:ng) {
        r1 <- ifelse(group == 1 | ng == 1, 1, cumsum(x$strata)[group - 
            1] + 1)
        r2 <- ifelse(ng == 1, length(x$time), cumsum(x$strata)[group])
        pSt0 <- matrix(0, nrow = nrow(x$pstate[r1:r2, ]), ncol = ne + 
            2)
        pSt0[, 1] <- x$pstate[r1:r2, 2]
        for (c in 2:ne) pSt0[, c] <- pSt0[, c - 1] + x$pstate[r1:r2, 
            c+1]
        pSt0[, ne + 1] <- 1
        pSt0[, ne + 2] <- x$time[r1:r2]
        pSt <- cbind(0, pSt0)
        pSt2 <- matrix(0, nrow = 2 * (dim(pSt)[1] - 1), ncol = dim(pSt)[2])
        pSt2[1, ne + 3] <- pSt[1, ne + 3]
        pSt2[nrow(pSt2), ne + 3] <- pSt[nrow(pSt), ne + 3]
        for (j in 1:(dim(pSt)[1] - 1)) {
            pSt2[2 * j - 1, ne + 3] <- pSt[j, ne + 3]
            pSt2[2 * j, ne + 3] <- pSt[j + 1, ne + 3]
            pSt2[2 * j - 1, 1:(ne + 2)] <- pSt[j, 1:(ne + 2)]
            pSt2[2 * j, 1:(ne + 2)] <- pSt[j, 1:(ne + 2)]
        }
        plot(as.numeric(pSt0[, ne + 2]), pSt0[, 2], type = "n", 
            ylim = ylim, yaxs = "i", xlab = xlab, ylab = ylab, 
            ...)
        for (i in 2:(ne + 1)) {
            polygon(c(pSt2[, ne + 3], rev(pSt2[, ne + 3])), c(pSt2[, 
                i], rev(pSt2[, i - 1])), col = fill[i - 1], border = NULL, 
                ...)
        }
        matlines( as.numeric(pSt0[, ne + 2]), pSt0[, 
            1:ne], type = "s", col = col)
        abline(v = max(x$time[r1:r2]), lwd = 2, col = "white")
        box()
    }
    else print(paste("Error: group indicator must be an integer from 1 to", 
        ng))
}
