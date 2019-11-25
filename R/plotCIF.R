plotCIF <- 
function (x,
      event = 1,
       xlab = "Time",
       ylab = "Cumulative incidence", 
       ylim = c(0, 1),
        lty = 1,
        col = "black",
        ... ) 
{
    ng <- length(x$n)
    co <- rep(col, ng)[1:ng]
    lt <- rep(lty, ng)[1:ng]
    ne <- dim(x$pstate)[2] - 1
    if (event %in% 1:ne) {
        time <- x$time
        gr <- rep(1, length(time))
        if (ng > 1) 
            for (g in 2:ng) for (t in (cumsum(x$strata)[g - 1] + 
                1):cumsum(x$strata)[g]) gr[t] <- g
        CI <- x$pstate[, 2]
        for (e in 1:ne) if (event == e) 
            CI <- x$pstate[, e+1]
        plot(c(0, time[gr == 1], max(time[gr == 1])), c(0, CI[gr == 
            1], max(CI[gr == 1])), type = "s", ylim = ylim, xlab = xlab, 
            ylab = ylab, col = co[1], lty = lt[1], ...)
        if (ng > 1) 
            for (g in 2:ng) lines(c(0, time[gr == g], max(time[gr == 
                g])), c(0, CI[gr == g], max(CI[gr == g])), type = "s", 
                lty = lt[g], col = co[g], ...)
    }
    else print(paste("Error: event must be an integer from 1 to", 
        ne))
}
