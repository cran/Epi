legendbox <-
 function(x, y,
      state = "State",
         py = "Person-time",
      begin = "no. begin",
        end = "no. end",
      trans = "Transitions",
      rates = "\n(Rate)",
       font = 1,
      right = !left,
       left = !right,
        ...)
{
btxt <- paste0(state, "\n",
                  py, "\n",
               begin, "          ", end)
ww <- strwidth (btxt) * 1.2
hh <- strheight(btxt) * 1.3
zz <- tbox(btxt,
           x = x, y = y,
           wd = ww,
           ht = hh,
         font = font,
         ...)
if (missing(right) & missing(left)) right = TRUE
if (right) text(x + ww / 1.8, y, paste0(trans, rates), adj = 0)
if (left)  text(x - ww / 1.8, y, paste0(trans, rates), adj = 1)
}
