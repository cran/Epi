fgrep <- function( pattern, x, ... )        x [grep( pattern,        x , ... )]
ngrep <- function( pattern, x, ... )  names(x)[grep( pattern,  names(x), ... )]
lgrep <- function( pattern, x, ... ) levels(x)[grep( pattern, levels(x), ... )]
