adtte2Lexis_usub <-
function(dat, start.state = "rand", Tr.ev, Term.st = NULL)
{
   # the crap to pass the check
   ADT <- ADT2 <- AGE <- CNSR <- Count <- EVNTDESC <- STARTDT <-
   USUBJID <- end.date <- end.state <- startdat <- NULL

  `%>%` <- magrittr::`%>%`

  dat <- as.data.frame(dat)
  dat2 <- dat %>%
    dplyr::group_by(USUBJID) %>%
    dplyr::mutate(startdat = min(STARTDT),
                      ADT2 = Epi::cal.yr(as.Date(ADT, format="%Y-%m-%d"))) %>%
    dplyr::select(USUBJID, ADT2, startdat, AGE, EVNTDESC, CNSR)

  exitstat <- dat2 %>%
    dplyr::group_by(USUBJID) %>%
    dplyr::arrange(ADT2) %>%
    dplyr::slice(dplyr::n())

  
  HH <- dat2[dat2$EVNTDESC %in% Tr.ev & dat2$CNSR!=1,]
  HH <- HH %>%
    dplyr::group_by(USUBJID, EVNTDESC) %>%
    dplyr::mutate( Count = dplyr::row_number(),
                   statecount = paste(EVNTDESC,  dplyr::row_number(), sep="_"),
                   lex.id = USUBJID)

  HHlast <- HH %>%
    dplyr::group_by(USUBJID) %>%
    dplyr::arrange(ADT2) %>%
    dplyr::slice(dplyr::n())

  exitstat$end.state <- HHlast$statecount
  exitstat$end.date <- exitstat$ADT2

  Lx <- Epi::Lexis(data = exitstat,
              entry = list(per = Epi::cal.yr(startdat, format="%Y-%m-%d"),
                           t.sE = 0,
                           age = AGE),
              exit  = list(per = end.date),
              id    = USUBJID,
              entry.status = factor(start.state),
              exit.status = end.state,
              notes = FALSE)

  HHc <- dat2[dat2$EVNTDESC %in% Tr.ev & dat2$CNSR==1,]
  HHc <- HHc %>%
    dplyr::group_by(USUBJID, CNSR) %>%
    dplyr::filter(dplyr::row_number(CNSR) == 1) %>%
    dplyr::mutate(Count = dplyr::row_number(),
                  statecount = "Censored",
                  lex.id = USUBJID)

  HH <- rbind(HH, HHc)

  HH1 <- subset(HH, Count == 1)
  HH2 <- subset(HH, Count != 1)

  if(nrow(HH1)!=0 & nrow(HH2)==0){
    HH1 <- HH1[, c("lex.id", "ADT2", "statecount")]
    colnames(HH1) <- c("lex.id", "cut", "new.state")

    }
  if(nrow(HH1)!=0 & nrow(HH2)!=0){
  HH1$stat <- HH1$statecount#
  HH1 <- HH1[, c("lex.id", "ADT2", "statecount")]
  HH2 <- HH2[, c("lex.id", "ADT2", "statecount")]

  colnames(HH1) <- colnames(HH2) <- c("lex.id", "cut", "new.state")
  }

  if(nrow(HH1)!=0 & nrow(HH2)==0){
   for(i in 1:length(unique(HH1$new.state))){
      Lx <- Epi::cutLexis(Lx,
                  cut= HH1[HH1$new.state == unique(HH1$new.state)[i],],
                  timescale = "per",
                  new.scale = paste(paste0("t.s", unique(HH1$new.state)[i])))
   }
  }
  if(nrow(HH1)!=0 & nrow(HH2)!=0){
     HH3 <- rbind(HH1, HH2)
     HH3 <- HH3[order(HH3$cut),]
     for(i in 1:length(unique(HH3$new.state))){
       Lx <- Epi::cutLexis(Lx,
                      cut= HH3[HH3$new.state == unique(HH3$new.state)[i],],
                      timescale = "per",
                      new.scale = ifelse(substr(unique(HH3$new.state)[i],
                                                nchar(unique(HH3$new.state)[i]),
                                                      nchar(unique(HH3$new.state)[i]))=="1",
                                         paste0("t.s", unique(HH1$new.state)[i]) ,FALSE))

     }
  }

 Lx2 <- subset(Lx, select = -c(EVNTDESC, AGE, CNSR, ADT2, startdat, end.state, end.date))
 if(!is.null(Term.st)){
 tsT <- noquote(paste("t.s" ,paste(Term.st, "_1", sep=""), sep=""))
 Lx2[,!(colnames(Lx2) %in% tsT)]
 }

  return(Lx2)
}

adtte2Lexis <- function(dat, start.state = "Rand", Tr.ev, Term.st = NULL)
{
  # the crap to pass the check
  USUBJID <- NULL
  `%>%` <- magrittr::`%>%`
    lst <- dat %>% dplyr::group_by(USUBJID) %>%
    dplyr::group_map(~.x, .keep=TRUE)
 do.call("rbind", lapply(lst, function(x)adtte2Lexis_usub(x,start.state, Tr.ev,
                                                          Term.st)))
}








