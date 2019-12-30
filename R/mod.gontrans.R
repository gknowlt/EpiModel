#' @title STI Transmission Module
#'
#' @description Stochastically simulates GC transmission given the current
#'              state of the edgelist.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
gon_trans <- function(dat, at) {
  
  # Parameters ----------------------------------------------------------
  
  # Acquisition probabilities given contact with infected man
  #first letter indicates the uninfected site, while second letter indicates infected site
  rugc.tprob <- dat$param$rugc.tprob
  rpgc.tprob <- dat$param$rpgc.tprob
  
  urgc.tprob <- dat$param$urgc.tprob
  upgc.tprob <- dat$param$upgc.tprob
  
  prgc.tprob <- dat$param$prgc.tprob
  pugc.tprob <- dat$param$pugc.tprob
  
  # Probability of symptoms given infection
  rgc.sympt.prob <- dat$param$rgc.sympt.prob
  ugc.sympt.prob <- dat$param$ugc.sympt.prob
  pgc.sympt.prob <- dat$param$pgc.sympt.prob
  
  # Relative risk of infection given condom use during act
  sti.cond.rr <- dat$param$sti.cond.rr
  
  # Cessation
  gc.prob.cease <- dat$param$gc.prob.cease
  
  # Attributes ----------------------------------------------------------
  
  # Current infection state
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  pGC <- dat$attr$pGC
  
  # n Times infected
  rGC.timesInf <- dat$attr$rGC.timesInf
  uGC.timesInf <- dat$attr$uGC.timesInf
  pGC.timesInf <- dat$attr$pGC.timesInf
  
  # Set disease status to 0 for new births
  newBirths <- which(dat$attr$arrival.time == at)
  rGC[newBirths] <- rGC.timesInf[newBirths] <- 0
  uGC[newBirths] <- uGC.timesInf[newBirths] <- 0
  pGC[newBirths] <- pGC.timesInf[newBirths] <- 0
  
  
  # Infection time
  rGC.infTime <- dat$attr$rGC.infTime
  uGC.infTime <- dat$attr$uGC.infTime
  pGC.infTime <- dat$attr$pGC.infTime
  
  
  # Infection symptoms (non-varying)
  rGC.sympt <- dat$attr$rGC.sympt
  uGC.sympt <- dat$attr$uGC.sympt
  pGC.sympt <- dat$attr$pGC.sympt
  
  # Men who cease sexual activity during symptomatic infection
  GC.cease <- dat$attr$GC.cease
  
  # Pull act list
  al <- dat$temp$al
  
  ###Need to figure out how to adapt this to 3 types of sex acts
  ###Code defining al (act list) is in mod.condoms and then mod.position 
  ## ins variable coding
  # ins = 0 : p2 is insertive
  # ins = 1 : p1 is insertive
  # ins = 2 : both p1 and p2 are insertive
  
  # Rectal GC ---------------------------------------------
  
  # Rectal infection from urethra -----------------------------------------------------------
  # Requires: uGC in insertive man, and no rGC in receptive man
  p1Inf_rugc <- which(uGC[al[, "p1"]] == 1 & uGC.infTime[al[, "p1"]] < at &
                       rGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(1, 2))
  p2Inf_rugc <- which(uGC[al[, "p2"]] == 1 & uGC.infTime[al[, "p2"]] < at &
                       rGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(0, 2))
  allActs_rugc <- c(p1Inf_rugc, p2Inf_rugc)
  
  # UAI modifier
  uai_rugc <- al[, "uai"][allActs_rugc]
  tprob_rugc <- rep(rugc.tprob, length(allActs_rugc))
  tprob_rugc[uai_rugc == 0] <- tprob_rugc[uai_rugc == 0] * sti.cond.rr
  
  # Stochastic transmission
  trans_rugc <- rbinom(length(allActs_rugc), 1, tprob_rugc)
  
  # Determine the infected partner
  idsInf_rugc <- NULL
  if (sum(trans_rugc) > 0) {
    transAL_rugc <- al[allActs_rugc[trans_rugc == 1], , drop = FALSE]
    idsInf_rugc <- unique(ifelse(uGC[transAL_rugc[, "p1"]] == 1,
                                transAL_rugc[, "p2"], transAL_rugc[, "p1"]))
  }
  
  # Rectal infection from pharynx -----------------------------------------------------------
  # Requires: pGC in insertive man, and no rGC in receptive man
  p1Inf_rpgc <- which(pGC[al[, "p1"]] == 1 & pGC.infTime[al[, "p1"]] < at &
                        rGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(1, 2))
  p2Inf_rpgc <- which(pGC[al[, "p2"]] == 1 & pGC.infTime[al[, "p2"]] < at &
                        rGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(0, 2))
  allActs_rpgc <- c(p1Inf_rpgc, p2Inf_rpgc)
  
  # UAI modifier
  uai_rpgc <- al[, "uai"][allActs_rpgc]
  tprob_rpgc <- rep(rpgc.tprob, length(allActs_rpgc))
  tprob_rpgc[uai_rpgc == 0] <- tprob_rpgc[uai_rpgc == 0] * sti.cond.rr
  
  # Stochastic transmission
  trans_rpgc <- rbinom(length(allActs_rpgc), 1, tprob_rpgc)
  
  # Determine the infected partner
  idsInf_rpgc <- NULL
  if (sum(trans_rpgc) > 0) {
    transAL_rpgc <- al[allActs_rpgc[trans_rpgc == 1], , drop = FALSE]
    idsInf_rpgc <- unique(ifelse(pGC[transAL_rpgc[, "p1"]] == 1,
                                 transAL_rpgc[, "p2"], transAL_rpgc[, "p1"]))
  }

  # Update attributes
  idsInf_rgc <- unique(c(idsInf_rugc, idsInf_rpgc))
  rGC[idsInf_rgc] <- 1
  rGC.infTime[idsInf_rgc] <- at
  rGC.sympt[idsInf_rugc] <- rbinom(length(idsInf_rugc), 1, rgc.sympt.prob)
  rGC.timesInf[idsInf_rgc] <- rGC.timesInf[idsInf_rgc] + 1
  
  
  # Urethral GC ---------------------------------------------------------
  
  # Urethral infection from rectum
  # Requires: rGC in receptive man, and no uGC in insertive man
  p1Inf_urgc <- which(rGC[al[, "p1"]] == 1 & rGC.infTime[al[, "p1"]] < at &
                       uGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(0, 2))
  p2Inf_urgc <- which(rGC[al[, "p2"]] == 1 & rGC.infTime[al[, "p2"]] < at &
                       uGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(1, 2))
  allActs_urgc <- c(p1Inf_urgc, p2Inf_urgc)
  
  # UAI modifier
  uai_urgc <- al[, "uai"][allActs_urgc]
  tprob_urgc <- rep(urgc.tprob, length(allActs_urgc))
  tprob_urgc[uai_urgc == 0] <- tprob_urgc[uai_urgc == 0] * sti.cond.rr
  
  # Stochastic transmission
  trans_urgc <- rbinom(length(allActs_urgc), 1, tprob_urgc)
  
  # Determine the newly infected partner
  idsInf_urgc <- NULL
  if (sum(trans_urgc) > 0) {
    transAL_urgc <- al[allActs_urgc[trans_urgc == 1],  , drop = FALSE]
    idsInf_urgc <- unique(ifelse(rGC[transAL_urgc[, "p1"]] == 1,
                                transAL_urgc[, "p2"], transAL_urgc[, "p1"]))
  }

  # Urethral infection from pharynx
  # Requires: pGC in receptive man, and no uGC in insertive man
  p1Inf_upgc <- which(pGC[al[, "p1"]] == 1 & pGC.infTime[al[, "p1"]] < at &
                       uGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(0, 2))
  p2Inf_upgc <- which(pGC[al[, "p2"]] == 1 & pGC.infTime[al[, "p2"]] < at &
                       uGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(1, 2))
  allActs_upgc <- c(p1Inf_upgc, p2Inf_upgc)
  
  # UAI modifier
  uai_upgc <- al[, "uai"][allActs_upgc]
  tprob_upgc <- rep(upgc.tprob, length(allActs_upgc))
  tprob_upgc[uai_upgc == 0] <- tprob_upgc[uai_upgc == 0] * sti.cond.rr
  
  # Stochastic transmission
  trans_upgc <- rbinom(length(allActs_upgc), 1, tprob_upgc)
  
  # Determine the newly infected partner
  idsInf_upgc <- NULL
  if (sum(trans_upgc) > 0) {
    transAL_upgc <- al[allActs_upgc[trans_upgc == 1],  , drop = FALSE]
    idsInf_upgc <- unique(ifelse(pGC[transAL_upgc[, "p1"]] == 1,
                                transAL_upgc[, "p2"], transAL_upgc[, "p1"]))
  }
    
  # Update attributes
  idsInf_ugc <- unique(c(idsInf_urgc, idsInf_upgc))
  uGC[idsInf_ugc] <- 1
  uGC.infTime[idsInf_ugc] <- at
  uGC.sympt[idsInf_ugc] <- rbinom(length(idsInf_ugc), 1, ugc.sympt.prob)
  uGC.timesInf[idsInf_ugc] <- uGC.timesInf[idsInf_ugc] + 1  
  
  # Pharyngeal GC ------------------------------------------------------
  
  # pharyngeal infection from urethra
  # Requires: uGC in insertive man, and no pGC in receptive man
  p1Inf_pugc <- which(uGC[al[, "p1"]] == 1 & uGC.infTime[al[, "p1"]] < at &
                       pGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(1, 2))
  p2Inf_pugc <- which(uGC[al[, "p2"]] == 1 & uGC.infTime[al[, "p2"]] < at &
                       pGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(0, 2))
  allActs_pugc <- c(p1Inf_pugc, p2Inf_pugc)
  
  # UAI modifier
  uai_pugc <- al[, "uai"][allActs_pugc]
  tprob_pugc <- rep(pugc.tprob, length(allActs_pugc))
  tprob_pugc[uai_pugc == 0] <- tprob_pugc[uai_pugc == 0] * sti.cond.rr
  
  # Stochastic transmission
  trans_pugc <- rbinom(length(allActs_pugc), 1, tprob_pugc)
  
  # Determine the newly infected partner
  idsInf_pugc <- NULL
  if (sum(trans_pugc) > 0) {
    transAL_pugc <- al[allActs_pugc[trans_pugc == 1],  , drop = FALSE]
    idsInf_pugc <- unique(ifelse(uGC[transAL_pugc[, "p1"]] == 1,
                                transAL_pugc[, "p2"], transAL_pugc[, "p1"]))
  }

  # pharyngeal infection from rectum
  # Requires: rGC in receptive man, and no pGC in insertive man
  p1Inf_prgc <- which(rGC[al[, "p1"]] == 1 & rGC.infTime[al[, "p1"]] < at &
                        pGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(0, 2))
  p2Inf_prgc <- which(rGC[al[, "p2"]] == 1 & rGC.infTime[al[, "p2"]] < at &
                        pGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(1, 2))
  allActs_prgc <- c(p1Inf_prgc, p2Inf_prgc)
  
  # UAI modifier
  uai_prgc <- al[, "uai"][allActs_prgc]
  tprob_prgc <- rep(prgc.tprob, length(allActs_prgc))
  tprob_prgc[uai_prgc == 0] <- tprob_prgc[uai_prgc == 0] * sti.cond.rr
  
  # Stochastic transmission
  trans_prgc <- rbinom(length(allActs_prgc), 1, tprob_prgc)
  
  # Determine the newly infected partner
  idsInf_prgc <- NULL
  if (sum(trans_prgc) > 0) {
    transAL_prgc <- al[allActs_prgc[trans_prgc == 1],  , drop = FALSE]
    idsInf_prgc <- unique(ifelse(rGC[transAL_prgc[, "p1"]] == 1,
                                 transAL_prgc[, "p2"], transAL_prgc[, "p1"]))
  }
  
  # Update attributes
  idsInf_pgc <- unique(c(idsInf_prgc, idsInf_pugc))
  pGC[idsInf_pgc] <- 1
  pGC.infTime[idsInf_pgc] <- at
  pGC.sympt[idsInf_pgc] <- rbinom(length(idsInf_pgc), 1, pgc.sympt.prob)
  pGC.timesInf[idsInf_pgc] <- pGC.timesInf[idsInf_pgc] + 1
  
  # Set activity cessation attribute for newly infected -----------------
  
  # Symptomatic GC
  GC.sympt <- which(is.na(GC.cease) & (rGC.sympt == 1 | uGC.sympt == 1 | pGC.sympt == 1))
  idsGC.cease <- GC.sympt[which(rbinom(length(GC.sympt),
                                       1, gc.prob.cease) == 1)]
  GC.cease[GC.sympt] <- 0
  GC.cease[idsGC.cease] <- 1
  
  # Output --------------------------------------------------------------
  
  # attributes
  dat$attr$rGC <- rGC
  dat$attr$uGC <- uGC
  dat$attr$pGC <- pGC
  
  dat$attr$rGC.infTime <- rGC.infTime
  dat$attr$uGC.infTime <- uGC.infTime
  dat$attr$pGC.infTime <- pGC.infTime
  
  dat$attr$rGC.timesInf <- rGC.timesInf
  dat$attr$uGC.timesInf <- uGC.timesInf
  dat$attr$pGC.timesInf <- pGC.timesInf

  
  dat$attr$rGC.sympt <- rGC.sympt
  dat$attr$uGC.sympt <- uGC.sympt
  dat$attr$pGC.sympt <- pGC.sympt
  
  dat$attr$GC.cease <- GC.cease
  
  # Summary stats
  dat$epi$incid.rgc[at] <- length(idsInf_rgc)
  dat$epi$incid.ugc[at] <- length(idsInf_ugc)
  dat$epi$incid.pgc[at] <- length(idsInf_pgc)
  dat$epi$incid.gc[at] <- length(idsInf_rgc) + length(idsInf_ugc) + length(idsInf_pgc)
  
  # Check all infected have all STI attributes
  stopifnot(all(!is.na(rGC.infTime[rGC == 1])),
            all(!is.na(rGC.sympt[rGC == 1])),
            all(!is.na(uGC.infTime[uGC == 1])),
            all(!is.na(uGC.sympt[uGC == 1])),
            all(!is.na(pGC.infTime[pGC == 1])),
            all(!is.na(pGC.sympt[pGC == 1])))
            
  if (is.null(dat$epi$times.rgc)) {
    dat$epi$times.rgc <- rep(NA, length(dat$epi$num))
    dat$epi$times.ugc <- rep(NA, length(dat$epi$num))
    dat$epi$times.pgc <- rep(NA, length(dat$epi$num))
    }
  dat$epi$times.rgc[at] <- mean(rGC.timesInf, na.rm = TRUE)
  dat$epi$times.ugc[at] <- mean(uGC.timesInf, na.rm = TRUE)
  dat$epi$times.pgc[at] <- mean(pGC.timesInf, na.rm = TRUE)
  return(dat)
}


#' @title STI Recovery Module
#'
#' @description Stochastically simulates GC/CT recovery.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
sti_recov <- function(dat, at) {
  
  # Parameters ----------------------------------------------------------
  
  rgc.asympt.int <- dat$param$rgc.asympt.int
  ugc.asympt.int <- dat$param$ugc.asympt.int
  pgc.asympt.int <- dat$param$pgc.asympt.int
  gc.tx.int <- dat$param$gc.tx.int
  gc.ntx.int <- dat$param$gc.ntx.int

  
  # GC Recovery ---------------------------------------------------------
  
  # Asymptomatic untreated
  idsRGC_asympt_ntx <- which(dat$attr$rGC == 1 &
                               dat$attr$rGC.infTime < at &
                               dat$attr$rGC.sympt == 0 &
                               (is.na(dat$attr$rGC.tx) | dat$attr$rGC.tx == 0))
  idsUGC_asympt_ntx <- which(dat$attr$uGC == 1 &
                               dat$attr$uGC.infTime < at &
                               dat$attr$uGC.sympt == 0 &
                               (is.na(dat$attr$uGC.tx) | dat$attr$uGC.tx == 0))
  idsPGC_asympt_ntx <- which(dat$attr$pGC == 1 &
                               dat$attr$pGC.infTime < at &
                               dat$attr$pGC.sympt == 0 &
                               (is.na(dat$attr$pGC.tx) | dat$attr$pGC.tx == 0))
  
  recovRGC_asympt_ntx <- idsRGC_asympt_ntx[which(rbinom(length(idsRGC_asympt_ntx), 1,
                                                        1/rgc.asympt.int) == 1)]
  recovUGC_asympt_ntx <- idsUGC_asympt_ntx[which(rbinom(length(idsUGC_asympt_ntx), 1,
                                                        1/ugc.asympt.int) == 1)]
  recovPGC_asympt_ntx <- idsPGC_asympt_ntx[which(rbinom(length(idsPGC_asympt_ntx), 1,
                                                        1/pgc.asympt.int) == 1)]
  
  # Symptomatic untreated
  idsRGC_sympt_ntx <- which(dat$attr$rGC == 1 &
                              dat$attr$rGC.infTime < at &
                              dat$attr$rGC.sympt == 1 &
                              (is.na(dat$attr$rGC.tx) | dat$attr$rGC.tx == 0))
  idsUGC_sympt_ntx <- which(dat$attr$uGC == 1 &
                              dat$attr$uGC.infTime < at &
                              dat$attr$uGC.sympt == 1 &
                              (is.na(dat$attr$uGC.tx) | dat$attr$uGC.tx == 0))
  idsPGC_sympt_ntx <- which(dat$attr$pGC == 1 &
                              dat$attr$pGC.infTime < at &
                              dat$attr$pGC.sympt == 1 &
                              (is.na(dat$attr$pGC.tx) | dat$attr$pGC.tx == 0))
  
  # If parameter is null, uses recovery rate of asytomatic untreated
  if (!is.na(gc.ntx.int)) {
    recovRGC_sympt_ntx <- idsRGC_sympt_ntx[which(rbinom(length(idsRGC_sympt_ntx), 1,
                                                        1/gc.ntx.int) == 1)]
    recovUGC_sympt_ntx <- idsUGC_sympt_ntx[which(rbinom(length(idsUGC_sympt_ntx), 1,
                                                        1/gc.ntx.int) == 1)]
    recovPGC_sympt_ntx <- idsPGC_sympt_ntx[which(rbinom(length(idsPGC_sympt_ntx), 1,
                                                        1/gc.ntx.int) == 1)]
  } else {
    recovRGC_sympt_ntx <- idsRGC_sympt_ntx[which(rbinom(length(idsRGC_sympt_ntx), 1,
                                                        1/rgc.asympt.int) == 1)]
    recovUGC_sympt_ntx <- idsUGC_sympt_ntx[which(rbinom(length(idsUGC_sympt_ntx), 1,
                                                        1/ugc.asympt.int) == 1)]
    recovPGC_sympt_ntx <- idsPGC_sympt_ntx[which(rbinom(length(idsPGC_sympt_ntx), 1,
                                                        1/pgc.asympt.int) == 1)]
  }
  ###stopped here
  # Treated (asymptomatic and symptomatic)
  idsRGC_tx <- which(dat$attr$rGC == 1 &
                       dat$attr$rGC.infTime < at &
                       (dat$attr$rGC.tx == 1 | dat$attr$rGC.tx.prep == 1))
  idsUGC_tx <- which(dat$attr$uGC == 1 &
                       dat$attr$uGC.infTime < at &
                       (dat$attr$uGC.tx == 1 | dat$attr$uGC.tx.prep == 1))
  
  recovRGC_tx <- idsRGC_tx[which(rbinom(length(idsRGC_tx), 1,
                                        1/gc.tx.int) == 1)]
  recovUGC_tx <- idsUGC_tx[which(rbinom(length(idsUGC_tx), 1,
                                        1/gc.tx.int) == 1)]
  
  recovRGC <- c(recovRGC_asympt_ntx, recovRGC_sympt_ntx, recovRGC_tx)
  recovUGC <- c(recovUGC_asympt_ntx, recovUGC_sympt_ntx, recovUGC_tx)
  
  dat$attr$rGC[recovRGC] <- 0
  dat$attr$rGC.sympt[recovRGC] <- NA
  dat$attr$rGC.infTime[recovRGC] <- NA
  dat$attr$rGC.tx[recovRGC] <- NA
  dat$attr$rGC.tx.prep[recovRGC] <- NA
  
  dat$attr$uGC[recovUGC] <- 0
  dat$attr$uGC.sympt[recovUGC] <- NA
  dat$attr$uGC.infTime[recovUGC] <- NA
  dat$attr$uGC.tx[recovUGC] <- NA
  dat$attr$uGC.tx.prep[recovUGC] <- NA
  
  dat$attr$GC.cease[c(recovRGC, recovUGC)] <- NA
  
  
  
  # CT Recovery ---------------------------------------------------------
  
  # Asymptomatic untreated
  idsRCT_asympt_ntx <- which(dat$attr$rCT == 1 &
                               dat$attr$rCT.infTime < at &
                               dat$attr$rCT.sympt == 0 &
                               (is.na(dat$attr$rCT.tx) | dat$attr$rCT.tx == 0) &
                               (is.na(dat$attr$rCT.tx.prep) | dat$attr$rCT.tx.prep == 0))
  idsUCT_asympt_ntx <- which(dat$attr$uCT == 1 &
                               dat$attr$uCT.infTime < at &
                               dat$attr$uCT.sympt == 0 &
                               (is.na(dat$attr$uCT.tx) | dat$attr$uCT.tx == 0) &
                               (is.na(dat$attr$uCT.tx.prep) | dat$attr$uCT.tx.prep == 0))
  
  recovRCT_asympt_ntx <- idsRCT_asympt_ntx[which(rbinom(length(idsRCT_asympt_ntx),
                                                        1, 1/rct.asympt.int) == 1)]
  recovUCT_asympt_ntx <- idsUCT_asympt_ntx[which(rbinom(length(idsUCT_asympt_ntx),
                                                        1, 1/uct.asympt.int) == 1)]
  
  # Symptomatic untreated
  idsRCT_sympt_ntx <- which(dat$attr$rCT == 1 &
                              dat$attr$rCT.infTime < at &
                              dat$attr$rCT.sympt == 1 &
                              (is.na(dat$attr$rCT.tx) | dat$attr$rCT.tx == 0) &
                              (is.na(dat$attr$rCT.tx.prep) | dat$attr$rCT.tx.prep == 0))
  idsUCT_sympt_ntx <- which(dat$attr$uCT == 1 &
                              dat$attr$uCT.infTime < at &
                              dat$attr$uCT.sympt == 1 &
                              (is.na(dat$attr$uCT.tx) | dat$attr$uCT.tx == 0) &
                              (is.na(dat$attr$uCT.tx.prep) | dat$attr$uCT.tx.prep == 0))
  
  if (!is.na(ct.ntx.int)) {
    recovRCT_sympt_ntx <- idsRCT_sympt_ntx[which(rbinom(length(idsRCT_sympt_ntx),
                                                        1, 1/ct.ntx.int) == 1)]
    recovUCT_sympt_ntx <- idsUCT_sympt_ntx[which(rbinom(length(idsUCT_sympt_ntx),
                                                        1, 1/ct.ntx.int) == 1)]
  } else {
    recovRCT_sympt_ntx <- idsRCT_sympt_ntx[which(rbinom(length(idsRCT_sympt_ntx),
                                                        1, 1/rct.asympt.int) == 1)]
    recovUCT_sympt_ntx <- idsUCT_sympt_ntx[which(rbinom(length(idsUCT_sympt_ntx),
                                                        1, 1/uct.asympt.int) == 1)]
  }
  
  # Treated (asymptomatic and symptomatic)
  idsRCT_tx <- which(dat$attr$rCT == 1 &
                       dat$attr$rCT.infTime < at &
                       (dat$attr$rCT.tx == 1 | dat$attr$rCT.tx.prep == 1))
  idsUCT_tx <- which(dat$attr$uCT == 1 &
                       dat$attr$uCT.infTime < at &
                       (dat$attr$uCT.tx == 1 | dat$attr$uCT.tx.prep == 1))
  
  recovRCT_tx <- idsRCT_tx[which(rbinom(length(idsRCT_tx),
                                        1, 1/ct.tx.int) == 1)]
  recovUCT_tx <- idsUCT_tx[which(rbinom(length(idsUCT_tx),
                                        1, 1/ct.tx.int) == 1)]
  
  
  recovRCT <- c(recovRCT_asympt_ntx, recovRCT_sympt_ntx, recovRCT_tx)
  recovUCT <- c(recovUCT_asympt_ntx, recovUCT_sympt_ntx, recovUCT_tx)
  
  dat$attr$rCT[recovRCT] <- 0
  dat$attr$rCT.sympt[recovRCT] <- NA
  dat$attr$rCT.infTime[recovRCT] <- NA
  dat$attr$rCT.tx[recovRCT] <- NA
  dat$attr$rCT.tx.prep[recovRCT] <- NA
  
  dat$attr$uCT[recovUCT] <- 0
  dat$attr$uCT.sympt[recovUCT] <- NA
  dat$attr$uCT.infTime[recovUCT] <- NA
  dat$attr$uCT.tx[recovUCT] <- NA
  dat$attr$uCT.tx.prep[recovUCT] <- NA
  
  dat$attr$CT.cease[c(recovRCT, recovUCT)] <- NA
  
  # Summary stats
  dat$epi$recov.rgc[at] <- length(unique(recovRGC))
  dat$epi$recov.ugc[at] <- length(unique(recovUGC))
  dat$epi$recov.rct[at] <- length(unique(recovRCT))
  dat$epi$recov.uct[at] <- length(unique(recovUCT))
  
  return(dat)
}


#' @title STI Treatment Module
#'
#' @description Stochastically simulates GC/CT diagnosis and treatment.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
sti_tx <- function(dat, at) {
  
  # Parameters
  gc.sympt.prob.tx <- dat$param$gc.sympt.prob.tx
  ct.sympt.prob.tx <- dat$param$ct.sympt.prob.tx
  
  gc.asympt.prob.tx <- dat$param$gc.asympt.prob.tx
  ct.asympt.prob.tx <- dat$param$ct.asympt.prob.tx
  
  prep.sti.screen.int <- dat$param$prep.sti.screen.int
  prep.sti.prob.tx <- dat$param$prep.sti.prob.tx
  
  prep.cont.stand.tx <- dat$param$prep.continue.stand.tx
  if (prep.cont.stand.tx == TRUE) {
    prep.stand.tx.grp <- 0:1
  } else {
    prep.stand.tx.grp <- 0
  }
  
  # symptomatic gc treatment
  idsRGC_tx_sympt <- which(dat$attr$rGC == 1 &
                             dat$attr$rGC.infTime < at &
                             dat$attr$rGC.sympt == 1 &
                             is.na(dat$attr$rGC.tx) &
                             dat$attr$prepStat %in% prep.stand.tx.grp)
  idsUGC_tx_sympt <- which(dat$attr$uGC == 1 &
                             dat$attr$uGC.infTime < at &
                             dat$attr$uGC.sympt == 1 &
                             is.na(dat$attr$uGC.tx) &
                             dat$attr$prepStat %in% prep.stand.tx.grp)
  idsGC_tx_sympt <- c(idsRGC_tx_sympt, idsUGC_tx_sympt)
  
  txGC_sympt <- idsGC_tx_sympt[which(rbinom(length(idsGC_tx_sympt), 1,
                                            gc.sympt.prob.tx) == 1)]
  txRGC_sympt <- intersect(idsRGC_tx_sympt, txGC_sympt)
  txUGC_sympt <- intersect(idsUGC_tx_sympt, txGC_sympt)
  
  # asymptomatic gc treatment
  idsRGC_tx_asympt <- which(dat$attr$rGC == 1 &
                              dat$attr$rGC.infTime < at &
                              dat$attr$rGC.sympt == 0 &
                              is.na(dat$attr$rGC.tx) &
                              dat$attr$prepStat %in% prep.stand.tx.grp)
  idsUGC_tx_asympt <- which(dat$attr$uGC == 1 &
                              dat$attr$uGC.infTime < at &
                              dat$attr$uGC.sympt == 0 &
                              is.na(dat$attr$uGC.tx) &
                              dat$attr$prepStat %in% prep.stand.tx.grp)
  idsGC_tx_asympt <- c(idsRGC_tx_asympt, idsUGC_tx_asympt)
  
  txGC_asympt <- idsGC_tx_asympt[which(rbinom(length(idsGC_tx_asympt), 1,
                                              gc.asympt.prob.tx) == 1)]
  txRGC_asympt <- intersect(idsRGC_tx_asympt, txGC_asympt)
  txUGC_asympt <- intersect(idsUGC_tx_asympt, txGC_asympt)
  
  # all treated GC
  txRGC <- union(txRGC_sympt, txRGC_asympt)
  txUGC <- union(txUGC_sympt, txUGC_asympt)
  
  idsRGC_tx <- union(idsRGC_tx_sympt, idsRGC_tx_asympt)
  idsUGC_tx <- union(idsUGC_tx_sympt, idsUGC_tx_asympt)
  
  
  # symptomatic ct treatment
  idsRCT_tx_sympt <- which(dat$attr$rCT == 1 &
                             dat$attr$rCT.infTime < at &
                             dat$attr$rCT.sympt == 1 &
                             is.na(dat$attr$rCT.tx) &
                             dat$attr$prepStat %in% prep.stand.tx.grp)
  idsUCT_tx_sympt <- which(dat$attr$uCT == 1 &
                             dat$attr$uCT.infTime < at &
                             dat$attr$uCT.sympt == 1 &
                             is.na(dat$attr$uCT.tx) &
                             dat$attr$prepStat %in% prep.stand.tx.grp)
  idsCT_tx_sympt <- c(idsRCT_tx_sympt, idsUCT_tx_sympt)
  
  txCT_sympt <- idsCT_tx_sympt[which(rbinom(length(idsCT_tx_sympt), 1,
                                            ct.sympt.prob.tx) == 1)]
  txRCT_sympt <- intersect(idsRCT_tx_sympt, txCT_sympt)
  txUCT_sympt <- intersect(idsUCT_tx_sympt, txCT_sympt)
  
  # asymptomatic ct treatment
  idsRCT_tx_asympt <- which(dat$attr$rCT == 1 &
                              dat$attr$rCT.infTime < at &
                              dat$attr$rCT.sympt == 0 &
                              is.na(dat$attr$rCT.tx) &
                              dat$attr$prepStat == 0)
  idsUCT_tx_asympt <- which(dat$attr$uCT == 1 &
                              dat$attr$uCT.infTime < at &
                              dat$attr$uCT.sympt == 0 &
                              is.na(dat$attr$uCT.tx) &
                              dat$attr$prepStat == 0)
  idsCT_tx_asympt <- c(idsRCT_tx_asympt, idsUCT_tx_asympt)
  
  txCT_asympt <- idsCT_tx_asympt[which(rbinom(length(idsCT_tx_asympt), 1,
                                              ct.asympt.prob.tx) == 1)]
  txRCT_asympt <- intersect(idsRCT_tx_asympt, txCT_asympt)
  txUCT_asympt <- intersect(idsUCT_tx_asympt, txCT_asympt)
  
  # all treated CT
  txRCT <- union(txRCT_sympt, txRCT_asympt)
  txUCT <- union(txUCT_sympt, txUCT_asympt)
  
  idsRCT_tx <- union(idsRCT_tx_sympt, idsRCT_tx_asympt)
  idsUCT_tx <- union(idsUCT_tx_sympt, idsUCT_tx_asympt)
  
  # Interval-based treatment for MSM on PrEP
  idsSTI_screen <- which(dat$attr$prepStartTime == at |
                           (at - dat$attr$prepLastStiScreen >= prep.sti.screen.int))
  
  dat$attr$prepLastStiScreen[idsSTI_screen] <- at
  
  
  idsRGC_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$rGC == 1 &
                                      dat$attr$rGC.infTime < at &
                                      is.na(dat$attr$rGC.tx.prep)))
  idsUGC_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$uGC == 1 &
                                      dat$attr$uGC.infTime < at &
                                      is.na(dat$attr$uGC.tx.prep)))
  idsRCT_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$rCT == 1 &
                                      dat$attr$rCT.infTime < at &
                                      is.na(dat$attr$rCT.tx.prep)))
  idsUCT_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$uCT == 1 &
                                      dat$attr$uCT.infTime < at &
                                      is.na(dat$attr$uCT.tx.prep)))
  
  txRGC_prep <- idsRGC_prep_tx[which(rbinom(length(idsRGC_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txUGC_prep <- idsUGC_prep_tx[which(rbinom(length(idsUGC_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txRCT_prep <- idsRCT_prep_tx[which(rbinom(length(idsRCT_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txUCT_prep <- idsUCT_prep_tx[which(rbinom(length(idsUCT_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  
  
  # update attributes
  dat$attr$rGC.tx[idsRGC_tx] <- 0
  dat$attr$rGC.tx[txRGC] <- 1
  
  dat$attr$uGC.tx[idsUGC_tx] <- 0
  dat$attr$uGC.tx[txUGC] <- 1
  
  dat$attr$rCT.tx[idsRCT_tx] <- 0
  dat$attr$rCT.tx[txRCT] <- 1
  
  dat$attr$uCT.tx[idsUCT_tx] <- 0
  dat$attr$uCT.tx[txUCT] <- 1
  
  dat$attr$rGC.tx.prep[idsRGC_prep_tx] <- 0
  dat$attr$rGC.tx.prep[txRGC_prep] <- 1
  
  dat$attr$uGC.tx.prep[idsUGC_prep_tx] <- 0
  dat$attr$uGC.tx.prep[txUGC_prep] <- 1
  
  dat$attr$rCT.tx.prep[idsRCT_prep_tx] <- 0
  dat$attr$rCT.tx.prep[txRCT_prep] <- 1
  
  dat$attr$uCT.tx.prep[idsUCT_prep_tx] <- 0
  dat$attr$uCT.tx.prep[txUCT_prep] <- 1
  
  
  # add tx at other site
  dat$attr$rGC.tx[which((dat$attr$uGC.tx == 1 | dat$attr$uGC.tx.prep == 1) & dat$attr$rGC == 1)] <- 1
  dat$attr$uGC.tx[which((dat$attr$rGC.tx == 1 | dat$attr$rGC.tx.prep == 1) & dat$attr$uGC == 1)] <- 1
  
  dat$attr$rCT.tx[which((dat$attr$uCT.tx == 1 | dat$attr$uCT.tx.prep == 1) & dat$attr$rCT == 1)] <- 1
  dat$attr$uCT.tx[which((dat$attr$rCT.tx == 1 | dat$attr$rCT.tx.prep == 1) & dat$attr$uCT == 1)] <- 1
  
  txRGC_all <- union(txRGC, txRGC_prep)
  txUGC_all <- union(txUGC, txUGC_prep)
  txRCT_all <- union(txRCT, txRCT_prep)
  txUCT_all <- union(txUCT, txUCT_prep)
  
  
  # summary stats
  if (is.null(dat$epi$num.asympt.tx)) {
    dat$epi$num.asympt.tx <- rep(NA, length(dat$epi$num))
    dat$epi$num.asympt.cases <- rep(NA, length(dat$epi$num))
    dat$epi$num.asympt.tx.prep <- rep(NA, length(dat$epi$num))
    dat$epi$num.asympt.cases.prep <- rep(NA, length(dat$epi$num))
    dat$epi$num.rect.tx <- rep(NA, length(dat$epi$num))
    dat$epi$num.rect.cases <- rep(NA, length(dat$epi$num))
    dat$epi$num.rect.tx.prep <- rep(NA, length(dat$epi$num))
    dat$epi$num.rect.cases.prep <- rep(NA, length(dat$epi$num))
  }
  
  asympt.tx <- c(intersect(txRGC_all, which(dat$attr$rGC.sympt == 0)),
                 intersect(txUGC_all, which(dat$attr$uGC.sympt == 0)),
                 intersect(txRCT_all, which(dat$attr$rCT.sympt == 0)),
                 intersect(txUCT_all, which(dat$attr$uCT.sympt == 0)))
  dat$epi$num.asympt.tx[at] <- length(unique(asympt.tx))
  asympt.cases <- c(idsRGC_tx_asympt, intersect(idsRGC_prep_tx, which(dat$attr$rGC.sympt == 0)),
                    idsUGC_tx_asympt, intersect(idsUGC_prep_tx, which(dat$attr$uGC.sympt == 0)),
                    idsRCT_tx_asympt, intersect(idsRCT_prep_tx, which(dat$attr$rCT.sympt == 0)),
                    idsUCT_tx_asympt, intersect(idsUCT_prep_tx, which(dat$attr$uCT.sympt == 0)))
  dat$epi$num.asympt.cases[at] <- length(unique(asympt.cases))
  
  
  asympt.tx.prep <- c(intersect(txRGC_prep, which(dat$attr$rGC.sympt == 0)),
                      intersect(txUGC_prep, which(dat$attr$uGC.sympt == 0)),
                      intersect(txRCT_prep, which(dat$attr$rCT.sympt == 0)),
                      intersect(txUCT_prep, which(dat$attr$uCT.sympt == 0)))
  dat$epi$num.asympt.tx.prep[at] <- length(unique(asympt.tx.prep))
  asympt.cases.prep <- c(intersect(idsRGC_prep_tx, which(dat$attr$rGC.sympt == 0)),
                         intersect(idsUGC_prep_tx, which(dat$attr$uGC.sympt == 0)),
                         intersect(idsRCT_prep_tx, which(dat$attr$rCT.sympt == 0)),
                         intersect(idsUCT_prep_tx, which(dat$attr$uCT.sympt == 0)))
  dat$epi$num.asympt.cases.prep[at] <- length(unique(asympt.cases.prep))
  
  
  rect.tx <- c(txRGC_all, txRCT_all)
  dat$epi$num.rect.tx[at] <- length(unique(rect.tx))
  rect.cases <- c(idsRGC_tx, idsRGC_prep_tx, idsRCT_tx, idsRCT_prep_tx)
  dat$epi$num.rect.cases[at] <- length(unique(rect.cases))
  
  rect.tx.prep <- c(txRGC_prep, txRCT_prep)
  dat$epi$num.rect.tx.prep[at] <- length(unique(rect.tx.prep))
  rect.cases.prep <- c(idsRGC_prep_tx, idsRCT_prep_tx)
  dat$epi$num.rect.cases.prep[at] <- length(unique(rect.cases.prep))
  
  return(dat)
}