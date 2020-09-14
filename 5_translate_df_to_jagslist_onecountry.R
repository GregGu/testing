
translate_df_to_jagslist_onecountry <- function(observed_data_preprocessed,
                                     meta_filtered
                                     # ref_date # was used in vignette, see its usage below, however it is not included in the model, instead year.ref is used from hyper
                                     ){
  #to be renamed below, easiest to just use old names within function for now
  meta <- meta_filtered
  datall <- observed_data_preprocessed
  jagsdata <- fixed_parameter_list

  #-------------------------------------
  # 1. model par and baselist


  ##aids parameters
  jagsdata$caids = meta$caids
  jagsdata$uaids = meta$uaids
  jagsdata$kaids = meta$kaids

  jagsdata$year.ref <- which(meta$year.t==fixed_parameter_list$referenceyear)

  X.cth <- meta$X.cth
  #meta$X.cth[1,,]
  B.ct <- meta$births.ct
  isssa.c <- ifelse(meta$mdg.c =="Sub-Saharan Africa",1,0)

  baselist <- list(
    X.cth = X.cth,
    isssa.c=isssa.c,
    a.ct  = meta$a.ct,
    E.ct = meta$deaths.ct, B.ct = B.ct,
    R = meta$R, getr.c = meta$getr.c, nyears = meta$nyears, C = meta$C,
    pi.constant = pi,
    input0.ct = matrix(0, meta$C, meta$nyears),# included for restriction R <1
    v.ct = meta$v.ct,
    Raids.ct = meta$v.ct*meta$uaids,
    muaids.ct = meta$v.ct*meta$uaids/B.ct*meta$a.ct*meta$deaths.ct
  )

  #----
  # add extra uncertainty for incomplete vr or study envelopes
  set.seed(123456)
  logbeta <- rnorm(1000000, 0, 1) # choose an sd here so that range of ratios are reasonable…

  #---------------
  # specialized studies with reported envelope, index is jinq
  isjinq.d <- (datall$type=="inq"
               & !is.na(datall$final_pm)
               &!is.na(datall$final_env)
              & datall$modelinclude
              & datall$completeness_inq > 0.95)
  if (sum(isjinq.d) == 0){
    datinq <- NULL
  } else {
    datainq <- datall[isjinq.d,]
    J <- dim(datainq)[1]
    getc.j = Getc.i(iso.i = paste(datainq$iso), iso.c = meta$iso.c)
    gett.j = Gett.i(years.i = floor(datainq$year), year.t = meta$year.t)
    # get the reference period
    start.j <- datainq$start
    end.j <- datainq$end
    gettstart.j <- Gett.i(years.i = floor(start.j), year.t = meta$year.t)
    gettend.j <- Gett.i(years.i = ceiling(end.j-1), year.t = meta$year.t)
    X.j <- gettend.j - gettstart.j +1
    partialtime.xj <- partialwhoenv.xj <- matrix(NA, max(X.j),J)
    for (j in 1:J){
      partialtime.xj[1:X.j[j],j] <- GetPartialTime(start = start.j[j], end = end.j[j], X = X.j[j])
      partialwhoenv.xj[1:X.j[j],j] <- partialtime.xj[1:X.j[j],j]*meta$deaths.ct[getc.j[j], gettstart.j[j]:gettend.j[j]]
    }
    #selogpm.jinq <- datainq$stoch_selogpm
    #selogpmobs.jinq <- ifelse(selogpm.jinq > 0.5, 0.5, selogpm.jinq)

    datinq <-  list(
      isjinq.d = isjinq.d,
      #tau.jinq  = 1/selogpmobs.jinq^2,
      #logpm.jinq = log(datainq$final_pm),
      Jinq = J,
      getc.jinq = getc.j, gett.jinq = gett.j,
      gettstart.jinq = gettstart.j,
      gettend.jinq = gettend.j,
      # for now, ceil end and floor mat
      env.jinq = ceiling(datainq$final_env),
      mat.jinq  = floor(datainq$final_pm*datainq$final_env),
      X.jinq = X.j,
      partialtime.xjinq  = partialtime.xj,
      partialwhoenv.xjinq  = partialwhoenv.xj
    )
  }


  # specialized studies with incomplete envelope, index is g
  isjinq.g <- (datall$type=="inq"
               & !is.na(datall$final_pm)
               &!is.na(datall$final_env)
               & datall$modelinclude
               & datall$completeness_inq <= 0.95)
  if (sum(isjinq.g) == 0){
    datinq_incomplete <- NULL
  } else {
    datainq <- datall[isjinq.g,] # same name
    J <- dim(datainq)[1]
    getc.j = Getc.i(iso.i = paste(datainq$iso), iso.c = meta$iso.c)
    gett.j = Gett.i(years.i = floor(datainq$year), year.t = meta$year.t)
    # get the reference period
    start.j <- datainq$start
    end.j <- datainq$end
    gettstart.j <- Gett.i(years.i = floor(start.j), year.t = meta$year.t)
    gettend.j <- Gett.i(years.i = ceiling(end.j-1), year.t = meta$year.t)
    X.j <- gettend.j - gettstart.j +1
    partialtime.xj <- partialwhoenv.xj <- matrix(NA, max(X.j),J)
    for (j in 1:J){
      partialtime.xj[1:X.j[j],j] <- GetPartialTime(start = start.j[j], end = end.j[j], X = X.j[j])
      partialwhoenv.xj[1:X.j[j],j] <- partialtime.xj[1:X.j[j],j]*meta$deaths.ct[getc.j[j], gettstart.j[j]:gettend.j[j]]
    }
    #selogpm.jinq <- datainq$stoch_selogpm
    #selogpmobs.jinq <- ifelse(selogpm.jinq > 0.5, 0.5, selogpm.jinq)

    # add completeness
    varfrombeta.g <- rep(NA, J)
    for (j in 1:J){
      varfrombeta.g[j] <- var(1/(datainq$completeness_inq[j] +
                                   (1-datainq$completeness_inq[j])*exp(logbeta)))
    }

    datinq_incomplete <-  list(
      isjinq.g = isjinq.g,
      #tau.jinq  = 1/selogpmobs.jinq^2,
      #logpm.jinq = log(datainq$final_pm),
      Ginq_incomplete = J,
      getc.g = getc.j, #gett.jinq = gett.j, #dont think we need this since we use gettstart and its dupplicate, same name as in complete inq data
      gettstart.g = gettstart.j,
      gettend.g = gettend.j,
      # for now, ceil end and floor mat
      env.g = ceiling(datainq$final_env),
      mat.g  = floor(datainq$final_pm*datainq$final_env),
      X.g = X.j,
      partialtime.xg  = partialtime.xj,
      partialwhoenv.xg  = partialwhoenv.xj,
      varfrombeta.g = varfrombeta.g
    )
  }

  #-------------------------------------
  # for VR data



  isj.d<-datall$type=="vr" & datall$modelinclude
  if (sum(isj.d)==0){
    datvr <- NULL
  } else {
    datavr <- datall[isj.d,]
    getj.k1 <- seq(1, dim(datavr)[1])
    K1 <- length(getj.k1)
    datvrmultiplier <- list(
      getj.k1 = getj.k1,
      K1 = K1
      )

    varfrombeta <- rep(NA, dim(datavr)[1])
    for (j in 1:dim(datavr)[1]){
      varfrombeta[j] <- var(1/(datavr$rhovrfinal_bmat[j] + (1-datavr$rhovrfinal_bmat[j])*exp(logbeta)))
    }

    getc.j <- Getc.i(iso.i = paste(datavr$iso), iso.c = meta$iso.c)
    gett.j <- Gett.i(years.i = floor(datavr$year), year.t = meta$year.t)
    datvr <- list(
      isj.d = isj.d, ##save it for later postprocessing
      J=dim(datavr)[1],
      getc.j = getc.j,
      gett.j = gett.j,
      pmobs.j = ifelse(datavr$final_pm > 0, datavr$final_pm, 0.0001),
      logpmobs.j = log(ifelse(datavr$final_pm > 0, datavr$final_pm, 0.0001)),
      mat.j  = round(datavr$final_pm*datavr$obs_env),
      env.j  = ceiling(datavr$final_env),
      var_sens.j = datavr$var_sens,
      var_spec.j = datavr$var_spec,
      cov_sesp.j = datavr$cov_sesp,
      sens.j = datavr$sens,
      spec.j = datavr$spec,
      sens_sq.j = datavr$sens_sq,
      oneminspec_sq.j = datavr$oneminspec_sq,
      crvs_completeness.j = datavr$rhovrfinal_bmat,
      input0forRstar.ct = matrix(0, meta$C, meta$nyears),
      input0forVR.ct = matrix(0, meta$C, meta$nyears),
      varfrombeta.j = varfrombeta
    )
  }

  #------------------------------------------------------------------------
  # other studies, index is jnew
  isjnew.d<- isj.d==FALSE & isjinq.d==FALSE & !is.na(datall$final_pm) & datall$modelinclude
  if (sum(isjnew.d) == 0){
    datnonvr <- NULL
  } else {
    datanonvr <- datall[isjnew.d, ]
    J <- dim(datanonvr)[1]
    getc.j = Getc.i(iso.i = paste(datanonvr$iso), iso.c = meta$iso.c)
    gett.j = Gett.i(years.i = floor(datanonvr$year), year.t = meta$year.t)
    start.j <- datanonvr$start
    end.j <- datanonvr$end
    gettstart.j <- Gett.i(years.i = floor(start.j), year.t = meta$year.t)
    gettend.j <- Gett.i(years.i = ceiling(end.j-1), year.t = meta$year.t)
    X.j <- gettend.j - gettstart.j +1
    partialtime.xj <- partialwhoenv.xj <- matrix(NA, max(X.j),J)
    for (j in 1:J){
      partialtime.xj[1:X.j[j],j] <- GetPartialTime(start = start.j[j], end = end.j[j], X = X.j[j])
      partialwhoenv.xj[1:X.j[j],j] <- partialtime.xj[1:X.j[j],j]*meta$deaths.ct[getc.j[j], gettstart.j[j]:gettend.j[j]]
    }
    aids.jnew <- rep(NA, J)
    for (j in 1:J){
      aids.jnew[j] <- baselist$v.ct[getc.j[j],gett.j[j]]*meta$a.ct[getc.j[j],gett.j[j]]
    }
    selogpm.jnew <- datanonvr$obs_selogpm
    selogpm.jnew[is.na(selogpm.jnew)] <- fixed_parameter_list$imputeSElogPM

    # data types are misc, with subtype dhs
    ismisc.jnew  <- datanonvr$type!="inq"
    isinq.jnew <- datanonvr$type=="inq"
    isdhs.jnew <- datanonvr$type=="dhs"

    getjnew.m <- which(ismisc.jnew)
    M <- length(getjnew.m)

    datnonvr <-  list(
      ismisc.jnew = ismisc.jnew,
      isinq.jnew = isinq.jnew,
      isdhs.jnew = isdhs.jnew,
      isjnew.d = isjnew.d,
      tausamp.jnew  = 1/selogpm.jnew^2,
      logpm.jnew = log(datanonvr$final_pm),
      isssa.jnew = isssa.c[getc.j],
      Jnew = J,
      getc.jnew = getc.j, gett.jnew = gett.j,
      ispreg.jnew = ifelse(datanonvr$definition=="pregn",1,0),
      gettstart.jnew = gettstart.j,
      gettend.jnew = gettend.j,
      X.jnew = X.j,
      partialtime.xjnew  = partialtime.xj,
      partialwhoenv.xjnew  = partialwhoenv.xj,
      aidsva.jnew = aids.jnew)
  }

  #-------------------------------------
  jagsdata <- c(jagsdata, baselist, datvr, datinq, datinq_incomplete, datnonvr)

  # for constraints
  jagsdata$input1.ct <- matrix(1, meta$C, meta$nyears)
  jagsdata$input1again.ct <- jagsdata$input1.ct
  # jagsdata$tref.c <- which(meta$year.t == ref_date)
  jagsdata_vr <- fixed_parameter_list %>% .[c(
                  "etaworld.k",
                  "rsigma.alpha",
                  "rsigma.beta",
                  "rrho.alphabeta",
                  "sigma.alpha",
                  "sigma.beta",
                  "rho.alphabeta",
                  "specmin",
                  "sensmin")]
  jagsdata_vr["K"] <- 2
  jagsdata <- c(jagsdata, jagsdata_vr)
  return(jagsdata)
} # end function

#----------------------
# The End!
