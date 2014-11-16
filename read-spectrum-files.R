
###################################################
####  function to read HIV projection outputs  ####
###################################################

read.hivproj.output.5.0 <- function(specdp.file){

  ## read .DP file
  dp <- read.csv(specdp.file, as.is=TRUE)


  ## projection parameters
  proj.tidx <- which(dp[,1] == "<General 3>")
  aids5.tidx <- which(dp[,1] == "<AIDS5>")
  yr.start <- as.numeric(dp[proj.tidx+6,4])
  yr.end <- as.numeric(dp[proj.tidx+8,4])
  proj.years <- yr.start:yr.end
  t0 <- as.numeric(dp[aids5.tidx+2,4])
  timedat.idx <- 4+1:length(proj.years)-1

  agegr.lab <- c(paste(0:15*5, 1:16*5, sep="-"), "80+")
  

  ## Number HIV+
  hivnum.tidx <- which(dp[,2] == "HIV")
  hivnum.age <- sapply(dp[hivnum.tidx+4:54, timedat.idx], as.numeric)
  rownames(hivnum.age) <- dp[hivnum.tidx+4:54, 3]
  hivnum.m <- hivnum.age[grep("Sex=1", rownames(hivnum.age)),]
  hivnum.f <- hivnum.age[grep("Sex=2", rownames(hivnum.age)),]
  dimnames(hivnum.m) <- dimnames(hivnum.f) <- list(agegr.lab, proj.years)

  ## Number new infections
  newinf.tidx <- which(dp[,2] == "New Infections")
  newinf.age <- sapply(dp[newinf.tidx+4:54, timedat.idx], as.numeric)
  rownames(newinf.age) <- dp[newinf.tidx+4:54, 3]
  newinf.m <- newinf.age[grep("Sex=1", rownames(newinf.age)),]
  newinf.f <- newinf.age[grep("Sex=2", rownames(newinf.age)),]
  dimnames(newinf.m) <- dimnames(newinf.f) <- list(agegr.lab, proj.years)
  
  ## Total population size
  totpop.tidx <- which(dp[,2] == "Total Population")
  totpop.m <- sapply(dp[totpop.tidx + 1:17*7 + 6, timedat.idx], as.numeric)
  totpop.f <- sapply(dp[totpop.tidx + 1:17*7 + 8, timedat.idx], as.numeric)
  dimnames(totpop.m) <- dimnames(totpop.f) <- list(agegr.lab, proj.years)
  
  ## ART need
  artneed.tidx <- which(dp[,2] == "Need FL")
  artneed.m <- sapply(dp[artneed.tidx+0:16*3+5, timedat.idx], as.numeric)
  artneed.f <- sapply(dp[artneed.tidx+0:16*3+6, timedat.idx], as.numeric)
  dimnames(artneed.m) <- dimnames(artneed.f) <- list(agegr.lab, proj.years)
  
  ## On ART
  artnum.tidx <- which(dp[,2] == "On FL")
  artnum.m <- sapply(dp[artnum.tidx+0:16*3+5, timedat.idx], as.numeric)
  artnum.f <- sapply(dp[artnum.tidx+0:16*3+6, timedat.idx], as.numeric)
  dimnames(artnum.m) <- dimnames(artnum.f) <- list(agegr.lab, proj.years)

  ## number deaths
  deathm.tidx <- which(dp[,2] == "Deaths - Male")
  deathf.tidx <- which(dp[,2] == "Deaths - Female")
  deaths.m <- sapply(dp[deathm.tidx+2:18, timedat.idx], as.numeric)
  deaths.f <- sapply(dp[deathf.tidx+2:18, timedat.idx], as.numeric)
  dimnames(deaths.m) <- dimnames(deaths.f) <- list(agegr.lab, proj.years)

  ## number AIDS deaths
  aidsdm.tidx <- which(dp[,2] == "AIDS deaths - Male")
  aidsdf.tidx <- which(dp[,2] == "AIDS deaths - Female")
  aidsdeaths.m <- sapply(dp[aidsdm.tidx+2:18, timedat.idx], as.numeric)
  aidsdeaths.f <- sapply(dp[aidsdf.tidx+2:18, timedat.idx], as.numeric)
  dimnames(aidsdeaths.m) <- dimnames(aidsdeaths.f) <- list(agegr.lab, proj.years)

  return(list("totpop.m" = totpop.m,
              "totpop.f" = totpop.f,
              "hivnum.m" = hivnum.m,
              "hivnum.f" = hivnum.f,
              "artnum.m" = artnum.m,
              "artnum.f" = artnum.f,
              "newinf.m" = newinf.m,
              "newinf.f" = newinf.f,
              "deaths.m" = deaths.m,
              "deaths.f" = deaths.f,
              "aidsdeaths.m" = aidsdeaths.m,
              "aidsdeaths.f" = aidsdeaths.f))
}


######################################################
####  function to read HIV projection parameters  ####
######################################################

read.hivproj.param.5.0 <- function(specdp.file){

  ## read .DP file
  dp <- read.csv(specdp.file, as.is=TRUE)

  ## find tag indexes (tidx)
  proj.tidx <- which(dp[,1] == "<General 3>")
  aids5.tidx <- which(dp[,1] == "<AIDS5>")
  nathist.tidx <- which(dp[,1] == "<AdultTransParam2>")
  cd4initdist.tidx <- which(dp[,1] == "<DistNewInfectionsCD4>")
  infectreduc.tidx <- which(dp[,1] == "<InfectReduc>")
  adult.artnumperc.tidx <- which(dp[,1] == "<HAARTBySexPerNum>")
  adult.art.tidx <- which(dp[,1] == "<HAARTBySex>")
  adult.arteligthresh.tidx <- which(dp[,1] == "<CD4ThreshHoldAdults>")

  ## state space dimensions
  NG <- 2
  AG <- 17
  DS <- 8
  TS <- 4

  ## projection parameters
  yr.start <- as.numeric(dp[proj.tidx+6,4])
  yr.end <- as.numeric(dp[proj.tidx+8,4])
  proj.years <- yr.start:yr.end
  t0 <- as.numeric(dp[aids5.tidx+2,4])
  timedat.idx <- 4+1:length(proj.years)-1

  ## scalar paramters
  relinfect.ART <- 1.0 - as.numeric(dp[infectreduc.tidx+1, 4])
  vert.trans <- 0  ## !! Spectrum vertical transmission not yet implemented
  fert.rat <- setNames(as.numeric(dp[aids5.tidx+185, 4+0:6]), seq(15, 45, 5))

  ## sex/age-specific incidence ratios (time varying)
  inc.sexrat <- setNames(as.numeric(dp[aids5.tidx+181,timedat.idx]), proj.years) # !!! Not sure what aids5.tidx+183 (Ratio of female to male prevalence)

  inc.agerat <- array(NA, c(AG, NG, length(proj.years)), list(0:(AG-1)*5, c("Male", "Female"), proj.years))
  inc.agerat[,"Male",] <- sapply(dp[aids5.tidx+281:297,timedat.idx], as.numeric)
  inc.agerat[,"Female",] <- sapply(dp[aids5.tidx+299:315,timedat.idx], as.numeric)

  ## hiv natural history
  cd4.initdist <- array(NA, c(DS-1, AG, NG), list(2:DS, 0:(AG-1)*5, c("Male", "Female")))
  cd4.initdist[,,"Male"] <- array(as.numeric(dp[cd4initdist.tidx+2, 4:31])/100, c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.initdist[,,"Female"] <- array(as.numeric(dp[cd4initdist.tidx+3, 4:31])/100, c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]

  cd4.prog <- array(NA, c(DS-2, AG, NG), list(2:(DS-1), 0:(AG-1)*5, c("Male", "Female")))
  cd4.prog[,,"Male"] <- array(1/as.numeric(dp[nathist.tidx+4, 4:27]), c(DS-2, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.prog[,,"Female"] <- array(1/as.numeric(dp[nathist.tidx+5, 4:27]), c(DS-2, 4))[,rep(1:4, c(5,2,2,8))]

  cd4.art.mort <- array(NA, c(TS, DS-1, AG, NG), list(1:TS, 2:DS, 0:(AG-1)*5, c("Male", "Female")))
  cd4.art.mort[1,,,"Male"] <- array(as.numeric(dp[nathist.tidx+7, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.art.mort[1,,,"Female"] <- array(as.numeric(dp[nathist.tidx+8, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.art.mort[2,,,"Male"] <- array(as.numeric(dp[nathist.tidx+10, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.art.mort[2,,,"Female"] <- array(as.numeric(dp[nathist.tidx+11, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.art.mort[3,,,"Male"] <- array(as.numeric(dp[nathist.tidx+13, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.art.mort[3,,,"Female"] <- array(as.numeric(dp[nathist.tidx+14, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.art.mort[4,,,"Male"] <- array(as.numeric(dp[nathist.tidx+16, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.art.mort[4,,,"Female"] <- array(as.numeric(dp[nathist.tidx+17, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]

  ## program parameters
  artnumperc.15plus <- sapply(dp[adult.artnumperc.tidx+3:4, timedat.idx], as.numeric)
  dimnames(artnumperc.15plus) <- list(c("Male", "Female"), proj.years)

  artnum.15plus <- sapply(dp[adult.art.tidx+3:4, timedat.idx], as.numeric)
  dimnames(artnum.15plus) <- list(c("Male", "Female"), proj.years)

  arteligthresh.15plus <- setNames(as.numeric(dp[adult.arteligthresh.tidx+2, timedat.idx]), proj.years)

  ## permute indices to match model state space (year, NG, AG, DS, TS)
  inc.agerat <- aperm(inc.agerat, 3:1)
  cd4.initdist <- aperm(cd4.initdist, 3:1)
  cd4.prog <- aperm(cd4.prog, 3:1)
  cd4.art.mort <- aperm(cd4.art.mort, 4:1)
  artnum.15plus <- aperm(artnum.15plus, 2:1)
  artnumperc.15plus <- aperm(artnumperc.15plus, 2:1)

  projp <- list("yr.start"=yr.start, "yr.end"=yr.end, "t0"=t0,
                "relinfect.ART"=relinfect.ART, "vert.trans"=vert.trans,
                "fert.rat"=fert.rat, "inc.sexrat"=inc.sexrat, "inc.agerat"=inc.agerat,
                "cd4.initdist"=cd4.initdist, "cd4.prog"=cd4.prog, "cd4.art.mort"=cd4.art.mort,
                "artnumperc.15plus"=artnumperc.15plus, "artnum.15plus"=artnum.15plus,
                "arteligthresh.15plus"=arteligthresh.15plus)
  return(projp)
}


###################################################################
####  function to read UN Population Division projection file  ####
###################################################################

## note: indices are hard coded -- better to parse files using <markers>

DEM.INIT.YEAR <- 1970

read.demog.param <- function(upd.file, age.intervals = 5){

  ## check age intervals and prepare age groups vector
  if(length(age.intervals) == 1){
    if(80 %% age.intervals != 0)
      stop("Invalid age interval interval")
    age.groups <- c(rep(seq(0, 79, age.intervals), each=age.intervals), 80)
    age.intervals <- c(rep(age.intervals, 80 / age.intervals), 1)
    ## } else if(length(age.groups) < 81 && age.groups[1] == 0 && tail(age.groups, 1) == 80){
    ##   stop("not defined yet") ## define thie one
  } else if(length(age.intervals) < 81){
    if(sum(age.intervals != 81))
      stop("Invalid vector of age intervals")
    age.groups <- rep(1:length(age.intervals), times=age.intervals)
  }

  ## Read and parse udp file
  upd <- read.csv(upd.file, header=FALSE, as.is=TRUE)
  
  bp.tidx <- which(upd[,1] == "<basepop>")
  lt.tidx <- which(upd[,1] == "<lfts>")
  pasfrs.tidx <- which(upd[,1] == "<pasfrs>")
  migration.tidx <- which(upd[,1] == "<migration>")
  tfr.tidx <- which(upd[,1] == "<tfr>")
  srb.tidx <- which(upd[,1] == "<srb>")
  
  bp <- setNames(data.frame(upd[bp.tidx+2:1459,1:4]), upd[bp.tidx+1,1:4])
  lt <- setNames(data.frame(upd[lt.tidx+2:13121,1:12]), upd[lt.tidx+1,1:12])
  pasfrs <- setNames(data.frame(upd[pasfrs.tidx+1+1:(80*35),1:3]), upd[pasfrs.tidx+1,1:3])
  migration <- setNames(data.frame(upd[migration.tidx+1+1:(80*2*81),1:4]), upd[migration.tidx+1,1:4])
  tfr <- setNames(as.numeric(upd[tfr.tidx+1+1:80,2]), upd[tfr.tidx+1+1:80,1])
  srb <- setNames(as.numeric(upd[srb.tidx+1+1:80,2]), upd[srb.tidx+1+1:80,1])


  ## Aggregate into specified age groups

  ## population size
  basepop <- array(as.numeric(bp$value), c(81, 2, 9))
  dimnames(basepop) <- list(0:80, c("Male", "Female"), seq(1970, 2010, 5))
  basepop <- apply(basepop, 2:3, tapply, age.groups, sum)

  ## mx
  Sx <- as.numeric(lt$Sx[-(1:(2*80)*82-1)]) # 80+ age group given twice
  dim(Sx) <- c(81, 2, 80)
  dimnames(Sx) <- list(0:80, c("Male", "Female"), 1970:2049)
  Sx <- apply(Sx, 2:3, tapply, age.groups, prod)
  mx <- -sweep(log(Sx), 1, age.intervals, "/")

  ## asfr
  asfd <- array(as.numeric(pasfrs$value), c(35, 80))
  dimnames(asfd) <- list(15:49, 1970:2049)
  asfr <- sweep(asfd, 2, tfr, "*")
  asfr <- apply(asfr, 2, tapply, age.groups[16:50], mean)

  ## migration
  netmigr <- array(as.numeric(migration$value), c(81, 2, 80))
  netmigr <- apply(netmigr, 2:3, tapply, age.groups, sum)

  ## permute indices to match model state space (year, NG, AG)
  basepop <- aperm(basepop, 3:1)
  mx <- aperm(mx, 3:1)
  Sx <- aperm(Sx, 3:1)
  asfr <- aperm(asfr, 2:1)
  netmigr <- aperm(netmigr, 3:1)

  demp <- list("basepop"=basepop, "mx"=mx, "Sx"=Sx, "asfr"=asfr, "srb"=srb, "netmigr"=netmigr)
  return(demp)
}


##################################################
####  Function to write outputs as CSV files  ####
##################################################

##  Outputs requested by Tim Hallett, 9 October 2014

write.outputs <- function(output.path, demp, projp, proj, out.yrs = as.character(1975:2013)){

  incid.m <- proj$newinf.m[,out.yrs]/(proj$totpop.m[,out.yrs] - proj$hivnum.m[,out.yrs])
  incid.f <- proj$newinf.f[,out.yrs]/(proj$totpop.f[,out.yrs] - proj$hivnum.f[,out.yrs])
  prev.m <- proj$hivnum.m[,out.yrs]/proj$totpop.m[,out.yrs]
  prev.f <- proj$hivnum.f[,out.yrs]/proj$totpop.f[,out.yrs]

  agegr.lab <- c(paste(0:15*5, 1:16*5, sep="-"), "80+")
  dimnames(demp$mx)[[3]] <- agegr.lab
  dimnames(demp$basepop)[[3]] <- agegr.lab
  dimnames(demp$asfr)[[2]] <- agegr.lab[4:10]

  write.csv(t(demp$basepop[min(out.yrs),,]), paste(output.path, "-basepop.csv", sep=""))
  write.csv(t(demp$asfr[out.yrs,]), paste(output.path, "-asfr.csv", sep=""))
  write.csv(t(demp$mx[out.yrs,"Male",]), paste(output.path, "-natmx-m.csv", sep=""))
  write.csv(t(demp$mx[out.yrs,"Female",]), paste(output.path, "-natmx-f.csv", sep=""))

  write.csv(proj$artnum.m[,out.yrs], paste(output.path, "-artnum-m.csv", sep=""))
  write.csv(proj$artnum.f[,out.yrs], paste(output.path, "-artnum-f.csv", sep=""))
  write.csv(incid.m[,out.yrs], paste(output.path, "-incid-m.csv", sep=""))
  write.csv(incid.f[,out.yrs], paste(output.path, "-incid-f.csv", sep=""))
  write.csv(prev.m[,out.yrs], paste(output.path, "-prev-m.csv", sep=""))
  write.csv(prev.f[,out.yrs], paste(output.path, "-prev-f.csv", sep=""))
  write.csv(proj$deaths.m[,out.yrs], paste(output.path, "-deaths-m.csv", sep=""))
  write.csv(proj$deaths.f[,out.yrs], paste(output.path, "-deaths-f.csv", sep=""))
  write.csv(proj$aidsdeaths.m[,out.yrs], paste(output.path, "-aidsdeaths-m.csv", sep=""))
  write.csv(proj$aidsdeaths.f[,out.yrs], paste(output.path, "-aidsdeaths-f.csv", sep=""))
  write.csv(proj$newinf.m[,out.yrs], paste(output.path, "-newinf-m.csv", sep=""))
  write.csv(proj$newinf.f[,out.yrs], paste(output.path, "-newinf-f.csv", sep=""))
  write.csv(proj$hivnum.m[,out.yrs], paste(output.path, "-hivnum-m.csv", sep=""))
  write.csv(proj$hivnum.f[,out.yrs], paste(output.path, "-hivnum-f.csv", sep=""))

  return(NULL)
}


#########################
####  run functions  ####
#########################

## Countries: South Africa, Botswana, Swaziland, Malawi

southafrica.file <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/South Africa 2014/SouthAfrica_p_2014June10-H-c.DP"
botswana.file <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Botswana 2014/Botswana 2014_Nat 19_06_14-c.DP"
malawi.file <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Malawi 2014/Malawi_2014_25June2014-c.DP"
swaziland.file <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Swaziland 2014/Swaziland Projections 31 Mar 14-c.DP"

southafrica.upd <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/unpop/South Africa_710.upd"
botswana.upd <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/unpop/Botswana_72.upd"
malawi.upd <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/unpop/Malawi_454.upd"
swaziland.upd <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/unpop/Swaziland_748.upd"

southafrica.demp <- read.demog.param(southafrica.upd)
southafrica.projp <- read.hivproj.param.5.0(southafrica.file)
southafrica.proj <- read.hivproj.output.5.0(southafrica.file)

botswana.demp <- read.demog.param(botswana.upd)
botswana.projp <- read.hivproj.param.5.0(botswana.file)
botswana.proj <- read.hivproj.output.5.0(botswana.file)

malawi.demp <- read.demog.param(malawi.upd)
malawi.projp <- read.hivproj.param.5.0(malawi.file)
malawi.proj <- read.hivproj.output.5.0(malawi.file)

swaziland.demp <- read.demog.param(swaziland.upd)
swaziland.projp <- read.hivproj.param.5.0(swaziland.file)
swaziland.proj <- read.hivproj.output.5.0(swaziland.file)


southafrica.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Spectrum-2014-outputs/southafrica"
botswana.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Spectrum-2014-outputs/botswana"
malawi.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Spectrum-2014-outputs/malawi"
swaziland.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Spectrum-2014-outputs/swaziland"


write.outputs(southafrica.path, southafrica.demp, southafrica.projp, southafrica.proj)
write.outputs(botswana.path, botswana.demp, botswana.projp, botswana.proj)
write.outputs(malawi.path, malawi.demp, malawi.projp, malawi.proj)
write.outputs(swaziland.path, swaziland.demp, swaziland.projp, swaziland.proj)


## Summaries for 15+ population
out.yrs <- as.character(1975:2013)

southafrica.15plus <- data.frame(newinf=colSums(southafrica.proj$newinf.m[4:17,out.yrs]) + colSums(southafrica.proj$newinf.f[4:17,out.yrs]),
           hivnum=colSums(southafrica.proj$hivnum.m[4:17,out.yrs]) + colSums(southafrica.proj$hivnum.f[4:17,out.yrs]),
           aidsdeaths=colSums(southafrica.proj$aidsdeaths.m[4:17,out.yrs]) + colSums(southafrica.proj$aidsdeaths.f[4:17,out.yrs]))

botswana.15plus <- data.frame(newinf=colSums(botswana.proj$newinf.m[4:17,out.yrs]) + colSums(botswana.proj$newinf.f[4:17,out.yrs]),
           hivnum=colSums(botswana.proj$hivnum.m[4:17,out.yrs]) + colSums(botswana.proj$hivnum.f[4:17,out.yrs]),
           aidsdeaths=colSums(botswana.proj$aidsdeaths.m[4:17,out.yrs]) + colSums(botswana.proj$aidsdeaths.f[4:17,out.yrs]))

malawi.15plus <- data.frame(newinf=colSums(malawi.proj$newinf.m[4:17,out.yrs]) + colSums(malawi.proj$newinf.f[4:17,out.yrs]),
           hivnum=colSums(malawi.proj$hivnum.m[4:17,out.yrs]) + colSums(malawi.proj$hivnum.f[4:17,out.yrs]),
           aidsdeaths=colSums(malawi.proj$aidsdeaths.m[4:17,out.yrs]) + colSums(malawi.proj$aidsdeaths.f[4:17,out.yrs]))

swaziland.15plus <- data.frame(newinf=colSums(swaziland.proj$newinf.m[4:17,out.yrs]) + colSums(swaziland.proj$newinf.f[4:17,out.yrs]),
           hivnum=colSums(swaziland.proj$hivnum.m[4:17,out.yrs]) + colSums(swaziland.proj$hivnum.f[4:17,out.yrs]),
           aidsdeaths=colSums(swaziland.proj$aidsdeaths.m[4:17,out.yrs]) + colSums(swaziland.proj$aidsdeaths.f[4:17,out.yrs]))

write.csv(southafrica.15plus, paste(southafrica.path, "-15plus-aggregate.csv", sep=""))
write.csv(botswana.15plus, paste(botswana.path, "-15plus-aggregate.csv", sep=""))
write.csv(malawi.15plus, paste(malawi.path, "-15plus-aggregate.csv", sep=""))
write.csv(swaziland.15plus, paste(swaziland.path, "-15plus-aggregate.csv", sep=""))
  


####  CD4 progression and mortality  ####


southafrica.projp$cd4.prog[1,,] == southafrica.projp$cd4.prog[2,,]
southafrica.projp$cd4.prog[1,,] == botswana.projp$cd4.prog[1,,]
southafrica.projp$cd4.prog[1,,] == malawi.projp$cd4.prog[1,,]
southafrica.projp$cd4.prog[1,,] == swaziland.projp$cd4.prog[1,,]

southafrica.projp$cd4.initdist[1,,] == southafrica.projp$cd4.initdist[2,,]
southafrica.projp$cd4.initdist[1,,] == botswana.projp$cd4.initdist[1,,]
southafrica.projp$cd4.initdist[1,,] == malawi.projp$cd4.initdist[1,,]
southafrica.projp$cd4.initdist[1,,] == swaziland.projp$cd4.initdist[1,,]

southafrica.projp$cd4.art.mort[1,,,1] == southafrica.projp$cd4.art.mort[2,,,1]
southafrica.projp$cd4.art.mort[1,,,1] == botswana.projp$cd4.art.mort[1,,,1]
southafrica.projp$cd4.art.mort[1,,,1] == malawi.projp$cd4.art.mort[1,,,1]
southafrica.projp$cd4.art.mort[1,,,1] == swaziland.projp$cd4.art.mort[1,,,1]

southafrica.projp$cd4.art.mort[1,,,2] == southafrica.projp$cd4.art.mort[2,,,2]
southafrica.projp$cd4.art.mort[1,,,2] == botswana.projp$cd4.art.mort[1,,,2]
southafrica.projp$cd4.art.mort[1,,,2] == malawi.projp$cd4.art.mort[1,,,2]
botswana.projp$cd4.art.mort[1,,,2] == swaziland.projp$cd4.art.mort[1,,,2]

mean(botswana.projp$cd4.art.mort == swaziland.projp$cd4.art.mort)
mean(southafrica.projp$cd4.art.mort == swaziland.projp$cd4.art.mort)
mean(malawi.projp$cd4.art.mort == swaziland.projp$cd4.art.mort)



agegr.lab <- c(paste(0:15*5, 1:16*5, sep="-"), "80+")
cd4.lab <- c(">500", "350-500", "250-350", "200-250", "100-200", "50-100", "<50")

cd4.prog <- southafrica.projp$cd4.prog[1,,]
cd4.initdist <- southafrica.projp$cd4.initdist[1,,]
hivmort.noart <- southafrica.projp$cd4.art.mort[1,,,1]

dimnames(cd4.prog) <- list(agegr.lab, cd4.lab[-7])
dimnames(cd4.initdist) <- list(agegr.lab, cd4.lab)
dimnames(hivmort.noart) <- list(agegr.lab, cd4.lab)

artmort.less6mo.southern.m <- southafrica.projp$cd4.art.mort[1,,,2]
artmort.less6mo.southern.f <- southafrica.projp$cd4.art.mort[2,,,2]
artmort.6to12mo.southern.m <- southafrica.projp$cd4.art.mort[1,,,3]
artmort.6to12mo.southern.f <- southafrica.projp$cd4.art.mort[2,,,3]
artmort.greater12mo.southern.m <- southafrica.projp$cd4.art.mort[1,,,4]
artmort.greater12mo.southern.f <- southafrica.projp$cd4.art.mort[2,,,4]

artmort.less6mo.eastern.m <- malawi.projp$cd4.art.mort[1,,,2]
artmort.less6mo.eastern.f <- malawi.projp$cd4.art.mort[2,,,2]
artmort.6to12mo.eastern.m <- malawi.projp$cd4.art.mort[1,,,3]
artmort.6to12mo.eastern.f <- malawi.projp$cd4.art.mort[2,,,3]
artmort.greater12mo.eastern.m <- malawi.projp$cd4.art.mort[1,,,4]
artmort.greater12mo.eastern.f <- malawi.projp$cd4.art.mort[2,,,4]

dimnames(artmort.less6mo.eastern.m) <- list(agegr.lab, cd4.lab)
dimnames(artmort.less6mo.eastern.f) <- list(agegr.lab, cd4.lab)
dimnames(artmort.6to12mo.eastern.m) <- list(agegr.lab, cd4.lab)
dimnames(artmort.6to12mo.eastern.f) <- list(agegr.lab, cd4.lab)
dimnames(artmort.greater12mo.eastern.m) <- list(agegr.lab, cd4.lab)
dimnames(artmort.greater12mo.eastern.f) <- list(agegr.lab, cd4.lab)

dimnames(artmort.less6mo.southern.m) <- list(agegr.lab, cd4.lab)
dimnames(artmort.less6mo.southern.f) <- list(agegr.lab, cd4.lab)
dimnames(artmort.6to12mo.southern.m) <- list(agegr.lab, cd4.lab)
dimnames(artmort.6to12mo.southern.f) <- list(agegr.lab, cd4.lab)
dimnames(artmort.greater12mo.southern.m) <- list(agegr.lab, cd4.lab)
dimnames(artmort.greater12mo.southern.f) <- list(agegr.lab, cd4.lab)



write.csv(cd4.prog, "progression-parameters/cd4-prog.csv")
write.csv(cd4.initdist, "progression-parameters/cd4-inditdist.csv")
write.csv(hivmort.noart, "progression-parameters/hivmort-noart.csv")

write.csv(artmort.less6mo.eastern.m, "progression-parameters/artmort-less6mo-eastern-m.csv")
write.csv(artmort.less6mo.eastern.f, "progression-parameters/artmort-less6mo-eastern-f.csv")
write.csv(artmort.6to12mo.eastern.m, "progression-parameters/artmort-6to12mo-eastern-m.csv")
write.csv(artmort.6to12mo.eastern.f, "progression-parameters/artmort-6to12mo-eastern-f.csv")
write.csv(artmort.greater12mo.eastern.m, "progression-parameters/artmort-greater12mo-eastern-m.csv")
write.csv(artmort.greater12mo.eastern.f, "progression-parameters/artmort-greater12mo-eastern-f.csv")

write.csv(artmort.less6mo.southern.m, "progression-parameters/artmort-less6mo-southern-m.csv")
write.csv(artmort.less6mo.southern.f, "progression-parameters/artmort-less6mo-southern-f.csv")
write.csv(artmort.6to12mo.southern.m, "progression-parameters/artmort-6to12mo-southern-m.csv")
write.csv(artmort.6to12mo.southern.f, "progression-parameters/artmort-6to12mo-southern-f.csv")
write.csv(artmort.greater12mo.southern.m, "progression-parameters/artmort-greater12mo-southern-m.csv")
write.csv(artmort.greater12mo.southern.f, "progression-parameters/artmort-greater12mo-southern-f.csv")
