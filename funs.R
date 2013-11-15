# models for fit
mod0 <- function() ~ 1
amod <- function() ~ age
famod <- function() ~ factor(age)
fymod <- function() ~ factor(year) 
faymod <- function() ~ factor(age) + factor(year)
bamod <- function(na) as.formula(paste("~bs(age,",na, ")"))
bymod <- function(ny) as.formula(paste("~bs(year,",ny, ")"))
baymod <- function(na, ny) as.formula(paste("~",as.character(bamod(na))[2],"+", as.character(bymod(ny))[2]))	
temod <- function(na, ny) as.formula(paste("~ te(age, year, bs=c(\"tp\", \"tp\"), k=c(", na,",",ny,"))"))

selDn <- function (params, data){
    a1 = FLQuant(1, dimnames = dimnames(data)) %*% params["a1"]
    s = FLQuant(1, dimnames = dimnames(data)) %*% params["sl"]
    sr = FLQuant(1, dimnames = dimnames(data)) %*% params["sr"]
    if (dims(data)$iter == 1 & dims(a1)$iter > 1) 
        data = propagate(data, dims(a1)$iter)
    s[data >= a1] = sr[data >= a1]
    2^(-((data - a1)/s * (data - a1)/s))
}

meanssq <- function(x){
	mean(x^2)
}

sim <- function(x){
	# BRP
	xx <- lh(gislasim(FLPar(linf=x$linf, k=x$k, s=x$s, v=x$v, a1=x$a1*x$a50, sr=x$sr, sl=x$sl, a50=x$a50, a=x$a, b=x$b)), sr=as.character(x$srMod), range=c(min=1, max=ceiling(x$amax*1.1), minfbar=1, maxfbar=ceiling(x$amax*1.1)-1, plusgroup=ceiling(x$amax*1.1)))
	# exploitation history
	Fc <- c(refpts(xx)["crash","harvest"]*0.8)
	Fmsy <- c(refpts(xx)["msy","harvest"])
	if(is.na(Fc)){
		Fc <- c(refpts(xx)["f0.1","harvest"]*5)
		Fmsy <- c(refpts(xx)["f0.1","harvest"])
	}
	Ftrg <- c(seq(0, Fc, len=19), rep(Fc, 20), seq(Fc, Fmsy, len=10))
	trg <- fwdControl(data.frame(year=c(2:50), quantity=rep('f', 49), val=Ftrg))
	stk <- as(xx, "FLStock")[,1:50]
	stk <- fwd(stk, ctrl=trg, sr=list(model=SRModelName(xx@model), params=xx@params))
	stk0 <- stk
	# Observation error
	ages <- range(stk)["min"]:min(10, range(stk)["max"])
	yrs <- range(stk)["minyear"]:range(stk)["maxyear"]
	range(stk)[c("maxfbar")] <- min(5, range(stk)[c("maxfbar")])
	stk <- setPlusGroup(stk, max(ages))
	N <- stock.n(stk)
	sel <- selDn(FLPar(a1=x$a1*(x$a50+0.5), sr=x$sr, sl=x$sl), data=FLQuant(ages+0.5,dimnames=list(age=ages)))
	logq <- FLQuant(log(x$oe.iq), dimnames=list(year=yrs))
	logq[] <- cumsum(c(logq))
	logq <- log(exp(logq[rep(1,length(ages))])*sel[,rep(1,length(yrs))]*0.01) # fixed absolute q value
	index <- FLIndex(index = N * exp(logq + rnorm(prod(dim(N)), 0, x$oe.icv)), index.q=exp(logq))
	dimnames(index.q(index)) <- dimnames(index(index))
	index@range[c("startf", "endf")] <- c(0.01, 1)		
	catch.n(stk) <- catch.n(stk) * exp(rnorm(prod(dim(N)), 0, x$oe.ccv))
	catch(stk) <- computeCatch(stk)
	lst <- list(scn=x, brp=xx, stock=stk0, 
				oe1=list(stk = window(stk, 5, 20), idx = window(index, 5, 20)),
				oe2=list(stk = window(stk, 15, 30), idx = window(index, 15, 30)),
				oe3=list(stk = window(stk, 25, 40), idx = window(index, 25, 40)),
				oe4=list(stk = window(stk, 35, 50), idx = window(index, 35, 50)),
				oe5=list(stk = stk, idx = index)
	)
	lst
}

getMods <- function(scn, na, ny){

	fmodel <- scn["fmodel"]
	if(fmodel=="fm1"){
		fmodel <- faymod()
	} else if(fmodel=="fm2"){
		fmodel <- baymod(na=na, ny=ny)
	} else {
		fmodel <- temod(na=na, ny=ny)
	}

	qmodel <- scn["qmodel"]
	if(qmodel=="qm0"){
		qmodel <- list(mod0())
	} else if(qmodel=="qm1"){
		qmodel <- list(amod())
	} else if(qmodel=="qm2"){
		qmodel <- list(famod())
	} else if(qmodel=="qm3"){
		qmodel <- list(bamod(na=na))
	} else {
		qmodel <- list(baymod(na=na, ny=ny))
	}

	rmodel <- scn["rmodel"]
	if(rmodel=="rm1"){
		rmodel <- fymod()
	} else {
		rmodel <- bymod(ny=ny)
	}

	list(fmodel=fmodel, qmodel=qmodel, rmodel=rmodel)
} 

doFits0 <- function(x){
		
	# build the models considering the age and year lengths
	ny <- 10	
	na <- min(range(x$oe5$stk)["max"], 4)

	# models
	mods <- getMods(x$scn, na=na, ny=ny)
	fmodel <- mods$fmodel
	qmodel <- mods$qmodel
	rmodel <- mods$rmodel

	fit <- try(a4aFit(fmodel=fmodel, qmodel=qmodel, rmodel=rmodel, stock=x$oe1$stk, indices=FLIndices(idx=x$oe1$idx)))
	if(class(fit)=="try-error") x$oe1$fit <- NA else x$oe1$fit <- fit
	rm(fit)  

	fit <- try(a4aFit(fmodel=fmodel, qmodel=qmodel, rmodel=rmodel, stock=x$oe2$stk, indices=FLIndices(idx=x$oe2$idx)))
	if(class(fit)=="try-error") x$oe2$fit <- NA else x$oe2$fit <- fit  
	rm(fit)  

	fit <- try(a4aFit(fmodel=fmodel, qmodel=qmodel, rmodel=rmodel, stock=x$oe3$stk, indices=FLIndices(idx=x$oe3$idx)))
	if(class(fit)=="try-error") x$oe3$fit <- NA else x$oe3$fit <- fit  
	rm(fit)  

	fit <- try(a4aFit(fmodel=fmodel, qmodel=qmodel, rmodel=rmodel, stock=x$oe4$stk, indices=FLIndices(idx=x$oe4$idx)))
	if(class(fit)=="try-error") x$oe4$fit <- NA else x$oe4$fit <- fit  
	rm(fit)  

	# for this last fit ny needs to change and models rebuild
	ny <- 30
	# models
	mods <- getMods(x$scn, na=na, ny=ny)
	fmodel <- mods$fmodel
	qmodel <- mods$qmodel
	rmodel <- mods$rmodel
	
	fit <- try(a4aFit(fmodel=fmodel, qmodel=qmodel, rmodel=rmodel, stock=window(x$oe5$stk, 5), indices=FLIndices(idx=window(x$oe5$idx, 5))))
	if(class(fit)=="try-error") x$oe5$fit <- NA else x$oe5$fit <- fit  
	rm(fit)  

	x
}

getStats <- function(fit, stk, ssb0, fbar0, rec0, cat0, cat0.wt, q0){
	yrs <- dimnames(stock.n(fit))[[2]]
	stk <- stk[,yrs]+fit
	fitssb <- ssb(stk)
	fitfbar <- fbar(stk)

	# mse
	q1 <- meanssq(fit@logq[[1]]-log(q0[,yrs]))
	ssb1 <- meanssq(fitssb-ssb0[,yrs])
	fbar1 <- meanssq(fitfbar-fbar0[,yrs])
	rec1 <- meanssq(fit@stock.n[1]-rec0[,yrs])
	cat1 <- meanssq(quantSums(exp(fit@catch.lhat) * cat0.wt[,yrs]) - cat0[,yrs])
	# rbias
	q2 <- mean(abs(fit@logq[[1]]/log(q0[,yrs])-1))
	ssb2 <- mean(abs(fitssb/ssb0[,yrs]-1))
	fbar2 <- mean(abs(fitfbar/fbar0[,yrs]-1))
	rec2 <- mean(abs(fit@stock.n[1]/rec0[,yrs]-1))
	cat2 <- mean(abs(quantSums(exp(fit@catch.lhat) * cat0.wt[,yrs]) / cat0[,yrs]-1))

	c(fit@fit.sum, catlvar=mean(fit@catch.lvar), idxvar=mean(fit@index.lvar[[1]]), ssbrbias=ssb2, ssbmse=ssb1, fbarrbias=fbar2, fbarmse=fbar1, recrbias=rec2, recmse=rec1, catrbias=cat2, catmse=cat1, qrbias=q2, qmse=q1)
}

getFitStats <- function(x){
	cat(x$scn$id, ",")
	stk <- x$stock	
	#stk <- x$oe5$stk	
	ssb0 <- ssb(stk)
	# need this tweak because range in a4aFit is not correct
	v <- x$oe5$stk@range[c("minfbar", "maxfbar")]
	fbar0 <- quantMeans(harvest(stk)[as.character(v[1]:v[2])])
	rec0 <- rec(stk)
	cat0 <- catch(stk)
	df0 <- data.frame(x$scn, fmsy=c(x$brp@refpts["msy","harvest"]), f0.1=c(x$brp@refpts["f0.1","harvest"]), m=c(mean(x$stock@m)), expl=NA, nopar=NA, nlogl=NA, maxgrad=NA, npar=NA, logDetHess=NA, catlvar=NA, idxlvar=NA, ssbrbias=NA, ssbmse=NA, fbarrbias=NA, fbarmse=NA, recrbias=NA, recmse=NA, catrbias=NA, catmse=NA, qrbias=NA, qmse=NA)
	df0 <- df0[rep(1,5),]
	df0$expl <- 1:5
	vars <- c("nopar", "nlogl", "maxgrad", "npar", "logDetHess", "catlvar", "idxlvar", "ssbrbias", "ssbmse", "fbarrbias", "fbarmse", "recrbias", "recmse", "catrbias", "catmse", "qrbias","qmse")
	if(!is.na(x$oe1$fit)) df0[1,vars] <- getStats(x$oe1$fit, x$oe1$stk, ssb0, fbar0, rec0, cat0, catch.wt(x$oe1$stk), index.q(x$oe1$idx)) 
	if(!is.na(x$oe2$fit)) df0[2,vars] <- getStats(x$oe2$fit, x$oe2$stk, ssb0, fbar0, rec0, cat0, catch.wt(x$oe2$stk), index.q(x$oe2$idx)) 
	if(!is.na(x$oe3$fit)) df0[3,vars] <- getStats(x$oe3$fit, x$oe3$stk, ssb0, fbar0, rec0, cat0, catch.wt(x$oe3$stk), index.q(x$oe3$idx)) 
	if(!is.na(x$oe4$fit)) df0[4,vars] <- getStats(x$oe4$fit, x$oe4$stk, ssb0, fbar0, rec0, cat0, catch.wt(x$oe4$stk), index.q(x$oe4$idx)) 
	if(!is.na(x$oe5$fit)) df0[5,vars] <- getStats(x$oe5$fit, x$oe5$stk, ssb0, fbar0, rec0, cat0, catch.wt(x$oe5$stk), index.q(x$oe5$idx)) 
	df0
}

getSumm <- function(fit, stk, ssb0, fbar0, rec0, cat0, cat0.wt, I0, B0, C0, qFa, qFa0, qR0, expl){
	yrs <- dimnames(stock.n(fit))[[2]]
	stk <- stk[,yrs]+fit
	stat <- c("S", "F", "R", "Y", "I", "B", "C", "qFa", "qR") 
	src <- c("obs", "hat")
	df0 <- data.frame(expand.grid(expl=expl, y=yrs, stat=stat, src=src), val=NA)	
	ssb1 <- ssb(stk)
	fbar1 <- fbar(stk)
	rec1 <- fit@stock.n[1]
	cat1 <- quantSums(exp(fit@catch.lhat) * cat0.wt[,yrs])

	I1 <- quantSums(exp(fit@index.lhat[[1]]))
	B1 <- quantSums(exp(fit@index.lhat[[1]]) * cat0.wt[,yrs])
	C1 <- quantSums(fit@catch.n)
	qFa1 <- fit@logq[[1]][qFa] 
	qR1 <- fit@logq[[1]][1]

	df0[,"val"] <- c(ssb0[,yrs], fbar0[,yrs], rec0[,yrs], cat0[,yrs], I0[,yrs], B0[,yrs], C0[,yrs], qFa0[,yrs], qR0[,yrs], ssb1, fbar1, rec1, cat1, I1, B1, C1, qFa1, qR1)

	df0
}

getFitSumm <- function(x){
	stk <- x$stock	
	idxq <- x$oe5$idx@index.q	
	ssb0 <- ssb(stk)
	# need this tweak because range in a4aFit is not correct
	a <- dimnames(idxq)[[1]]
	y <- dimnames(idxq)[[2]]
	v <- x$oe5$stk@range[c("minfbar", "maxfbar")]
	fbar0 <- quantMeans(harvest(stk)[as.character(v[1]:v[2])])
	rec0 <- rec(stk)
	cat0 <- catch(stk)
	I0 <- quantSums(stock.n(stk)[a,y]*idxq)
	B0 <- quantSums(stock.n(stk)[a,y]*idxq*stock.wt(stk)[a,y])
	C0 <- quantSums(catch.n(stk))
	qFa <- as.character(floor(x$scn["a1"] * (x$scn["a50"]+0.5)))
	if(!(qFa %in% a)) qFa <- a[ceiling(length(a)/2)]
	qFa0 <- log(idxq[qFa]) 
	qR0 <- log(idxq[1])
	df0 <- data.frame(expl=NA, y=NA, stat=NA, src=NA, val=NA)	
	if(!is.na(x$oe1$fit)) df0 <- rbind(df0, getSumm(x$oe1$fit, x$oe1$stk, ssb0, fbar0, rec0, cat0, catch.wt(x$oe1$stk), I0, B0, C0, qFa, qFa0, qR0, "1")) 
	if(!is.na(x$oe2$fit)) df0 <- rbind(df0, getSumm(x$oe2$fit, x$oe2$stk, ssb0, fbar0, rec0, cat0, catch.wt(x$oe2$stk), I0, B0, C0, qFa, qFa0, qR0, "2")) 
	if(!is.na(x$oe3$fit)) df0 <- rbind(df0, getSumm(x$oe3$fit, x$oe3$stk, ssb0, fbar0, rec0, cat0, catch.wt(x$oe3$stk), I0, B0, C0, qFa, qFa0, qR0, "3")) 
	if(!is.na(x$oe4$fit)) df0 <- rbind(df0, getSumm(x$oe4$fit, x$oe4$stk, ssb0, fbar0, rec0, cat0, catch.wt(x$oe4$stk), I0, B0, C0, qFa, qFa0, qR0, "4")) 
	if(!is.na(x$oe5$fit)) df0 <- rbind(df0, getSumm(x$oe5$fit, x$oe5$stk, ssb0, fbar0, rec0, cat0, catch.wt(x$oe5$stk), I0, B0, C0, qFa, qFa0, qR0, "5")) 
	attr(df0, "scn") <- x$scn
	df0[-1,]
}



## HAD TO FIX A SMALL BUG ##
#### Life History Generator ####################################################
lh=function(par,
            growth       =vonB,
            fnM          =function(par,len) exp(par["M1"]+par["M2"]*log(len)),
#            fnM          =function(par,len,T=290,a=FLPar(c(a=-2.1104327,b=-1).7023068,c=1.5067827,d=0.9664798,e=763.5074169),iter=dims(par)$iter))
#                                    exp(a[1]+a[2]*log(len) + a[3]*log(par["linf"]) + a[4]*log(par["k"]) + a[5]/T),
            fnMat        =logistic,
            fnSel        =dnormal,
            sr           ="bevholt",
            range        =c(min=1,max=40,minfbar=1,maxfbar=40,plusgroup=40),
            spwn         = 0,
            fish = 0.5, # proportion of year when fishing happens
            units=if("units" %in% names(attributes(par))) attributes(par)$units else NULL,
            ...){

  # Check that m.spwn and harvest.spwn are 0 - 1
  if (spwn > 1 | spwn < 0 | fish > 1 | fish < 0)
    stop("spwn and fish must be in the range 0 to 1\n")
  
  if (("m.spwn" %in% names(args)))
     m.spwn =args[["m.spwn"]]
  else
    m.spwn=FLQuant(spwn, dimnames=list(age=range["min"]:range["max"]))

  if (("harvest.spwn" %in% names(args)))
    harvest.spwn =args[["harvest.spwn"]]
  else
    harvest.spwn=FLQuant(spwn, dimnames=list(age=range["min"]:range["max"]))

  age=propagate(FLQuant(range["min"]:range["max"],dimnames=list(age=range["min"]:range["max"])),length(dimnames(par)$iter))

   # Get the lengths through different times of the year
   stocklen   <- growth(par[c("linf","t0","k")],age+m.spwn)    # stocklen is length at spawning time
   catchlen   <- growth(par[c("linf","t0","k")],age+fish) # catchlen is length when fishing happens
   midyearlen <- growth(par[c("linf","t0","k")],age+0.5) # midyear length used for natural mortality

   # Corresponding weights
   swt=par["a"]*stocklen^par["b"]
   cwt=par["a"]*catchlen^par["b"]
   if ("bg" %in% dimnames(par)$param)  
      swt=par["a"]*stocklen^par["bg"]
  
   args<-list(...)

   m.   =fnM(  par=par,len=midyearlen) # natural mortality is always based on mid year length
   mat. =fnMat(par,age + m.spwn) # maturity is biological therefore + m.spwn
   sel. =fnSel(par,age + fish) # selectivty is fishery  based therefore + fish

   ## create a FLBRP object to   calculate expected equilibrium values and ref pts
   dms=dimnames(m.)
   res=FLBRP(stock.wt       =swt,
             landings.wt    =cwt,
             discards.wt    =cwt,
             bycatch.wt     =cwt,
             m              =m.,
             mat            =FLQuant(mat., dimnames=dimnames(m.)),
             landings.sel   =FLQuant(sel., dimnames=dimnames(m.)),
             discards.sel   =FLQuant(0,    dimnames=dimnames(m.)),
             bycatch.harvest=FLQuant(0,    dimnames=dimnames(m.)),
             harvest.spwn   =FLQuant(harvest.spwn,    dimnames=dimnames(m.)),
             m.spwn         =FLQuant(m.spwn,    dimnames=dimnames(m.)),
             availability   =FLQuant(1,    dimnames=dimnames(m.)),
             range          =range)

   ## FApex
   #if (!("range" %in% names(args))) range(res,c("minfbar","maxfbar"))[]<-as.numeric(dimnames(landings.sel(res)[landings.sel(res)==max(landings.sel(res))][1])$age)

   ## replace any slot passed in as an arg
   for (slt in names(args)[names(args) %in% names(getSlots("FLBRP"))[names(getSlots("FLBRP"))!="fbar"]])
     slot(res, slt)<-args[[slt]]
   params(res)=propagate(params(res),dims(res)$iter)
   ## Stock recruitment relationship
   model(res) =do.call(sr,list())$model

  if (dims(par)$iter>1) {
     warning("Scarab, iters dont work for SRR:sv/ab etc")
  
     params(res)=FLPar(c(a=NA,b=NA),iter=dims(par)$iter)
       
     for (i in seq(dims(par)$iter))
       params(res)[,i][]=unlist(c(ab(par[c("s","v"),i],sr,spr0=iter(spr0(res),i))[c("a","b")]))

      warning("iter(params(res),i)=ab(par[c(s,v),i],sr,spr0=iter(spr0(res),i))[c(a,b)] assignment doesnt work")
      warning("iter(FLBRP,i) doesnt work")
  }else
    params(res)=ab(par[c("s","v")],sr,spr0=spr0(res))[c("a","b")]

   dimnames(refpts(res))$refpt[5]="crash"

   res=brp(res)
   
   if ("fbar" %in% names(args)) 
      fbar(res)<-args[["fbar"]] else 
   if (any((!is.nan(refpts(res)["crash","harvest"])))) 
	  # the line below has a bug, the quant is not "age" 
	  # which later creates problems. Fixed on the next line 
      #fbar(res)<-FLQuant(seq(0,1,length.out=101))*refpts(res)["crash","harvest"]
      fbar(res)<-FLQuant(seq(0,1,length.out=101), quant="age")*refpts(res)["crash","harvest"]
  
   res=brp(res)

   if (!("units" %in% names(attributes(par))))  return(res)

   res <- setUnits(res, par)

  return(res)}


setMethod('sv', signature(x='FLPar', model='character'),
  function(x, model, spr0=NA){
 
   a=x["a"]
   b=x["b"]
   s=FLPar(a=1,dimnames=dimnames(a))  
   v=FLPar(b=1,dimnames=dimnames(a))  
   spr0=FLPar(spr0,dimnames=dimnames(a))  

   if ("spr0" %in% dimnames(x)$params)
     spr0=x["spr0"] 

   c=FLPar(c=1,dimnames=dimnames(a))  
   d=FLPar(d=1,dimnames=dimnames(a))  
   if (("c" %in% dimnames(x)$params))  c=x["c"]
   if (("d" %in% dimnames(x)$params))  d=x["d"]

   v <- v*spr2v(model, spr0, a, b, c, d)
   s <- s*srr2s(model, ssb=v*.2, a=a, b=b, c=c, d=d) / srr2s(model, ssb=v, a=a, b=b, c=c, d=d)
  
   res=rbind(s, v, spr0)
 
   if ("c" %in% dimnames(x)$params)
     res=rbind(res, c)
 
   if ("d" %in% dimnames(x)$params)
     res=rbind(res, d)
 
   res=rbind(res, spr0)
 
   return(res)})

abPars. <- function(x,spr0=NA,model){
  s=x["s"]
  v=x["v"]
  if ("c" %in% names(x))
     c=x["c"]
  if ("d" %in% names(x))
     d=x["d"]
  if ("spr0" %in% names(x))
     spr0=x["spr0"]
  # converts a & b parameterisation into steepness & virgin biomass (s & v)
  switch(model,
    "bevholt"   ={a=(v%+%(v%-%s%*%v)%/%(5%*%s%-%1))%/%spr0; b=(v%-%s%*%v)%/%(5%*%s%-%1)},
    "bevholtSV" ={a=(v+(v-s*v)/(5*s-1))/spr0; b=(v-s*v)/(5*s-1)},
    "ricker"    ={b=log(5*s)/(v*0.8); a=exp(v*b)/spr0},
    "rickerSV"  ={b=log(5*s)/(v*0.8); a=exp(v*b)/spr0},
    "cushing"   ={b=log(s)/log(0.2); a=(v^(1-b))/(spr0)},
    "cushingSV" ={b=log(s)/log(0.2); a=(v^(1-b))/(spr0)},
    "shepherd"  ={b=v*(((0.2-s)/(s*0.2^c-0.2))^-(1/c)); a=((v/b)^c+1)/spr0},
    "shepherdSV"={b=v*(((0.2-s)/(s*0.2^c-0.2))^-(1/c)); a=((v/b)^c+1)/spr0},
    "mean"      ={a=v/spr0;b=NULL},
    "meanSV"    ={a=v/spr0;b=NULL},
    "segreg"    ={a=5*s/spr0; b=v/(a*spr0)},
    "segregSV"  ={a=5*s/spr0; b=v/(a*spr0)},
    {stop("model name not recognized")})

  res <- c(a=a, b=b)
  return(res[!is.null(res)])} 


# setMethod('ab', signature(x='FLPar', model='character'),
#   function(x, model, spr0=NA){
#  
#    s=x["a"]
#    v=x["b"]
#    a=FLPar(a=1,dimnames=dimnames(s))  
#    b=FLPar(b=1,dimnames=dimnames(v)) 
#    
#    if ("spr0" %in% dimnames(x)$params)
#       spr0=x["spr0"]  else 
#       spr0=FLPar(spr0,dimnames=dimnames(a)) 
# 
#    c=FLPar(c=1,dimnames=dimnames(a))  
#    d=FLPar(d=1,dimnames=dimnames(a))  
#    if (("c" %in% dimnames(x)$params))  c=x["c"]
#    if (("d" %in% dimnames(x)$params))  d=x["d"]
# 
#    v <- v*spr2v(model, spr0, a, b, c, d)
#    s <- s*srr2s(model, ssb=v*.2, a=a, b=b, c=c, d=d) / srr2s(model, ssb=v, a=a, b=b, c=c, d=d)
#   
#    res=rbind(s, v, spr0)
#  
#    if ("c" %in% dimnames(x)$params)
#      res=rbind(res, c)
#  
#    if ("d" %in% dimnames(x)$params)
#      res=rbind(res, d)
#  
#    res=rbind(res, spr0)
#  
#    return(res)})

gislasim=function(par,t0=-0.1,a=0.00001,b=3,ato95=1,sl=2,sr=5000,s=0.9,v=1000){
  
  names(dimnames(par)) <- tolower(names(dimnames(par)))
  
  if (!("t0"    %in% dimnames(par)$params)) par=rbind(par,FLPar("t0"    =t0, iter=dims(par)$iter))
  if (!("a"     %in% dimnames(par)$params)) par=rbind(par,FLPar("a"     =a,  iter=dims(par)$iter))
  if (!("b"     %in% dimnames(par)$params)) par=rbind(par,FLPar("b"     =b,  iter=dims(par)$iter))
  if (!("asym"  %in% dimnames(par)$params)) par=rbind(par,FLPar("asym"  =1,  iter=dims(par)$iter))
  if (!("bg"    %in% dimnames(par)$params)) par=rbind(par,FLPar("bg"    =b,  iter=dims(par)$iter))
  if (!("sl"    %in% dimnames(par)$params)) par=rbind(par,FLPar("sl"    =sl, iter=dims(par)$iter))
  if (!("sr"    %in% dimnames(par)$params)) par=rbind(par,FLPar("sr"    =sr, iter=dims(par)$iter))
  if (!("s"     %in% dimnames(par)$params)) par=rbind(par,FLPar("s"     =s,  iter=dims(par)$iter))
  if (!("v"     %in% dimnames(par)$params)) par=rbind(par,FLPar("v"     =v,  iter=dims(par)$iter))

  ## growth parameters
  if (!("k"     %in% dimnames(par)$params)) par=rbind(par,FLPar("k"=3.15*par["linf"]^(-0.64), iter=dims(par)$iter)) # From Gislason et al 2008, all species combined
  
  # Natural mortality parameters from Model 2, Table 1 Gislason 2010
  par=rbind(par,FLPar(M1=0.55+1.44*log(par["linf"])+log(par["k"]), iter=dims(par)$iter),
                FLPar(M2=-1.61                                   , iter=dims(par)$iter))

  if (!("ato95" %in% dimnames(par)$params)) par=rbind(par,FLPar("ato95" =ato95, iter=dims(par)$iter))
  if (!("sl"    %in% dimnames(par)$params)) par=rbind(par,FLPar("sl"    =sl,    iter=dims(par)$iter))
  if (!("sr"    %in% dimnames(par)$params)) par=rbind(par,FLPar("sr"    =sr,    iter=dims(par)$iter))
 
  ## maturity parameters from http://www.fishbase.org/manual/FishbaseThe_MATURITY_Table.htm
  if (!("asym"    %in% dimnames(par)$params)) par=rbind(par,FLPar("asym"    =asym, iter=dims(par)$iter))

  if (!("a50" %in% dimnames(par)$params)){
    par=rbind(par,FLPar(a50=0.72*par["linf"]^0.93, iter=dims(par)$iter))
    par["a50"]=invVonB(par,c(par["a50"]))
    }

  ## selectivity guestimate
  a1=par["a50"]
 
  dimnames(a1)$params="a1"
 
  par=rbind(par,a1)
  
  attributes(par)$units=c("cm","kg","1000s")
  
  return(par)}

setUnits=function(res, par){

    units=attributes(par)$units
    #browser()
    allUnits=list("params"=      "",          
               "refpts"=         "",            
               "fbar"=           "",        
               "fbar.obs"=       "",    
               "landings.obs"=   paste(units[2],units[3]),
               "discards.obs"=   paste(units[2],units[3]),
               "rec.obs"=        units[3],         
               "ssb.obs"=        paste(units[2],units[3]),
               "stock.obs"=      paste(units[2],units[3]),
               "profit.obs"=     "",     
 #              "revenue.obs"=    "",    
               "landings.sel"=   "",    
               "discards.sel"=   "", 
               "bycatch.harvest"="",        
               "stock.wt"=       units[2],     
               "landings.wt"=    units[2],     
               "discards.wt"=    units[2],      
               "bycatch.wt"=     units[2],               
               "m"=              "",             
               "mat"=            "proportion", 
               "harvest.spwn"=   "proportion",          
               "m.spwn"=         "proportion",    
               "availability"=   "proportion",           
               "price"=          "",           
               "vcost"=          "",           
               "fcost"=          "")            

    
    units(res)[names(allUnits)]=allUnits
    
    return(res)}


