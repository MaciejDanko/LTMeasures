require(Matrix)
require(splines)

compmat<-function(x, 
                  y, 
                  x.div=10, 
                  zero.class.add=T, 
                  zero.class.end=NULL,
                  zero.class.frac=0.2,
                  exposures=NULL, 
                  bs.use=(x.div*length(x))>=bs.df.max, 
                  bs.method=c('MortalitySmooth','bs'), 
                  bs.df=c('max','thumb',floor(length(x)/4)),
                  bs.df.max=150,
                  bs.deg=3){ 
  
  #Construct composition matrix object and automatically recalibrate age classes
  
  #x = age vector
  #y = vector with counts, e.g. DX. It must have the same length as x
  
  #x.div = number of small age classes within original one unit of age
  #zero.class.add = Adding open interval. 
  #zero.class.end = anticipated end of the open interval = age with expected zero counts.
  #     - if set to NULL and zero.class.add=TRUE then it will be anticipated automatically
  #zero.class.frac = the fraction of total range of x added as the last zero interval if zero.class.end = NULL
  #exposures = optionally for mortality rates. A vector of the same length as x and y
  #bs.use = use B-spline basis to speed-up computations
  #bs.method = Basis for B-spline:
  #     -'MortalitySmooth' = "P-splines" basis based on MortalitySmooth:::MortSmooth_bbase() (preferable), 
  #     -'bs' = basic B-splines basis based on splines:::bs()
  #bs.df = B-spline degree of freedom (number of inner knots) or a way to its calculation: 
  #     -'maxprec' = equal to the number of small age classes up to maximum defined in bs.df.max (preferable), 
  #     -'thumb' = rule of thumb: one knot for the B-spline basis each 4-5 observations, up to maximum defined in bs.df.max
  #     - or any value
  #bs.df.max = maximum number of knots (df) for B-spline basis.
  #bs.deg = degree of the piecewise polynomial
  
  require(Matrix)
  require(splines)
  
  .MortSmooth_Bbase<-function (x, xl=min(x), xr=max(x), df=floor(length(x)/4), deg=bs.deg) {
    #Modified version of MortalitySmooth:::MortSmooth_bbase
    ndx=df-deg
    dx <- (xr - xl)/ndx
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    P <- outer(x, knots, function (x, t)  (x - t)^deg * (x > t))
    D <- diff(diag(dim(P)[2]), diff = deg + 1)/(gamma(deg + 1) * dx^deg)
    B <- (-1)^(deg + 1) * P %*% t(D)
    B
  }
  
  warn.list=list()
  
  if (any(abs(as.integer(x)-x)>1e-6)) {
    WARN='Fractional values in age vector. Values were rounded'
    warning(WARN) 
    warn.list=c(warn.list,WARN)
  }
  x=as.integer(round(x))
  
  if ((length(x)<3+bs.deg) && (bs.use)) {
    WARN='Not enough data classes to use B-spline basis. Exact method was used'
    warning(WARN)
    warn.list=c(warn.list,WARN)
    bs.use=F
  }
  
  if (zero.class.add) {
    if (length(zero.class.end)==0){
      mstep=min(diff(x))
      r=diff(range(x))
      zero.class.end=ceiling(ceiling(r*zero.class.frac)/mstep)*mstep+max(x)
    }
    y=c(y,0)
    open.int.len=zero.class.end-x[length(x)]
    x=c(x,zero.class.end)
  } else open.int.len=NA
  
  segm=round(x*x.div)
  C2=approx(x=segm+1,y=1:length(segm),method='constant',xout=(min(segm):(max(segm)))+1,rule=2) 
  SparMat=Matrix:::sparseMatrix(C2$y,C2$x)
  C.=1*Matrix(SparMat,sparse=F)
  if (bs.use) {
    if (length(exposures)>0) stop('exposures cannot be used together with B-Splines')
    bs.method=bs.method[1]
    bs.df=bs.df[1]
    if (bs.df=='max') bs.df=min(bs.df.max,length(x)*x.div-2) 
    else if (tolower(bs.df)=='thumb') bs.df=min(bs.df.max,(length(x)*x.div-2)/4) 
    else if(is.na(as.numeric(bs.df))) stop(paste('Wrong bs.df parameter:',bs.df))
    bs.df=max(length(segm)-2,min(bs.df,length(segm)*x.div-2))
    if (bs.method=='bs'){
      X=splines:::bs(C2$x,df=bs.df) 
    } else if (tolower(bs.method)=='mortalitysmooth'){
      X=.MortSmooth_Bbase(C2$x,df=bs.df) 
    } else stop('Unknown B-spline basis method')
  } else{ 
    if (length(exposures)>0){
      expo=rep(exposures,dim(C.)[1])
      dim(expo)=rev(dim(C.))
      expo=t(expo)
      C.=C.*expo 
    } else expo=NULL
    bs.df=NA
    X=diag(dim(C.)[2])
  }
  list(C=as.matrix(C.),X=X,x=x,y=y,bs.df=bs.df,step=1/x.div,bs.use=bs.use,open.int.len=open.int.len,warn.list=warn.list)
}

#The extended version of PCLM proceure described in
#Rizzi S, Gampe J, Eilers PHC. Efficient estimation of smooth distributions
#from coarsely grouped data. Am J Epidemiol. 2015;182:138?47.
pclm <- function(CompositionMatrix,  lambda=1, deg = 2, max.iter=100,
                 lsfit.tol=.Machine$double.eps^0.5,pclm.tol=.Machine$double.eps^0.5){
  
  # Fit a PCLM (estimate b in ) E(y) = C %*% exp(X %*% b)
  # CompositionMatrix$y = the vector of observed counts of length i
  # CompositionMatrix$C = the composition matrix of dimension IxJ
  # CompositionMatrix$X = the identity matrix of dimension JxJ; or B-spline basis
  # lambda = smoothing parameter
  # deg = order of differences of the components of b
  # show = indicates whether iteration details should be shown
  # Fit the penalized composite link model
  
  warn.list=list()
  
  C= CompositionMatrix$C
  X= CompositionMatrix$X
  y= CompositionMatrix$y
  y=as.matrix(as.vector(y))
  nx <- dim(X)[2] #number of small classes
  D <- base:::diff(diag(nx), diff=deg) 
  bstart <- log(sum(y) / nx);
  b <- rep(bstart, nx);
  was.break=F
  for (it in 1:max.iter) {
    b0 <- b
    eta <- X %*% b
    gam <- exp(eta)
    mu <- C %*% gam
    w <- c(1 / mu, rep(lambda, nx - deg))
    Gam <- gam %*% rep(1, nx)
    Q <- C %*% (Gam * X)
    z <- c(y - mu + Q %*% b, rep(0, nx - deg))
    ls.x=rbind(Q, D)
    Fit <- lsfit(ls.x, z, wt = w, intercept = F,tolerance=lsfit.tol)
    b <- Fit$coef
    db <- max(abs(b - b0))
    if (db < pclm.tol) {
      was.break=T
      break
    }
  }
  if (!was.break) {
    WARN='Maximum iteration reached without convergence. Try to increase max.iter parameter.'
    warning(WARN)
    warn.list=c(warn.list,WARN)
  }
  R <- t(Q) %*% diag(c(1 / mu)) %*% Q
  H <- solve(R + lambda * t(D) %*% D) %*% R
  fit <- list()
  fit$trace <- sum(diag(H))
  ok <- y > 0 & mu > 0
  fit$dev <- 2 * sum(y[ok] * log(y[ok] / mu[ok]))
  fit$gamma <- gam
  fit$aic <- fit$dev + 2 * fit$trace
  #fit$aicc <- AICc(fit$aic,n=length(y),k=fit$trace)
  fit$bic <- fit$dev + fit$trace*log(length(y))
  fit$mu <- mu
  fit$warn.list <- warn.list
  fit
}

#optimize the smooth parameter lambda for pclm method
opt.pclm<-function(CompositionMatrix, opt.method=c('BIC','AIC'), opt.tol=.Machine$double.eps^0.5,...){
  
  #INPUT
  #CompositionMatrix = Composition Matrix calculated by compmat()
  #opt.method = one of the maximization criteria for lambda optimization: AIC and BIC (preferable)
  #opt.tol = Optimization tolerance
  #... options passed to the pclm method
  #OUTPUT
  # $X = ungrouped age classes
  # $Y = ungrouped counts
  # $lambda = optimal smooth parameter
  # $fit = fit object
  # $CompositionMatrix = CompositionMatrix object
  
  opt.method=opt.method[1]
  tryme<-function(G,Altern=1e200) suppressWarnings(if (class(try(G,silent=T))=="try-error") Altern else try(G,silent=T))
  if (toupper(opt.method)=='AIC') opty<-function (log10lam) tryme(pclm(CompositionMatrix, lambda = 10^log10lam,...)$aic) else
    if (toupper(opt.method)=='BIC') opty<-function (log10lam) tryme(pclm(CompositionMatrix, lambda = 10^log10lam,...)$bic) else 
      stop('Unknown method of lambda optimization')
  res.opt=optimize(f=opty,interval=c(-12,22),tol = opt.tol)$minimum
  lambda = 10^res.opt
  fit=pclm(CompositionMatrix, lambda = lambda, ...)
  X=seq(min(CompositionMatrix$x),max(CompositionMatrix$x),CompositionMatrix$step)
  Y=fit$gamma
  list(X=X,Y=Y,lambda=lambda,fit=fit,CompositionMatrix=CompositionMatrix)
}

#Calculate life table from the fit object
get.smooth.life.table<-function(fit, out.step=NULL, count.type=c('DX','LX')){
  #INPUT
  #fit = pclm fit object
  #out.step = age interval length in output life table, if out.step = NULL then the original data vector x is used
  #OUTPUT
  #grouped smoothed life table and raw smoothed life table tables 
  dig.prec=floor(-log10(.Machine$double.eps^0.5))
  count.type=count.type[1]
  p<-function(x) round(x,dig.prec) #to avoid numerical mistakes
  warn.list=fit$warn.list
  Y=fit$Y
  X=fit$X
  n=c(diff(X),diff(X)[length(X)-1])
  if (length(out.step)==0) x=p(fit$CompositionMatrix$x) else {
    if (out.step<diff(X)[1]) {
      WARN='Output age interval length (out.step) too small. Now equals the smallest age class'
      warning(WARN)
      warn.list=c(warn.list,WARN)
      out.step=diff(X)[1]
    }
    tmp=round(out.step/diff(X)[1])*diff(X)[1]
    if (p(tmp)!=p(out.step)) {
      WARN='Output age interval length (out.step) rounded to an integer multiple of the smallest age class'
      warning(WARN)
      warn.list=c(warn.list,WARN)
    }
    out.step=tmp
    x=p(unique(c(seq(X[1],X[length(X)],out.step),X[length(X)])))
  }  
  
  if (toupper(count.type)=='DX'){
    ax=rep(NA,length(x)-1)
    Dx=rep(NA,length(x)-1)
    for (j in 1:(length(x)-1)){
      ind=(x[j]<=X)&(x[j+1]>X)
      #print(X[ind])
      sDx=sum(Y[ind])
      mX=sum(Y[ind]*(X[ind]-X[ind][1]+n[j]/2))
      ax[j]=mX/sDx
      Dx[j]=sDx
    }
    grouped=data.frame(x=x[-length(x)],lx=sum(Dx)-c(0,cumsum(Dx)[-length(Dx)]),
                       dx=Dx,ax=ax,n=diff(x),Ax=ax/diff(x))
    raw=data.frame(x=X,lx=sum(Y)-c(0,cumsum(Y)[-length(Y)]),dx=Y,
                   ax=c(0.5*diff(X),0.5*diff(X)[length(X)-1]),
                   n=n,Ax=0.5)
  } else if (toupper(count.type)=='LX'){
    Lx=rep(NA,length(x)-1)
    for (j in 1:(length(x)-1)){
      ind=(x[j]<=X)&(x[j+1]>X)
      sLx=sum(Y[ind])
      Lx[j]=sLx
    }
    grouped=data.frame(x=x[-length(x)],Lx=Lx,n=diff(x))
    raw=data.frame(x=X,Lx=Y,n=n)
  } else stop('Unknown life-table type')
  object=list(grouped=grouped,
              raw=raw,
              fit=fit,
              warn.list=warn.list)
  class(object)='subPCLM'
  object
}

#fractional part 
frac<-function(x,digits=floor(-log10(.Machine$double.eps^0.5))) round(x-floor(round(x,digits)),digits)

#Main procedure to calculate PCLM with automated step
#More detailed description of prameters can be found in compmat and pclm function
PCLM<-function(x,y,      #DATA
               x.div=10, #Increase for higher precision
               x.max.ext=20, #internal parameter, maximal multiple of age interval
               x.auto.trans=T, #Automatically multiple age intervals to remove fractions?
               count.type=c('DX','LX'),
               zero.class.add=T, #T: It typically should be added
               zero.class.end=NULL,  #put here manually the end of open interval if necessary
               zero.class.frac=0.2,  #internal parameter
               exposures=NULL, #this method is equivalent to unsmoothing LX and DX separately to calcualte mu
               bs.use=(x.div*length(x))>=bs.df.max, #internal parameter, F force extensive but precise computations
               bs.method=c('MortalitySmooth','bs'), #internal parameter, use first method only
               bs.df=c('max','thumb',floor(length(x)/4)), #internal parameter, first method preferable
               bs.df.max=100, #internal parameter, but decrease it if computations are too slow
               bs.deg=3, #internal parameter, preferable value is 3
               out.step=min(diff(x),1), #age interval length in output file, automatically corrected
               opt.method=c('BIC','AIC'), #internal, BIC is more preferable - more smooth output
               opt.tol=.Machine$double.eps^0.5, #internal parameter
               pclm.deg = 2, #internal parameter, higher values may not work
               pclm.max.iter=100, #internal parameter, but increase on warning
               pclm.lsfit.tol=.Machine$double.eps^0.5, #internal parameter
               pclm.tol=.Machine$double.eps^0.5){ #internal parameter
  
  
  WARN=list() 
  if (all(order(x)!=1:length(x))) stop ('Age classes are not ordered')
  if (x.auto.trans) {
    
    extend.age<-function(x) {
      for (j in 1:(x.max.ext+1)){
        y=frac(j*x)
        #print(y)
        if (all(y==0)) break
      }
      j
    }
    
    m=extend.age(x)
    if (m>x.max.ext){
      WARN='Too small age interval found. Age classes were rounded.'
      x=round(x*x.max.ext)/x.max.ext
      m=extend.age(x)
      warning(WARN)
    } else {
      if (m>1) WARN='Age vector automaticaly transformed.'
    }
    x=as.integer(round(x*m)) #transform x
  } else m=1
  
  CM=compmat(x,y,x.div=x.div,
             zero.class.add=zero.class.add, 
             zero.class.end=zero.class.end,
             zero.class.frac=zero.class.frac,
             exposures=exposures, 
             bs.use=bs.use, 
             bs.method=bs.method[1], 
             bs.df=bs.df[1],
             bs.df.max=bs.df.max,
             bs.deg=bs.deg) 
  
  fit=opt.pclm(CM,
               opt.method=opt.method[1], 
               opt.tol=opt.tol,
               deg=pclm.deg,
               max.iter=pclm.max.iter,
               lsfit.tol=pclm.lsfit.tol,
               pclm.tol=pclm.tol)

  fit$CompositionMatrix$x=fit$CompositionMatrix$x/m
  fit$X=fit$X/m
  fit$CompositionMatrix$open.int.len=fit$CompositionMatrix$open.int.len/m
  
  GLT=get.smooth.life.table(fit,out.step=out.step,count.type=count.type[1])
  warn.list=c(CM$warn.list,GLT$warn.list,WARN)
  GLT$m=m
  GLT$x.div=x.div
  GLT$warn.list=warn.list
  class(GLT)='PCLM'
  GLT
}

#Plotting PCLM object
PCLM.plot<-function(object,ylab='Counts / interval length',xlab='Age'){
 if (class(object)!='PCLM') stop ('Object of class PCLM needed')

 n1=diff(object$fit$CompositionMatrix$x)
 n1=c(n1,n1[length(n1)])
 z0=n1[1]/2
 n2=diff(object$fit$X)
 n2=c(n2,n2[length(n2)])

 xlim=range(c(object$fit$X,object$fit$CompositionMatrix$x))
 ylim=range(c(object$fit$Y/n2,object$fit$CompositionMatrix$y/n1))
 ylim[1]=0
 ylim[2]=ylim[2]+diff(ylim)*0.2
 
 tmp.lwd=par('lwd'); par(lwd=2,xaxs='i',yaxs='i')
 barplot(width=n1,space=0,height=object$fit$CompositionMatrix$y/n1,xlab=xlab,ylab=ylab,
         col='gold2',border='white',xlim=xlim,ylim=ylim)
 par(lwd=tmp.lwd)
 
 lines(object$fit$CompositionMatrix$x,object$fit$CompositionMatrix$y/n1,type='s')
 axis(1)
 
 lines(object$fit$X, object$fit$Y/n2,type='s',col=2,lwd=2)
 
 if (!is.na(object$fit$CompositionMatrix$open.int.len)) ind=-(length(n1):(length(n1)-1)) else ind=1:length(n1)
 legend('topleft',c(paste('PCLM total classes =',length(n2)),
                    paste('Number of smoothing parameters =',object$fit$CompositionMatrix$bs.df),
                    paste('PCLM per class divider=',object$x.div),
                    paste('PCLM int. length =',round(min(n2),3)),
                    paste('Original min. int. length =',round(min(n1[ind]),3)),
                    paste('Original max. int. length =',round(max(n1[ind]),3)),
                    paste('Open int. length =',round(object$fit$CompositionMatrix$open.int.len,3)),
                    paste('PCLM classes per min. int. =',round(min(n1[ind])/min(n2),3)),
                    paste('PCLM classes per max. int. =',round(max(n1[ind])/max(n2),3)),
                    paste('Age class multiple =',object$m)),bty='n')
 box()
}
