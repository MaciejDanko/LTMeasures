#Some basic functions
tryme<-function(G,Altern) if (class(try(G,silent=T))=="try-error") Altern else try(G,silent=T)
rm.inf<-function(Vec,minn=1e100) {
  ind=(abs(Vec)==(Inf))
  ind[is.na(ind)]=F
  Vec[ind]=sign(Vec[ind])*minn
  return(Vec)
}

rm.nan<-function(Vec,altern=0) {
  ind=is.nan(Vec)
  ind[is.na(ind)]=F
  Vec[ind]=altern
  return(Vec)
}

#####################################
#Constant mortality model
Inv.C.surv<-function(Sx,a) -log(Sx)/a
C.surv<-function(x,a) exp(-a*x)
C.ax<-function(x,a,step=c(diff(x),diff(x)[length(x)-1]))
  sapply(1:length(x),function(i){
    x.=x[i]
    Lx=integrate(f=function(y) C.surv(y,a),lower=x.,upper=x.+step[i], subdivisions=1e4)$value
    lx=C.surv(x.,a)
    lx1=C.surv(x.+step[i],a)
    ((Lx-step[i]*lx1)/(lx-lx1))/step[i]
  })

C.mort<-function(x,a) a
C.e0<-function(a) 1/a
C.sq.e0<-function(a) 1/(2*a)
C.eR<-function(R,a) exp(-a*R)/a *a
C.dager<-function(a) 1/a
C.gini<-function(a) 1-C.sq.e0(a)/C.e0(a)
C.mode<-function(a) 0
C.LAR<-function(x,a) rep(0,length(x))
C.max.x<-function(a) min(1e12,max(c(10,Inv.C.surv(1e-12*a,a))))
C.var<-function(a,max.x=C.max.x(a)) integrate(f=function (x) rm.nan(((x-1/a)^2)*exp(-a*x)*a),
                                              lower=0,upper=max.x,subdivisions=1e4)$value
C.CV<-function(a,max.x=max(c(1e3,Inv.C.surv(1e-12*a,a)))) sqrt(C.var(a,max.x))/C.e0(a)

#####################################
#Gompertz
Inv.G.surv<-function(Sx,a,b) if (b<1e-12) {
  return(Inv.C.surv(Sx,a))
} else return(log((1-log(Sx)*(b/a)))/b)
G.surv<-function(x,a,b) if (b!=0) exp(-a*(exp(b*x)-1)/b) else exp(-a*x)
G.mort<-function(x,a,b) a*exp(b*x)
G.max.x<-function(a,b) min(1e12,max(c(10,Inv.G.surv(1e-12*a,a,b))))
G.ax<-function(x,a,b,step=c(diff(x),diff(x)[length(x)-1]))
  sapply(1:length(x),function(i){
    x.=x[i]
    Lx=integrate(f=function(y) G.surv(y,a,b),lower=x.,upper=x.+step[i], subdivisions=1e4)$value
    lx=G.surv(x.,a,b)
    lx1=G.surv(x.+step[i],a,b)
    ((Lx-step[i]*lx1)/(lx-lx1))/step[i]
  })

G.e0<-function(a,b,max.x=G.max.x(a,b)) if (b<1e-12) {
  return(C.e0(a))
} else return(integrate(f=function (x) (G.surv(x,a,b)),lower=0,upper=max.x,subdivisions=1e4)$value)
G.sq.e0<-function(a,b,max.x=G.max.x(a,b)) if (b<1e-12) {
  return(C.sq.e0(a))
} else return(integrate(f=function (x) rm.nan(G.surv(x,a,b)^2),lower=0,upper=max.x,subdivisions=1e4)$value)
G.eR<-function(R,a,b,max.x=G.max.x(a,b)) if (b<1e-12) {
  return(C.eR(R,a))
} else if (G.surv(R,a,b)<=0) return(0) else return(integrate(f=function (x) (G.surv(x,a,b)/G.surv(R,a,b)),lower=R,upper=max.x,subdivisions=1e4)$value)

G.var<-function(a,b,max.x=G.max.x(a,b)) {
  e0=G.e0(a,b,max.x)
  if (b==0) {
    return(C.var(a,b)) 
  } else return(integrate(f=function (x) rm.nan(((x-e0)^2)*G.surv(x,a,b)*G.mort(x,a,b)),
                          lower=0,upper=max.x,subdivisions=1e4)$value)
}

G.CV<-function(a,b,max.x=G.max.x(a,b)) sqrt(G.var(a,b,max.x))/G.e0(a,b,max.x)

G.dager<-function(a,b,max.x=G.max.x(a,b)) if (b<1e-12) {
  return(C.dager(a))
} else return(integrate(f=function (x) rm.nan(-G.surv(x,a,b)*log(G.surv(x,a,b))),lower=0,upper=max.x,subdivisions=1e4)$value)

G.gini<-function(a,b,max.x=G.max.x(a,b)) 1-G.sq.e0(a,b,max.x)/G.e0(a,b,max.x)

G.mode<-function(a,b) if (b<1e-12) {
  return(C.mode(a))
} else return((-log(a/b)/b)*((-log(a/b)/b)>0))
G.LAR<-function(x,a,b) b*rep(1,length(x))

#####################################
#Gompertz and Gompertz-Makeham
Inv.GM.surv<-function(Sx,a,b,c) if (b<1e-12){
  return(Inv.C.surv(Sx,a+c))
} else {
  z=(a*exp(-(b*log(Sx) - a)/c))/c #detect numerical problems, they are close to Gompertz
  Index=((c<=0) | (z>1e280))
  InvSur=Sx*NA
  InvSur[Index]=log((1-log(Sx[Index])*(b/a)))/b #Gompertz
  InvSur[!Index]= -(log(Sx[!Index]) + (c*gsl:::lambert_W0( z[!Index]))/b - (a/b))/c
  return(InvSur)
}
GM.max.x<-function(a,b,c) min(1e12,max(c(10,Inv.GM.surv(1e-12*(a+c),a,b,c))))
GM.surv<-function(x,a,b,c) if(b!=0) exp(-c*x -a*(exp(b*x)-1)/b) else exp(-c*x-a*x) 
GM.mort<-function(x,a,b,c) a*exp(b*x)+c

GM.ax<-function(x,a,b,c,step=c(diff(x),diff(x)[length(x)-1]))
  sapply(1:length(x),function(i){
    x.=x[i]
    Lx=integrate(f=function(y) GM.surv(y,a,b,c),lower=x.,upper=x.+step[i], subdivisions=1e4)$value
    lx=GM.surv(x.,a,b,c)
    lx1=GM.surv(x.+step[i],a,b,c)
    ((Lx-step[i]*lx1)/(lx-lx1))/step[i]
  })

GM.e0<-function(a,b,c,max.x=GM.max.x(a,b,c)) if (b<1e-12) {
  return(C.e0(a+c)) 
} else if (c==0) {
  return(G.e0(a,b))
} else return(integrate(f=function (x) (GM.surv(x,a,b,c)),lower=0,upper=max.x,
                        subdivisions=1e4)$value)

GM.var<-function(a,b,c,max.x=GM.max.x(a,b,c)) {
  e0=GM.e0(a,b,c,max.x)
  if (c==0) {
    return(G.var(a,b,c)) 
  } else return(integrate(f=function (x) rm.nan(((x-e0)^2)*GM.surv(x,a,b,c)*GM.mort(x,a,b,c)),
                          lower=0,upper=max.x,subdivisions=1e4)$value)
}

GM.CV<-function(a,b,c,max.x=GM.max.x(a,b,c)) sqrt(GM.var(a,b,c,max.x))/GM.e0(a,b,c,max.x)

GM.sq.e0<-function(a,b,c,max.x=GM.max.x(a,b,c)) if (b<1e-12) {
  return(C.sq.e0(a+c)) 
} else if (c==0) {
  return(G.sq.e0(a,b))
} else return(integrate(f=function (x) rm.nan(GM.surv(x,a,b,c)^2),lower=0,upper=max.x,
                        subdivisions=1e4)$value)
GM.eR<-function(R,a,b,c,max.x=GM.max.x(a,b,c)) if (b<1e-12) {
  return(C.eR(R,a+c)) 
} else if (c==0) {
  return(G.eR(R,a,b))
} else if (GM.surv(R,a,b,c)<=0) return(0) else return(integrate(f=function (x) (GM.surv(x,a,b,c)/
                                                                                  GM.surv(R,a,b,c)),lower=R,upper=max.x,subdivisions=1e4)$value)

GM.dager<-function(a,b,c,max.x=GM.max.x(a,b,c)) if (b<1e-12) {
  return(C.dager(a+c)) 
} else if (c==0) {
  return(G.dager(a,b))
} else return(integrate(f=function (x) rm.nan(-GM.surv(x,a,b,c)*log(GM.surv(x,a,b,c))),
                        lower=0,upper=max.x,subdivisions=1e4)$value)

GM.gini<-function(a,b,c,max.x=GM.max.x(a,b,c)) 1-GM.sq.e0(a,b,c,max.x)/GM.e0(a,b,c,max.x)

GM.mode<-function(a,b,c) if (b<1e-12) {
  return(C.mode(a+c))
} else if (c==0) {
  return(G.mode(a,b))
} else return(optim(fn=function(x) (-GM.surv(x,a,b,c)*GM.mort(x,a,b,c)),par=GM.e0(a,b,c),method='Brent',
                    lower=0,upper=min(1000,Inv.GM.surv(1e-10,a,b,c)), 
                    control=list(abstol=1e-10,maxit=2000))$par)

#####################################
#Gamma-Gompertz
Inv.GG.surv<-function(Sx,a,b,s) if (s<1e-12) {
  return(Inv.G.surv(Sx,a,b))
} else if (b<1e-12) {
  return(Inv.C.surv(Sx,a))
} else return(log((Sx^(-s)-1)*b/a/s+1)/b)
GG.max.x<-function(a,b,s) min(1e12,max(c(10,Inv.GG.surv(1e-12*(a),a,b,s))))
GG.surv<-function(x,a,b,s) if (s!=0) (1+a*s*(exp(b*x)-1)/b)^(-1/s) else G.surv(x,a,b)
GG.mort<-function(x,a,b,s) a*exp(b*x)/(1+a*s*(exp(b*x)-1)/b)

GG.ax<-function(x,a,b,s,step=c(diff(x),diff(x)[length(x)-1]))
  sapply(1:length(x),function(i){
    x.=x[i]
    Lx=integrate(f=function(y) GG.surv(y,a,b,s),lower=x.,upper=x.+step[i], subdivisions=1e4)$value
    lx=GG.surv(x.,a,b,s)
    lx1=GG.surv(x.+step[i],a,b,s)
    ((Lx-step[i]*lx1)/(lx-lx1))/step[i]
  })

GG.e0<-function(a,b,s,max.x=GG.max.x(a,b,s)) if (s<1e-12) {
  return(G.e0(a,b))
} else if (b<1e-12) {
  return(C.e0(a))
} else return(integrate(f=function (x) (GG.surv(x,a,b,s)),lower=0,upper=max.x,
                        subdivisions=1e4)$value) 
GG.sq.e0<-function(a,b,s,max.x=GG.max.x(a,b,s)) if (s<1e-12) {
  return(G.sq.e0(a,b))
} else if (b<1e-12) {
  return(C.sq.e0(a))
} else return(integrate(f=function (x) rm.nan(GG.surv(x,a,b,s)^2),lower=0,upper=max.x,
                        subdivisions=1e4)$value) 
GG.eR<-function(R,a,b,s,max.x=GG.max.x(a,b,s)) if (s<1e-12) {
  return(G.eR(R,a,b))
} else if (b<1e-12) {
  return(C.eR(R,a))
} else if (GG.surv(R,a,b,s)<=0) return(0) else return(integrate(f=function (x) 
  (GG.surv(x,a,b,s)/GG.surv(R,a,b,s)),lower=R,upper=max.x,subdivisions=1e4)$value) 

GG.var<-function(a,b,s,max.x=GG.max.x(a,b,s)) {
  e0=GG.e0(a,b,s,max.x)
  if ((s<1e-12) ) { 
    return(G.var(a,b)) 
  } else return(integrate(f=function (x) rm.nan(((x-e0)^2)*GG.surv(x,a,b,s)*GG.mort(x,a,b,s)),
                          lower=0,upper=max.x,subdivisions=1e4)$value)
}

GG.CV<-function(a,b,s,max.x=GG.max.x(a,b,s)) sqrt(GG.var(a,b,s,max.x))/GG.e0(a,b,s,max.x)

GG.dager<-function(a,b,s,max.x=GG.max.x(a,b,s)) if (s<1e-12) {
  return(G.dager(a,b))
} else if (b<1e-12) {
  return(C.dager(a))
} else return(integrate(f=function (x) rm.nan(-GG.surv(x,a,b,s)*log(GG.surv(x,a,b,s))),
                        lower=0,upper=max.x,subdivisions=1e4)$value) 

GG.gini<-function(a,b,s,max.x=GG.max.x(a,b,s)) 1-GG.sq.e0(a,b,s,max.x)/GG.e0(a,b,s,max.x)

#####################################
#Gamma-Gompertz-Makeham
Inv.GGM.surv<-function(Sx,a,b,s,c,upper=1e7) if (b<1e-12) {
  return(Inv.C.surv(Sx,a+c)) 
} else if ((s<1e-12) && (c==0)) { 
  return(Inv.G.surv(Sx,a,b)) 
} else if (s<1e-12) {
  return(Inv.GM.surv(Sx,a,b,c)) 
} else if (c==0) {
  return(Inv.GG.surv(Sx,a,b,s)) 
} else return(sapply(Sx, function (y) uniroot(f=function(x) 
  (GGM.surv(x,a,b,s,c)-y), lower=0, upper=upper,extendInt='no')$root))

GGM.max.x<-function(a,b,s,c) min(1e7,max(c(10,Inv.GGM.surv(1e-12*(a+c),a,b,s,c))))

GGM.surv<-function(x,a,b,s,c) if ((s!=0)&&(b!=0)) exp(-c*x)*(1+a*s*(exp(b*x)-1)/b)^(-1/s) else GM.surv(x,a,b,c)
GGM.mort<-function(x,a,b,s,c) a*exp(b*x)/(1+a*s*(exp(b*x)-1)/b)+c

GGM.Lx<-function(x,a,b,s,c,step=c(diff(x),diff(x)[length(x)-1])) sapply(1:length(x),function(i) integrate(f=function(y) GGM.surv(y,a,b,s,c),lower=x[i],upper=x[i]+step[i], subdivisions=1e4)$value)

GGM.ax<-function(x,a,b,s,c,step=c(diff(x),diff(x)[length(x)-1]))
  sapply(1:length(x),function(i){
    x.=x[i]
    Lx=integrate(f=function(y) GGM.surv(y,a,b,s,c),lower=x.,upper=x.+step[i], subdivisions=1e4)$value
    lx=GGM.surv(x.,a,b,s,c)
    lx1=GGM.surv(x.+step[i],a,b,s,c)
    ((Lx-step[i]*lx1)/(lx-lx1))/step[i]
  })

GGM.ax(x=0:10, a=0.06,b=0.14,s=0,c=0.00)
GGM.e0<-function(a,b,s,c,max.x=GGM.max.x(a,b,s,c)) if ((s<1e-12) && (c==0)) { 
  return(G.e0(a,b)) 
} else if (s<1e-12) {
  return(GM.e0(a,b,c)) 
} else if (c==0) {
  return(GG.e0(a,b,s)) 
} else return(integrate(f=function (x) (GGM.surv(x,a,b,s,c)),lower=0,upper=max.x,
                        subdivisions=1e4)$value)
GGM.sq.e0<-function(a,b,s,c,max.x=GGM.max.x(a,b,s,c)) if ((s<1e-12) && (c==0)) { 
  return(G.sq.e0(a,b)) 
} else if (s<1e-12) {
  return(GM.sq.e0(a,b,c)) 
} else if (c==0) {
  return(GG.sq.e0(a,b,s)) 
} else return(integrate(f=function (x) rm.nan(GGM.surv(x,a,b,s,c)^2),lower=0,upper=max.x,subdivisions=1e4)$value)
GGM.eR<-function(R,a,b,s,c,max.x=GGM.max.x(a,b,s,c)) if ((s<1e-12) && (c==0)) { 
  return(G.eR(R,a,b)) 
} else if (s<1e-12) {
  return(GM.eR(R,a,b,c)) 
} else if (c==0) {
  return(GG.eR(R,a,b,s)) 
} else if (GGM.surv(R,a,b,s,c)<=0) return(0) else return(
  integrate(f=function (x) (GGM.surv(x,a,b,s,c)/GGM.surv(R,a,b,s,c)),lower=R,upper=max.x,subdivisions=1e4)$value)

GGM.dager<-function(a,b,s,c,max.x=GGM.max.x(a,b,s,c)) if ((s<1e-12) && (c==0)) { 
  return(G.dager(a,b)) 
} else if (s<1e-12) {
  return(GM.dager(a,b,c)) 
} else if (c==0) {
  return(GG.dager(a,b,s)) 
} else return(integrate(f=function (x) rm.nan(-GGM.surv(x,a,b,s,c)*log(GGM.surv(x,a,b,s,c))),
                        lower=0,upper=max.x,subdivisions=1e4)$value)


GGM.var<-function(a,b,s,c,max.x=GGM.max.x(a,b,s,c)) {
  e0=GGM.e0(a,b,s,c,max.x)
  if ((s<1e-12) && (c==0)) { 
    return(G.var(a,b)) 
  } else if (s<1e-12) {
    return(GM.var(a,b,c)) 
  } else if (c==0) {
    return(GG.var(a,b,s)) 
  } else return(integrate(f=function (x) rm.nan(((x-e0)^2)*GGM.surv(x,a,b,s,c)*GGM.mort(x,a,b,s,c)),
                          lower=0,upper=max.x,subdivisions=1e4)$value)
}

GGM.CV<-function(a,b,s,c,max.x=GGM.max.x(a,b,s,c)) sqrt(GGM.var(a,b,s,c,max.x))/GGM.e0(a,b,s,c,max.x)

GGM.mode<-function(a,b,s,c) if ((s<1e-12) && (c==0)) { 
  return(G.mode(a,b)) 
} else if (s<1e-12) {
  return(GM.mode(a,b,c)) 
} else if (c==0) {
  return(GG.mode(a,b,s)) 
} else return(optim(fn=function(x) (-GGM.surv(x,a,b,s,c)*GGM.mort(x,a,b,s,c)),
                    par=GGM.e0(a,b,s,c),method='Brent',lower=0,upper=min(10000,Inv.GGM.surv(1e-10,a,b,s,c)), 
                    control=list(abstol=1e-10,maxit=2500))$par)

GGM.gini<-function(a,b,s,c,max.x=GGM.max.x(a,b,s,c)) 1-GGM.sq.e0(a,b,s,c,max.x)/GGM.e0(a,b,s,c,max.x)
