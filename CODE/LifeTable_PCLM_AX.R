#This is only a small part of the main program 
#The program test the effect of interval length on calculation some life-table measures, like: 
#Keyfitz entropy, Gini, Coefficient of variation, etc.
#Please note that most computationally advanced routines can be found in attached LIB files (especially in LIB_PCLM)

require(magicaxis)
require(gplots)

Dir=paste0(getwd(),'/')
Lib_Dir='./CODE/'

#Here is the source code for some functions used
#They will be placed in my local r-package, currently attached as "source" r-files

source(paste0(Lib_Dir,'LIB_PASH.r'))
source(paste0(Lib_Dir,'LIB_PCLM.r'))
source(paste0(Lib_Dir,'LIB_LT.r'))
source(paste0(Lib_Dir,'LIB_THEOR.r'))

#Creating a custom palette
#n = number of colors
#alpha = transparency
#rgb = rgb format T/F
#colors = color segments
my.pal<-function(n,alpha=1,rgb=F,colors=c('#000064FF','darkblue','blue3','blue2','blue','white','red','red2','red3','darkred','#640000FF')){
  z=col2rgb(colors)
  x=seq(1,3,length.out=length(colors))
  r=approx(x=x,z[1,],n=n)$y/255
  g=approx(x=x,z[2,],n=n)$y/255
  b=approx(x=x,z[3,],n=n)$y/255
  rgb.m <- matrix(c(r, g, b), ncol = 3, dimnames = list(NULL,c("red", "green", "blue")))
  col <- mapply(rgb, r, g, b, alpha)
  if (rgb) 
    attr(col, "rgb") <- cbind(rgb.m, alpha)
  return(col)
}

#Modified version of filled.contour {graphics} function
#The function extends original routines by adding multiple panel plotting, custom labeling
#Additional parameters:
#First.Plot = flag denoting if this is the first panel
#ScaleText = text label for the scale
my.filled.contour<-function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1,length.out = ncol(z)), 
                             z, xlim = range(x, finite = TRUE),my.mfrow=c(1,1) ,First.Plot=T, neg.z=F,ScaleText='',
                             ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), main='',sub='[%]',
                             levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, cex.main,cex.lab,
                             col = color.palette(length(levels) - 1), plot.title, plot.axes, col.main='darkblue',
                             key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
                             axes = TRUE, frame.plot = axes, XXX=0.3,...) {
  if (missing(cex.main)) cex.main=par('cex.main')
  if (missing(cex.lab)) cex.lab=par('cex.lab')
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  #on.exit(par(par.orig))
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  Ele=my.mfrow[1]*my.mfrow[2]
  if (First.Plot){
    new.order=as.vector(sapply(seq(0,Ele*2-2,2),function(k) (k+c(2, 1))))
    if (length(XXX)==0) XXX=lcm(w)
    layout(matrix(new.order, ncol = 2*my.mfrow[2], nrow=1*my.mfrow[1],byrow=T), widths = rep(c(1-XXX, XXX),Ele))
    layout.show(8)
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 1
    par(mar = mar)
  } else {
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 1
    par(mar = mar)
    
  }
  par(las = las)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
              yaxs = "i")
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
  title(main = ScaleText,line=1,cex.main=cex.lab)
  title(xlab='',ylab='')  
  if (missing(key.axes)) {
    if (axes) {
      G=axis(4,labels=F)
      if (neg.z) axis(4,at=G,labels=-G) else axis(4)
    }
  }
  else key.axes
  box()
  if (!missing(key.title)) 
    key.title
  
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  
  
  plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  graphics:::.filled.contour(x, y, z, levels, col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "", sub="")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) {
    title(sub='',main='',cex.lab=cex.lab,cex.main=cex.main,col.main=col.main,...)
    title(main=main,line=1,col.main=col.main,cex.main=cex.main)
  }
  else plot.title
  invisible()
}

#Additional information on the plot
ptext<-function(r.,m.,x=0.05,y=0.1,adj=0,cex=1.1,col='black'){
  m.=round(m.,2)
  r.=round(r.,2)
  txt=paste('Mean absolute difference: ',m.,' %\nRange from ',r.[1],'% to ',r.[2],' %',sep='')
  text(x=x,y=y,col=col,adj=adj,cex=cex,labels=txt)
}

#Performs computations for a space of Gamma-Gompertz-Makeham parameters
#x - age vector
#a,b - vectors for a and b Gompertz parameters
#c,s - scalars for Makeham and heterogeneity parameters
Perform.Analysis<-function(x,a,b,c=0,s=0,trace=T){
  
  #Allocating memory
  Empty=matrix(NA,length(a),length(b))
  Mat.Theo.e0<-Mat.Theo.K<-Mat.Theo.G<-Mat.Theo.CV<-  
    Mat.Empi.e0a<-Mat.Empi.e0b<-Mat.Empi.e0c<-Mat.Empi.e0d<-
    Mat.Empi.K1<- Mat.Empi.K2<- Mat.Empi.K3<- Mat.Empi.K4<-
    Mat.Empi.G1<- Mat.Empi.G2<- Mat.Empi.G3<- Mat.Empi.G4<-
    Mat.Empi.CV1<-Mat.Empi.CV2<-Mat.Empi.CV3<-Mat.Empi.CV4<-Empty
  
  for (ai in 1:length(a))
    for (bi in 1:length(b)){
      
      pos=round(100*((ai-1)*length(b)+bi)/(length(a)*length(b)),2)
      cat(paste('\r\rCompleteted: ',pos,'%',paste(rep(' ',20),sep=''),sep=''))
      a_=a[ai]
      b_=b[bi]
      
      #Creating "infinite" population for a set of Gamma-Gompertz-Makeham parameters
      lx=(1e12*GGM.surv(x=x,a=a_,b=b_,s,c))
      nlx=round(lx)
      x_=x[nlx>0] #tailoring x
      lx=nlx[nlx>0] #tailoring lx
      dx=lx_2_dx(lx)
      
      #fittin PCLM model
      pclmres=PCLM(x_,dx,out.step=NULL,bs.df.max=50)
      #ploting PCLM model
      if (trace){
        graphics.off()
        PCLM.plot(pclmres)
      }
      #Calculating life-table measures
      HT1=PASH(lx=pclmres$grouped$lx,x=pclmres$grouped$x,iax=pclmres$grouped$Ax)$M
      HT2=PASH(lx=pclmres$raw$lx,x=pclmres$raw$x,iax=pclmres$raw$Ax)$M
      HT3=PASH(lx=lx,x=x_,iax=pclmres$grouped$Ax[1:length(x_)])$M
      HT4=PASH(lx=lx,x=x_,iax=0.5)$M
      
      #preparing output
      #theoretical measures
      Mat.Theo.e0[ai,bi]=GGM.e0(a_,b_,s,c); Mat.Theo.K[ai,bi]=(GGM.dager(a_,b_,s,c)/Mat.Theo.e0[ai,bi]);
      Mat.Theo.G[ai,bi]=GGM.gini(a_,b_,s,c); Mat.Theo.CV[ai,bi]=GGM.CV(a_,b_,s,c)
      #empirical measures                   
      Mat.Empi.e0a[ai,bi]=HT1[1]; Mat.Empi.e0b[ai,bi]=HT2[1]
      Mat.Empi.e0c[ai,bi]=HT3[1]; Mat.Empi.e0d[ai,bi]=HT4[1]
      Mat.Empi.G1[ai,bi]=HT1[3]; Mat.Empi.G2[ai,bi]=HT2[3]
      Mat.Empi.G3[ai,bi]=HT3[3]; Mat.Empi.G4[ai,bi]=HT4[3]
      Mat.Empi.CV1[ai,bi]=HT1[4]; Mat.Empi.CV2[ai,bi]=HT2[4]
      Mat.Empi.CV3[ai,bi]=HT3[4]; Mat.Empi.CV4[ai,bi]=HT4[4]
      Mat.Empi.K1[ai,bi]=HT1[2]; Mat.Empi.K2[ai,bi]=HT2[2]
      Mat.Empi.K3[ai,bi]=HT3[2]; Mat.Empi.K4[ai,bi]=HT4[2]
    }
  
  cat('\n')
  A=round(a,5)
  B=round(b,5)
  colnames(Mat.Empi.e0a)=B; rownames(Mat.Empi.e0a)=A
  colnames(Mat.Empi.e0b)=B; rownames(Mat.Empi.e0b)=A
  colnames(Mat.Empi.e0c)=B; rownames(Mat.Empi.e0c)=A
  colnames(Mat.Empi.e0d)=B; rownames(Mat.Empi.e0d)=A
  colnames(Mat.Empi.K1)=B; rownames(Mat.Empi.K1)=A
  colnames(Mat.Empi.K2)=B; rownames(Mat.Empi.K2)=A
  colnames(Mat.Empi.K3)=B; rownames(Mat.Empi.K3)=A
  colnames(Mat.Empi.K4)=B; rownames(Mat.Empi.K4)=A
  colnames(Mat.Empi.G1)=B; rownames(Mat.Empi.G1)=A
  colnames(Mat.Empi.G2)=B; rownames(Mat.Empi.G2)=A
  colnames(Mat.Empi.G3)=B; rownames(Mat.Empi.G3)=A
  colnames(Mat.Empi.G4)=B; rownames(Mat.Empi.G4)=A
  colnames(Mat.Empi.CV1)=B; rownames(Mat.Empi.CV1)=A
  colnames(Mat.Empi.CV2)=B; rownames(Mat.Empi.CV2)=A
  colnames(Mat.Empi.CV3)=B; rownames(Mat.Empi.CV3)=A
  colnames(Mat.Empi.CV4)=B; rownames(Mat.Empi.CV4)=A
  colnames(Mat.Theo.G)=B; rownames(Mat.Theo.G)=A
  colnames(Mat.Theo.CV)=B; rownames(Mat.Theo.CV)=A
  colnames(Mat.Theo.e0)=B; rownames(Mat.Theo.e0)=A
  colnames(Mat.Theo.K)=B; rownames(Mat.Theo.K)=A
  
  list(Mat.Theo.K=Mat.Theo.K, Mat.Theo.e0=Mat.Theo.e0,
       Mat.Theo.G=Mat.Theo.G, Mat.Theo.CV=Mat.Theo.CV,
       Mat.Empi.K1=Mat.Empi.K1, Mat.Empi.K2=Mat.Empi.K2,
       Mat.Empi.K3=Mat.Empi.K3, Mat.Empi.K4=Mat.Empi.K4,
       Mat.Empi.e0a=Mat.Empi.e0a, Mat.Empi.e0b=Mat.Empi.e0b,
       Mat.Empi.e0c=Mat.Empi.e0c, Mat.Empi.e0d=Mat.Empi.e0d,
       Mat.Empi.G1=Mat.Empi.G1, Mat.Empi.G2=Mat.Empi.G2,
       Mat.Empi.G3=Mat.Empi.G3, Mat.Empi.G4=Mat.Empi.G4,
       Mat.Empi.CV1=Mat.Empi.CV1, Mat.Empi.CV2=Mat.Empi.CV2,
       Mat.Empi.CV3=Mat.Empi.CV3, Mat.Empi.CV4=Mat.Empi.CV4,
       a=a,b=b,c=c,s=s, x=x,Age.Interval=diff(x)[1])
}

#Perform plotting based on previouse calculations
#DDD - Results to be presented
#Panel.id.pos - position of panel letter
#Panel.id.col - color of panel letter
#Stat.info.col - Color of stat. text in plot
#PAL1, PAL2 - color segments for different kinds of plot
Perform.Plotting<-function(DDD,Dir=paste0(getwd(),'/'),
                           Panel.id.pos=0.025,Panel.id.col='red3',Stat.info.col='black',
                           ColSegm1=c('blue','white','red'),
                           ColSegm2=c('#000064FF','darkblue','blue3','blue2','blue',
                                      'white','red','red2','red3','darkred','#640000FF')) {
  
  PAL1=function(n,alpha=1,rgb=F) my.pal(n=n,alpha=alpha,rgb=rgb,colors=ColSegm1)
  PAL2=function(n,alpha=1,rgb=F) my.pal(n=n,alpha=alpha,rgb=rgb,colors=ColSegm2)
  
  #Age interval
  S=DDD$Age.Interval
  
  #Vectors with Gompertz parameters
  a.=DDD$a
  b.=DDD$b
  
  #Basic statistics
  RelDiff<-function(Theo,Empi) 100*((Theo-Empi)/Theo)[Theo!=0]
  
  e0ar=range(RelDiff(DDD$Mat.Theo.e0,DDD$Mat.Empi.e0a))
  e0am=mean(abs(RelDiff(DDD$Mat.Theo.e0,DDD$Mat.Empi.e0a)))
  e0br=range(RelDiff(DDD$Mat.Theo.e0,DDD$Mat.Empi.e0b))
  e0bm=mean(abs(RelDiff(DDD$Mat.Theo.e0,DDD$Mat.Empi.e0b)))
  e0cr=range(RelDiff(DDD$Mat.Theo.e0,DDD$Mat.Empi.e0c))
  e0cm=mean(abs(RelDiff(DDD$Mat.Theo.e0,DDD$Mat.Empi.e0c)))
  e0dr=range(RelDiff(DDD$Mat.Theo.e0,DDD$Mat.Empi.e0d))
  e0dm=mean(abs(RelDiff(DDD$Mat.Theo.e0,DDD$Mat.Empi.e0d)))
  
  K1r=range(RelDiff(DDD$Mat.Theo.K,DDD$Mat.Empi.K1))
  K1m=mean(abs(RelDiff(DDD$Mat.Theo.K,DDD$Mat.Empi.K1)))
  K2r=range(RelDiff(DDD$Mat.Theo.K,DDD$Mat.Empi.K2))
  K2m=mean(abs(RelDiff(DDD$Mat.Theo.K,DDD$Mat.Empi.K2)))
  K3r=range(RelDiff(DDD$Mat.Theo.K,DDD$Mat.Empi.K3))
  K3m=mean(abs(RelDiff(DDD$Mat.Theo.K,DDD$Mat.Empi.K3)))
  K4r=range(RelDiff(DDD$Mat.Theo.K,DDD$Mat.Empi.K4))
  K4m=mean(abs(RelDiff(DDD$Mat.Theo.K,DDD$Mat.Empi.K4))) 
  
  G1r=range(RelDiff(DDD$Mat.Theo.G,DDD$Mat.Empi.G1))
  G1m=mean(abs(RelDiff(DDD$Mat.Theo.G,DDD$Mat.Empi.G1)))
  G2r=range(RelDiff(DDD$Mat.Theo.G,DDD$Mat.Empi.G2))
  G2m=mean(abs(RelDiff(DDD$Mat.Theo.G,DDD$Mat.Empi.G2)))
  G3r=range(RelDiff(DDD$Mat.Theo.G,DDD$Mat.Empi.G3))
  G3m=mean(abs(RelDiff(DDD$Mat.Theo.G,DDD$Mat.Empi.G3)))
  G4r=range(RelDiff(DDD$Mat.Theo.G,DDD$Mat.Empi.G4))
  G4m=mean(abs(RelDiff(DDD$Mat.Theo.G,DDD$Mat.Empi.G4)))
  
  CV1r=range(RelDiff(DDD$Mat.Theo.CV,DDD$Mat.Empi.CV1))
  CV1m=mean(abs(RelDiff(DDD$Mat.Theo.CV,DDD$Mat.Empi.CV1)))
  CV2r=range(RelDiff(DDD$Mat.Theo.CV,DDD$Mat.Empi.CV2))
  CV2m=mean(abs(RelDiff(DDD$Mat.Theo.CV,DDD$Mat.Empi.CV2)))
  CV3r=range(RelDiff(DDD$Mat.Theo.CV,DDD$Mat.Empi.CV3))
  CV3m=mean(abs(RelDiff(DDD$Mat.Theo.CV,DDD$Mat.Empi.CV3)))
  CV4r=range(RelDiff(DDD$Mat.Theo.CV,DDD$Mat.Empi.CV4))
  CV4m=mean(abs(RelDiff(DDD$Mat.Theo.CV,DDD$Mat.Empi.CV4)))
  
  ###### BEGINING OF FIRST PLOT #######
  
  #z axis calculation
  R=ceiling(max(abs(K1r),abs(K2r),abs(K3r))*10)/10
  
  Xsize=7;Ysize=6;FName=paste('AX_G_K_',S,'.tiff',sep='');Fac=333
  graphics.off()
  tiff(paste(Dir,FName,sep=''), width=Xsize*Fac, height=Ysize*Fac, res=Fac,compression="lzw",pointsize=10)
  par(mar=c(5.1,4.1,2.1,2.1))
  par(cex=0.9);par(cex.axis=1.2);
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.K-DDD$Mat.Empi.K1)/DDD$Mat.Theo.K,nlevels=32,
                    color.palette=PAL1,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    XXX=0.225,my.mfrow=c(2,2),First.Plot=T,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    main=expression(bold('PCLM, grouped x, Ax computed')),
                    ScaleText=expression('[%]'),
                    las=0)
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('A')),cex=1.5,col=Panel.id.col)
  ptext(K1r,K1m,col=Stat.info.col)
  
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.K-DDD$Mat.Empi.K2)/DDD$Mat.Theo.K,nlevels=32,color.palette=PAL1,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    my.mfrow=c(2,2),First.Plot=F,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    main=expression(bold('PCLM, not grouped x, Ax = 0.5')),
                    ScaleText=expression('[%]'),
                    las=0)
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('B')),cex=1.5,col=Panel.id.col)
  ptext(K2r,K2m,col=Stat.info.col)
  
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.K-DDD$Mat.Empi.K3)/DDD$Mat.Theo.K,nlevels=32,color.palette=PAL1,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    my.mfrow=c(2,2),First.Plot=F,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    sub=expression(paste(100,'*',(K[Theoretical]-K[Empirical])/K[Theoretical])),
                    main=expression(bold('No PCLM, Ax computed from PCLM')),
                    ScaleText=expression('[%]'),
                    las=0)
  
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('C')),cex=1.5,col=Panel.id.col)
  ptext(K3r,K3m,col=Stat.info.col)
  
  #z axis calculation
  R=ceiling(max(abs(K4r))*10)/10
  
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.K-DDD$Mat.Empi.K4)/DDD$Mat.Theo.K,nlevels=32,
                    color.palette=PAL2,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    my.mfrow=c(2,2),First.Plot=F,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    sub=expression(paste(100,'*',(K[Theoretical]-K[Empirical])/K[Theoretical])),
                    main=expression(bold('No PCLM, Ax = 0.5')),
                    ScaleText=expression('[%]'),
                    las=0)
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('D')),cex=1.5,col=Panel.id.col)
  ptext(K4r,K4m,col=Stat.info.col)
  dev.off()
  ###### END OF FIRST PLOT #######
  
  
  ###### BEGINING OF SECOND PLOT #######
  Xsize=7;Ysize=6;FName=paste('AX_G_G_',S,'.tiff',sep='');Fac=333
  graphics.off()
  tiff(paste(Dir,FName,sep=''), width=Xsize*Fac, height=Ysize*Fac, res=Fac,compression="lzw",pointsize=10)
  par(mar=c(5.1,4.1,2.1,2.1))
  par(cex=0.9);par(cex.axis=1.2);
  #z axis calcualtion
  R=ceiling(max(abs(G1r),abs(G2r),abs(G3r),abs(G4r))*10)/10
  
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.G-DDD$Mat.Empi.G1)/DDD$Mat.Theo.G,nlevels=32,
                    color.palette=PAL1,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    XXX=0.225,my.mfrow=c(2,2),First.Plot=T,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    main=expression(bold('PCLM, grouped x, Ax computed')),
                    ScaleText=expression('[%]'),
                    las=0)
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('A')),cex=1.5,col=Panel.id.col)
  ptext(G1r,G1m,col=Stat.info.col)
  
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.G-DDD$Mat.Empi.G2)/DDD$Mat.Theo.G,nlevels=32,color.palette=PAL1,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    my.mfrow=c(2,2),First.Plot=F,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    main=expression(bold('PCLM, not grouped x, Ax = 0.5')),
                    ScaleText=expression('[%]'),
                    las=0)
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('B')),cex=1.5,col=Panel.id.col)
  ptext(G2r,G2m,col=Stat.info.col)
  
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.G-DDD$Mat.Empi.G3)/DDD$Mat.Theo.G,nlevels=32,color.palette=PAL1,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    my.mfrow=c(2,2),First.Plot=F,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    sub=expression(paste(100,'*',(K[Theoretical]-K[Empirical])/K[Theoretical])),
                    main=expression(bold('No PCLM, Ax computed from PCLM')),
                    ScaleText=expression('[%]'),
                    las=0)
  
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('C')),cex=1.5,col=Panel.id.col)
  ptext(G3r,G3m,col=Stat.info.col)
  
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.K-DDD$Mat.Empi.K4)/DDD$Mat.Theo.K,nlevels=32,
                    color.palette=PAL1,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    my.mfrow=c(2,2),First.Plot=F,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    sub=expression(paste(100,'*',(K[Theoretical]-K[Empirical])/K[Theoretical])),
                    main=expression(bold('No PCLM, Ax = 0.5')),
                    ScaleText=expression('[%]'),
                    las=0)
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('D')),cex=1.5,col=Panel.id.col)
  ptext(G4r,G4m,col=Stat.info.col)
  dev.off()
  ###### END OF SECOND PLOT #######
  
  
  ###### BEGINING OF THIRD PLOT #######
  Xsize=7;Ysize=6;FName=paste('AX_G_CV_',S,'.tiff',sep='');Fac=333
  graphics.off()
  #z axis calcualtion
  R=ceiling(max(abs(CV1r),abs(CV2r),abs(CV3r))*10)/10
  
  tiff(paste(Dir,FName,sep=''), width=Xsize*Fac, height=Ysize*Fac, res=Fac,compression="lzw",pointsize=10)
  par(mar=c(5.1,4.1,2.1,2.1))
  par(cex=0.9);par(cex.axis=1.2);
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.CV-DDD$Mat.Empi.CV1)/DDD$Mat.Theo.CV,nlevels=32,
                    color.palette=PAL1,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    XXX=0.225,my.mfrow=c(2,2),First.Plot=T,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    main=expression(bold('PCLM, grouped x, Ax computed')),
                    ScaleText=expression('[%]'),
                    las=0)
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('A')),cex=1.5,col=Panel.id.col)
  ptext(CV1r,CV1m,col=Stat.info.col)
  
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.CV-DDD$Mat.Empi.CV2)/DDD$Mat.Theo.CV,nlevels=32,color.palette=PAL1,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    my.mfrow=c(2,2),First.Plot=F,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    main=expression(bold('PCLM, not grouped x, Ax = 0.5')),
                    ScaleText=expression('[%]'),
                    las=0)
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('B')),cex=1.5,col=Panel.id.col)
  ptext(CV2r,CV2m,col=Stat.info.col)
  
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.CV-DDD$Mat.Empi.CV3)/DDD$Mat.Theo.CV,nlevels=32,color.palette=PAL1,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    my.mfrow=c(2,2),First.Plot=F,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    main=expression(bold('No PCLM, Ax computed from PCLM')),
                    ScaleText=expression('[%]'),
                    las=0)
  
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('C')),cex=1.5,col=Panel.id.col)
  ptext(CV3r,CV3m,col=Stat.info.col)
  
  #z axis calculation
  R=ceiling(max(abs(CV4r))*10)/10
  
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.CV-DDD$Mat.Empi.CV4)/DDD$Mat.Theo.CV,nlevels=32,
                    color.palette=PAL2,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    my.mfrow=c(2,2),First.Plot=F,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    main=expression(bold('No PCLM, Ax = 0.5')),
                    ScaleText=expression('[%]'),
                    las=0)
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('D')),cex=1.5,col=Panel.id.col)
  ptext(CV4r,CV4m,col=Stat.info.col)
  dev.off()
  ###### END OF THIRD PLOT #######
  
  ###### BEGINING OF FOURTH PLOT #######
  Xsize=7;Ysize=6;FName=paste('AX_G_e0_',S,'.tiff',sep='');Fac=333
  graphics.off()
  R=ceiling(max(abs(e0ar),abs(e0br),abs(e0cr))*10)/10
  tiff(paste(Dir,FName,sep=''), width=Xsize*Fac, height=Ysize*Fac, res=Fac,compression="lzw",pointsize=10)
  par(mar=c(5.1,4.1,2.1,2.1))
  par(cex=0.9);par(cex.axis=1.2);
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.e0-DDD$Mat.Empi.e0a)/DDD$Mat.Theo.e0,nlevels=32,
                    color.palette=PAL1,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    XXX=0.225,my.mfrow=c(2,2),First.Plot=T,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    main=expression(bold('PCLM, grouped x, Ax computed')),
                    ScaleText=expression('[%]'),
                    las=0)
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('A')),cex=1.5,col=Panel.id.col)
  ptext(e0ar,e0am,col=Stat.info.col)
  
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.e0-DDD$Mat.Empi.e0b)/DDD$Mat.Theo.e0,nlevels=32,color.palette=PAL1,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    my.mfrow=c(2,2),First.Plot=F,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    main=expression(bold('PCLM, not grouped x, Ax = 0.5')),
                    ScaleText=expression('[%]'),
                    las=0)
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('B')),cex=1.5,col=Panel.id.col)
  ptext(e0br,e0bm,col=Stat.info.col)
  
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.e0-DDD$Mat.Empi.e0c)/DDD$Mat.Theo.e0,nlevels=32,color.palette=PAL1,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    my.mfrow=c(2,2),First.Plot=F,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    main=expression(bold('No PCLM, Ax computed from PCLM')),
                    ScaleText=expression('[%]'),
                    las=0)
  
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('C')),cex=1.5,col=Panel.id.col)
  ptext(e0cr,e0cm,col=Stat.info.col)
  
  #z axis calculation
  R=ceiling(max(abs(e0dr))*10)/10
  my.filled.contour(a.,b.,100*(DDD$Mat.Theo.e0-DDD$Mat.Empi.e0d)/DDD$Mat.Theo.e0,nlevels=91,
                    color.palette=PAL2,xlab="Gompertz's a",ylab="Gompertz's b",cex.lab=1.5,
                    my.mfrow=c(2,2),First.Plot=F,neg.z=F,cex.main=1.3,zlim=c(-R,R),
                    main=expression(bold('No PCLM, Ax = 0.5')),
                    ScaleText=expression('[%]'),
                    las=0)
  mtext(at=Panel.id.pos,line=-2.5,expression(bold('D')),cex=1.5,col=Panel.id.col)
  ptext(e0dr,e0dm,col=Stat.info.col)
  dev.off()
  ###### END OF FOURTH PLOT #######
  
}


###############################################################################################
# MAIN PROGRAM
###############################################################################################

a.=seq(0.0125,0.25,0.0125)
b.=seq(0.025,0.3,0.0125)
c=0
s=0

load.flag=T # Set to F if you want re-run computations

if (load.flag) {
  load(file=paste(Dir,'AX3_5.dat',sep=''))
  load(file=paste(Dir,'AX3_2.dat',sep=''))
  load(file=paste(Dir,'AX3_1.dat',sep='')) 
} else {
  #Quite long computations!!!
  DDD_5=Perform.Analysis(x=seq(0,1e4,5),a=a.,b=b.,c=c,s=s)
  DDD_2=Perform.Analysis(x=seq(0,1e4,2),a=a.,b=b.,c=c,s=s)
  DDD_1=Perform.Analysis(x=seq(0,1e4,1),a=a.,b=b.,c=c,s=s)  
  save(DDD_5,file=paste(Dir,'AX3_5.dat',sep=''))
  save(DDD_2,file=paste(Dir,'AX3_2.dat',sep=''))
  save(DDD_1,file=paste(Dir,'AX3_1.dat',sep=''))
}

Perform.Plotting(DDD=DDD_1)
Perform.Plotting(DDD=DDD_2)
Perform.Plotting(DDD=DDD_5)
