#Calculating life-table and its measures
#lx - survivorship or number of surviving
#x - age
#iax - life-table Ax
#open.int - account for open interval
#pop - population size
#...Method - different methods of approximation
PASH<-function(lx,x=(1:length(lx))-1,iax=0.5,open.int=F,pop=lx[1],
               LXMethod=c('Preston','Mortality'),
               GiniMethod=c('PASH','PASCARIU2016','DANKO2016','D2016',
                            'SKHOLNIKOVANDREEW2014','SA2014','DANKO2016','ALL'),
               CVMethod=c('PASH','P2016','PASCARIU2016','SKHOLNIKOVANDREEW2014',
                          'SA2014','WRYCZA2014','DW2014','ALL')){
                 
  is.mono<-function(V) all(diff(V)>=0) | all(diff(V)<=0)           
                 
  CVMethod=toupper(CVMethod)[1]
  GiniMethod=toupper(GiniMethod)[1]
  LXMethod=toupper(LXMethod)[1]
  
  if (!(GiniMethod %in% c('PASH','PASCARIU2016','DANKO2016','D2016',
                          'SKHOLNIKOVANDREEW2014','SA2014','DANKO2016','ALL'))) 
    stop('Unknown method for Gini coeficient')
  if (!(CVMethod %in% c('PASH','P2016','PASCARIU2016','SKHOLNIKOVANDREEW2014',
                        'SA2014','WRYCZA2014','DW2014','ALL')))
    stop('Unknown method for coeficient of variation')
  
  lx=pop*lx/lx[1]
  MaxLifeSpan=x[round(lx,1)>=0.5]
  MaxLifeSpan=MaxLifeSpan[length(MaxLifeSpan)]
  
  #test for some obvious errors
  if (length(MaxLifeSpan)==0) {
    cat('lx: ',lx,'\n')
    cat('x: ',x,'\n')
    cat('Pop: ',pop,'\n')
    MaxLifeSpan=NA
  }
  lx[is.na(lx)]=0
  
  #test if lx is monotonic
  if (!is.mono(lx)) {
    cat('lx is not monotonic.\n')
    return(EmptyShapes)
  }
  
  #Basic life-table computation
  step=diff(x)
  Last=length(x)
  step=c(step,step[Last-1]) #Preston's n
  dx=c(-diff(lx),lx[Last]) #assuming a cohort lifetable
  qx=dx/lx
  ax=iax*step #equivalent to n_a_x = n * a_x
  if (toupper(LXMethod)=='MORTALITY') {
    Lx=dx/(-log(1-qx)/step) #
  } else if (toupper(LXMethod)=='PRESTON') {
    Lx=step*(lx-dx)+ax*dx #or Lx=step*lx-(1-ax)*dx 
  } else stop('Unknown Lx method')
  mx2=-log(1-qx) #in half of interval
  if (open.int) Lx[Last]=lx[Last]*ax[Last]  
  mx=dx/Lx
  Tx=(sum(Lx)-c(0,cumsum(Lx)))[1:Last]
  if (open.int) Tx[Last]=Lx[Last]
  ex=Tx/lx
  ex[lx==0]=0
  e0=ex[1]
  Ax=ax/step
  
  #Life-table measures
  
  #############
  #Keyfitz Entropy
  tmp=dx*(c(ex[2:Last],0)*Ax+ex*(1-Ax))
  tmp=c(sum(tmp))
  e.dager=(1/lx[1])*tmp
  K=(e.dager/ex[1])
  K[lx==0]=0
  K=K[1]
  
  #############
  #Coeficient of variation
  if ((CVMethod=='SKHOLNIKOVANDREEW2014')||(CVMethod==4)||(CVMethod=='SA2014')||(toupper(CVMethod)=='ALL')) {
    #Coef of var #4 classic #1
    y=x+ax
    tmp=y*y*dx
    tmp=sum(tmp)
    STD=suppressWarnings(sqrt((tmp/lx[1])-(ex+x)[1]^2))
    CV4=STD/(e0+x[1]) 
    CV=CV4
  } 
  if ((CVMethod=='WRYCZA2014')||(CVMethod==1)||(CVMethod=='dW2014')||(toupper(CVMethod)=='ALL')) {
    #Coef of var #1 - discrette form of Wrycza 2014, Trapezoidal rule
    tmp=dx*((c(ex[2:Last],0)^2)*Ax+(ex^2)*(1-Ax))/lx[1]
    CV1=sqrt(sum(tmp))/e0
    CV=CV1
  } 
  if ((CVMethod=='PASCARIU2016')||(CVMethod==2)||(CVMethod=='PASH')||(CVMethod=='P2016')||(toupper(CVMethod)=='ALL')) {
    #Coef of var #2 - MP 2016
    tmp=dx*((c(ex[2:Last],0)*Ax+ex*(1-Ax))^2)/lx[1]
    CV2=sqrt(sum(tmp))/e0
    CV=CV2
  } 
  if ((CVMethod=='DANKO2016')||(CVMethod==3)||(toupper(CVMethod)=='ALL')) {
    #Coef of var #3 - Danko 2016, 
    tmp=dx*((x+step-e0)^3-(x-e0)^3)/(3*step)/lx[1]
    CV3=sqrt(sum(tmp))/e0
    CV=CV3
  }
  
  if (toupper(CVMethod)=='ALL') CV=c(CV1,CV2,CV3,CV4)
  
  #############
  #Gini coeficient
  if ((GiniMethod=='SKHOLNIKOVANDREEW2014')||(GiniMethod==1)||(GiniMethod=='SA2014')||(toupper(GiniMethod)=='ALL')) {
    #Gini coeficient #1 - Classic 1, Shkolnikov and Andreev 2014
    tmp=step*(c(lx[2:Last],0)^2+Ax*(lx^2-c(lx[2:Last],0)^2))
    tmp=sum(tmp)
    G=(1-(1/(lx[1]*lx[1]*e0))*tmp)[1]
    G1=G
  }
  
  if ((GiniMethod=='PASCARIU2016')||(GiniMethod==2)||(GiniMethod=='PASH')||(GiniMethod=='P2016')||(toupper(GiniMethod)=='ALL')) {
    #Gini - MP 2016
    tmp=step*(lx*(1-Ax)+Ax*c(lx[2:Last],0))^2
    tmp=sum(tmp)
    G=(1-(1/(lx[1]*lx[1]*e0))*tmp)[1]
    G2=G
  }  
  
  if ((GiniMethod=='DANKO2016')||(GiniMethod==3)||(GiniMethod=='D2016')||(toupper(GiniMethod)=='ALL')) {
    #Gini - Danko 2016
    Z=dx*dx/6
    tmp=step*(lx^2*(1-Ax)+Ax*c(lx[2:Last],0)^2-Z)
    tmp=sum(tmp)
    G=(1-(1/(lx[1]*lx[1]*e0))*tmp)[1]
    G3=G
  }
  
  if (toupper(GiniMethod)=='ALL') G=c(G1,G2,G3)
  
  #LifeTable
  LT=as.data.frame(cbind(x=x,lx=lx,dx=dx,qx=qx,Ax=iax,ax=ax,Lx=Lx,Tx=Tx,mx=mx,ex=ex))
  
  Shape.G=1/G-1
  Shape.H=1-K[1]
  Shape.CV=1-CV[1]
  Shape.Crude=MaxLifeSpan/e0
  
  #Output
  Measures=c(e0=e0,
             Entropy=K[1],
             Gini=G,
             CV=CV,
             MaxLifeSpan=MaxLifeSpan,
             Shape.H=Shape.H,
             Shape.G=Shape.G,
             Shape.CV=Shape.CV,
             Shape.Crude=Shape.Crude
  )
  list(M=Measures,LT=LT)
}

#to be consistent with older version of code
CalculateShapeMeasures=PASH
CalculateShapeMeasures.lx=PASH