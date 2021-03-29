
	program covid

	implicit none
c*********************************************************

	integer seed1,nsum,sigma(-1:1),ifinal
	integer i,j,k,l,m,imin,ini,ndayt,kloop
	real Br(-100:5000),laplBr(-100:5000)
	real muj(-100:5000),ran2,a,ap,xx,yp,y
	real*8  LL,LLp,LL1,LL2,LL3,dLL3,S(0:5000),mu,dmu,dBr
	real*8 f(0:5000),g,dtau,xt,gnorm,taup
	real n(-100:5000),alfa(0:5000)
	character finput*40,foutput*40

c*********************************************************


C****************************************************************
C******************** MAIN CODE *********************************
C****************************************************************


c*****  INPUT PARAMETERS **********
	
	print*,'Please, enter the name of the input file'
	read(*,'(a)')finput
	print*,'Please, enter the name of the output file'
	read(*,'(a)')foutput
	print*,'Please, enter the Gamma distribution parameter tau'
	read(*,*)taup 
	print*,'Please, enter the Gamma distribution parameter a'
	read(*,*)ap
		
c*****  READING INPUT FILE *******
	
	open(55,file=finput)
	do i=1,10000
	   read(55,*,end=99)xx,n(i) !days and n(i)
	   nsum=nsum+n(i)
	enddo
 99	ndayt=i-1 !number of total days
	ini=2	!Initial day for LL evaluation
		
c*****  Initial setting	**********

	sigma(-1)=1
	sigma(0)=-2
	sigma(1)=1
 
	do i=imin,ndayt+1
	   laplBr(i)=0  !
	   alfa(i)=(1/3.5)*sqrt(nsum*1./5E5) !alfa=1/V in Eq.3
	enddo

c***** Integral of Gamma distribution *******	
	
	i=1
	g=0d0
	dtau=taup/100000d0 !integration step		 
	do j=1,int(taup/dtau)*300
	   xt=j*dtau
	   g=g+dtau*(xt**(ap-1d0)*exp(-xt*1d0/taup))
	   if(xt.ge.i)then
	      f(i)=g  
	      i=i+1
	   endif
	enddo
	ifinal=i-1
	gnorm=1d0/g   !normalization
	f(0)=0
	do i=1,ndayt+100
	   if(i.le.ifinal)then
	      f(i)=f(i)*gnorm
	   else
	      f(i)=f(ifinal)
	   endif
	enddo
	
c*****  Initial evaluation for the first 5 days of Br(i) ****
c*****  by means of the Walling-Teunis algorithm ****
	
	do i=2,5
	   Br(i)=0  !Br is R_s(i)
	   do j=0,ndayt-i
	      s(i)=0d0
	      do k=0,i+j
		 s(i)=s(i)+n(i+j-k)*(f(k)-f(k-1))
	      enddo
	      if(s(i).gt.0)Br(i)=Br(i)+n(i+j)*(f(j)-f(j-1))/s(i)
	   enddo
	enddo
	Br(1)=Br(2)
	mu=0
	
c******* Initial setting of Br(i) for the subsequent days ****
c******* Initial setting of immigrants rate muj(i) *****
	
	do i=5,ndayt+1
	   Br(i)=Br(2)-0.05*(i-2)
	   Br(i)=max(Br(i),0.01)
	   muj(i)=0
	enddo
	do i=ini,ndayt
	   laplBr(i)=(Br(i-1)+Br(i+1)-2*Br(i))
	   mu=mu+muj(i)      
	enddo
	
c****   Initial estimate of LL *******
	
	LL2=0d0
	LL3=0d0
	do i=2,ndayt
	   S(i)=0
	   do j=1,i-1
	      S(i)=S(i)+n(j)*Br(j)*(f(i-j)-f(i-j-1))
	   enddo
	   LL3=LL3+alfa(i)*laplBr(i)**2
	enddo
	do i=1,ndayt
	   LL2=LL2+n(i)*Br(i)*f(ndayt-i)
	enddo
	LL1=0d0
	do i=2,ndayt
	   LL1=LL1+n(i)*dlog(S(i)+muj(i))
	enddo
	LL=LL1-mu-LL2-LL3 !Initial LL

c********************************************
c***  Markov-Chain-Monte-Carlo method *******
c********************************************
	
	do kloop=1,300	!loop on Monte Carlo steps
	   do l=ini,ndayt*20
	      j=int(ran2(seed1)*(ndayt-2))+2   !chose a random day
 41		 dBr=max(Br(j),0.2)*0.02*(ran2(seed1)-.5) !trial dBr 
		 if(Br(j)+dBr.le.0.or.Br(j)+dBr.gt.9)goto 41 !lower and upper bounds for Br(i)
!                If j is an internal day
		 if(j.ne.ndayt-1)then   
		    dLL3=0
		    do m=-1,1
		       dLL3=dll3+alfa(j+m)*((laplBr(j+m)+sigma(m)*dBr)**2-
     &		            (laplBr(j+m))**2)
		    enddo
		    LLp=0d0
		    do i=j+1,ndayt
		       LLp=LLp+n(i)*dlog(S(i)+muj(i)+dBr*n(j)*(f(i-j)-f(i-j-1)))
		    enddo
		    do i=ini,j
		       LLp=LLp+n(i)*dlog(S(i)+muj(i))
		    enddo
		    LLp=LLp-mu-LL2-dBr*n(j)*f(ndayt-j)-(LL3+dLL3) !trial LL
		    if(LLp.gt.LL)then !acceptance of the trial
		       do i=j+1,ndayt
			  S(i)=S(i)+dBr*n(j)*(f(i-j)-f(i-j-1))
		       enddo
		       LL2=LL2+dBr*n(j)*f(ndayt-j)
		       LL3=LL3+dLL3
		       LL=LLp
		       Br(j)=Br(j)+dBr
		       do m=-1,1
			  laplBr(j+m)=laplBr(j+m)+sigma(m)*dBr
		       enddo
		    endif
		 endif
 !              If j is the last day		 
		 if(j.eq.ndayt-1)then		    
		    dLL3=alfa(j-1)*((laplBr(j-1)+dBr)**2-
     &	         	 (laplBr(j-1))**2)
		    dLL3=dll3+alfa(j)*((laplBr(j)-dBr)**2-
     &		         (laplBr(j))**2)
		    LLp=0d0
		    do i=j+1,ndayt
		       LLp=LLp+n(i)*dlog(S(i)+muj(i)+dBr*n(j)*(f(i-j)-f(i-j-1)))
		    enddo
		    do i=ini,j
		       LLp=LLp+n(i)*dlog(S(i)+muj(i))
		    enddo
		    LLp=LLp-mu-LL2-dBr*n(j)*f(ndayt-j)-(LL3+dLL3) !trial LL
		    if(LLp.gt.LL)then !acceptance of the trial
		       do i=j+1,ndayt
			  S(i)=S(i)+dBr*n(j)*(f(i-j)-f(i-j-1))
		       enddo
		       LL2=LL2+dBr*n(j)*f(ndayt-j)
		       LL3=LL3+dLL3
		       LL=LLp
		       Br(j)=Br(j)+dBr
		       laplBr(j-1)=laplBr(j-1)+dBr
		       laplBr(j)=laplBr(j)-dBr
		       Br(j+1)=Br(j)
		    endif
		 endif
c************************************************************		 
	      enddo   !loop over 20 MCs
c****** Loop on trials for muj(i) ******* 
	      do l=ini,ndayt
		 j=int(ran2(seed1)*(ndayt-ini))+ini
 57		 dmu=.2*sign(1.,ran2(seed1)-.5) !trial dmu
		 if(muj(j)+dmu.le.0)goto 57
		 LLp=LL+n(j)*(dlog(S(j)+muj(j)+dmu)-dlog(S(j)+muj(j)))-dmu
		 if(LLp.ge.LL)then   !acceptance of the trial
		    muj(j)=muj(j)+dmu
		    mu=mu+dmu
		    LL=LLp
		 endif
	      enddo
c**********************************************************	      
	   enddo  !end loop over all Monte Carlo steps

c*********************************************************
c*********************WRITING ****************************
	   
	   open(200,file=foutput)	!output file
	   write(200,*)'# tau=',taup,'a=',ap,'logL=',LL
	   do i=1,ndayt-7
	      write(200,*)i,Br(i),muj(i)
	   enddo
	   call flush(200)
	   close(200)
	   call flush(200)

	STOP
	end

	
C****************************************************************
C******************** SUBROUTINES AND FUNCTIONS *****************
C****************************************************************

	
	FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END




