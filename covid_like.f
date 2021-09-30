ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     EVALUATION OF REPRODUCTION CASE NUMBER
c     from
c     assume generation time distribution is a Gamma function w(z)=(tau**(-a)/Gamma(a))*(xt**(-1)*exp(-xt*1/tau))
c     with two input parameters tau and a
c     read from input file with 2 columns containing day and daily number of infected
c     write output file with 3 columns containing day, reproduction number Rc and daily number of immigrants
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	program rc_evaluation

	implicit none
c*********************************************************
c     VARIABLE DECLARATION
c*******************************************************
	integer nsum,sigma(-1:1),ifinal
	integer i,j,k,l,m,imin,ini,ndayt,kloop
	real Br(-100:5000),laplBr(-100:5000)
	real muj(-100:5000),a,ap,xx,yp,y
	real*8  LL,LLp,LL1,LL2,LL3,dLL3,S(0:5000),mu,dmu,dBr
	real*8 f(0:5000),g,dtau,xt,gnorm,taup
	real n(-100:5000),alfa(0:5000)
	character finput*40,foutput*40
c*********************************************************

cccccccccccccccccccccccccccccccccccc
c*****   READ INPUT PARAMETERS **********
        
	print*,'Please, enter the name of the input file'
        print*,'File input format: day, number of infected'
	read(*,'(a)')finput 
	print*,'Please, enter the name of the output file'
	read(*,'(a)')foutput
	print*,'Please, enter the Gamma distribution parameter tau'
        print*,'Suggested value tau \in (0.05,0.5)'
	read(*,*)taup 
	print*,'Please, enter the Gamma distribution parameter a'
        print*,'Suggested value a \simeq 6/taup'
	read(*,*)ap
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
c*****  READING INPUT FILE *******
	
	open(55,file=finput)
	do i=1,10000
	   read(55,*,end=99)xx,n(i) !days and n(i)=daily number of infected
	   nsum=nsum+n(i)
	enddo
 99	ndayt=i-1 !total number of days
	ini=2     !Initial day for LL evaluation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
c*****  Initial setting	**********

	sigma(-1)=1
	sigma(0)=-2
	sigma(1)=1
 
	do i=imin,ndayt+1
	   laplBr(i)=0  !
	   alfa(i)=0.3*sqrt(nsum*1.)  !alfa=1/V in Eq.3 form ref. Lippiello et. al
	enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       NUMERICAL INTEGRAL OF GAMMA DISTRIBUTION  cc
	
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
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*****  Initial evaluation for the first 5 days of the case reproduction number Rc ****
c*****  by means of the Walling-Teunis algorithm ****
c*****  Rc is indicated by Br(i) at day i ****	
	
	do i=2,5
	   Br(i)=0  
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
c******* Initial setting of Br(i) for the subsequent days>5 ****
c******* Initial setting of immigrants rate muj *****
	
	do i=5,ndayt+1
	   Br(i)=Br(2)-0.05*(i-2)
	   Br(i)=max(Br(i),0.01)
	   muj(i)=0
	enddo
	
	do i=ini,ndayt
	   laplBr(i)=(Br(i-1)+Br(i+1)-2*Br(i))
	   mu=mu+muj(i)      
	enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
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
	      j=int(rand()*(ndayt-2))+2   !chose a random day
 41		 dBr=max(Br(j),0.2)*0.02*(rand()-.5) !trial dBr 
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
		 j=int(rand()*(ndayt-ini))+ini
 57		 dmu=.2*sign(1.,rand()-.5) !trial dmu
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
c*********************WRITING ***************************
c       output file: first column=day, second column=Rc, third column=\mu
	   
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


