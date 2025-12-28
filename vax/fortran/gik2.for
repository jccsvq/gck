 
	program gik

	character*250 file
	integer lib$init_timer,lib$stat_timer
     1  ,lib$free_timer,status
	type *,char(27)//'[2J'//CHAR(27)//'[2;1H'
	type *,char(27)//'[2j'
	type *,char(27)//'#3'//char(27)//'#4'//char(27)//'#6'
     1  //char(27)//'[5m'//char(27)//'[1mCALCULO DE GI-K'
	type *,char(27)//'[0m'

	status=lib$init_timer(iatime)
	if(.not.status) call lib$signal(%val(status))

	ist=lib$get_command(file,'Fichero de datos:',if)
	if(.not.ist) call lib$signal(%val(ist))
	open(unit=5,file=file(1:if),status='old',err=100)
	goto 103
100	type *,'Error al abrir el fichero'
	stop
103	ntot=0
101	read (5,*,end=102)
	ntot=ntot+1
	goto 101
102	rewind(5)
	type *,'Numero de puntos :',ntot
	ist=lib$get_vm(4*ntot,iax)
	if(.not.ist) call lib$signal(%val(ist))
	ist=lib$get_vm(4*ntot,iay)
	if(.not.ist) call lib$signal(%val(ist))
	ist=lib$get_vm(4*ntot,iaz)
	if(.not.ist) call lib$signal(%val(ist))
	ist=lib$get_vm(4*ntot,ianeig)
	if(.not.ist) call lib$signal(%val(ist))
	ist=lib$get_vm(4*ntot,iadist)
	if(.not.ist) call lib$signal(%val(ist))
	ist=lib$get_vm(4*ntot,iaoct)
	if(.not.ist) call lib$signal(%val(ist))
	call lee_datos(ntot,%val(iax),%val(iay),%val(iaz))
	close(unit=5)



c
c	Todo a una subrutina mogollonica
c
	call megalin(ntot,%val(iax),%val(iay),%val(iaz),
     1  %val(ianeig),%val(iadist),%val(iaoct))


	status=lib$stat_timer(2,item,iatime)
	if(.not.status) call lib$signal(%val(status))
	tiempo=.01*item
	status=lib$free_timer(iatime)
	if(.not.status) call lib$signal(%val(status))
        type *,tiempo,' Segundos de CPU'

	stop
	end

	subroutine lee_datos(ntot,x,y,z)

	real x(ntot),y(ntot),z(ntot)

	do i=1,ntot
	read (5,*) a,z(i),x(i),y(i)
	end do

	return
	end

	subroutine megalin(ntot,x,y,z,neigh,dist,oct)


	
	dimension nmon(3)
	data nmon/1,3,6/
	dimension x(ntot),y(ntot),z(ntot),neigh(ntot),dist(ntot)
	integer oct(ntot),octan
	type *,'Entre el orden k (0,1,2):'
	accept *,nork
	nork=nork+1
1	type *,'Entre numero de vecinos a usar:'
	accept *,nvec
	if(nvec.gt.32) then
	type *,'Debe ser menor o igual a 32'
	goto 1
	end if

	ist=lib$get_vm(4*(nvec+nmon(nork))**2,iamat)
	if(.not.ist) call lib$signal(%val(ist))

	ist=lib$get_vm(4*(nvec+nmon(nork)),iavec)
	if(.not.ist) call lib$signal(%val(ist))

	ist=lib$get_vm(4*(nvec+nmon(nork)),iavec1)
	if(.not.ist) call lib$signal(%val(ist))



	do ix=1,ntot
	xx=x(ix)
	yy=y(ix)

	do i=1,ntot
	dist(i)=sqrt((xx-x(i))**2+(yy-y(i))**2)
	oct(i)=octan(x(i)-xx,Y(i)-yy)
	end do

	call busca(ntot,dist,neigh,nvec+1,oct,ifracaso)
	if (ifracaso.ne.1) then
	call quitame_aca_ese_punto(ntot,dist,neigh,nvec)

	call krigea_gik(ntot,nvec,x,y,z,dist,neigh,zestim,
     1	sigma,%val(iamat),%val(iavec),%val(iavec1),nork,nmon,xx,yy)

	if(zestim.gt.0.) then
	CALL ESCRIBE_GKI(NTOT,NEIGH,NVEC,%VAL(IAVEC),NORK,IX)
	end if
	end if

	end do

	return
	end



c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



	subroutine krigea_gik(ntot,nvec,x,y,z,dist,neigh,zestim,
     1	sigma,a,b,bb,nork,nmon,xx,yy)

	dimension nmon(3)
	dimension x(ntot),y(ntot),z(ntot),neigh(ntot),dist(ntot)
	dimension a(nvec+nmon(nork),nvec+nmon(nork)),b(nvec+nmon(nork))
	dimension bb(nvec+nmon(nork))
	real monomio

c
c	establecemos las ecuaciones de krigea_gikge
c

c	Vector

	do i=1,nvec
	b(i)=estructura_gki(dist(i))
	bb(i)=b(i)
	end do
	do i=1,nmon(nork)
	b(nvec+i)=-monomio(xx,yy,i)
	bb(nvec+i)=b(nvec+i)
	end do

c	diagonal principal

	c0=estructura_gki(0.)
	do i=1,nvec
	a(i,i)=c0
	end do

c	caja nula

	do i=1,nmon(nork)
	do j=1,nmon(nork)
	a(nvec+i,nvec+j)=0.
	end do
	end do


c	parte simetrica 1a. caja

	do i=1,nvec-1
	ii=neigh(i)
	do j=i+1,nvec
	jj=neigh(j)
	a(i,j)=estructura_gki(sqrt((x(ii)-x(jj))**2+(y(ii)-y(jj))**2))
	a(j,i)=a(i,j)
	end do
	end do

c	cajas laterales

	do i=1,nvec
	ii=neigh(i)
	do j=1,nmon(nork)
	a(i,nvec+j)=-monomio(x(ii),y(ii),j)
	a(nvec+j,i)=a(i,nvec+j)
	end do
	end do

c	resolucion de las ecuaciones 

	call simq(a,b,nvec+nmon(nork),k1)
	if(k1.eq.0)then
	zestim=0.
	sigma=0.
	do i=1,nvec
	zestim=zestim+b(i)*z(neigh(i))
	end do
	do i=1,nvec+nmon(nork)
	sigma=sigma+b(i)*bb(i)
	end do
	sigma=sigma-c0
	if(sigma.gt.0.) then
	sigma=sqrt(sigma)
	else
	zestim=-9.e-29
	sigma=9.e-29
	end if
	else
	sigma=-9.e-29
	zestim=-.9e-30
	end if
	return
	end


	REAL FUNCTION estructura_gki(H)

	estructura_gki=H

	return
	end


	SUBROUTINE QUITAME_ACA_ESE_PUNTO(NTOT,DIST,NEIGH,NVEC)

	REAL DIST(NTOT)
	INTEGER NEIGH(NTOT)

	DO I=1,NVEC
	NEIGH(I)=NEIGH(I+1)
	DIST(I)=DIST(I+1)
	END DO

	RETURN
	END




	real function monomio(x,y,i)


	goto (100,200,300,400,500,600),i

100	monomio=1.
	return

200	monomio=x
	return

300	monomio=y
	return

400	monomio=x**2
	return

500	monomio=y**2
	return

600	monomio=x*y
	return

	end


	SUBROUTINE ESCRIBE_GKI(NTOT,NEIGH,NVEC,PONDE,NORK,IX)

	REAL PONDE(NTOT)
	INTEGER NEIGH(NTOT)

C	CALCULAMOS LOS PONDERADORES:

	SQ=1.

	DO I=NVEC+1,2,-1
	PONDE(I)=-PONDE(I-1)
	SQ=SQ+PONDE(I)**2
	END DO

	PONDE(1)=1.
	SQ=SQRT(SQ)

	DO I=1,NVEC+1
	PONDE(I)=PONDE(I)/SQ
	END DO


	WRITE(6,*)NORK,NVEC+1
	WRITE(6,*)IX,(NEIGH(I),I=1,NVEC)
	WRITE(6,*)(PONDE(I),I=1,NVEC+1)

	RETURN
	END
	SUBROUTINE SIMQ(A,B,N,KS)
	DIMENSION A(1),B(1)
	TOL=0.0
	KS=0
	JJ=-N
	DO 65 J=1,N
	JY=J+1
	JJ=JJ+N+1
	BIGA=0
	IT=JJ-J
	DO 30 I=J,N
	IJ=IT+I
	IF(ABS(BIGA)-ABS(A(IJ))) 20,30,30
20	BIGA=A(IJ)
	IMAX=I
30	CONTINUE
	IF(ABS(BIGA)-TOL) 35,35,40
35	KS=1
	RETURN
40	I1=J+N*(J-2)
	IT=IMAX-J
	DO 50 K=J,N
	I1=I1+N
	I2=I1+IT
	SAVE=A(I1)
	A(I1)=A(I2)
	A(I2)=SAVE
50	A(I1)=A(I1)/BIGA
	SAVE=B(IMAX)
	B(IMAX)=B(J)
	B(J)=SAVE/BIGA
	IF(J-N) 55,70,55
55	IQS=N*(J-1)
	DO 65 IX=JY,N
	IXJ=IQS+IX
	IT=J-IX
	DO 60 JX=JY,N
	IXJX=N*(JX-1)+IX
	JJX=IXJX+IT
60	A(IXJX)=A(IXJX)-(A(IXJ)*A(JJX))
65	B(IX)=B(IX)-(B(J)*A(IXJ))
70	NY=N-1
	IT=N*N
	DO 80 J=1,NY
	IA=IT-J
	IB=N-J
	IC=N
	DO 80 K=1,J
	B(IB)=B(IB)-A(IA)*B(IC)
	IA=IA-N
80	IC=IC-1
	RETURN
	END

	SUBROUTINE BUSCA (ntot,dist,neigh,nvec,oct,ifrac)
C
C Esta rutina ordena el vector dist de menor a mayor
C
	REAL        dist(ntot)
	INTEGER   neigh(ntot),oct(ntot)

	ifrac=0
C
	DO I=1,ntot
	neigh(I) = I
	END DO
C
	DO I=1,ntot-1
	DO J=I+1,ntot
	IF (dist(I).LE.dist(J)) GOTO 1
	VV    = dist(I)
	dist(I)  = dist(J)
	dist(J)  = VV
	K     = neigh(I)
	neigh(I) = neigh(J)
	neigh(J) = K
	K     = oct(I)
	oct(I) = oct(J)
	oct(J) = K
  1	CONTINUE
	END DO
	END DO
C
	call octavino(ntot,dist,neigh,nvec,oct,ifrac)
	RETURN
	END

	subroutine octavino(ntot,dist,neigh,nvec,oct,ifrac)
	dimension n(8),na(7,8),b(7,8),dist(ntot),neigh(ntot)
	integer oct(ntot)
	do i=1,8
	n(i)=0
	end do
	do i=1,7
	do j=1,8
	na(i,j)=0
	b(i,j)=0.
	end do
	end do
	do i=1,ntot
	if(n(oct(i)).lt.7) then
	n(oct(i))=n(oct(i))+1
	na(n(oct(i)),oct(i))=neigh(i)
	b(n(oct(i)),oct(i))=dist(i)
	end if
	end do
	nonul=0
	do i=1,8
	if(n(i).ne.0) nonul=nonul+1
	end do
	if(nonul.le.5) then
	ifrac=1
	return
	end if
	nnn=0
	do i=1,7
	do j=1,8
	if(na(i,j).ne.0) then
	nnn=nnn+1
	neigh(nnn)=na(i,j)
	dist(nnn)=b(i,j)
	end if
	end do
	end do
	if(nnn.ge.nvec) then
	ifrac=0
	return
	else
	ifrac=1
	return
	end if
	end

	integer function octan(x,y)

	if(x.lt.0.) then

	 if(y.lt.0.) then
	  if(-x.ge.-y) then
	   octan=5
	  else
	   octan=6
	  end if
	 else
	  if(-x.ge.y) then
	   octan=4
	  else
	   octan=3
	  end if
	 end if

	else

	 if(y.lt.0.) then
	  if(x.ge.-y) then
	   octan=8
	  else
	   octan=7
	  end if
	 else
	  if(x.ge.y) then
	   octan=1
	  else
	   octan=2
	  end if
	 end if

	end if
	return
	end

 
