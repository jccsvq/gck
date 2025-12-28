 
	program geko

	character*250 file
	integer lib$init_timer,lib$stat_timer
     1  ,lib$free_timer,status
	status=lib$init_timer(iatime)
	if(.not.status) call lib$signal(%val(status))
	type *,char(27)//'[2J'//CHAR(27)//'[2;1H'
	type *,char(27)//'[2j'
	type *,char(27)//'#3'//char(27)//'#4'//char(27)//'#6'
     1  //char(27)//'[5m'//char(27)//'[1mESTIMACION COVARIANZA
     1 GENERALIZADA'
	type *,char(27)//'[0m'


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
	call lee_datos(ntot,%val(iax),%val(iay),%val(iaz))
	close(unit=5)


	ist=lib$get_command(file,'Fichero de GIK:',if)
	if(.not.ist) call lib$signal(%val(ist))
	open(unit=5,file=file(1:if),status='old',err=200)
	goto 203
200	type *,'Error al abrir el fichero'
	stop
203	nlam=0
	nvec1=0
201	read (5,*,end=202)nork,nvec
	read (5,*,end=202) (n,i=1,nvec)
	read (5,*,end=202) (a,i=1,nvec)
	nvec=max(nvec,nvec1)
	nvec1=nvec
	nlam=nlam+1
	goto 201
202	rewind(5)
	type *,'Numero de GIK :',nlam
	type *,'Numero max de vecinos :',nvec
	type *,'El orden es k=',nork-1
	ist=lib$get_vm(4*nlam*nvec,ialam)
	if(.not.ist) call lib$signal(%val(ist))
	ist=lib$get_vm(4*nlam*nvec,ianeig)
	if(.not.ist) call lib$signal(%val(ist))
	ist=lib$get_vm(4*nlam,ianneig)
	if(.not.ist) call lib$signal(%val(ist))
	call lee_gik(nlam,nvec,%val(ialam),%val(ianeig),%val(ianneig))
	close(unit=5)

c
c	Todo a una subrutina mogollonica
c
	call pandemonium(ntot,nlam,nvec,%val(iax),%val(iay),%val(iaz),
     1  %val(ialam),%val(ianeig),%val(ianneig),nork-1)


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

	subroutine lee_gik(nlam,nvec,lambda,neigh,nneigh)

	real lambda(nlam,nvec)
	integer neigh(nlam,nvec),nneigh(nlam)

	do i=1,nlam

	read(5,*) nork,nneigh(i)
	read(5,*) (neigh(i,j),j=1,nneigh(i))
	read(5,*) (lambda(i,j),j=1,nneigh(i))

	end do

	return
	end


	subroutine pandemonium(ntot,nlam,nvec,x,y,z,
     1  lambda,neigh,nneigh,nork)


	real x(ntot),y(ntot),z(ntot),lambda(nlam,nvec)
	real qq(5,5),aa(5)
	integer neigh(nlam,nvec),nneigh(nlam)

	nmod=4

	ist=lib$get_vm(4*nlam,iayy)
	if(.not.ist) call lib$signal(%val(ist))
	ist=lib$get_vm(4*nlam,iaw)
	if(.not.ist) call lib$signal(%val(ist))
	ist=lib$get_vm(4*nlam*(nmod+1),iaxx)
	if(.not.ist) call lib$signal(%val(ist))

	call matri1(ntot,nlam,nvec,nmod,x,y,z,lambda,neigh,nneigh,
     1  %val(iayy),%val(iaxx),%val(iaw))

	call matri2(ntot,nlam,nvec,nmod,
     1  %val(iayy),%val(iaxx),%val(iaw),qq,aa)

c	ist=lib$free_vm(4*nlam,iayy)
c	if(.not.ist) call lib$signal(%val(ist))
c	ist=lib$free_vm(4*nlam,iaw)
c	if(.not.ist) call lib$signal(%val(ist))
c	ist=lib$free_vm(4*nlam*(nmod+1),iaxx)
c	if(.not.ist) call lib$signal(%val(ist))

	call modelitos(qq,aa,nork,%val(iayy),%val(iaxx),%val(iaw),
     1  nlam,nmod)

	 return
	 end

	subroutine modelitos(qq,aa,nork,yy,xx,w,nlam,nmod)

	real qq(5,5),aa(5),q(5,5),a(5)
	real xx(nlam,nmod+1),yy(nlam),w(nlam)
	integer iq(5)
	character*250 file

	q0=0.
	z0=0.
	do i=1,nlam
	z0=z0+yy(i)
	q0=q0+w(i)*(yy(i))**2
	end do

	ist=lib$get_command(file,'Fichero de modelos >',if)
	if(.not.ist) call lib$signal(%val(ist))
	open(unit=5,file=file(1:if),status='old',err=100)
	goto 101
100	type *,'Imposible abrir fichero de modelos'
	stop
101	read(5,*,end=300) (iq(i),i=1,5)
	do i=1,5
	a(i)=aa(i)*iq(i)
	do j=1,5
	q(i,j)=qq(i,j)*iq(i)*iq(j)
	q(j,i)=q(i,j)
	end do
	if(iq(i).eq.0) q(i,i)=1
	end do
	call simq(q,a,5,ks)
	if(ks.eq.0) then
	if(((a(5).eq.0.).and.(a(1).ge.0.).and.(a(2).le.0.).and.(a(4).le.0.)
     1  .and.(a(3).ge.(-3.33333*sqrt(a(2)*a(4)) )) ).or.
     2  ((a(4).eq.0.).and.(a(1).ge.0.).and.(a(2).le.0.).and.(a(3).ge.0.)
     3  .and.(a(5).ge.(-1.5*sqrt(-a(3)*a(2)))))) then
	zb=0.
	qb=0.
	do i=1,nlam
	qb=qb+w(i)*(yy(i)-a(1)*xx(i,1)-a(2)*xx(i,2)
     1  -a(3)*xx(i,3)-a(4)*xx(i,4)-a(5)*xx(i,5))**2
	zb=zb+a(1)*xx(i,1)+a(2)*xx(i,2)+a(3)*xx(i,3)
     1  +a(4)*xx(i,4)+a(5)*xx(i,5)
	end do
	ratio =qb/q0
	ro=zb/z0
	write(6,*) '(QS/Q0)= ',ratio,'ro= ',ro
	write(6,*) (a(i),i=1,5)
	end if
	else
	write(6,*) 'Mal condicionado.'
	end if
	goto 101
300	close(unit=5)
	return
	end

	subroutine matri2(ntot,nlam,nvec,nmod,
     1  yy,xx,w,qq,aa)

	real xx(nlam,nmod+1),yy(nlam),w(nlam),qq(5,5),aa(5)

	do i=1,5
	aa(i)=0.
	do j=1,nlam
	aa(i)=aa(i)+yy(j)*xx(j,i)*w(j)
	end do
	end do

	do i=1,5
	do j=i,5
	qq(i,j)=0.	
	do k=1,nlam
	qq(i,j)=qq(i,j)+w(k)*xx(k,i)*xx(k,j)
	end do
	if(i.ne.j) qq(j,i)=qq(i,j)
	end do
	end do

	return
	end
	 
	subroutine matri1(ntot,nlam,nvec,nmod,x,y,z,
     1  lambda,neigh,nneigh,yy,xx,w)


	real x(ntot),y(ntot),z(ntot),lambda(nlam,nvec)
	integer neigh(nlam,nvec),nneigh(nlam)
	real xx(nlam,nmod+1),yy(nlam),w(nlam)

	do i=1,nlam
	yy(i)=0.
	do j=1,nneigh(i)
	yy(i)=yy(i)+lambda(i,j)*z(neigh(i,j))
	end do
	yy(i)=yy(i)**2
	end do

	do i=1,nlam
	xx(i,1)=1
	end do

	do i=1,nlam
	do j=2,nmod+1
	xx(i,j)=0.
	do k1=1,nneigh(i)
	do k2=1,nneigh(i)
	h=sqrt((x(neigh(i,k1))-x(neigh(i,k2)))**2
     1  +(y(neigh(i,k1))-y(neigh(i,k2)))**2)
	xx(i,j)=xx(i,j)+lambda(i,k1)*lambda(i,k2)*covmodel(h,j-1)
	end do
	end do
	end do
	end do

	do i=1,nlam
	w(i)=1./(xx(i,2)**2)
	end do

	return
	end


	real function covmodel(h,i)

	goto (100,200,300,400),i

100	covmodel=h
	return

200	covmodel=h**3
	return

300	covmodel=h**5
	return

400	if(h.le.0.) then
	covmodel=0.
	else
	covmodel=(h**2)*log(h)
	end if
	return
	end
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
 
