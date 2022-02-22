! Two island model Temporal evolution-this program generates data for Figs 4c, 5 and S5. It also generates data for the movies 
!Diversity patterns and speciation processes ina two-island system with continuous migration
!Debora Princepe, Simone Czarnobai, Thiago M. Pradella, Rodrigo A. Caetano,Flavia M. D. Marquitti, Marcus A.M. de Aguiar, and Sabrina B. L. Araujo
program lake
	implicit none 
	integer n, num, nger, m, idum, i, j, nfilhos, NG
	integer acada,medesp,mesp,ti,tf
	real*8 c,q, tmut, deltac, cmin, cmax, conect
	parameter (nger=500, m=1000, q=0.05d0, tmut=0.001d0, c=0.02d0)! the parameter m means genome size and c migration rate 
	parameter (num=400,acada=1, ti=300,tf=500)
	integer contl1,contl2,ne,ns, l,ll
	integer pai, mae, dif, nachou, ntent, intq,m1, m2, nl, auxiliar(m+2)
	integer  contmae, contpai, ncontesp,ncontesp1,ncontesp2, npais1,npais2,npais12 , pai0,novo , lagoa, cont !SA
	real*8 ran2, ram, mut
	character*6 nomeB, nomec, nomer
	Integer, ALLOCATABLE ::  P(:,:), Pn(:,:),Pmeio(:,:)
	intq=nint(m*q) 
	nl=m+2
	idum=1  ! seed for random number
	write(nomeB,'(I6)') m
	OPEN(UNIT=60,FILE='Temporal.txt',STATUS='UNKNOWN')
	write(60,*) 'pop ', 'B ','mig ','generation', 'esp ', 'abund1 ',  'abund2 '
	!close(60)
	OPEN(UNIT=65,FILE='Network.txt',STATUS='UNKNOWN')
	write(65,*) 'generation ', 'ind1 ', 'ind2 ', 'H/B'
	OPEN(UNIT=66,FILE='Ind_information.txt',STATUS='UNKNOWN')
	write(66,*) 'generation ', 'ind ', 'esp ', 'island ', 'migrate_yes/no'
	!		Close(65)
	ALLOCATE( P(num, m+2), Pn(num, m+2), Pmeio(num/2,m+2))
			medesp=0
			
			P=0
			Pn=0
			! initial condition
			P(1:num/2,m+2)=1 
			P(num/2+1:num,m+2)=2 
			P(:,m+1)=1
			ns=1
			mesp=1
			do i=1, nger	
				! REPRODUCTION
					Pn=0
					contmae=0
					nfilhos=0 
					npais1=count(P(:,m+2).eq.1)
					npais2=count(P(:,m+2).eq.2)
					do while (nfilhos<num)
						if (nfilhos<num/2) then
							npais12=npais1
							pai0=1	
						else
							npais12=npais2
							pai0=npais1+1
						endif
						mae=pai0+Int(ran2(idum)*(npais12)) 
						nachou=0 
						ntent=npais12
						do contpai=1, ntent	
								pai=pai0+Int(ran2(idum)*(npais12)) 
								if (mae.ne.pai) then
									dif=0 
									do j=1, m  
										dif=dif+abs(P(mae,j)-P(pai,j))
									enddo
									if (dif<=intq) then
										nachou=1 
										nfilhos=nfilhos+1
										do j=1, m 
											if(ran2(idum)<=0.5) then
												Pn(nfilhos,j)=P(mae,j) ! inherits the mother's locus
												else
												Pn(nfilhos,j)=P(pai,j) !inherits the father's locus
											endif
											mut=(ran2(idum))
											if (mut<=tmut) Pn(nfilhos,j)= abs(Pn(nfilhos,j)-1)!mutation
										enddo
										Pn(nfilhos,m+2)= P(mae, m+2)!same island as parents
										Pn(nfilhos,m+1)= P(mae,m+1)!species of the parents
									endif
								endif
								if (nachou==1) exit
						enddo
					enddo
					P=Pn
					if(i.ge.ti.and.i.le.tf) then ! selecting the data to the outpt: from time ti to tf
						if(mod(i,acada).eq.0) then
							call species (P,num,nl,intq,ncontesp,mesp)!call the subroutine that identify species
							!OPEN(UNIT=60,FILE='TemporalB'//trim(adjustl(nomeB))//'.txt',STATUS='old',Access = 'append')
							!OPEN(UNIT=65,FILE='Network_B'//trim(adjustl(nomeB))//'.txt',STATUS='old',Access = 'append')
							Do l=1, num
								ne=P(l,m+1)
								Do ll=(l+1), num
									if(P(ll,m+1).eq.ne) then
										 dif=0 
										do j=1, m  
											dif=dif+abs(P(l,j)-P(ll,j)) 
										enddo
										if(dif.le.intq) then
											write(65,*) i,l,ll, Real(dif*100)/Real(m)
										endif
									endif
								Enddo
							Enddo
							Do ne=1, mesp
							contl1=0
							contl2=0
								Do j=1, num
									if(P(j,m+2).eq.1.and.P(j,m+1).eq.ne) then
										contl1=contl1+1
									endif
									if(P(j,m+2).eq.2.and.P(j,m+1).eq.ne) then
										contl2=contl2+1
									endif
								Enddo
								if(contl1.ne.0.or.contl2.ne.0)then
								write(60,*) num,m,c,i, ne,contl1,contl2
								endif
							Enddo	
							!close(60)
						endif
					endif	
					!MIGRATION
					do j=1, num
						conect=ran2(idum)
						novo=0
						if (conect<c) then
							if (P(j,m+2).eq.1) then
								novo=2
							else
								novo=1
							endif	
							P(j,m+2)=novo
						endif
					if(i.ge.ti.and.i.le.tf) then
						if(mod(i,acada).eq.0) then
						write(66,*) i, j,P(j,m+1), P(j,m+2), novo
						endif
					endif	
					enddo
					!ordenando
					!ordenamento dos parasitos
					Pn=0
					cont=1
					Do lagoa=1, 2
						Do j=1, num
							if (P(j,m+2).eq.lagoa)then
							Pn(cont,:)=P(j,:)
							cont=cont+1
							endif
						enddo
					enddo	
					P=Pn
				enddo
			close(60)
			close(65)
			close(66)
100 format(f8.6,1X,I4,1X,f8.2)
	end
			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!))!!!!!
	subroutine species(MaV,np,nl,nG,ncontesp,mesp) !(P,num,nl,intq,ncontesp,mesp)
	IMPLICIT REAL*8 (A-H,O-Z)
	Dimension::MaV(np,nl)
	Integer, ALLOCATABLE :: NES(:),Nes2(:),Nfreq(:,:),Nlo(:),Nv(:)
	Allocate(Nlo(np), Nv(np), NES(mesp))
	nfim=0
	nGes=NG
	ngrupo=0
	Nv=0
	Nlo=0
	Do while(nfim.eq.0)
		nfim=1
		Do l=1,np 
			if(Nv(l).eq.0)then 
				nfim=0
				ngrupo=ngrupo+1
				Nv(l)=ngrupo ! identifying groups of individuals from the same species
				go to 22	! to verify which are the individuals that will belong to the same species
			endif
		Enddo
22  	continue
		nbus=0
		Do while (nbus.eq.0) 
			nbus=1
			Do i=1,np
			if(Nv(i).eq.ngrupo.and.Nlo(i).eq.0) then	
					Nlo(i)=1 
					Do ii=1,np
						if(Nv(ii).eq.0) then ! if the individual does not yet have a group
							ndiff=0
							Do j=1,nl-2 !genoma
							 n1=(MaV(i,j))
							 n2=(MaV(ii,j))
							 ndiff=ndiff+abs(n1-n2)
							Enddo
							if(ndiff.le.NG)then !if the diff is less than G they are the same species
								nbus=0
								Nv(ii)=ngrupo
							else
								if(ndiff.le.NGes) then
									if(MaV(i,nl-1).eq.MaV(ii,nl-1)) then
										nbus=0
										Nv(ii)=ngrupo
									endif
								endif
							endif
						endif
					enddo
				endif
			Enddo
		Enddo
	Enddo
	
	!!!!!!!!!!
	ALLOCATE(Nfreq(ngrupo,2))
	Nes=0
	nmax=mesp
	ncontesp=0 
	Do n=1,nmax
		Nfreq=0
		Do i=1,np
			mva=(MaV(i,nl-1))! name of the original species
			if(mva.eq.n)then 
				Nfreq(nv(i),2)=Nfreq(nv(i),2)+1	
				Nfreq(nv(i),1)=nv(i)	
			endif
		Enddo
		call sortc(Nfreq,ngrupo)! identify the most abundant species
		ngesco=Nfreq(1,1) !most abundant
		if(Nes(n).eq.0)then
			Do i=1,np
				IF(Nv(i).eq.ngesco)then
					MaV(i,nl-1)=n ! same name as original species (to save the colors in the plot)
					Nes(n)=1 
				endif
			Enddo
		endif
		if(Nes(n).gt.0)	ncontesp=ncontesp+1 
		if(ngrupo.gt.1) then
			do nn=2,ngrupo
					if(Nfreq(nn,2).ge.1)then !SI 2020 1 indivíduo define uma espécie
						mesp=mesp+1
						ncontesp=ncontesp+1 !SA
						Do i=1,np
							if(Nv(i).eq.Nfreq(nn,1))then
								MaV(i,nl-1)=mesp
							endif
						enddo
					Endif
			Enddo
		endif
	enddo
endsubroutine species
!!!!!!!!!!!!!!!!!
subroutine sortc(a,n)
    implicit none
    INTEGER :: n, i, j, x, x2, a
    Dimension :: a(n,2)
    DO i = 2, n
        x = a(i,2)
		x2 = a(i,1)
        j = i - 1
        DO WHILE (j.ge.1)
            if (a(j,2).ge.x) exit
            a(j + 1, 2) = a(j,2)
			a(j + 1, 1) = a(j,1)
			j = j - 1
        end do
        a(j + 1, 2) = x
		a(j + 1, 1) = x2
	ENDDO	
	end subroutine sortc

!##############################################
FUNCTION ran2(idum)
IMPLICIT NONE
INTEGER, PARAMETER :: im1=2147483563,im2=2147483399,imm1=im1-1,&
     ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,&
     ir2=3791,ntab=32,ndiv=1+imm1/ntab
DOUBLE PRECISION , PARAMETER ::   am=1.d0/im1,eps=1.d-14,rnmx=1.d0-eps
DOUBLE PRECISION :: ran2
INTEGER, DIMENSION(ntab) :: iv
INTEGER :: idum,idum2,j,k,iy

save iv,iy,idum2
data idum2/123456789/,iv /ntab*0/,iy /0/
      
if(idum.le.0) then
  idum=max(-idum,1)
  idum2=idum
  do j=ntab+8,1,-1
    k=idum/iq1
    idum=ia1*(idum-k*iq1)-ir1*k
    if(idum.lt.0) idum=idum+im1
    if(j.le.ntab) iv(j)=idum
   end do
   iy=iv(1)
endif

k=idum/iq1
idum=ia1*(idum-k*iq1)-ir1*k
if(idum.lt.0) idum=idum+im1
k=idum2/iq2
idum2=ia2*(idum2-k*iq2)-ir2*k

if(idum2.lt.0) idum2=idum2+im2

j=1+iy/ndiv
iy=iv(j)-idum2
iv(j)=idum
if (iy.lt.1)iy=iy+imm1
ran2=min(am*iy,rnmx)

END FUNCTION
