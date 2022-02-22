! infinite_model.f90
! Runs the infinite genes DH model on two islands (Manzo & Peliti)
! First island has M1 individuals
! Second island has M2 individuals
! Use M2 >= M1 always.
! use a seed file to run the random number generator

! Random choice of second parent
! Compute species before migration

! Only eps1 is given, since eps2 = M1*eps1/M2

! Compute q1, q2 and p

! at t=0 the overlap matrix is q(i,j) = 1

! only non-diagonal elements are saved in vector x(k)
! position of q(i,j) goes into x(k) with k=(2*nc-i)*(i-1)/2 + (j-i)
! Careful: j > i always

! Outputs:
! number.dat - time, sp.total, sp.island1, sp.island2
! numbertot.dat -  3 lines for each time step:
!                  time, abundances of species in island1
!                  time, abundances of species in island2
!                  time, abundances of all species in the system

! main program

PROGRAM dh2d
IMPLICIT REAL*8(A-H,O-Z)
INTEGER, ALLOCATABLE :: order(:)
INTEGER, ALLOCATABLE :: ispv(:),ispv1(:),ispv2(:),ispecies(:,:),ispecies1(:,:),ispecies2(:,:)
INTEGER, ALLOCATABLE :: list1(:),list2(:)
INTEGER, ALLOCATABLE :: fat(:),mot(:),v(:),vk1(:),vk2(:)
INTEGER :: iseed(33),deltat
REAL*8, ALLOCATABLE :: nspstotal(:),nspilha1(:),nspilha2(:)
REAL*8, ALLOCATABLE :: x(:),x1(:),x2(:),xl(:)
REAL*8 :: nu1,nu2,mut,qmin,mutexp


eps1 = 0.01   ! migration probability
qmin = 0.9    ! mating restriction
ntime = 2000  ! total time
m1 = 200      ! pop island 1
m2 = 200      ! pop island 2
mut = 0.001   ! mutation rate per locus

! parameters for averages
ncut = ntime/2
nstep = 50
nsample = ncut/nstep
ALLOCATE(nspstotal(nsample),nspilha1(nsample),nspilha2(nsample))

! internal parameters
deltat = 5
maxmate1 = m1
maxmate2 = m2


! dynamical parameters
ntot = m1 + m2
am1 = dfloat(m1)
am2 = dfloat(m2)
antot = am1 + am2
r = am1/am2
eps2 = eps1*am1/am2
mutexp = 0.25D0*exp(-4.0D0*mut)
nu1 = 4.0D0*am1*mut
nu2 = 4.0D0*am2*mut
sigma = eps1*am1
aux1 = nu1+2.0d0*sigma+1.d0
aux2 = nu2+2.0d0*sigma+1.d0
cx1 = nu1+(1.d0+r)*sigma
dux2 = aux1/aux2
deno = aux1*cx1-2.0d0*sigma**2*(dux2+r)

qbar1 = cx1/deno
qbar2 = cx1*dux2/deno
pbar = (r+dux2)*sigma/deno
tau = 0.0d0
m1t = m1
m2t = m2
m1ta = int(m1*1.2)  ! allocate vectors with 20% tolerance for fluctuations
m2ta = int(m2*1.2)
ms1t = m1t*(m1t-1)/2
ms2t = m2t*(m2t-1)/2
ms1ta = m1ta*(m1ta-1)/2
ms2ta = m2ta*(m2ta-1)/2
mstot = ntot*(ntot-1)/2

! allocate vectors
ALLOCATE (x(mstot),xl(mstot),x1(ms1ta),x2(ms2ta),order(ntot))
ALLOCATE (ispv(ntot),ispv1(m1),ispv2(m2))
ALLOCATE (ispecies(ntot,ntot),ispecies1(m1ta,m1ta),ispecies2(m2ta,m2ta))
ALLOCATE (fat(ntot),mot(ntot),v(ntot),vk1(ntot),vk2(ntot))
ALLOCATE (list1(m1),list2(m2))
x = 1.0D0

! initialize random number generator
OPEN(UNIT=50,FILE='seed.in',STATUS='OLD')
READ(50,*) iseed
CLOSE(50)
CALL RANDOM_SEED(put=iseed)


OPEN(UNIT=15,FILE="number.dat",STATUS='unknown')

OPEN(UNIT=17,FILE="numbertot.dat",STATUS='unknown')


WRITE(6,202) 'qmin = ',qmin,'qbar = ',qbar,'tau = ',int(tau),'total number of individuals = ',m1+m2
WRITE(6,201) 'total time = ',ntime
WRITE(6,*)
201 FORMAT(A,I6,5x,A,I5)
202 FORMAT(A,F6.4,5x,A,F6.4,5x,A,I5,5x,A,I6)


! Time evolution: mating and mutation in each deme + migration
DO jtime=1,ntime
    DO i=1,m1t
        vk1(i) = i
    END DO
    DO i=1,m2t
        vk2(i) = m1t + i
    END DO
    !Mating - Deme 1
    k = 0
    DO WHILE(k < m1)                  ! if random individual kmother cannot find
        CALL RANDOM_NUMBER(aux)       ! a compatible mate after maxmate trials, a
        kmother = int(aux*m1t)+1      ! systematic search is performed.
        icount = 0
        imate = 0
        DO WHILE (imate == 0 .and. icount < maxmate1)
            DO
                CALL RANDOM_NUMBER(aux)
                kmate = int(aux*m1t)+1  ! get a kmate /= kmother
                IF(kmate /= kmother) EXIT
            END DO
            icount = icount + 1
            ip = kmother
            jp = kmate
            IF(ip > jp) THEN
                iaux = ip
                ip = jp
                jp = iaux
            END IF
            kp = (2*ntot-ip)*(ip-1)/2 + (jp-ip)
            IF(x(kp) >= qmin) THEN
                k = k + 1
                imate = 1
                mot(k) = kmother
                fat(k) = kmate
            END IF
        END DO
    END DO

    !Mating - Deme 2
    k = m1
    DO WHILE(k < ntot)                  ! if random individual kmother cannot find
        CALL RANDOM_NUMBER(aux)         ! a compatible mate after maxmate trials, a
        kmother = int(aux*m2t)+1+m1t    ! systematic search is performed.
        icount = 0
        imate = 0
        DO WHILE (imate == 0 .and. icount < maxmate2)
            DO
                CALL RANDOM_NUMBER(aux)
                kmate = int(aux*m2t)+1+m1t  ! get a kmate /= kmother
                IF(kmate /= kmother) EXIT
            END DO
            icount = icount + 1
            ip = kmother
            jp = kmate
            IF(ip > jp) THEN
                iaux = ip
                ip = jp
                jp = iaux
            END IF
            kp = (2*ntot-ip)*(ip-1)/2 + (jp-ip)
            IF(x(kp) >= qmin) THEN
                k = k + 1
                imate = 1
                mot(k) = kmother
                fat(k) = kmate
            END IF
        END DO
    END DO


    ! update overlap matrix before migration
    DO i=1,ntot
        ifat = fat(i)
        imot = mot(i)
            DO j=i+1,ntot
                jfat = fat(j)
                jmot = mot(j)
                k = (2*ntot-i)*(i-1)/2 + (j-i)

                IF(ifat /= jfat) THEN  ! first contribution
                    ip = ifat
                    jp = jfat
                    IF(ip > jp) THEN
                        iaux = ip
                        ip = jp
                        jp = iaux
                    END IF
                    kp = (2*ntot-ip)*(ip-1)/2 + (jp-ip)
                    xl(k) = mutexp*x(kp)
                ELSE
                    xl(k) = mutexp
                END IF

                IF(ifat /= jmot) THEN  ! second contribution
                    ip = ifat
                    jp = jmot
                    IF(ip > jp) THEN
                        iaux = ip
                        ip = jp
                        jp = iaux
                    END IF
                    kp = (2*ntot-ip)*(ip-1)/2 + (jp-ip)
                    xl(k) = xl(k) + mutexp*x(kp)
                ELSE
                    xl(k) = xl(k) + mutexp
                END IF

                IF(imot /= jfat) THEN  ! third contribution
                    ip = imot
                    jp = jfat
                    IF(ip > jp) THEN
                        iaux = ip
                        ip = jp
                        jp = iaux
                    END IF
                    kp = (2*ntot-ip)*(ip-1)/2 + (jp-ip)
                    xl(k) = xl(k) + mutexp*x(kp)
                ELSE
                    xl(k) = xl(k) + mutexp
                END IF

                IF(imot /= jmot) THEN  ! forth contribution
                    ip = imot
                    jp = jmot
                    IF(ip > jp) THEN
                        iaux = ip
                        ip = jp
                        jp = iaux
                    END IF
                    kp = (2*ntot-ip)*(ip-1)/2 + (jp-ip)
                    xl(k) = xl(k) + mutexp*x(kp)
                ELSE
                    xl(k) = xl(k) + mutexp
                END IF
        END DO
    END DO

    ! update matrices before migration to calculate species
    x = xl

    ! for deme 1
    x1 = 0.0
    DO i=1,m1
        DO j=i+1,m1
            k = (2*ntot-i)*(i-1)/2 + (j-i)
            k1 = (2*m1-i)*(i-1)/2 + (j-i)
            x1(k1) = x(k)
        END DO
    END DO

    ! for deme 2
    x2 = 0.0
    DO i=1,m2
        ii = i + m1
        DO j=i+1,m2
            jj = j + m1
            kk = (2*ntot-ii)*(ii-1)/2 + (jj-ii)
            k2 = (2*m2-i)*(i-1)/2 + (j-i)
            x2(k2) = x(kk)
        END DO
    END DO

    ! search for species every deltat generations
    IF (MOD(jtime,deltat) == 0) THEN
        CALL FINDSPECIES(x1,ispv1,ispecies1,igt1,m1,ms1ta,qmin)
        CALL FINDSPECIES(x2,ispv2,ispecies2,igt2,m2,ms2ta,qmin)
        CALL FINDSPECIES(x,ispv,ispecies,igt,ntot,mstot,qmin)
        N1 = 0
        N2 = 0 
        DO i=1,igt
             n1aux = 0
             n2aux = 0
             DO j=1,ispv(i)
                 if(ispecies(i,j) <= m1) n1aux = 1
                 if(ispecies(i,j) > m1) n2aux = 1  
             END DO
             if(n1aux == 1) N1 = N1 + 1
             if(n2aux == 1) N2 = N2 + 1
        END DO
        
        write(6,*) jtime,igt,igt1,igt2,N1,N2
        write(15,*) jtime,igt,igt1,igt2
        write(17,906) jtime,(ispv(k),k=1,igt)
        write(17,906) jtime,(ispv1(k),k=1,igt1)
        write(17,906) jtime,(ispv2(k),k=1,igt2)
    END IF

    IF(jtime > ncut) THEN
        IF(mod(jtime,nstep) == 0) THEN
            jc = (jtime-ncut)/nstep
            nspstotal(jc) = igt
            nspilha1(jc) = igt1
            nspilha2(jc) = igt2
        END IF
    END IF
 
    ! migration: migrants from 1->2 might be different than from to 2->1
    mig12 = 0
    DO k=1,m1  ! list of migrants from 1 -> 2
        CALL RANDOM_NUMBER(aux)
        IF(aux < eps1) THEN
            mig12 = mig12 + 1
            list1(mig12) = k
        END IF
    END DO

    mig21 = 0
    DO k=m1+1,ntot  ! list of migrants from 2 -> 1
        CALL RANDOM_NUMBER(aux)
        IF(aux < eps2) THEN
            mig21 = mig21 + 1
            list2(mig21) = k
        END IF
    END DO

    ! order migrants
    l=0
    ! first island
    DO i=1,m1   ! first the non-migrants
        inlist = 0
        DO j=1,mig12
            IF(i == list1(j)) THEN
                inlist = 1
                EXIT
            END IF
        END DO
        IF(inlist == 0) THEN
            l = l + 1
            order(l) = i
        END IF
    END DO
    DO i=1,mig21  ! now the migrants from 2 -> 1
        l = l + 1
        order(l) = list2(i)
    END DO
    m1t = l
    ms1t = m1t*(m1t-1)/2

    ! second island
    DO i=m1+1,ntot   ! first the non-migrants
        inlist = 0
        DO j=1,mig21
            IF(i == list2(j)) THEN
                inlist = 1
                EXIT
            END IF
        END DO
        IF(inlist == 0) THEN
            l = l + 1
            order(l) = i
        END IF
    END DO
    DO i=1,mig12  ! now the migrants from 1 -> 2
        l = l + 1
        order(l) = list1(i)
    END DO
    m2t = l - m1t
    ms2t = m2t*(m2t-1)/2
    
!    WRITE(6,*) m1t,m2t,m1t+m2t,l
!    PAUSE
!    WRITE(6,*) (list2(i),i=1,mig21)
!    PAUSE
!    WRITE(6,*) order
!    PAUSE

    ! exchange migrants
    DO i=1,ntot
        ii = order(i)
        DO j=i+1,ntot
            jj = order(j)
            k = (2*ntot-i)*(i-1)/2 + (j-i)
            ip = ii
            jp = jj
            IF(ip > jp) THEN
                iaux = ip
                ip = jp
                jp = iaux
            END IF
            kk = (2*ntot-ip)*(ip-1)/2 + (jp-ip)
            x(k) = xl(kk)
        END DO
    END DO
    
    ! for deme 1
    x1 = 0.0
    DO i=1,m1t
        DO j=i+1,m1t
            k = (2*ntot-i)*(i-1)/2 + (j-i)
            k1 = (2*m1t-i)*(i-1)/2 + (j-i)
            x1(k1) = x(k)
        END DO
    END DO

    ! for deme 2
    x2 = 0.0
    DO i=1,m2t
        ii = i + m1t
        DO j=i+1,m2t
            jj = j + m1t
            kk = (2*ntot-ii)*(ii-1)/2 + (jj-ii)
            k2 = (2*m2t-i)*(i-1)/2 + (j-i)
            x2(k2) = x(kk)
        END DO
    END DO

    write(6,*) jtime


END DO
CLOSE(17)
CLOSE(15)


CALL RANDOM_SEED(get=iseed)
OPEN(UNIT=50,FILE='seed.in',STATUS='OLD', POSITION='REWIND')
WRITE (50,*) iseed
close(50)



CLOSE(15)
CLOSE(17)

906 FORMAT(100(1x,i5))

END PROGRAM dh2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find species                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDSPECIES(x,ispv,ispecies,igt,nc,no,qc)
IMPLICIT REAL*8(A-H,O-Z)
INTEGER, ALLOCATABLE :: species(:),pop1(:),pop2(:)
INTEGER :: ispecies(nc,nc),ispv(nc)
REAL*8 :: x(no),qc

ALLOCATE (species(nc),pop1(nc),pop2(nc))

itot = 0  ! count total population
igt = 0   ! count number of species
i2 = nc
DO i=1,i2  ! initialize pop2 with the list of all individuals
    pop2(i) = i
END DO

DO WHILE (itot < nc)

    icr = pop2(1)   !take first individual and find its species
    ispold = 1
    i1 = 0
    pop1 = 0
    isp = 1
    species(1) = icr
    DO i=2,i2
        ind2 = pop2(i)
        ki = ind2
        kj = icr
        IF (ki > kj) THEN
            kj = ind2
            ki = icr
        END IF
        k=(2*nc-ki)*(ki-1)/2 + (kj-ki)
        IF(x(k) < qc) THEN
            i1 = i1 + 1
            pop1(i1) = ind2      !put creatures with similarity < qc into pop1
        ELSE
            isp = isp + 1
            species(isp) = ind2   !collect individuals with similarity >= qc from icr
        END IF
    END DO

    !check if individuals in aux1 have to be included; put the rest in pop2
    itest = 1
    DO WHILE(itest /= 0)
        i2 = 0
        pop2 = 0
        itest = 0
        isp0 = isp
        loop1: DO i=1,i1
            ind1 = pop1(i)
            jtest = 0
            DO j=ispold+1,isp0
                indsp = species(j)
                ki = ind1
                kj = indsp
                IF (ki > kj) THEN
                    kj = ind1
                    ki = indsp
                END IF
                k=(2*nc-ki)*(ki-1)/2 + (kj-ki)
                IF (x(k) >= qc) THEN
                    isp = isp + 1
                    species(isp) = pop1(i)   !collect individuals with similarity >= qc from indsp
                    itest = 1
                    jtest = 1
                    CYCLE loop1
                END IF
            END DO
            IF (jtest == 0) THEN
                i2 = i2 + 1
                pop2(i2) = pop1(i)      !put creatures with similarity < qc into pop2
            END IF
        END DO loop1
        pop1 = pop2   ! pop1 contains the individuals that are part of the species
        i1 = i2
        ispold = isp0
    END DO

    itot = itot + isp    !total number of individuals classified into species
    igt = igt + 1        !number of species
    !WRITE(6,*) igt,isp

    ! save species info
    DO i=1,isp
        ispecies(igt,i) = species(i)
    END DO
    ispv(igt) = isp          ! number of individuals in species

END DO

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Order                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE RANDVEC(n,vk,v)
INTEGER v(n),vk(n)

v = vk

! Randomize vector v(i)

CALL RANDOM_NUMBER(aux)
DO i=1,2*n
    CALL RANDOM_NUMBER(aux)
    i1 = int(aux*n) + 1
    CALL RANDOM_NUMBER(aux)
    i2 = int(aux*n) + 1
    k = v(i1)
    v(i1) = v(i2)
    v(i2) = k
END DO

RETURN
END
