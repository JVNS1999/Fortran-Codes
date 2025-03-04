!******************************************
! PROGRAMA PARA O MODIFIED SQUARE SPIN ICE
!******************************************

! Usa condicoes de contorno abertas
! Visita os sitios sequencialmente
! delta e = 0.5
! imprime o histograma de energias e da fracao de vertices com carga nao nula.

!+++++++++++++++++++!
MODULE integer_kind
  IMPLICIT NONE

  INTEGER, PARAMETER :: K15=selected_int_kind(15)
  INTEGER, PARAMETER :: K4B=selected_int_kind(9)
END MODULE integer_kind
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!Variables for the pseudo-random numbers generator
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE raset1
  USE integer_kind
  IMPLICIT NONE

  REAL(8), SAVE :: u(97), cc, cd, cm
  INTEGER(K15), SAVE :: i97, j97
END MODULE raset1

!------------------------------------------------------------------------------!
!       MODULO DE VARIAVEIS
!------------------------------------------------------------------------------!

MODULE variaveis

USE integer_kind
IMPLICIT NONE

REAL(8)::t0,dt,tf,beta,t,r,sq2,e,de,prob,qs2,e_fin,qt,ns1,mx,my,m,nss1

REAL(8),ALLOCATABLE::rx(:),ry(:),ex(:),ey(:),aij(:,:),vx(:),vy(:),ener(:)
REAL(8),ALLOCATABLE::ener_n(:),qtot(:),hist_m(:)

REAL(4)::tend,tbegin

INTEGER(4)::i,j,k,L,nmcs,imc,nterm,iterm,it,nt,rseed,nspins,cont,cont2
INTEGER(4)::nvert,ei,ef,seq,L4,eint,nmult
INTEGER(K15)::iseed,iseed2

INTEGER(4),ALLOCATABLE::s(:),q(:),he(:),sn(:)

CHARACTER(7)::dado
CHARACTER(70)::arquivo

CONTAINS

! Gerador que o julio mandou que o grupo do michael estah usando.
!
! Instrucoes
!
! Escolher duas raizes
!iseed1 = 0 -> 31328
!iseed2 = 0 -> 30081
!
!Chamar a rotinha de inicializacao
!CALL rmarin(iseed1,iseed2)
!
!depois para r= 0 -> 1, so chamar
!r = ranmar()
!


!------------------------------------------------------------------------------!
!Pseudorandom numbers generator
!------------------------------------------------------------------------------!
SUBROUTINE rmarin(ij, kl)
!  This subroutine and the next function generate random numbers. See
!  the comments for SA for more information. The only changes from the
!  orginal code is that (1) the test to make sure that RMARIN runs first
!  was taken out since SA assures that this is done (this test didn't
!  compile under IBM's VS Fortran) and (2) typing ivec as integer was
!  taken out since ivec isn't used. With these exceptions, all following
!  lines are original.

! This is the initialization routine for the random number generator
!     RANMAR()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081
USE integer_kind
USE raset1
IMPLICIT NONE
INTEGER(K15), INTENT(IN) :: ij, kl
INTEGER(K15) :: i, j, k, l, ii, jj, m
DOUBLE PRECISION :: s, t

IF( ij < 0  .OR.  ij > 31328  .OR. kl < 0  .OR.  kl > 30081 ) THEN
  WRITE(*, '(A)') ' The first random number seed must have a value ',  &
               'between 0 AND 31328'
  WRITE(*, '(A)') ' The second seed must have a value between 0 and 30081'
  STOP
END IF

i = MOD(ij/177, 177) + 2
j = MOD(ij, 177) + 2
k = MOD(kl/169, 178) + 1
l = MOD(kl, 169)
DO ii = 1, 97
  s = 0.0D0
  t = 0.5D0
  DO jj = 1, 24
    m = MOD(MOD(i*j, 179)*k, 179)
    i = j
    j = k
    k = m
    l = MOD(53*l + 1, 169)
    IF (MOD(l*m, 64) >= 32) THEN
      s = s + t
    END IF
    t = 0.5D0*t
  END DO
  u(ii) = s
END DO
cc = 362436.0D0/16777216.0D0
cd = 7654321.0D0/16777216.0D0
cm = 16777213.0D0/16777216.0D0
i97 = 97
j97 = 33

RETURN
END SUBROUTINE rmarin
!------------------------------------------------------------------------------!
! in Florida State University Report: FSU-SCRI-87-50
!------------------------------------------------------------------------------!
FUNCTION ranmar() RESULT(fn_val)
USE raset1
IMPLICIT NONE
DOUBLE PRECISION :: fn_val
! Local variable
DOUBLE PRECISION:: uni

  uni = u(i97) - u(j97)
  IF( uni < 0.0D0 ) uni = uni + 1.0D0
  u(i97) = uni
  i97 = i97 - 1
  IF(i97 == 0) i97 = 97
  j97 = j97 - 1
  IF(j97 == 0) j97 = 97
  cc = cc - cd
  IF( cc < 0.0D0 ) cc = cc + cm
  uni = uni - cc
  IF( uni < 0.0D0 ) uni = uni + 1.0D0
!  IF( uni == 0.0D0 ) uni = 2.0D-38
  fn_val = uni

  RETURN
END FUNCTION ranmar

!------------------------------------------------------------------------------!
!The same generator as above but for a vector of n random numbers 
!------------------------------------------------------------------------------!
SUBROUTINE vranmar(n, fn_val)
USE integer_kind
USE raset1
IMPLICIT NONE
  INTEGER(K4B), INTENT(IN) :: n
  DOUBLE PRECISION, DIMENSION(1:n),INTENT(OUT) :: fn_val

  DOUBLE PRECISION :: uni
  INTEGER(K4B) :: i

  DO i = 1, n
    uni = u(i97) - u(j97)
    IF( uni < 0.0D0 ) uni = uni + 1.0D0
    u(i97) = uni
    i97 = i97 - 1
    IF(i97 == 0) i97 = 97
    j97 = j97 - 1
    IF(j97 == 0) j97 = 97
    cc = cc - cd
    IF( cc < 0.0D0 ) cc = cc + cm
    uni = uni - cc
    IF( uni < 0.0D0 ) uni = uni + 1.0D0
!    IF( uni == 0.0D0 ) uni = 2.0D-38
    fn_val(i) = uni
  END DO

  RETURN
END SUBROUTINE vranmar
!------------------------------------------------------------------------------!


!************************************************
!********* Gerador de numeros aleatorios ********
!************************************************

FUNCTION ale(idum)
IMPLICIT NONE
INTEGER, PARAMETER :: K4B=selected_int_kind(9)
INTEGER(K4B), INTENT(INOUT) :: idum
REAL :: ale
INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
REAL, SAVE :: am
INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
if (idum <= 0 .or. iy < 0) then 
am=nearest(1.0,-1.0)/IM
iy=ior(ieor(888889999,abs(idum)),1)
ix=ieor(777755555,abs(idum))
idum=abs(idum)+1 
end if
ix=ieor(ix,ishft(ix,13)) 
ix=ieor(ix,ishft(ix,-17))
ix=ieor(ix,ishft(ix,5))
k=iy/IQ 
iy=IA*(iy-k*IQ)-IR*k
if (iy < 0) iy=iy+IM
ale=am*ior(iand(IM,ieor(ix,iy)),1) 
END FUNCTION ale 

subroutine init_random_seed()

integer :: values(1:8), k
integer, dimension(:), allocatable :: seed

call date_and_time(values=values)

call random_seed(size=k)
allocate(seed(1:k))
seed(:) = values(8)
call random_seed(put=seed)

end subroutine init_random_seed


END MODULE VARIAVEIS

!====================
! INICIO DO PROGRAMA
!====================

PROGRAM massi

USE variaveis
IMPLICIT NONE

CALL dados_ini
CALL rede
CALL energia
CALL temperatura

END PROGRAM massi

!=======================
! INICIO DAS SUBROTINAS
!=======================

SUBROUTINE temperatura

USE variaveis
IMPLICIT NONE


CALL print_conf

   DO it=1,nt
   
    t=t+dt
    beta=-1.d0/t
    PRINT*,'temperatura',t
   
    DO iterm=1,nterm
      CALL mcstep_seq
      CALL mcstep_seq
      CALL mcstep_seq
      CALL mcstep_seq
      CALL mcstep_seq
      CALL mcstep_mult
    END DO
   
    print*,'Termalizou!'
   
    call flush()
   
    he=0
    qtot=0.d0
    hist_m=0.d0
    DO imc=1,10
     DO iterm=1,nmcs/10
      CALL mcstep_seq
      CALL mcstep_seq
      CALL mcstep_seq
      CALL mcstep_seq
      CALL mcstep_seq
      CALL mcstep_mult
!histograma
      eint=NINT(2*e)
      he(eint)=he(eint)+1
!carga nos vertices
      DO i=1,nvert-1
        q(i)=S(2*i-1)-S(2*i+1)
      END DO
      DO i=2*L,nvert,2*L
        q(i)=S(2*i-1)-S(2*(i-1)+1)
      END DO
      qt=sum(abs(q))*ns1
      qtot(eint)=qtot(eint)+qt
!magnetizacao
      mx=sum(ex*S)
      my=sum(ey*S)
      m=sqrt(mx*mx+my*my)*nss1
      hist_m(eint)=hist_m(eint)+m
     END DO
     Print*,'Progresso',imc
     call print_conf
     call flush()
    END DO
    
    write(arquivo,"('hist_L=',I2.2,'_t=',F6.4,'.dat')")L,t
    open(1,file=arquivo)
    DO i=ei,ef
     write(1,*)0.5d0*i,he(i),qtot(i),hist_m(i)
    END DO
    close(1)
    call flush()
   
    t0=t
   
   END DO
   
   
   CALL print_conf
   
RETURN
END SUBROUTINE temperatura

!------------------------------

SUBROUTINE dados_ini

USE variaveis
IMPLICIT NONE

call cpu_time(tbegin)

OPEN(1,FILE='dados_ini.dat')

532 FORMAT(A7,I10)
235 FORMAT(A7,E15.7)

READ(1,532)dado,rseed
print*,dado,rseed
READ(1,532)dado,iseed
print*,dado,iseed
READ(1,532)dado,iseed2
print*,dado,iseed2
READ(1,*)dado
print*,dado
READ(1,532)dado,L
print*,dado,L
READ(1,*)dado
print*,dado
READ(1,235)dado,t0
print*,dado,t0
READ(1,235)dado,tf
print*,dado,tf
READ(1,235)dado,dt
print*,dado,dt
READ(1,*)dado
print*,dado
READ(1,532)dado,nterm
print*,dado,nterm
READ(1,532)dado,nmcs
print*,dado,nmcs
READ(1,532)dado,nmult
print*,dado,nmult


CLOSE(1)
call flush()

L4=4*L
nspins=4*L*L
nvert=2*L*L
sq2=dsqrt(2.d0)*0.5d0
ns1=0.5d0/dble(nspins)
nss1=1.d0/dble(nspins)

ei=-3*nspins
ef=nspins

t=t0-dt
nt=NINT(ABS((t-tf)/dt))
print*,'numero de passos de temperatura',nt

SELECT CASE(rseed)

CASE (1)

print*,'**************************************************'
print*,'------ ATENCAO!!! USOU SEMENTES ALEATORIAS -------'
print*,'**************************************************'

call init_random_seed()
call random_number(r)

iseed=INT(r*31328)
call random_number(r)
iseed2=INT(r*30081)

print*,'iseed = ',iseed
print*,'iseed2= ',iseed2

CALL rmarin(iseed,iseed2)


CASE (0)

print*,'------ ATENCAO!!! USOU SEMENTES FORNECIDAS -------'

print*,'iseed = ',iseed
print*,'iseed2= ',iseed2

CALL rmarin(iseed,iseed2)

END SELECT

END SUBROUTINE dados_ini

!------------------------------

SUBROUTINE print_conf

USE variaveis
IMPLICIT NONE

100 FORMAT(1x,' Al ',3F20.7,' atom_vector',3F20.7)
200 format(1x,' Cu ',3F20.7)
300 format(1x,' Mg ',3F20.7)
400 format(1x,' Fe ',3F20.7)

write(arquivo,"('config_L=',I3.3,'_t=',F7.5,'.xyz')")L,t
open(1,file=arquivo,position='append')

write(1,*)nspins+nvert
write(1,*)'L=20',L,'t=',t,'imc=',imc
DO i=1,nspins
 write(1,100)rx(i)-0.5d0*ex(i)*S(i),ry(i)-0.5d0*ey(i)*S(i),0.d0,ex(i)*S(i),ey(i)*S(i),0.d0
END DO
DO i=1,nvert
 IF (q(i)==-2) THEN
  write(1,200)vx(i),vy(i),0.d0
 ELSE IF (q(i)==2) THEN
  write(1,300)vx(i),vy(i),0.d0
 ELSE
  write(1,400)vx(i),vy(i),0.d0
 END IF
END DO

RETURN
END SUBROUTINE print_conf

!------------------------------

SUBROUTINE rede

USE variaveis
IMPLICIT NONE

REAL(8)::r1x,r2x,r3x,r4x,r1y,r2y,r3y,r4y,a
REAL(8)::esc1,esc2,esc3,rijx,rijy,rij,r3,r5
INTEGER(4)::nx,ny

ALLOCATE(S(nspins),rx(nspins),ry(nspins),ex(nspins),ey(nspins),vx(nvert),vy(nvert),&
   &aij(nspins,nspins),ener(nspins),q(0:nvert),he(ei:ef),sn(nspins),ener_n(nspins),&
   &qtot(ei:ef),hist_m(ei:ef))

!tamanho da celula unitaria
a=2.d0*(1.d0+sq2)

!vetores de base da celula unitaria
r1x= -0.5d0*(1.d0+sq2); r1y=1.d0+sq2*0.5d0
r2x=0.d0; r2y=0.d0
r3x=0.5d0*(1.d0+qs2); r3y=1.d0+sq2*0.5d0
r4x=1.d0+sq2; r4y=2.d0+sq2

! vetores da rede
cont=0
cont2=0
DO j=1,L
 DO i=1,L
  cont=cont+1
  rx(cont)=i*a+r1x
  ry(cont)=j*a+r1y
  ex(cont)=sq2
  ey(cont)=-sq2

  cont=cont+1
  rx(cont)=i*a+r2x
  ry(cont)=j*a+r2y
  ex(cont)=1.d0
  ey(cont)=0.d0

  cont=cont+1
  rx(cont)=i*a+r3x
  ry(cont)=j*a+r3y
  ex(cont)=sq2
  ey(cont)=sq2

  cont=cont+1
  rx(cont)=i*a+r4x
  ry(cont)=j*a+r4y
  ex(cont)=1.d0
  ey(cont)=0.d0

  cont2=cont2+1
  vx(cont2)=i*a
  vy(cont2)=j*a+0.5d0

  cont2=cont2+1
  vx(cont2)=i*a+1.d0+sq2
  vy(cont2)=j*a+1.5d0+sq2
 END DO
END DO

! Matriz de interacao PBC 500 copias

aij=0.d0
!percorre nx negativo e todos ny
DO nx=-10,10
 DO ny=-10,10
    DO i=1,nspins
     DO j=1,nspins
      rijx=rx(i)-rx(j)+dble(nx*L*a)
      rijy=ry(i)-ry(j)+dble(ny*L*a)
      rij=dsqrt(rijx*rijx+rijy*rijy)
      r3=rij*rij*rij
      r5=r3*rij*rij
      esc1=ex(i)*ex(j)+ey(i)*ey(j)
      esc2=ex(i)*rijx+ey(i)*rijy
      esc3=ex(j)*rijx+ey(j)*rijy
      aij(i,j)=aij(i,j)+esc1/r3-3.d0*esc2*esc3/r5
     END DO
    END DO
 END DO
END DO

print*,aij

print*,'========================================================='

DO i=1,nspins
 aij(i,i)=0.d0
END DO

print*,aij
pause

! !percorre nx positivo e todos ny
! DO nx=1,10
!  DO ny=-10,10
!     DO i=1,nspins
!      DO j=1,nspins
!       rijx=rx(i)-rx(j)+dble(nx*L*a)
!       rijy=ry(i)-ry(j)+dble(ny*L*a)
!       rij=dsqrt(rijx*rijx+rijy*rijy)
!       r3=rij*rij*rij
!       r5=r3*rij*rij
!       esc1=ex(i)*ex(j)+ey(i)*ey(j)
!       esc2=ex(i)*rijx+ey(i)*rijy
!       esc3=ex(j)*rijx+ey(j)*rijy
!       aij(i,j)=aij(i,j)+esc1/r3-3.d0*esc2*esc3/r5
!      END DO
!     END DO
!  END DO 
! END DO
! !percorre nx=0 e ny negativo
! nx=0
!  DO ny=-10,-1
!     DO i=1,nspins
!      DO j=1,nspins
!       rijx=rx(i)-rx(j)+dble(nx*L*a)
!       rijy=ry(i)-ry(j)+dble(ny*L*a)
!       rij=dsqrt(rijx*rijx+rijy*rijy)
!       r3=rij*rij*rij
!       r5=r3*rij*rij
!       esc1=ex(i)*ex(j)+ey(i)*ey(j)
!       esc2=ex(i)*rijx+ey(i)*rijy
!       esc3=ex(j)*rijx+ey(j)*rijy
!       aij(i,j)=aij(i,j)+esc1/r3-3.d0*esc2*esc3/r5
!      END DO
!     END DO
!  END DO 
! !percorre nx=0 e ny positivo
!  DO ny=1,10
!     DO i=1,nspins
!      DO j=1,nspins
!       rijx=rx(i)-rx(j)+dble(nx*L*a)
!       rijy=ry(i)-ry(j)+dble(ny*L*a)
!       rij=dsqrt(rijx*rijx+rijy*rijy)
!       r3=rij*rij*rij
!       r5=r3*rij*rij
!       esc1=ex(i)*ex(j)+ey(i)*ey(j)
!       esc2=ex(i)*rijx+ey(i)*rijy
!       esc3=ex(j)*rijx+ey(j)*rijy
!       aij(i,j)=aij(i,j)+esc1/r3-3.d0*esc2*esc3/r5
!      END DO
!     END DO
!  END DO 
! !percorre nx e ny iguais a zero, excluindo i=j
!  ny=0
!  DO i=1,nspins
!      DO j=1,i-1
!       rijx=rx(i)-rx(j)
!       rijy=ry(i)-ry(j)
!       rij=dsqrt(rijx*rijx+rijy*rijy)
!       r3=rij*rij*rij
!       r5=r3*rij*rij
!       esc1=ex(i)*ex(j)+ey(i)*ey(j)
!       esc2=ex(i)*rijx+ey(i)*rijy
!       esc3=ex(j)*rijx+ey(j)*rijy
!       aij(i,j)=aij(i,j)+esc1/r3-3.d0*esc2*esc3/r5
!      END DO
!      DO j=i+1,nspins
!       rijx=rx(i)-rx(j)
!       rijy=ry(i)-ry(j)
!       rij=dsqrt(rijx*rijx+rijy*rijy)
!       r3=rij*rij*rij
!       r5=r3*rij*rij
!       esc1=ex(i)*ex(j)+ey(i)*ey(j)
!       esc2=ex(i)*rijx+ey(i)*rijy
!       esc3=ex(j)*rijx+ey(j)*rijy
!       aij(i,j)=aij(i,j)+esc1/r3-3.d0*esc2*esc3/r5
!      END DO
!     END DO
!     
! configuração inicial aleatória
DO i=1,nspins
 IF (ranmar()<0.5d0) THEN
  S(i)=1
 Else
  S(i)=-1
 END IF
END DO

RETURN
END SUBROUTINE rede

!------------------------------

SUBROUTINE energia

USE variaveis
IMPLICIT NONE

ener=0.d0
e=0.d0

DO i=1,nspins
 DO j=1,nspins
  ener(i)=ener(i)+aij(i,j)*S(j)
 END DO
 e=e+ener(i)*S(i)
END DO
e=0.50d0*e

print*,e

RETURN
END SUBROUTINE energia

!------------------------------

SUBROUTINE energia_n

USE variaveis
IMPLICIT NONE

ener_n=0.d0
e_fin=0.d0

DO i=1,nspins
 DO j=1,nspins
  ener_n(i)=ener_n(i)+aij(i,j)*Sn(j)
 END DO
 e_fin=e_fin+ener_n(i)*Sn(i)
END DO
e_fin=0.50d0*e_fin

RETURN
END SUBROUTINE energia_n


!------------------------------

SUBROUTINE mcstep_seq

USE variaveis
IMPLICIT NONE

INTEGER(4)::i1,k1,j1

DO i1=1,nspins

 k1=i1
 de=-2.d0*S(k1)*ener(k1)

 IF (de<0.d0) THEN
  
  e=e+de
  S(k1)=-S(k1)
  DO j1=1,nspins
   ener(j1)=ener(j1)+2.d0*aij(j1,k1)*S(k1)
  END DO

 ELSE

  prob=dexp(beta*de)
  IF (ranmar()<prob) THEN

   e=e+de
   S(k1)=-S(k1)
   DO j1=1,nspins
    ener(j1)=ener(j1)+2.d0*aij(j1,k1)*S(k1)
   END DO

  END IF
 END IF

END DO

print*,e
call energia

pause 

RETURN
END SUBROUTINE mcstep_seq

!------------------------------

SUBROUTINE mcstep_nseq

USE variaveis
IMPLICIT NONE

INTEGER(4)::i1,k1,j1

DO i1=1,nspins

 k1=INT(ranmar()*nspins)+1
 de=-2.d0*S(k1)*ener(k1)

 IF (de<0.d0) THEN
  
  e=e+de
  S(k1)=-S(k1)
  DO j1=1,nspins
   ener(j1)=ener(j1)+2.d0*aij(j1,k1)*S(k1)
  END DO

 ELSE

  prob=dexp(beta*de)
  IF (ranmar()<prob) THEN

   e=e+de
   S(k1)=-S(k1)
   DO j1=1,nspins
    ener(j1)=ener(j1)+2.d0*aij(j1,k1)*S(k1)
   END DO

  END IF
 END IF

END DO

RETURN
END SUBROUTINE mcstep_nseq


! 
! PASSOS MULTIPLOS
!
! FLIP DE FAIXAS
!   para selecionar uma faixa basta pegar um numero ci entre 1 e L. Os spins da faixa estao entre
!    i1=(ci-1)*L*4+1 e i2=ci*L*4
! 
! FLIP DE LINHA
!   linha de base
!       para selecionar uma linha de base basta pegar um numero ci entre 1 e L. Os spins da linha de base
!       estao entre i1=(ci-1)*L*4+2 e i2=ci*L*4-2 de 4 em 4
!   linha de topo
!       para selecionar uma linha de topo basta pegar um numero ci entre 1 e L. Os spins da linha de topo
!       estao entre i1=(ci-1)*L*4+4 e i2=ci*L*4 de 4 em 4
!
! FLIP DE LINHA DIAGONAL
!
! UP
!  CONDICAO DE CONTORNO PERIODICA
!    para selecionar a linha basta selecionar a celula inicial ci=1,2,...,L o numero de spins a ser flipado 
!    e L. O spin inicial e i1=(ci-1)*L*4+3, o passo e ds=(L+1)*4 e se i> i2=L*L*4-1, faz-se i=i-L*L*4
!
!  CONDICAO DE CONTORNO ABERTA
!     Da para usar o mesmo algoritmo. Basta ou começar depois do if ou terminar antes
!
! DOWN
!   CONDICAO DE CONTORNO PERIODICA
!    para selecionar a linha basta selecionar a celula inicial ci=1,2,...,L o numero de spins a ser flipado
!    e L. O spin inicial e i1=(ci-1)*L*4+1, o passo e ds=-(L-1)*4 e se i<0 faz-se i=i+L*L*4.
!
!  CONDICAO DE CONTORNO ABERTA
!     Da para usar o mesmo algoritmo. Basta ou começar depois do if ou terminar antes
!
! LOOPS

!------------------------------

SUBROUTINE mcstep_mult

USE variaveis
IMPLICIT NONE

INTEGER(4)::i1,i2,ci,ds,passo,ii,ili
REAL(8)::e_ini

DO ili=1,nmult
  
  passo=INT(ranmar()*5)+1
  ci=INT(ranmar()*L)+1
  ! INICIA CORRETAMENTE AS VARIAVEIS USADAS
  e_ini=e
  sn=s
  
  SELECT CASE(passo)
  
  CASE(1)     !'faixas'
   
    i1=(ci-1)*L4+1
    i2=ci*L4
    do i=i1,i2
     sn(i)=-sn(i)
    end do
   
  
  CASE(2) ! LINHA DE BASE
  
    i1=(ci-1)*L4+2
    i2=ci*L4-2
    ds=4
    do i=i1,i2,ds
      sn(i)=-sn(i)
    end do
   
  CASE(3) ! LINHA NO TOPO
  
    i1=(ci-1)*L4+4
    i2=ci*L4
    ds=4
    do i=i1,i2,ds
      sn(i)=-sn(i)
    end do
   
  CASE(4) ! Diagonal UP 
   
    i1=(ci-1)*L4+3
    I2=L*L4-1
    ds=L4+4
    i=i1-ds
  ! condicao de contorno aberta. Escolhe entre uma das duas metadades da rede
    IF (ranmar()<0.5d0) THEN
      DO ii=1,L
        i=i+ds
        IF (i>I2) EXIT
        sn(i)=-sn(i)
      END DO
    ELSE
      j=0
      DO ii=1,L
        i=i+ds
        IF (i>I2) THEN
                i=i-L*L4
                j=10
        END IF
        IF (j>5) Sn(i)=-Sn(i)
      END DO
    END IF
  
  
  CASE(5) ! Diagonal DOWN 
   
    i1=(ci-1)*L4+1
    ds=-(L-1)*4
    i=i1-ds
  
    ! condicao de contorno aberta. Escolhe entre uma das duas metadades da rede
    IF (ranmar()<0.5d0) THEN
      DO ii=1,L
        i=i+ds
        IF (i<0) Exit
        sn(i)=-sn(i)
      END DO
    ELSE 
      j=0
      DO ii=1,L
        i=i+ds
        IF (i<0) THEN
                i=i+L*L4
                j=10
        END IF
        IF (j>5) sn(i)=-sn(i)
      END DO
    END IF
  
  
  END SELECT
  
  
  
  call energia_n
    
  de=e_fin-e_ini
  IF (de<0) THEN
        e=e_fin
        S=Sn
        ener=ener_n
  
  ELSE
        prob=dexp(beta*de)
        IF (ranmar()<prob) THEN
              e=e_fin
              S=Sn
              ener=ener_n
        END IF
  END IF
  
END DO
 
RETURN
END SUBROUTINE mcstep_mult

