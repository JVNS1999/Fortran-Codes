!****************************************************
! PROGRAMA PARA DINAMICA COM CAMPO NA REDE QUADRADA
!****************************************************


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

! -----------------------------
! modulo de variaveis
! -----------------------------

MODULE variaveis

        USE integer_kind
        IMPLICIT NONE

INTEGER(K15)::iseed,iseed2

CHARACTER(7)::dado
CHARACTER(70)::arquivo
REAL(4)::tend,tbegin

REAL(8),ALLOCATABLE::rx(:),ry(:),rxv(:),ryv(:),ex(:),ey(:),nx(:),ny(:),r3(:)
REAL(8),ALLOCATABLE::bx(:),by(:),thre(:)
INTEGER(4),ALLOCATABLE::nviz(:),viz(:),s(:),posq1p(:),posq1m(:)
INTEGER(4),ALLOCATABLE::posq2p(:),posq2m(:)
LOGICAL,ALLOCATABLE::tested(:),mask(:)
LOGICAL::pri

REAL(8)::hc,h0,dh,sigh,thetah,hx,hy,dx,dy,dist,esc,esc1,esc2
REAL(8)::mx,my,nss1,nvs1,alcover,ms(100,30),qs(100,30)
INTEGER(4)::i,j,k,cont,cont1,cont2,ac,L,nspins,nvert,nh
INTEGER(4)::q1,q2,q1p,q1m,q2p,q2m,rseed,conf,sample

CONTAINS

! Gerador que o julio mandou que o grupo do michael estah usando.
!
! Instrucoes
!
! Escolher duas raizes
!iseed1 = 0 -> 31328
!iseed2 = 0 -> 30081
!
!Chamar a rotina de inicializacao
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
       

subroutine init_random_seed()

integer :: values(1:8), k
integer, dimension(:), allocatable :: seed

call date_and_time(values=values)

call random_seed(size=k)
allocate(seed(1:k))
seed(:) = values(8)
call random_seed(put=seed)

end subroutine init_random_seed



END MODULE variaveis

!====================
! INICIO DO PROGRAMA
!====================

PROGRAM din

USE variaveis
IMPLICIT NONE

CALL dados_ini
CALL rede

CALL varredura

END PROGRAM din

!=======================
! INICIO DAS SUBROTINAS
!=======================

!------------------------------

SUBROUTINE varredura

USE variaveis
IMPLICIT NONE

INTEGER(4)::iloop


!write(arquivo,*)"saida.dat"

pri=.false. 

DO sample=1,100

   CALL conf_ini
   call campo
   hx=0.d0
   
   DO iloop=1,nh
    hx=hx+dh
    CALL varre
    CALL cargas
    ms(sample,iloop)=mx/cont
    qs(sample,iloop)=dble(q1)/cont
   END DO
   print*,sample
   
END DO


open(2,file='saida.dat')
open(3,file='medias.dat')
DO sample=1,100
 write(2,*)sample
 hx=0.d0
 DO iloop=1,nh
  hx=hx+dh
  write(2,*)hx,ms(sample,iloop),qs(sample,iloop)
 END DO
END DO
hx=0.d0
DO iloop=1,nh
  hx=hx+dh
  write(3,*)hx,sum(ms(:,iloop))*0.01d0,sum(qs(:,iloop))*0.01d0
END DO

CLOSE(3)
CLOSE(2)

123 FORMAT(3E15.7)

END SUBROUTINE varredura

!------------------------------

SUBROUTINE dados_ini

USE variaveis
IMPLICIT NONE

INTEGER(4)::nvizdip

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
READ(1,532)dado,conf
print*,dado,conf
READ(1,*)dado
print*,dado
READ(1,235)dado,hc
print*,dado,hc
READ(1,235)dado,h0
print*,dado,h0
READ(1,235)dado,dh
print*,dado,dh
READ(1,235)dado,sigh
print*,dado,sigh


CLOSE(1)
call flush()

nspins=2*(L*L-L)
nvizdip=nspins*nspins-nspins

nss1=1.d0/dble(nspins)
nvs1=1.d0/dble((L-2)*(L-2))

hx=0.d0
nh=NINT(h0/dh)+1
print*,'numero de passos de temperatura',nh

ALLOCATE(rx(nspins),ry(nspins),rxv(L*L),ryv(L*L),ex(nspins),ey(nspins)&
        &,nx(nvizdip),ny(nvizdip),r3(nvizdip),bx(nspins),by(nspins),thre(nspins),&
        &nviz(0:nspins),viz(nvizdip),s(nspins),posq1p(L*L),posq1m(L*L),posq2p(L*L),&
        &posq2m(L*L),tested(nspins),mask(nspins))


SELECT CASE(rseed)

CASE (1)

print*,'**************************************************'
print*,'------ ATENCAO!!! USOU SEMENTES ALEATORIAS -------'
print*,'**************************************************'

call init_random_seed()
call random_number(esc)

iseed=INT(esc*31328)
call random_number(esc)
iseed2=INT(esc*30081)

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


!-------------------------------------
!--- Gera uma configuracao inicial ---
!-------------------------------------

SUBROUTINE conf_ini

USE variaveis
IMPLICIT NONE
REAL(8)::zpi
zpi=8.d0*datan(1.d0)


SELECT CASE(conf)

 CASE(1)
  ! GERA UMA CONFIGURACAO ALEATORIA
  DO i=1,nspins
   esc=ranmar()
   IF (esc<0.50d0) then
    S(i)=-1
   ELSE
    S(i)=1
   END IF
  END DO

 CASE(2)
  ! GERA O ESTADO ORDENADO
  S=-1
! distribuição dos campos de reversao
!===================================================================
! MUITA ATENCAO NESSE FATOR MULTIPLICATIVO COLOCADO APENAS NO CAMPO
!   E APENAS NA HORA DE TESTAR SE FLIPA. O CALCULO DO CAMPO NAO FOI
!   ALTERADO!!!
!===================================================================
! Escolhi 0.3 Oe.  Pelas contas que fiz, basicamente, Ms~8 kA/m, 
! vol=2,4X10^(-20) m^3, mu=20x10^(-15) Am^2, a=4x10^(-6) m
! que da para D/mu (contante que multiplica o b) 3x10^(-5) T = 0.3 Oe
! 
! Quando reduzido em 5% vai para 0.285
!
616 alcover=0.3d0     

DO i=1,nspins,2
 IF (ranmar()<0.8) THEN
  hc=147.d0
  sigh=0.075d0*hc
  esc1=ranmar()
  esc2=ranmar()
  esc=sigh*dsqrt(-2.d0*log(esc1))
  thre(i)=hc+esc*cos(zpi*esc2)
  thre(i+1)=hc+esc*sin(zpi*esc2)
  IF (thre(i)<0.d0) goto 616
  IF (thre(i+1)<0.d0) goto 616
 ELSE
  hc=70.d0
  sigh=0.20d0*hc
  esc1=ranmar()
  esc2=ranmar()
  esc=sigh*dsqrt(-2.d0*log(esc1))
  thre(i)=hc+esc*cos(zpi*esc2)
  thre(i+1)=hc+esc*sin(zpi*esc2)
  IF (thre(i)<0.d0) goto 616
  IF (thre(i+1)<0.d0) goto 616
 END IF 
END DO

thre=-thre

 CASE(3)
!  ! GERA O ESTADO FUNDAMENTAL DO MOLLER
!  cont=0
!  DO j=1,L
!   DO i=1,L
!    cont=cont+1
!    if (mod(j,2)==0) then
!     sh(cont)=-1
!    else
!     sh(cont)=1
!    end if
!    if (mod(i,2)==0) then
!     sv(cont)=-1
!    else
!     sv(cont)=1
!    end if
!   END DO
!  END DO
!
 CASE(4)
  ! GERA O ESTADO FUNDAMENTAL DO SPIN ICE NORMAL
  s(1)=1
  DO i=2,nspins
   s(i)=-s(i-1)
  END DO

 CASE(5)
  ! GERA O ESTADO DE MAIOR ENERGIA (SO CARGA DUPLA)
!  sh(0)=1
!  sv(0)=1
!  cont=0
!  DO j=1,L
!   DO i=1,L
!    cont=cont+1
!    sh(cont)=-sh(cont-1)
!    sv(cont)=-sv(cont-1)
!    IF (i==1) THEN
!     sh(cont)=sh(cont-1)
!     sv(cont)=sv(cont-1)
!    END IF
!  END DO
! END DO
!
END SELECT


RETURN
END SUBROUTINE conf_ini

!----------------------------------
!--- Monta a rede e os vizinhos ---
!----------------------------------

SUBROUTINE rede

USE variaveis
IMPLICIT NONE

! Define as posicoes dos vertices
cont=0
DO j=1,L
 DO i=1,L
  cont=cont+1
  rxv(cont)=DBLE(i)
  ryv(cont)=DBLE(j)
 END DO
END DO

! Define as posicoes dos spins 
! na primeira linha, temos que
! spins de 1 a L-1 sao horizontais
! spins de L a 2L-1 sao verticais
! o raciocinio se repete nas demais linhas
! os vetores ex e ey fornecem a direcao dos spins

! total de vertices na rede LxL
! total de spins 2LL -2L
! numero de vertices completos (com os 4 spins) - (L-2)x(L-2)
! indice do vertice na coluna i, linha j - (L*(j-1))+i
! indice spin vert baixo (i,j) - (j-1)*(2*L-1)+i-L
! indice spin vert cima (i,j) - j*(2*L-1)+i-L
! indice spin hor esq (i,j) - (j-1)*(2L-1)+i-1
! indice spin hor dir (i,j) - (j-1)*(2L-1)+i


cont=0
DO j=1,L-1
 DO i=1,L-1 
  cont=cont+1
  rx(cont)=DBLE(i)+0.50d0
  ry(cont)=DBLE(j)
  ex(cont)=1.d0
  ey(cont)=0.d0
 END DO
 DO i=1,L
  cont=cont+1
  rx(cont)=DBLE(i)
  ry(cont)=DBLE(j)+0.50d0
  ex(cont)=0.d0
  ey(cont)=1.d0
 END DO
END DO
j=L
DO i=1,L-1
 cont=cont+1
 rx(cont)=DBLE(i)+0.50d0
 ry(cont)=DBLE(j)
 ex(cont)=1.d0
 ey(cont)=0.d0
END DO

! vetores para calculo do campo
cont=0
DO i=1,nspins
 DO j=1,i-1
  dx=rx(i)-rx(j)
  dy=ry(i)-ry(j)
  dist=dsqrt(dx*dx+dy*dy)
  cont=cont+1
  viz(cont)=j
  nx(cont)=dx/dist
  ny(cont)=dy/dist
  r3(cont)=1.d0/dist/dist/dist
 END DO
 DO j=i+1,nspins
  dx=rx(i)-rx(j)
  dy=ry(i)-ry(j)
  dist=dsqrt(dx*dx+dy*dy)
  cont=cont+1
  viz(cont)=j
  nx(cont)=dx/dist
  ny(cont)=dy/dist
  r3(cont)=1.d0/dist/dist/dist
 END DO
 nviz(i)=cont
END DO 
nviz(0)=0 

RETURN
END SUBROUTINE rede

!-------------------------
!--- Calcula o campo ---
!-------------------------

SUBROUTINE campo

USE variaveis
IMPLICIT NONE

bx=0.d0
by=0.d0

DO i=1,nspins
 DO j=nviz(i-1)+1,nviz(i)
  esc1 =ex(viz(j))*nx(j)+ey(viz(j))*ny(j)
  esc2 = 3.d0*nx(j)*esc1*r3(j)-ex(viz(j))*r3(j)
  bx(i)=bx(i)+s(viz(j))*esc2
  esc2 = 3.d0*ny(j)*esc1*r3(j)-ey(viz(j))*r3(j)
  by(i)=by(i)+s(viz(j))*esc2
 END DO
END DO

RETURN
END SUBROUTINE campo


!-------------------------
!--- Varre a rede ---
!-------------------------

SUBROUTINE varre

USE variaveis
IMPLICIT NONE

outer: DO   ! testa ate que nao tenha ninguem pra ser flipado

 ac=0
 
 tested=.false.
 do1: DO i=1,INT(0.75*nspins)
  k=INT(ranmar()*nspins)+1                      ! escolhe aleatoriamente
  esc=s(k)*ex(k)*(alcover*bx(k)+hx)+s(k)*ey(k)*(alcover*by(k)+hy)
  tested(k)=.true.
  IF (esc .le. thre(k)) THEN                    ! testa o escolhido
   DO j=nviz(k-1)+1,nviz(k)                     ! atualiza os vizinhos caso flipe
    esc1 =ex(k)*nx(j)+ey(k)*ny(j)
    esc2=3.d0*nx(j)*esc1*r3(j)-ex(k)*r3(j)
    bx(viz(j))=bx(viz(j))-2.d0*s(k)*esc2
    esc2=3.d0*ny(j)*esc1*r3(j)-ey(k)*r3(j)
    by(viz(j))=by(viz(j))-2.d0*s(k)*esc2
   END DO
   s(k)=-s(k)                                   !flipa
   ac=ac+1
   pri=.true.
   exit do1                                     ! se flipou, volta a escolher aleatorio do zero...
  END IF
 END DO do1

 IF (ac==0) THEN                                ! testa o resto da rede se ninguem tiver sido flipado
  DO k=1,nspins
   IF (.not. tested(k)) THEN                          ! nao testa um mesmo spin 2 vezes
    esc=s(k)*ex(k)*(alcover*bx(k)+hx)+s(k)*ey(k)*(alcover*by(k)+hy)
    IF (esc .le. thre(k)) THEN
     DO j=nviz(k-1)+1,nviz(k)
      esc1 =ex(k)*nx(j)+ey(k)*ny(j)
      esc2=3.d0*nx(j)*esc1*r3(j)-ex(k)*r3(j)
      bx(viz(j))=bx(viz(j))-2.d0*s(k)*esc2
      esc2=3.d0*ny(j)*esc1*r3(j)-ey(k)*r3(j)
      by(viz(j))=by(viz(j))-2.d0*s(k)*esc2
     END DO
     s(k)=-s(k)
     ac=ac+1
    END IF
   END IF
  END DO
 END IF
 
 IF (ac==0) exit outer                           ! se ninguem flipou apos testar todos, sai pro proximo campo

END DO outer

RETURN
END SUBROUTINE varre


!-------------------------
!--- Varre a rede ---
!-------------------------
!
!SUBROUTINE varresh
!
!USE variaveis
!IMPLICIT NONE
!
!outer: DO   ! testa ate que nao tenha ninguem pra ser flipado
!
! ac=0
! 
! tested=.false.
! do1: DO i=1,INT(0.75*nspins)
!  k=INT(ranmar()*nspins)+1                      ! escolhe aleatoriamente
!  esc=s(k)*ex(k)*(alcover*bx(k))+s(k)*ey(k)*(alcover*by(k))
!  tested(k)=.true.
!  IF (esc .le. thre(k)) THEN            ! testa o escolhido
!          pause
!   DO j=nviz(k-1)+1,nviz(k)                     ! atualiza os vizinhos caso flipe
!    esc1 =ex(k)*nx(j)+ey(k)*ny(j)
!    esc2=3.d0*nx(j)*esc1*r3(j)-ex(k)*r3(j)
!    bx(viz(j))=bx(viz(j))-2.d0*s(k)*esc2
!    esc2=3.d0*ny(j)*esc1*r3(j)-ey(k)*r3(j)
!    by(viz(j))=by(viz(j))-2.d0*s(k)*esc2
!   END DO
!   s(k)=-s(k)                                   !flipa
!   ac=ac+1
!   pri=.true.
!   exit do1                                     ! se flipou, volta a escolher aleatorio do zero...
!  END IF
! END DO do1
!
! IF (ac==0) THEN                                ! testa o resto da rede se ninguem tiver sido flipado
!  DO k=1,nspins
!   IF (.not. tested(k)) THEN                          ! nao testa um mesmo spin 2 vezes
!    esc=s(k)*ex(k)*(alcover*bx(k))+s(k)*ey(k)*(alcover*by(k))
!    IF (esc .le. thre(k)) THEN
!            pause
!     DO j=nviz(k-1)+1,nviz(k)
!      esc1 =ex(k)*nx(j)+ey(k)*ny(j)
!      esc2=3.d0*nx(j)*esc1*r3(j)-ex(k)*r3(j)
!      bx(viz(j))=bx(viz(j))-2.d0*s(k)*esc2
!      esc2=3.d0*ny(j)*esc1*r3(j)-ey(k)*r3(j)
!      by(viz(j))=by(viz(j))-2.d0*s(k)*esc2
!     END DO
!     s(k)=-s(k)
!     ac=ac+1
!    END IF
!   END IF
!  END DO
! END IF
! 
! IF (ac==0) exit outer                           ! se ninguem flipou apos testar todos, sai pro proximo campo
!
!END DO outer
!
!RETURN
!END SUBROUTINE varresh
!
!-------------------------------------------------
!--- Calcula a quantidade de cargas magneticas ---
!-------------------------------------------------

SUBROUTINE cargas

USE variaveis
IMPLICIT NONE

INTEGER(4)::hd,he,vc,vb,prov_int

q1=0
q2=0
q1p=0
q1m=0
q2p=0
q2m=0
posq1p=0
posq1m=0
posq2p=0
posq2m=0
mx=0.d0

cont=0

DO j=6,L-6      !coluna
 DO i=6,L-6     !linha
  hd= (j-1)*(2*L-1)+i      !horizontal a direita 
  he= (j-1)*(2*L-1)+i-1    !horizontal a esquerda
  vc= j*(2*L-1)+i-L        !vertical acima
  vb= (j-1)*(2*L-1)+i-L    !vertical abaixo
  prov_int=s(hd)+s(vc)+s(vb)+s(he)
  SELECT CASE(prov_int)
   CASE(0)
    IF (s(hd)*s(vc)>0) THEN
     q2=q2+1                            !carga dupla
     IF (s(hd)*s(vc)*s(vb)>0) THEN
      q2m=q2m+1                         !carga dupla negativa				
      posq2m(q2m)=L*(j-1)+i
     ELSE
      q2p=q2p+1                         !carga dupla positiva
      posq2p(q2p)=L*(j-1)+i
     END IF
    END IF
   CASE(2)
    q1=q1+1
    IF (s(hd)*s(vc)>0) THEN
     q1p=q1p+1                          !carga simples positiva
     posq1p(q1p)=L*(j-1)+i
    ELSE
     q1m=q1m+1                          !carga simples negativa
     posq1m(q1m)=L*(j-1)+i
    END IF
   CASE(-2)
    q1=q1+1
    IF (s(hd)*s(vc)>0) THEN
     q1m=q1m+1                          !carga simples negativa
     posq1m(q1m)=L*(j-1)+i
    ELSE
     q1p=q1p+1                          !carga simples positiva
     posq1p(q1p)=L*(j-1)+i
    END IF
  END SELECT
  IF (ex(hd)>0) THEN
    mx=mx+dble(s(hd))
  END IF
  cont=cont+1
 END DO
END DO

mx=mx

RETURN

END SUBROUTINE cargas


!****************************************************
!********** IMPRIME UMA CONFIGURACAO ****************
!****************************************************

!----> ABRE O ARQUIVO COM NOME NA VARIAVEL ARQUIVO

SUBROUTINE figura

USE variaveis

!WRITE(arquivo,*)"configuracoes.xyz"
OPEN(1,file='configuracoes.xyz',POSITION='APPEND')

WRITE(1,*)nspins+q1p+q1m+q2p+q2m
WRITE(1,*)'H=',dsqrt(hx*hx+hy*hy)
DO i=1,nspins
 WRITE(1,100)rx(i)-S(i)*ex(i)*0.4,ry(i)-S(i)*ey(i)*0.4,0.d0,s(i)*ex(i),s(i)*ey(i),0.0d0
END DO
DO i=1,q1p
 write(1,200)rxv(posq1p(i)),ryv(posq1p(i)),0.0d0
END DO
DO i=1,q1m
 write(1,300)rxv(posq1m(i)),ryv(posq1m(i)),0.0d0
END DO
DO i=1,q2p
 write(1,400)rxv(posq2p(i)),ryv(posq2p(i)),0.0d0
END DO
DO i=1,q2m
 write(1,500)rxv(posq2m(i)),ryv(posq2m(i)),0.0d0
END DO

CLOSE(1,status='keep')
100 FORMAT(1x,' Al ',3F20.7,' atom_vector',3F20.7)
200 format(1x,' Cu ',3F20.7)
300 format(1x,' Mg ',3F20.7)
400 format(1x,' Fe ',3F20.7)
500 format(1x,' Mn ',3F20.7)
666 FORMAT(A12,I2,A3,f7.4,A3,I1,A4)





RETURN
END SUBROUTINE figura




