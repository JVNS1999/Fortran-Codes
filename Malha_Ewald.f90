!Este programa visa estudar as propriedades físicas 
!de materiais magnéticos na forma de malha, onde cada
!elemento da malha corresponde a uma nanoilha alongada.

!As "ilhas" interagem entre si via energia dipolar e
!implementamos a soma de Ewald com intuito de calcular
!a energia para um sistema submetido a condições de 
!contorno periódicas.







!########DEFINIÇÃO DE VARIAVÉIS
MODULE Variaveis 
implicit none
  integer(4), parameter:: L=70				!tamanho da rede
  integer(4):: s(1:L*L), ex(1:L*L), ey(1:L*L)		!sigma, componentes da direção do spin
  real(8):: rx(1:L*L), ry(1:L*L)			!componentes da posição de cada sítio
  real(8), parameter:: Ti=1.81d0, Tf=0.01d0, dT=-0.1d0!temperatura inicial, final e passo
  real(8):: T!, dT					!temperatura
  real(8):: M, Mx, My					!magnetização total do sistema
  real(8):: E, E_real, E_fourier, E_self, E1_fourier	!Energia total, da parte real, recíproca e auto-energia
  real(8):: dist					!distância entre sítios
  real(8):: x, y					!variaveis aleatorias
  integer(4), parameter:: MC=20e3, terma_step=10e3	!passo máximo de monte carlo, passo da termalização do sistema
  integer(4), parameter::num_max=10, config_ini=1	!número máximo de amostras, configuração inicial do sistema
  integer(4):: n					!contador do número de amostras
  integer(4), allocatable::viz(:), Nviz(:)		!vizinhos e número de vizinhos
  real(8), allocatable:: A(:), B(:)			!upgrade da energia
  real(8), parameter:: rc=18.0d0				!raio de corte
  integer(4):: nx, ny					!componentes das direções das cópias do sistema
  real(8), allocatable:: dx(:), dy(:)			!componentes da distância entre sítios
  integer(4), parameter:: nc=10				!número de cópias para cada direção
  real(8), parameter:: Pi=4.0d0*datan(1.0d0)		!constante pi
  real(8), parameter:: alpha=0.11d0				!parâmetro de convergência
  integer(4), parameter::nfour=5				!raio de corte no espaço de fourier
  real(8), allocatable:: B_real(:), C_real(:)		!coeficientes da energia no espaço real
  real(8), allocatable:: Dij(:)			!upgrade da energia no espaço real
  real(8), allocatable:: Gx(:), Gy(:)			!vetor de onda no espaço recíproco
  integer(4):: Gmax					!vetor de onda máximo
  real(8), allocatable:: seno(:), coseno(:)		!vetores seno e cosseno
  real(8), allocatable:: h1(:), Fs(:), Fc(:)		!coeficientes da energia no espaço de fourier
  real(8), allocatable:: F1s(:), F1c(:)		!coeficientes auxiliares para o cálculo da variação de energia
ENDMODULE




!########ESCOPO PRINCIPAL DO PROGRAMA
PROGRAM Principal

use Variaveis
implicit none
  character(50)::arquivo1, arquivo2, arquivo3, arquivo4	!índice do aquivo
  integer(4):: step						!contadores
  real(4):: t_end, t_begin					!instantes de tempo inicial e final do programa
 
 
 
200 FORMAT(A15,I2,A4)
300 FORMAT(A20,I2,A4)
400 FORMAT(A15, F4.2, A4)

    call CPU_TIME(t_begin)
    call Rede
    call Coeficientes
    !call Upgrade
    
    !#####################################################
!     call TESTANDO
!     call CPU_TIME(t_end)
!   
!     Print*, (t_end-t_begin), 'tempo de programa'
    !#####################################################
    
    
    !dT=-0.1d0
    Open(1, file='parametros.dat', form='unformatted')
      write(1) Ti, Tf, dT, terma_step, MC, num_max, L
    Close(1)

    write(arquivo2, 200) 'Energia_x_MC_L=',L,'.dat'
    Open(2, file=arquivo2)

    write(arquivo3, 300)'Magnetizacao_x_MC_L=',L,'.dat'
    Open(3, file=arquivo3)

  !loop de amostras
  do n=1, num_max
    call Rede
    call Energia
    call Magnetizacao

!     call Imagem; 
    !T=Ti; dT=-0.1d0
    !loop de temperatura
    !do while(nint(100.0d0*T)>=nint(100.0d0*Tf)) 
    do T=Ti, Tf, dT

    !condições sobre o dT
!     if(T<=0.9d0 .and. T>0.7d0) then
!       dT=-0.01d0
!     else if(T<=0.7d0 .and. T>0.1d0) then
!       dT=-0.1d0
!     else if(T<=0.1d0) then
!       dT=-0.01d0
!     end if

    write(arquivo4, 400) 'termalizacao_T=',T,'.dat'
    Open(4, file=arquivo4, form='unformatted', position='append')
	
	!loop de passos de monte carlo
	do step=1, MC
	  !call Energia
	  !call Magnetizacao
	  if(nint(100.0d0*T)==nint(100.0d0*Ti)) then
	    write(2,*) step, E
	    write(3,*) step, dsqrt(Mx*Mx+My*My)
	  endif
	  if(step>=terma_step) then
	    write(4) E, dsqrt(Mx*Mx+My*My)
	  endif
	  call Metropolis
	enddo
	!fim do loop de passos de monte carlo

	if(nint(100.0d0*T)==nint(100.0d0*Ti)) then
	    write(2,*)
	    write(3,*)
	endif
 
    Close(4)

    
	if (n==1) then
	  call Video
	endif
	if(nint(100.0d0*T)==nint(100.0d0*Tf)) then
	  call Imagem
	endif
	!T=T+dT
      enddo
      !fim do loop de temperatura
  enddo
  !fim do loop de amostras

    Close(2)
    Close(3)

  call CPU_TIME(t_end)
  Print*, (t_end-t_begin), E/dble(L*L)
  
  call system('./dadosNOrho')
  
ENDPROGRAM



!########LISTA DE VIZINHOS E ETC...
SUBROUTINE Upgrade

  use Variaveis
  implicit none
  real(8):: termo1, termo2
  integer(4):: i, j, cont				!contadores

  
  allocate(Nviz(0:L*L), viz(1:L*L*nint(3.15*rc*rc)), A(1:L*L*nint(3.15*rc*rc)))
  allocate(B(1:L*L))
  Nviz(0)=0; viz=0; A=0.0d0

  !Lista de vizinhos e calculo do vetor A
  cont=0
  do i=1, L*L
    do nx=-nc, nc !direções
    do ny=-nc, nc !das cópias
      do j=1, L*L
	dist=dsqrt((rx(j)-rx(i)+dble(nx*L))*(rx(j)-rx(i)+dble(nx*L)) + (ry(j)-ry(i)+dble(ny*L))*(ry(j)-ry(i)+dble(ny*L)))
	if(dist<rc .and. dist>0.01d0) then
	  cont=cont+1
	  viz(cont)=j

	  termo1= dble(ex(i))*dble(ex(j)) + dble(ey(i))*dble(ey(j))
	  termo2=( dble(ex(i))*(rx(j)-rx(i)+dble(nx*L)) + dble(ey(i))*(ry(j)-ry(i)+dble(ny*L)) )*( dble(ex(j))*(rx(j)-rx(i)+dble(nx*L)) + dble(ey(j))*(ry(j)-ry(i)+dble(ny*L)) )
	  A(cont)= (termo1/dist**(3.0d0)) - ((3.0d0*termo2)/dist**(5.0d0))
	endif
      enddo
    enddo
    enddo
      Nviz(i)=cont
  enddo
ENDSUBROUTINE



!########MONTAGEM DA REDE(MALHA)
SUBROUTINE Rede

use Variaveis
implicit none
  integer(4):: i, j, cont					!contadores

  select case(config_ini)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    !configuração aleatória
    case(1)
      !loop das posições e dos spins
      cont=0
      do i=0, L-1
	do j=0, L-1

	  cont=cont+1
	  rx(cont)=1.0d0*i
	  ry(cont)=1.0d0*j

	  call random_number(x)
	  if(mod(i+j,2)==0) then	!ilhas verticais
	    s(cont) =-1 + 2*nint(x)
	    ex(cont)= 1
	    ey(cont)= 0
	  else				!ilhas horizontais
	    s(cont) =-1 + 2*nint(x)
	    ex(cont)= 0
	    ey(cont)= 1
	  endif

	enddo
      enddo
      
!-------------------------------------------------------------------------------------------  

    !configuração de vórtice (estado fundamental)
    case(2)
      !loop das posições e dos spins
      cont=0
      do i=0, L-1
	do j=0, L-1

	  cont=cont+1
	  rx(cont)=1.0d0*i
	  ry(cont)=1.0d0*j

	  if(mod(i+j,2)==0) then	!ilhas verticais
	    if( ry(cont)>dble(L-1)/2.0d0 ) then
	      s(cont)=1
	    else 
	      s(cont)=-1
	    endif
	    ex(cont)=1
	    ey(cont)=0
	  else				!ilhas horizontais
	    if( rx(cont)>dble(L-1)/2.0d0 ) then
	      s(cont)=-1
	    else
	      s(cont)=1
	    endif
	    ex(cont)=0
	    ey(cont)=1
	  endif

	enddo
      enddo
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  end select

ENDSUBROUTINE



!########CALCULO DOS COEFICIENTES DA ENERGIA TOTAL
SUBROUTINE Coeficientes

use variaveis
implicit none

  integer(4):: i, j 					!contadores da posição dos sítios da rede
  integer(4):: cont					!contador
  integer(4):: kx, ky					!contadores do espaço de fourier
  real(8):: termo1, termo2				!termos dos coeficientes reais
  real(8):: arg, G					!termos dos coeficientes do espaço de fourier
  
  
  Gmax=(2*nfour+1)*(2*nfour+1)-1
  
  allocate(Nviz(0:L*L), viz(1:L*L*nint(3.15*rc*rc)))
  allocate(B_real(1:L*L*nint(3.15*rc*rc)), C_real(1:L*L*nint(3.15*rc*rc)), dx(1:L*L*nint(3.15*rc*rc)), dy(1:L*L*nint(3.15*rc*rc)))
  allocate(h1(1:Gmax), Gx(1:Gmax), Gy(1:Gmax), seno(1:Gmax*L*L), coseno(1:Gmax*L*L), F1s(1:Gmax), F1c(1:Gmax))
  allocate(Dij(1:L*L), Fs(1:Gmax), Fc(1:Gmax))
  
  Nviz(0)=0; viz=0; dx=0.0d0; dy=0.0d0
  B_real=0.0d0; C_real=0.d0
  h1=0.0d0; Gx=0.0d0; Gy=0.0d0; seno=0.0d0; coseno=0.0d0
  F1s=0.0d0; F1c=0.0d0
  
  
  
  cont=0
  !coeficientes da parte real
  do i=1, L*L
    do nx=-nc, nc
    do ny=-nc, nc
      do j=1, L*L
	dist=dsqrt((rx(j)-rx(i)+dble(nx*L))**(2.0d0) + (ry(j)-ry(i)+dble(ny*L))**(2.0d0))
	if(dist>0.1d0 .and. dist<rc) then
	  cont=cont+1
	  viz(cont)=j
	  dx(cont)=rx(j)-rx(i)+dble(nx*L)
	  dy(cont)=ry(j)-ry(i)+dble(ny*L)
	  termo1=erfc(alpha*dist)
	  termo2=(2.0d0*alpha/dsqrt(Pi))*dexp(-(alpha*dist)**(2.0d0))
	  B_real(cont)=-termo1/dist**(3.0d0) - termo2/dist**(2.0d0)
	  C_real(cont)=3.0d0*termo1/dist**(5.0d0) + termo2*(3.0d0/(dist*dist) + 2.0d0*alpha*alpha)/dist**(2.0d0)
	endif
      enddo
    enddo
    enddo
	  Nviz(i)=cont
  enddo

  
  cont=0
  !coeficientes da parte recíproca (espaço de fourier)
    !loop dos vetores de onda do espaço recíproco
    do kx=-nfour, nfour
      do ky=-nfour, nfour
	if(kx/=0 .or. ky/=0) then
	  cont=cont+1
	  Gx(cont)=2.0d0*Pi*dble(kx)/dble(L)
	  Gy(cont)=2.0d0*Pi*dble(ky)/dble(L)
	endif
      enddo
    enddo

    cont=0
    !loop do coeficiente h1 e do argumento do seno e coseno
    do i=1, Gmax
      G=dsqrt(Gx(i)*Gx(i)+Gy(i)*Gy(i))
      h1(i)=erfc(0.5d0*G/alpha)/G
      do j=1, L*L
	cont=cont+1
	arg=Gx(i)*rx(j)+Gy(i)*ry(j)
	seno(cont)=dsin(arg)
	coseno(cont)=dcos(arg)
      enddo
    enddo
  
ENDSUBROUTINE



!########CALCULO DA ENERGIA DIPOLAR
SUBROUTINE Energia

use Variaveis
implicit none
  integer(4):: i, j, cont				!contadores
  real(8):: termo1, termo2, termo3, termo4
  
  
  E=0.0d0; E_real=0.0d0;  cont=0; E_fourier=0.0d0; Fs=0.0d0; Fc=0.0d0; Dij=0.0d0;
  
  !Cálculo da auto-energia
  E_self=-2.0d0*alpha*alpha*alpha*dble(L*L)/(3.0d0*dsqrt(Pi))
    
    
  !Calculo da parte real da energia
  do i=1, L*L
    do j=Nviz(i-1)+1, Nviz(i)
      termo1=dble(ex(i))*dble(ex(viz(j))) + dble(ey(i))*dble(ey(viz(j)))
      termo2=dble(ex(i))*dx(j) + dble(ey(i))*dy(j)
      termo3=dble(ex(viz(j)))*dx(j) + dble(ey(viz(j)))*dy(j)
      Dij(i)=Dij(i) + (B_real(j)*termo1 + C_real(j)*termo2*termo3)*dble(s(viz(j)))
    enddo
    E_real=E_real+dble(s(i))*Dij(i)
  enddo
    E_real=-0.5d0*E_real
    
    
  !Cálculo da parte recíproca (espaço de Fourier) da energia
  do i=1, Gmax
    do j=1, L*L
      cont=cont+1
      termo4=Gx(i)*dble(s(j)*ex(j)) + Gy(i)*dble(s(j)*ey(j))
      Fs(i)=Fs(i)+termo4*seno(cont)
      Fc(i)=Fc(i)+termo4*coseno(cont)
    enddo
    E_fourier=E_fourier + h1(i)*(Fs(i)*Fs(i)+Fc(i)*Fc(i))
  enddo
    E_fourier=Pi*E_fourier/(L*L)
    E1_fourier=E_fourier		!valor provisório da energia de fourier
    
    
  !Energia total
  E=E_real + E_fourier + E_self
ENDSUBROUTINE



!########CALCULO DA MAGNETIZAÇÃO
SUBROUTINE Magnetizacao

use Variaveis
implicit none
  integer(4):: k					!contador

  Mx=0.0d0
  My=0.0d0
  M=0.0d0
  do k=1, L*L
    Mx=Mx+dble(s(k))*dble(ex(k))
    My=My+dble(s(k))*dble(ey(k))
  enddo
  M=dsqrt(Mx*Mx+My*My)
ENDSUBROUTINE



!########CÁLCULO DA VARIAÇÃO DE ENERGIA
FUNCTION Delta_E(sitio)

use Variaveis
implicit none
  integer(4), intent(in):: sitio			!fornece a localização de uma nanoilha
  real(8):: Delta_E					!variação da energia total
  real(8):: DE_real					!variação da energia real
  real(8):: termo1
  integer(4):: i, cont					!contadores
  
  E_fourier=0.0d0
  cont=0
  
  !variação da energia real
  DE_real=2.0d0*dble(s(sitio))*Dij(sitio)
  
  !variação da parte recíproca da energia
  do i=1, Gmax
    cont=(i-1)*L*L+sitio
    termo1=Gx(i)*dble(s(sitio)*ex(sitio)) + Gy(i)*dble(s(sitio)*ey(sitio))
    F1s(i)=Fs(i)-termo1*seno(cont)
    F1c(i)=Fc(i)-termo1*coseno(cont)
    
    termo1=Gx(i)*dble(-s(sitio)*ex(sitio)) + Gy(i)*dble(-s(sitio)*ey(sitio))		!termo com o spin invertido
    F1s(i)=F1s(i)+termo1*seno(cont)
    F1c(i)=F1c(i)+termo1*coseno(cont)
    
    E_fourier=E_fourier + h1(i)*( F1s(i)*F1s(i) + F1c(i)*F1c(i) )
  enddo
    E_fourier=Pi*E_fourier/(L*L)
  
  !variação da energia total
  Delta_E=DE_real + E_fourier - E1_fourier
  
RETURN  
ENDFUNCTION



!########ALGORITMO DE METROPOLIS
SUBROUTINE Metropolis

use Variaveis
implicit none
  real(8)::DE, Delta_E, Ef, Ei				!variação da energia total
  integer(4):: p					!sítio aleatorio
  integer(4):: k, i, j					!contador
  real(8):: termo1, termo2, termo3			!constantes

  do k=1, L*L

    !sorteio de sítio aleatorio
    call random_number(x)
    p=nint(1.0d0+(L*L-1)*x)
 
    Ei=E
    !calcular variação de energia
    DE=Delta_E(p)
     
    
    
    if(DE<=0.0d0) then
      
      !atualização da energia e magnetização
      E=E+DE
      Mx=Mx-2.0d0*dble(s(p))*dble(ex(p))
      My=My-2.0d0*dble(s(p))*dble(ey(p))
    
      !atualização dos coeficientes Dij
      do i=Nviz(p-1)+1, Nviz(p)
	termo1=dble(ex(p))*dble(ex(viz(i))) + dble(ey(p))*dble(ey(viz(i)))
	termo2=dble(ex(p))*dx(i) + dble(ey(p))*dy(i)
	termo3=dble(ex(viz(i)))*dx(i) + dble(ey(viz(i)))*dy(i)
	Dij(viz(i))=Dij(viz(i))-2.0d0*(B_real(i)*termo1 + C_real(i)*termo2*termo3)*s(p)
      enddo
      
      !atualizações importantes para o cálculo da energia
      s(p)=-1*s(p)!;pause 1; call Energia; Ef=E; print*, DE - (Ef-Ei)
      E1_fourier=E_fourier
      Fs=F1s
      Fc=F1c
      
    else
    
      call random_number(y)
      if(dexp(-DE/T)>y) then
      
	!atualização da energia e magnetização
	E=E+DE
	Mx=Mx-2.0d0*dble(s(p))*dble(ex(p))
	My=My-2.0d0*dble(s(p))*dble(ey(p))
	
	!atualização dos coeficientes Dij
	do i=Nviz(p-1)+1, Nviz(p)
	  termo1=dble(ex(p))*dble(ex(viz(i))) + dble(ey(p))*dble(ey(viz(i)))
	  termo2=dble(ex(p))*dx(i) + dble(ey(p))*dy(i)
	  termo3=dble(ex(viz(i)))*dx(i) + dble(ey(viz(i)))*dy(i)
	  Dij(viz(i))=Dij(viz(i))-2.0d0*(B_real(i)*termo1 + C_real(i)*termo2*termo3)*s(p)
	enddo
	
	!atualizações importantes para o cálculo da energia
	s(p)=-1*s(p)!;pause 2; call Energia; Ef=E; print*, DE - (Ef-Ei)
	E1_fourier=E_fourier
	Fs=F1s
	Fc=F1c
	
      endif
      
    endif

  enddo

ENDSUBROUTINE



!########VIDEO DA TERMALIZAÇÃO
SUBROUTINE Video

use Variaveis
implicit none
  character(50)::arquivo
  integer(4):: k					!contador

500 FORMAT(A8,I2,A4,F4.2,A4,F4.2,A4)
  !escrevendo a rede
  write(arquivo,500)'video_L=',L,'_Ti=',Ti,'_Tf=',Tf,'.xyz'
  Open(5, file=arquivo, position='append')
    501 FORMAT(1x,'H',3E20.7,'atom_vector',3E20.7)
    write(5,*) L*L
    write(5,*) T
    do k=1, L*L
      write(5,501) rx(k), ry(k), 0.0d0, dble(s(k))*dble(ex(k)), dble(s(k))*dble(ey(k)), 0.0d0
    enddo
  Close(5)

ENDSUBROUTINE



!########IMAGENS DOS "GROUND-STATES"
SUBROUTINE Imagem

use Variaveis
implicit none
  character(50)::arquivo
  integer(4):: k					!contador

600 FORMAT(A10,I2,A4,F4.2,A8,I2,A4)
  !escrevendo a rede
  write(arquivo,600)'Imagens_L=',L,'_Tf=',Tf,'_num_max',num_max,'.xyz'
  Open(6, file=arquivo, position='append')
    501 FORMAT(1x,'H',3E20.7,'atom_vector',3E20.7)
    write(6,*) L*L
    write(6,*) n
    do k=1, L*L
      write(6,501) rx(k), ry(k), 0.0d0, dble(s(k))*dble(ex(k)), dble(s(k))*dble(ey(k)), 0.0d0
    enddo
  Close(6)

ENDSUBROUTINE

!Open(3, file='rede_teste_L10.xyz', position='append')
!Lembrete: O comando position='append' permite continuar escrevendo
!no mesmo arquivo sem altera-lo.



SUBROUTINE Energia_dipolar

use Variaveis
implicit none
  integer(4):: i, j					!contadores

  
  E=0.0d0; B=0.0d0
  do i=1, L*L
    do j=Nviz(i-1)+1, Nviz(i)
      B(i)=B(i)+A(j)*dble(s(viz(j)))
    enddo
    E=E+dble(s(i))*B(i)
  enddo
    E=E*0.5d0

ENDSUBROUTINE


! 
! SUBROUTINE TESTANDO
! 
! use Variaveis
! implicit none
!   integer(4):: i, j, k, q
!   real(4)::z
!   real(4)::t_begin, t_end
!   real(8)::Ei_ex, Ajuste_linear
!   
!   999 FORMAT(5(E15.7, 3X))
!   
!   
!   call CPU_TIME(t_begin)
!   Open(100, file='E_dipolar.dat')
!   Open(200, file='E_ewald.dat')
!   !do i=1, 1	!loop das configurações
!     call Rede
!     do k=10, 40		!loop do raio de corte
!       nc=k; rc=dsqrt(2.0d0)*L*nc
!       call Upgrade
!       call Energia_dipolar
!       deallocate(Nviz, viz, A, B)
!       !Print*, 1.0d0/dble(nc), E
!       write(100,*) 1.0d0/dble(nc), E
!     enddo
!   Close(100)
!       
!       !REGRESSÃO LINEAR
!       !################################################################################################################
!       Ei_ex=Ajuste_linear()
!       
!       !################################################################################################################
!       
!       
!     do j=5, 20	!raio de corte no espaço real
!       do q=4, 12 !raio de corte no espaço de fourier
! 	z=0.05
! 	do while(z<=0.6)
! 	  rc=dble(j); nfour=q; alpha=z
! 	  call Coeficientes
! 	  call Energia
! 	  call CPU_TIME(t_end)
! 	  deallocate(Nviz, viz, B_real, C_real, dx, dy, h1, Gx, Gy, seno, coseno, F1s, F1c, Dij, Fs, Fc)
! 	  if((abs(E-Ei_ex)/Ei_ex)<0.0005) then
! 	  write(200,999) rc, dble(nfour), alpha, E, dble(t_end-t_begin)
! 	  endif
! 	  z=z+0.01
! 	enddo
!       enddo
!     enddo
!       
!   !enddo
!   Close(200)
! 
! ENDSUBROUTINE
! 
! 
! 
! FUNCTION Ajuste_linear
! 
!       real(8)::Ajuste_linear
!       Real(8):: a0, a1 !coeficientes
!       Real(8):: x, x2, y, xy
!       Integer(4):: N !quantidade de dados
!       Real(8), allocatable:: Xi(:), Yi(:)
!       
!       Open(100, file='E_dipolar.dat')      
!       N=31
!       allocate(Xi(N), Yi(N))!Designando o tamanho do vetor
!       !leitura dos dados do vetor
!       Do i=1,N
! 	Read(100,*) Xi(i), Yi(i)
!       End do      
!       Close(100)     
!       x=0.0; y=0.0; x2=0.0; xy=0.0
!       Do i=1,N
! 	x=x+Xi(i)
! 	y=y+Yi(i)
! 	x2=x2+Xi(i)*Xi(i)
! 	xy=xy+Xi(i)*Yi(i)
!       End do
! 	a1=(xy - x*y/N)/(x2 - x*x/N)	!angular
! 	a0=y/N - a1*x/N 		!linear
!       print*, a0, a1
! 
!       Ajuste_linear=a0
!       return 
! 
! 
! ENDFUNCTION
