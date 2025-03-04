!Este programa visa estudar as propriedades físicas 
!de materias magnéticos na forma de malha, onde cada
!elemento da malha corresponde a uma nanoilha alongada.

!As "ilhas" interagem entre si via energia dipolar.







!########DEFINIÇÃO DE VARIAVÉIS
MODULE Variaveis 
implicit none
  integer(4), parameter:: L=40				!tamanho da rede
  integer(4):: s(1:L*L), ex(1:L*L), ey(1:L*L)		!sigma, componentes da direção do spin
  real(8):: rx(1:L*L), ry(1:L*L)			!componentes da posição de cada sítio
  real(8), parameter:: Ti=0.01d0, Tf=0.01d0, dT=-0.02d0!temperatura inicial, final e passo
  real(8):: T!, dT					!temperatura
  real(8):: E, M, Mx, My				!energia e magnetização total do sistema
  real(8):: dist					!distância entre sítios
  real(8):: x, y					!variaveis aleatorias
  integer(4), parameter:: MC=20e3, terma_step=5e3	!passo máximo de monte carlo, passo da termalização do sistema
  integer(4), parameter::num_max=1, config_ini=2	!número máximo de amostras, configuração inicial do sistema
  integer(4):: n					!contador do número de amostras
  integer(4), allocatable::viz(:), Nviz(:)		!vizinhos e número de vizinhos
  real(8), allocatable:: A(:), B(:)			!upgrade da energia, coeficientes da energia
  real(8), parameter:: rc=50d0				!raio de corte
  real(8):: Hx, Hy, dH, Hxmin, Hxmax			!campo magnético externo
ENDMODULE



!########ESCOPO PRINCIPAL DO PROGRAMA
PROGRAM Principal

use Variaveis
implicit none
  character(50)::arquivo1, arquivo2, arquivo3, arquivo4	!índice do aquivo
  integer(4):: step						!contadores
 
200 FORMAT(A15,I2,A4)
300 FORMAT(A20,I2,A4)
400 FORMAT(A12,F4.2,A3,I2,A4)

    call Rede
    call Upgrade
    
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
    T=Ti
    dH=0.01d0; Hxmin=0.0d0; Hxmax=0.6d0; Hx=Hxmin

    do while(Hx<=Hxmax .or. Hx>=Hxmin)
       Hy=0.0d0

    write(arquivo4, 400) 'histerese_T=',T,'_L=',L,'.dat'
    Open(4, file=arquivo4, position='append')
	
	!loop de passos de monte carlo
	do step=1, MC
	  !call Energia
	  !call Magnetizacao
	  if(n==1) then
	    write(2,*) step, E
	    write(3,*) step, dsqrt(Mx*Mx+My*My)
	  endif
	  
	  if(step>=terma_step) then
	    write(4) step, Hx, Mx
	  endif
	  call Metropolis  
	enddo
	!fim do loop de passos de monte carlo

	if(n==1) then
	    write(2,*)
	    write(3,*)
	endif
 
    Close(4)
	Hx=Hx+dH
	if(nint(100*Hx)>=nint(100*Hxmax)) then
	dH=-dH; pause 0; print*, Hx, dH
	endif
	if(nint(100*Hx)==0 .and. dH<0.0d0) then
	dH=0.0d0; pause 1
	endif
	
	if (n==1) then
	  call Video
	endif
! 	if (nint(100.0d0*T)==nint(100.0d0*Tf)) then
! 	  call Imagem
! 	endif
	enddo !fim loop do campo magnético
	
  enddo
  !fim do loop de amostras

    Close(2)
    Close(3)
  
END PROGRAM



!########LISTA DE VIZINHOS E ETC...
SUBROUTINE Upgrade

  use Variaveis
  implicit none
  real(8):: termo1, termo2
  integer(4):: i, j, cont				!contadores

  
  allocate(Nviz(0:L*L), viz(1:nint(3.15d0*rc*rc)*L*L), A(1:nint(3.15d0*rc*rc)*L*L))
  allocate(B(1:L*L))
  Nviz(0)=0; viz=0; A=0.0d0

  !Lista de vizinhos e calculo do vetor A
  cont=0
  do i=1, L*L
    do j=1, L*L
      dist=dsqrt((rx(i)-rx(j))*(rx(i)-rx(j)) + (ry(i)-ry(j))*(ry(i)-ry(j)))
      if(dist<rc .and. dist>0.1d0) then
	cont=cont+1
	viz(cont)=j

	termo1= dble(ex(i))*dble(ex(j)) + dble(ey(i))*dble(ey(j))
	termo2=( dble(ex(i))*(rx(j)-rx(i)) + dble(ey(i))*(ry(j)-ry(i)) )*( dble(ex(j))*(rx(j)-rx(i)) + dble(ey(j))*(ry(j)-ry(i)) )
	A(cont)= (termo1/dist**(3.0d0)) - ((3.0d0*termo2)/dist**(5.0d0))
      endif
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
      !vacância
      s(28)=0
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



!########CALCULO DA ENERGIA DIPOLAR com o termo do campo externo
SUBROUTINE Energia

use Variaveis
implicit none
  integer(4):: i, j					!contadores
  
  E=0.0d0; B=0.0d0
  do i=1, L*L
    do j=Nviz(i-1)+1, Nviz(i)
      B(i)=B(i)+A(j)*dble(s(viz(j)))
    enddo
    E=E+dble(s(i))*B(i) - s(i)*(dble(ex(i))*Hx+dble(ey(i))*Hy)
  enddo
    E=E*0.5d0
    
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



!########ALGORITMO DE METROPOLIS
SUBROUTINE Metropolis

use Variaveis
implicit none
  real(8)::Delta_E, Ei, Ef				!Variação de energia
  integer(4):: p					!sítio aleatorio
  integer(4):: k, i, j					!contador

  do k=1, L*L

    !sorteio de sítio aleatorio
    call random_number(x)
    p=nint(1.0d0+(L*L-1)*x)
 
 
    !calcular variação de energia
    Delta_E=-2.0d0*dble(s(p))*B(p)-2.0d0*s(p)*(dble(ex(p))*Hx+dble(ey(p))*Hy)    
    
    
    if(Delta_E<0.0d0) then
    
      E=E+Delta_E
      Mx=Mx-2.0d0*dble(s(p))*dble(ex(p))
      My=My-2.0d0*dble(s(p))*dble(ey(p))
      !atualização dos coeficientes B(p)
      do i=Nviz(p-1)+1, Nviz(p)
	B(viz(i))=B(viz(i))-2.0d0*A(i)*s(p)
      enddo
      s(p)=-1*s(p)
      
    else
    
      call random_number(y)
      if(dexp(-Delta_E/T)>y) then
	E=E+Delta_E
	Mx=Mx-2.0d0*dble(s(p))*dble(ex(p))
	My=My-2.0d0*dble(s(p))*dble(ey(p))
	!atualização dos coeficientes B(p)
	do i=Nviz(p-1)+1, Nviz(p)
	  B(viz(i))=B(viz(i))-2.0d0*A(i)*s(p)
	enddo
	s(p)=-1*s(p)
      endif
      
    endif

  enddo

ENDSUBROUTINE



!########ALGORITMO DE METROPOLIS MODIFICADO
SUBROUTINE Metropolis_mod

use Variaveis
implicit none
  real(8):: Ei, Ef, p1					!Energia anterior e posterior a flipagem
  integer(4):: p,p0					!sítio aleatorio
  integer(4):: i, j					!contadores


    !sorteio de sítio aleatorio
    call random_number(x)
    p=nint(1.0d0+(L*L-1)*x)

    !sorteio linha ou coluna
    call random_number(p1)
    
    !energia inicial
    call Energia
    Ei=E

    p0=-1
    if(p1>0.5d0) then
      !flipar coluna
      do i=0,L-1
	if(p+i<=L*L) then
	  s(p+i)=-1*s(p+i)
	else
	  s(p+i-L*L)=-1*s(p+i-L*L)
	endif
      enddo
    else
      !flipar linha
      do j=0, L-1
	if(p+j*L<=L*L) then
	  s(p+j*L)=-1*s(p+j*L)
	else
	  if(mod(p+j*L,L)==0) then
	    p0=0
	  endif
	  p0=p0+1
	  s(mod(p+j*L,L)+4*p0)=-1*s(mod(p+j*L,L)+4*p0)
	endif
      enddo
    endif

    !energia final
    call Energia
    Ef=E

    if((Ef-Ei)>0.0d0) then
      call random_number(y)
      if(dexp(-(Ef-Ei)/T)<y) then

	p0=-1
	if(p1>0.5d0) then
	  !flipar coluna
	  do i=0,L-1
	    if(p+i<=L*L) then
	      s(p+i)=-1*s(p+i)
	    else
	      s(p+i-L*L)=-1*s(p+i-L*L)
	    endif
	  enddo
	else
	  !flipar linha
	  do j=0, L-1
	    if(p+j*L<=L*L) then
	      s(p+j*L)=-1*s(p+j*L)
	    else
	      if(mod(p+j*L,L)==0) then
		p0=0
	      endif
	      p0=p0+1
	      s(mod(p+j*L,L)+4*p0)=-1*s(mod(p+j*L,L)+4*p0)
	    endif
	  enddo
	endif

      endif
    endif

ENDSUBROUTINE



!########VIDEO DA TERMALIZAÇÃO
SUBROUTINE Video

use Variaveis
implicit none
  character(50)::arquivo
  integer(4):: k					!contador

500 FORMAT(A8,I2,A7,F4.2,A7,F4.2,A4)
  !escrevendo a rede
  write(arquivo,500)'video_L=',L,'_Hxmin=',Hxmin,'_Hxmax=',Hxmax,'.xyz'
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

