!Este programa visa estudar as propriedades físicas 
!de materias magnéticos na forma de malha, onde cada
!malha corresponde a uma nanoilha alongada




!########DEFINIÇÃO DE VARIAVÉIS
MODULE Variaveis 
implicit none
  integer(4), parameter:: L=4				!tamanho da rede
  real(8):: Sx(1:L*L), Sy(1:L*L)			!componentes do spin de cada sítio
  real(8):: rx(1:L*L), ry(1:L*L)			!componentes da posição de cada sítio
  real(8), parameter:: Ti=2.0d0, Tf=0.1d0, dT=-0.1d0	!temperatura inicial, final e passo
  real(8):: T						!temperatura
  real(8):: E, M					!energia e magnetização total do sistema
  real(8):: dist					!distância entre sítios
  real(8):: x, y					!variaveis aleatorias
  integer(4), parameter:: MC=10e3, terma_step=3e3	!passo máximo de monte carlo, passo da termalização do sistema
  integer(4):: i, j, k, n, step, cont			!contadores
  integer(4), parameter::num_max=20, config_ini=1	!número máximo de amostras, configuração inicial do sistema
ENDMODULE



!########ESCOPO PRINCIPAL DO PROGRAMA
PROGRAM Principal

use Variaveis
implicit none
character(50)::arquivo1, arquivo2, arquivo3, arquivo4	!índice do aquivo
 
100 FORMAT(A14)
200 FORMAT(A15, F4.2, A4)
300 FORMAT(A15,I2,A4)
400 FORMAT(A20,I2,A4)

    
    Open(1, file='parametros.dat', form='unformatted')
      write(1) Ti, Tf, dT, terma_step, MC, num_max, L
    Close(1)

    write(arquivo3, 300) 'Energia_x_MC_L=',L,'.dat'
    Open(3, file=arquivo3)
    write(arquivo4, 400)'Magnetizacao_x_MC_L=',L,'.dat'
    Open(4, file=arquivo4)

  !loop de temperatura
  do T=Ti, Tf, dT

    write(arquivo2, 200) 'termalizacao_T=',T,'.dat'
    Open(2, file=arquivo2, form='unformatted')

      !loop de amostras
      do n=1, num_max
	call Rede
	
	!loop de passos de monte carlo
	do step=1, MC
	  call Energia
	  call Magnetizacao
	  if(n==1) then
	    write(3,*) step, E
	    write(4,*) step, M
	  endif
	  if(step>=terma_step) then
	    write(2) step, E, M
	  endif
	  call Metropolis  
	enddo
	!fim do loop de passos de monte carlo

	if (n==1) then
	  call Video
	endif
      enddo
      !fim do loop de amostras

	    write(3,*)
	    write(4,*)
 
    Close(2)

  enddo
  !fim do loop de temperatura

    Close(3)
    Close(4)
  
ENDPROGRAM



!########MONTAGEM DA REDE(MALHA)
SUBROUTINE Rede

use Variaveis

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
	    Sx(cont)=-1.0d0 + 2.0d0*nint(x)
	    Sy(cont)=0.0d0

	  else				!ilhas horizontais
	    Sx(cont)=0.0d0
	    Sy(cont)=-1.0d0 + 2.0d0*nint(x)
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
	    Sy(cont)=0.0d0
	    if( ry(cont)>dble(L-1)/2.0d0 ) then
	      Sx(cont)=1.0d0
	    else 
	      Sx(cont)=-1.0d0
	    endif

	  else				!ilhas horizontais
	    Sx(cont)=0.0d0
	    if( rx(cont)>dble(L-1)/2.0d0 ) then
	      Sy(cont)=-1.0d0
	    else
	      Sy(cont)=1.0d0
	    endif
	  endif

	enddo
      enddo
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  end select

ENDSUBROUTINE



!########CALCULO DA ENERGIA DIPOLAR
SUBROUTINE Energia

use Variaveis
implicit none
  real(8):: termo1, termo2

  E=0.0d0
  do i=1, L*L
    do j=i+1, L*L
      
      !if(j==i) CYCLE
      dist = dsqrt( (rx(j)-rx(i))**(2.0d0) + (ry(j)-ry(i))**(2.0d0) )
      termo1 = Sx(i)*Sx(j) + Sy(i)*Sy(j)
      termo2 = ( Sx(i)*(rx(j)-rx(i)) + Sy(i)*(ry(j)-ry(i)) )*( Sx(j)*(rx(j)-rx(i)) + Sy(j)*(ry(j)-ry(i)) )

      E = E + (termo1)/(dist**(3.0d0)) - 3.0d0*(termo2)/(dist**(5.0d0))
      
    enddo
  enddo
!Print*, E
ENDSUBROUTINE



!########CALCULO DA MAGNETIZAÇÃO
SUBROUTINE Magnetizacao

  use Variaveis
  implicit none
  real(8):: Mx, My					!componentes x e y da magnetização

  cont=0
  Mx=0.0d0
  My=0.0d0
  M=0.0d0
  do 
    cont=cont+1
    Mx=Mx+Sx(cont)
    My=My+Sy(cont)
    if(cont==L*L) EXIT
  enddo
  M=dsqrt(Mx*Mx+My*My)
ENDSUBROUTINE



!########ALGORITMO DE METROPOLIS
SUBROUTINE Metropolis

use Variaveis
implicit none
  real(8):: Ei, Ef					!Energia anterior e posterior a flipagem
  integer(4):: p					!sítio aleatorio

  do k=1, L*L

    if(mod(k,5)==0) then
      call Metropolis_mod
    endif

    !sorteio de sítio aleatorio
    call random_number(x)
    p=nint(1.0d0+(L*L-1)*x)

    !energia inicial
    call Energia
    Ei=E

    !flipar sítio
    Sx(p)=-1.0d0*Sx(p)
    Sy(p)=-1.0d0*Sy(p)

    !energia final
    call Energia
    Ef=E

    if((Ef-Ei)>0.0d0) then
      call random_number(y)
      if(dexp(-(Ef-Ei)/T)<y) then
	Sx(p)=-1.0d0*Sx(p)
	Sy(p)=-1.0d0*Sy(p)
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

!    do k=1, L*L

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
	  Sx(p+i)=-1.0d0*Sx(p+i)
	  Sy(p+i)=-1.0d0*Sy(p+i)
	else
	  Sx(p+i-L*L)=-1.0d0*Sx(p+i-L*L)
	  Sy(p+i-L*L)=-1.0d0*Sy(p+i-L*L)
	endif
      enddo
    else
      !flipar linha
      do j=0, L-1
	if(p+j*L<=L*L) then
	  Sx(p+j*L)=-1.0d0*Sx(p+j*L)
	  Sy(p+j*L)=-1.0d0*Sy(p+j*L)
	else
	  if(mod(p+j*L,L)==0) then
	    p0=0
	  endif
	  p0=p0+1
	  Sx(mod(p+j*L,L)+4*p0)=-1.0d0*Sx(mod(p+j*L,L)+4*p0)
	  Sy(mod(p+j*L,L)+4*p0)=-1.0d0*Sy(mod(p+j*L,L)+4*p0)
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
	      Sx(p+i)=-1.0d0*Sx(p+i)
	      Sy(p+i)=-1.0d0*Sy(p+i)
	    else
	      Sx(p+i-L*L)=-1.0d0*Sx(p+i-L*L)
	      Sy(p+i-L*L)=-1.0d0*Sy(p+i-L*L)
	    endif
	  enddo
	else
	  !flipar linha
	  do j=0, L-1
	    if(p+j*L<=L*L) then
	      Sx(p+j*L)=-1.0d0*Sx(p+j*L)
	      Sy(p+j*L)=-1.0d0*Sy(p+j*L)
	    else
	      if(mod(p+j*L,L)==0) then
		p0=0
	      endif
	      p0=p0+1
	      Sx(mod(p+j*L,L)+4*p0)=-1.0d0*Sx(mod(p+j*L,L)+4*p0)
	      Sy(mod(p+j*L,L)+4*p0)=-1.0d0*Sy(mod(p+j*L,L)+4*p0)
	    endif
	  enddo
	endif

      endif
    endif

!    enddo
ENDSUBROUTINE



!########DADOS DO VÍDEO
SUBROUTINE Video

use Variaveis
implicit none
character(50)::arquivo

500 FORMAT(A13,I2,A3,F4.2,A4)
  !escrevendo a rede
  write(arquivo,500)'rede_final_L=',L,'_T=',T,'.xyz'
  Open(5, file=arquivo)
    501 FORMAT(1x,'H',3E20.7,'atom_vector',3E20.7)
    write(5,*) L*L
    write(5,*) 'Malha'
    do k=1, L*L
      write(5,501) rx(k), ry(k), 0.0d0, Sx(k), Sy(k), 0.0d0
    enddo
  Close(5)


ENDSUBROUTINE

!Open(3, file='rede_teste_L10.xyz', position='append')
!Lembrete: O comando position='append' permite continuar escrevendo
!no mesmo arquivo sem altera-lo.

