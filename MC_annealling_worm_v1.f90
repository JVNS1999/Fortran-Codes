
module parametros
 integer, parameter :: dp = kind(1.0d0)
 integer :: Ns
 integer :: N_mc
 integer :: N_temp
 integer :: N_viz
 integer :: N_metro
 integer, dimension(:), allocatable :: jviz
 integer, dimension(:), allocatable :: Nviz
 real(dp), dimension(:), allocatable :: Aij
 real(dp), dimension(:), allocatable :: rx, ry, sx, sy
 real(dp) :: t, t0, tf, dt, beta
end module parametros

module variaveis
 use parametros, only : dp
 integer, dimension(:), allocatable :: Spin
 real(dp), dimension(:), allocatable :: Bij
 real(dp) :: En_total
 real(dp) :: Mx, My, M_total
end module variaveis

module loop
 use parametros, only : dp
 integer :: Nvkt
 integer :: try,acc
 integer, dimension(:), allocatable :: Vk, Sk
 integer, dimension(:), allocatable :: Vt, St
 real(dp), dimension(:), allocatable :: Rk, Rt
end module loop

subroutine inicial
 use parametros
 use variaveis
 use loop
 implicit none
 integer,parameter :: seed=45126
 integer :: i,j,k
 integer :: ii,jj
 real(dp) :: E1,E2
 
 call srand(seed)
 
 open(10,file='input_ini.dat',status='old',action='read')
 read(10,*) Ns
 read(10,*) N_viz
 read(10,*) N_mc
 read(10,*) N_temp
 read(10,*) t0
 read(10,*) tf
 read(10,*) Nvkt
 close(10)
 dt = (tf-t0)/real((N_temp-1),dp)
 
 N_metro = Ns
 
 allocate(Aij(N_viz),jviz(N_viz))
 allocate(Nviz(0:Ns),Spin(Ns))
 allocate(rx(Ns),ry(Ns),sx(Ns),sy(Ns))
 allocate(Bij(Ns))
 allocate(Vk(Nvkt),Sk(Nvkt),Rk(Nvkt))
 allocate(Vt(Nvkt),St(Nvkt),Rt(Nvkt))
 
 open(20,file='Aij.dat',status='old',action='read')
 do i = 1,N_viz
  read(20,*) jviz(i),Aij(i)
 end do
 close(20)

 open(30,file='Nviz.dat',status='old',action='read')
 open(40,file='config0.xyz',status='old',action='read')
 read(40,*) i
 read(40,*) 
 read(30,*) Nviz(0)

 do i = 1,Ns
  read(30,*) Nviz(i)
  read(40,*) j,rx(i),ry(i),sx(i),sy(i)
 end do
 close(30)
 close(40)

 Spin = 1

 Bij = 0.0_dp
 En_total = 0.0_dp
 do i = 1,Ns
  do k = Nviz(i-1)+1,Nviz(i)
   j = jviz(k)
   Bij(i) = Bij(i) + Spin(j)*Aij(k)
  end do
  En_total = En_total + Spin(i)*Bij(i)
 end do
 
 En_total = 0.5d0*En_total
 
 open(50,file='S_k.dat',status='old',action='read')
 open(60,file='S_t.dat',status='old',action='read')
 open(70,file='V_k.dat',status='old',action='read')
 open(80,file='V_t.dat',status='old',action='read')
 
 do i = 1,Nvkt
  read(50,*) j,Sk(i)
  read(60,*) j,St(i)
  read(70,*) j,Vk(i),Rk(i)
  read(80,*) j,Vt(i),Rt(i)
 end do

 close(50); close(60); close(70); close(80)
 
 t = t0
 
 Mx = sum(Spin(:)*sx(:))
 My = sum(Spin(:)*sy(:))
 M_total = sqrt(Mx*Mx + My*My)
 
 open(100,file='ener_test.dat')
 open(200,file='config_test.xyz')
 open(300,file='config_worm.xyz')
 open(400,file='dE_worm.dat')
 write(100,*) N_mc
 write(100,*) t
 
 try = 0
 acc = 0

 !call Teste_energia
 
 !stop
 
 do i = 1,N_mc/2
   call metropolis
   !if (rand()<0.5d0) then
    call worm_K
   !else
   ! call worm_T
   !end if
   call salva_dados(i)
 end do
 close(100)
 close(200)
 close(300)
 close(400)
 
 print*, acc,try,real(acc,8)/real(try,8)
  

 return
end subroutine inicial

subroutine metropolis
 use parametros, only : Ns, N_metro, t
 use variaveis, only : dp, En_total, Bij, Spin
 implicit none
 integer :: i_metro,i
 real(dp) :: dE,b
 
 b = 1.0d0/t
 
 do i_metro = 1,N_metro
  i = int(Ns*rand()) + 1
  dE = -2.0_dp*Spin(i)*Bij(i)
  if ( (dE < 0.0_dp) .or. (rand() < dexp(-dE/t)) ) then
   En_total = En_total + dE
   call aceita(i)
  end if
 end do
 
 return
end subroutine metropolis

subroutine aceita(o)
 use parametros, only : dp, Nviz, jviz, Aij
 use variaveis, only : Spin, Bij
 implicit none
 integer, intent(in) :: o
 integer :: j,k
 real(dp) :: dBij
 
 do k = Nviz(o-1)+1,Nviz(o)
  j = jviz(k)
  dBij = -2.0_dp*Spin(o)*Aij(k)
  Bij(j) = Bij(j) + dBij
 end do
 Spin(o) = -Spin(o)
 
 return
end subroutine aceita
 
subroutine avevar(data,ave,var)
 use parametros, only : dp
 implicit none
 real(dp), dimension(:), intent(in) :: data
 real(dp), intent(out) :: ave,var
 integer :: N
 real(dp), dimension(size(data)) :: s
 
 N = size(data)
 
 ave = sum(data(:))/real(N,dp)
 s(:) = data(:) - ave
 var = dot_product(s,s)
 var = (var - sum(s)**2/real(N,dp))/real(N-1,dp)
 
 return
end subroutine avevar

subroutine worm_K
 use parametros, only : dp,Ns
 use variaveis, only : Spin
 use loop, only : Nvkt, Vk, Sk, Rk, try, acc
 implicit none
 integer :: i,j
 integer :: iworm, cont
 integer :: v_0, ivk, isk
 integer :: v_worm(Nvkt), s_worm(Ns)
 integer :: v_sequ(2*Nvkt), s_sequ(2*Nvkt)
 integer :: pool(3)
 
 v_worm(:) = 0
 s_worm(:) = 0
 v_sequ(:) = 0
 s_sequ(:) = 0
 
 v_0 = int(Nvkt*rand()/3.0_dp) + 1
 ivk = v_0
 v_worm(ivk) = 1
 v_sequ(1) = ivk
 
 iworm = 0
 
 try = try + 1
 
 do while (iworm <= 2*Nvkt)
  iworm = iworm + 1
  cont = 0
  pool(:) = 0
  
  do i = 3*(ivk-1)+1,3*(ivk-1)+3
   j = Vk(i)
   if (Rk(i)*Spin(j) < 0.0_dp) then
    cont = cont + 1
    pool(cont) = j
   end if
  end do
  
  if (cont .ne. 0) then
   isk = pool(int(cont*rand())+1)
   s_worm(isk) = 1
   s_sequ(iworm) = isk
  else
   return
  end if
  
  if (Sk(2*(isk-1)+1) .ne. ivk) then
   ivk = Sk(2*(isk-1)+1)
  else
   ivk = Sk(2*(isk-1)+2)
  end if
  
  if (v_worm(ivk) .eq. 1) then
   if (ivk .eq. v_0) then
    call metropolis_worm(s_worm)
    return
   else
    call corta_rabo(Ns,2*Nvkt,ivk,s_worm,s_sequ,v_sequ)
    call metropolis_worm(s_worm)
    return
   end if
  else
   v_worm(ivk) = 1
   v_sequ(iworm+1) = ivk
  end if
 end do
 
 return
 
end subroutine worm_K

subroutine corta_rabo(Ns,N,ivk,s_worm,s_sequ,v_sequ)
 implicit none
 integer, intent(in) :: Ns, N, ivk
 integer, intent(in) :: s_sequ(N), v_sequ(N)
 integer, intent(inout) :: s_worm(Ns)
 integer :: i,j
 
 do i = 1,N
  if (v_sequ(i) .ne. ivk) then
   j = s_sequ(i)
   s_worm(j) = 0
  else if (v_sequ(i) .eq. ivk) then
   return
  end if
 end do
 
 return
 
end subroutine corta_rabo

subroutine metropolis_worm(s_worm)
 use parametros, only : dp, Ns, Nviz, jviz, Aij, t, rx, ry, sx, sy
 use variaveis, only : Spin, Bij, En_total
 use loop, only : acc,try
 implicit none
 integer, dimension(Ns), intent(in) :: s_worm
 integer :: i,j,k
 integer, dimension(Ns) :: Sn
 real(dp) :: E0, En, dE
 real(dp), dimension(Ns) :: Bij_n
 
 Sn = Spin
 do i = 1,Ns
  if (s_worm(i) == 1) then
   Sn(i) = -Sn(i)
  end if
 end do
 
 E0 = En_total
 En = 0.0_dp
 
 Bij_n(:) = 0.0_dp
 do i = 1,Ns
  do k = Nviz(i-1)+1,Nviz(i)
   j = jviz(k)
   Bij_n(i) = Bij_n(i) + Sn(j)*Aij(k)
  end do
  En = En + Sn(i)*Bij_n(i)
 end do
 
 En = 0.5d0*En
 
 dE = En - E0
 
 write(400,*) try, dE, -dE/t
 
 if (dE < 0.0_dp) then
  Spin = Sn
  En_total = En
  Bij = Bij_n
  acc = acc + 1
  
  write(300,*) Ns
  write(300,*) ' '
  do i = 1,Ns
   write(300,*) rx(i),ry(i),Spin(i)*sx(i),Spin(i)*sy(i),s_worm(i)
  end do
  
 else if ((rand()<dexp(-dE/t))) then
  Spin = Sn
  En_total = En
  Bij = Bij_n
  acc = acc + 1
  
  write(300,*) Ns
  write(300,*) ' '
  do i = 1,Ns
   write(300,*) rx(i),ry(i),Spin(i)*sx(i),Spin(i)*sy(i),2*s_worm(i)
  end do
 end if
 
 !write(300,*) Ns
 !write(300,*) ' '
 !do i = 1,Ns
 ! write(300,*) rx(i),ry(i),Spin(i)*sx(i),Spin(i)*sy(i),s_worm(i)
 !end do
 
 return
end subroutine metropolis_worm
    
subroutine salva_dados(i)
 use parametros, only : dp, Ns, rx, ry, sx, sy
 use variaveis, only : Spin, En_total
 implicit none
 integer,intent(in)::i
 integer :: j
 real(dp) :: Mx, My, M_total
 
 
 Mx = sum(Spin*sx)/real(Ns,dp)
 My = sum(Spin*sy)/real(Ns,dp)
 M_total = sqrt(Mx*Mx + My*My)
 
 write(100,*) i,En_total/real(Ns,dp),Mx,My,M_total
 
 if (mod(i,1) == 0) then
  write(200,*) Ns
  write(200,*) ' '
  do j = 1,Ns
   write(200,*) rx(j),ry(j),Spin(j)*sx(j),Spin(j)*sy(j)
  end do
 end if

end subroutine salva_dados

subroutine worm_T
 use parametros, only : dp,Ns
 use variaveis, only : Spin
 use loop, only : Nvkt, Vt, St, Rt, try, acc
 implicit none
 integer :: i,j
 integer :: iworm, cont
 integer :: v_0, ivk, isk
 integer :: v_worm(Nvkt), s_worm(Ns)
 integer :: v_sequ(2*Nvkt), s_sequ(2*Nvkt)
 integer :: pool(6)
 
 v_worm(:) = 0
 s_worm(:) = 0
 v_sequ(:) = 0
 s_sequ(:) = 0
 
 v_0 = int(Nvkt*rand()/6.0_dp) + 1
 ivk = v_0
 v_worm(ivk) = 1
 v_sequ(1) = ivk
 
 iworm = 0
 
 try = try + 1
 
 do while (iworm <= 2*Nvkt)
  iworm = iworm + 1
  cont = 0
  pool(:) = 0
  
  do i = 6*(ivk-1)+1,6*(ivk-1)+6
   j = Vt(i)
   if (Rt(i)*Spin(j) < 0.0_dp) then
    cont = cont + 1
    pool(cont) = j
   end if
  end do
  
  if (cont .ne. 0) then
   isk = pool(int(cont*rand())+1)
   s_worm(isk) = 1
   s_sequ(iworm) = isk
  else
   return
  end if
  
  if (St(2*(isk-1)+1) .ne. ivk) then
   ivk = St(2*(isk-1)+1)
  else
   ivk = St(2*(isk-1)+2)
  end if
  
  if (v_worm(ivk) .eq. 1) then
   if (ivk .eq. v_0) then
    call metropolis_worm(s_worm)
    return
   else
    call corta_rabo(Ns,2*Nvkt,ivk,s_worm,s_sequ,v_sequ)
    call metropolis_worm(s_worm)
    return
   end if
  else
   v_worm(ivk) = 1
   v_sequ(iworm+1) = ivk
  end if
 end do
 
 return
 
end subroutine worm_T


subroutine Teste_energia
 use parametros, only : dp, Ns, Nviz, jviz, Aij
 use variaveis, only : Spin, Bij, En_total
 implicit none
 integer :: ii,jj
 integer :: i,j,k
 real(dp) :: E0, En
 real(dp), dimension(Ns) :: Bij_n
 
 open(600,file='TESTE_energia_erro.dat')
 
 do ii = 1,1000
  call metropolis
 end do
  
 do ii = 1,100
  do jj = 1,1000
   call metropolis
  end do
 
  E0 = En_total
  En = 0.0_dp
  
  Bij_n(:) = 0.0_dp
  do i = 1,Ns
   do k = Nviz(i-1)+1,Nviz(i)
    j = jviz(k)
    Bij_n(i) = Bij_n(i) + Spin(j)*Aij(k)
   end do
   En = En + Spin(i)*Bij_n(i)
  end do
  
  En = 0.5d0*En
  
  print*, ii*jj, E0-En
 end do
 
 return
end subroutine Teste_energia
  
  



program teste
 call inicial
end program teste
 
 
 
 
 
 
 
 
 
 
 
 
  
