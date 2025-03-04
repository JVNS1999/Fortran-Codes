module parametros
 integer::Ns = 5
 integer::N_mc
 integer::N_temp
 integer::N_viz
 integer::N_metro
 integer,dimension(:),allocatable::jviz
 integer,dimension(:),allocatable::Nviz
 real(8),dimension(:),allocatable::Aij
 real(8)::t,t0,tf,dt
end module parametros

module variaveis
 integer,dimension(:),allocatable::Spin
 real(8),dimension(:),allocatable::Bij,E,M
 real(8),dimension(:),allocatable::rx,ry,sx,sy
 real(8)::En_total,Mx,My
 real(8)::time0,timef
end module variaveis

program Annealing
 use variaveis,only:time0,timef
 implicit none
 
 call cpu_time(time0)
 call inicial
 call MC_termal
 call cpu_time(timef)
 
 print*, timef-time0
 
end program Annealing

subroutine inicial
 use parametros
 use variaveis
 implicit none
 integer,parameter::seed=4567
 integer::i,j,k
 
 call srand(seed)
 
 open(1,file='input_ini.dat',status='old',action='read')
 read(1,*) Ns
 read(1,*) N_viz
 read(1,*) N_mc
 read(1,*) N_temp
 read(1,*) t0
 read(1,*) tf
 close(1)
 dt=(tf-t0)/real((N_temp-1),8)
 
 N_metro = Ns
 
 allocate(Aij(N_viz),jviz(N_viz))
 open(2,file='Aij.dat',status='old',action='read')
 do i=1,N_viz
  read(2,*) jviz(i),Aij(i)
 end do
 close(2)
 
 allocate(Nviz(0:Ns),Spin(Ns))
 open(3,file='Nviz.dat',status='old',action='read')
 do i=0,Ns
  read(3,*) Nviz(i)
 end do
 close(3)
 
 allocate(rx(Ns),ry(Ns),sx(Ns),sy(Ns))
 open(4,file='config0.xyz',status='old',action='read')
 read(4,*) i
 read(4,*) 
 do i = 1,Ns
  read(4,*) rx(i),ry(i),sx(i),sy(i)
 end do
 close(4)
 
 Spin = 1
 do i = 1,Ns
  if (rand()<0.5d0) Spin(i) = -1
 end do
 
 Mx = sum(Spin(:)*sx(:))
 My = sum(Spin(:)*sy(:))
 
 !-----Cálculo do campo local e da energia total inicial-----!
 allocate(Bij(Ns))
 Bij = 0.0d0
 En_total = 0.0d0
 do i=1,Ns
  do k=Nviz(i-1)+1,Nviz(i)
   j=jviz(k)
   Bij(i) = Bij(i) + Spin(j)*Aij(k)
  end do
  En_total = En_total + Spin(i)*Bij(i)
 end do
  
 En_total = 0.5d0*En_total
 print*, En_total
 
 !-----------------------------------------------------------!
end subroutine inicial
  
subroutine MC_termal
 use parametros
 use variaveis,only:En_total,E,Mx,My,M,Spin,rx,ry,sx,sy
 implicit none
 integer,dimension(Ns)::S_med
 integer::i,N_terma
 integer::i_mc,i_temp,i_terma
 real(8)::E1,E2,sig,A
 real(8)::M1,M2,sig2
 
 N_terma=1.d0*Ns
 allocate(E(N_mc),M(N_mc))
 A = 1.0d0/real(N_mc-1,8)
 
 !open(10,file='Annealing.dat')
 open(11,file='Config_annealing.xyz')
 open(12,file='Energia.dat')
 open(13,file='Magnetizacao.dat')
 open(14,file='Histograma.dat')
 !print*, dt,N_temp
 dt = 0.0d0
 do i_temp = 1,N_temp
  t = t0 + real(i_temp-1,8)*dt
  !E = 0.0d0
  !M = 0.0d0
  !write(10,*) t
  do i_terma = 1,N_terma
   call metropolis
  end do
  !S_med = 0
  do i_mc = 1,N_mc
   call metropolis
   !S_med = S_med + Spin
   write(14,*) En_total
   !E(i_mc) = En_total
   !M(i_mc) = sqrt(Mx**2 + My**2)
   if (mod(i_mc,5).eq.0) then
    call config
   end if
   !write(10,*) i_mc*i_temp,En_total
  end do
  !E = E/real(Ns,8)
  !M = M/real(Ns,8)
  !call avevar(N_mc,E,E1,E2)
  !call avevar(N_mc,M,M1,M2)
  !sig = sqrt(A*(E2 - E1*E1))
  !sig2 = sqrt(A*(M2-M1*M1))
  !write(12,*) t,E1,sig,(E2-E1*E1)/(t*t)
  !write(13,*) t,M1,sig2,(M2-M1*M1)/t
  !call flush()
  
 end do
 !close(10)
 !open(11,file='Config_annealing.xyz')
 !do i = 1,Ns
 ! write(11,*) Spin(i)
 !end do
 close(11)
 close(12)
 close(13)
 close(14)
 
 open(15,file='Spin_med.xyz')
 write(15,*) Ns
 write(15,*) ' '
 do i = 1,Ns
  write(15,*) rx(i),ry(i),sx(i)*S_med(i)/real(N_mc,8),sy(i)*S_med(i)/real(N_mc,8)
 end do
 close(15)


 deallocate(E)
 
 return
end subroutine MC_termal

subroutine config
 use variaveis
 use parametros
 implicit none
 integer::i
 
 write(11,*) Ns
 write(11,*) t
 do i = 1,Ns
  write(11,*) rx(i),ry(i),sx(i)*Spin(i),sy(i)*Spin(i),Spin(i)*Bij(i)
 end do
 return
end subroutine config

subroutine avevar(N,data,A1,A2)
 implicit none
 integer,intent(in)::N
 real(8),dimension(N),intent(in)::data
 real(8),intent(out)::A1,A2
 
 A1 = sum(data)/real(N,8)
 A2 = sum(data*data)/real(N,8)
 
 return
end subroutine avevar

subroutine metropolis
 use variaveis,only:En_total,Bij,Spin,Mx,My,sx,sy
 use parametros,only:Ns,N_metro,t
 implicit none
 integer::i_metro,i
 real(8)::dE
 
 do i_metro = 1,N_metro
  i = int(Ns*rand())+1
  !call dEn(i,dE)
  dE = -2.0d0*Spin(i)*Bij(i)
  if (dE<0.0d0 .or. (rand()<exp(-dE/t))) then
   En_total = En_total + dE
   Mx = Mx - 2.0d0*Spin(i)*sx(i)
   My = My - 2.0d0*Spin(i)*sy(i)
   call aceita(i)
  end if
 end do
 
 return
end subroutine metropolis

subroutine aceita(i)
 use parametros
 use variaveis
 implicit none
 integer,intent(in)::i
 integer::j,k
 real(8)::dBij
 
 do k = Nviz(i-1)+1,Nviz(i)
  j = jviz(k)
  dBij = -2.0d0*Spin(i)*Aij(k)
  Bij(j) = Bij(j) + dBij
 end do
 Spin(i) = -Spin(i)
 
 return
end subroutine aceita

