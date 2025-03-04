module parametros
 integer(4),parameter :: prec = selected_real_kind(8)
 
 real(prec),parameter :: pi = 4.0_prec*atan(1.0_prec), pi2 = 2.0_prec*pi
 real(prec),parameter :: itiny = 1.0d-30
end module parametros

module variaveis
 use parametros

!======================================================================!
!=========================Parâmetros de Entrada========================!
!======================================================================!
 integer(4) :: N_sp, N_ress, N_pass, Lxy, N_samp, deltaN, NN
 real(prec) :: iL, iA, iW, alpha, h, damp, min_torque
 real(prec) :: iTroca, iDip, iAnis, iZee
 real(prec) :: tpc, pulso_width, angle, Bmax,omega
!======================================================================!
!======================================================================!
!======================================================================!

 integer(4) :: N_sl, N_st, N_dip, P, N_pulso,simu_i
 integer,parameter :: iuni_e = 10, iuni_mag = 11, iuni_xyz = 12, iuni_pul = 13
 integer,parameter :: iuni_mgSp = 14, iuni_mgSl = 15, iuni_hist = 16, iuni_fft = 17
 real(prec) :: iSp, iSl, time, idamp, L2, Ms
 real(prec) :: Bextx, Bexty, Bextz, tpmin
 real(prec),dimension(:,:),allocatable :: posixy
 
 real(prec),dimension(:),allocatable :: Sx,Sy,Sz,Snx,Sny,Snz,anisx,anisy !====!
 real(prec),dimension(:),allocatable :: rijx,rijy,distij !====!
 
 real(prec),dimension(:),allocatable :: B_pulso
 
 character(500) :: dir1,dir2,dir3,dir4,dir5,arqui

end module variaveis

module ABM
 use parametros, only : prec
 real(prec),dimension(:),allocatable :: F0x,F1x,F2x,F3x,F4x,EffFieldx   !====!
 real(prec),dimension(:),allocatable :: F0y,F1y,F2y,F3y,F4y,EffFieldy   !====!
 real(prec),dimension(:),allocatable :: F0z,F1z,F2z,F3z,F4z,EffFieldz   !====!
 real(prec) :: h24,Erro_passo,Erro_Total,Max_Torque
end module ABM

subroutine inicial()
 use parametros
 use variaveis
 use ABM
 implicit none
 integer(4) :: i,ierro
 real(prec) :: tp,B_pul
 logical :: direx
 
 !=====================================================================!
 !========================Ler dados iniciais===========================!
 open(unit=1,file='input',status='old',iostat=ierro)
 
 if (ierro .eq. 0) then
  read(1,*) Lxy
  read(1,*) iL, iA, iW
  read(1,*) N_samp, deltaN
  read(1,*) h, damp, min_torque
  read(1,*) iTroca, iDip, iAnis, iZee
  read(1,*) tpc, pulso_width, Bmax, angle, omega
 else
  write(*,*) 'Erro na leitura do arquivo de entrada'
  stop
 end if
 
 close(1)

!======================================================================!
!======================================================================!
 
!======================================================================!
!      Determina os parâmetros gerais da rede e das nano-ilhas         !
!======================================================================!

 iL = iL/iA; iW = iW/iA; iA = 1.0_prec
 L2 = 0.5_prec*real(Lxy,prec) + 1.0_prec
 
 iSp = iW*(iL - iW)
 iSl = (pi*iW**2)/8.0_prec
 
 P = 2*Lxy + 1
 N_sp = 2*Lxy*(Lxy + 1)
 N_sl = 2*N_sp
 N_st = N_sp + N_sl
 N_dip = N_st*(N_st-3)
 
 alpha = 0.5_prec*(iL - iW) + (2.0_prec*iW)/(3.0_prec*pi)
 idamp = 1.0_prec/(1.0_prec + damp**2)
 
 time = 0.0_prec
 Ms = 1.0_prec/(real(N_Sp,prec)*(iSp+2.0_prec*iSl))
 
 N_samp = 2**N_samp
 NN = N_samp

!======================================================================!
!  Define o intervalo de tempo em que o pulso é aplicado               !
!======================================================================!
 tp = 0.0_prec
 tpmin = 0.0_prec
 N_pulso = INT(2.0_prec*tpc/h)
 do i = 0,N_pulso
  tp = i*h
  B_pul = Bmax*exp(-(((tp - tpc))/(pulso_width))**2)
  if (B_pul > 10.0d-5) then
   tpmin = tp
   exit
  end if
 end do
 N_pulso = INT(2.0_prec*(tpc-tpmin)/h)
 angle = angle*pi/180.0_prec
 allocate(B_pulso(1:N_pulso))
 B_pulso = 0.0_prec
 do i = 1,N_pulso
  tp = i*h + tpmin
  B_pulso(i) = Bmax*exp(-(((tp - tpc))/(pulso_width))**2)
 end do

!======================================================================!

 !call rede()
 call rede1()

!======================================================================!
!=========================Preparar diretórios==========================!
!======================================================================!

 write(arqui,"('paramtros_',I4.4,'_.dat')") simu_i
 open(110,file=arqui)
 write(110,*) '==================================================='
 write(110,*) 'Lxy = ',Lxy
 write(110,*) ' '
 write(110,*) 'iL, iA, iW = ',iL,iA,iW
 write(110,*) '==================================================='
 write(110,*) 'h, damp, min_torque = ', h, damp, min_torque
 write(110,*) ' '
 write(110,*) 'iTroca, iDip, iAnis, iZee = ', iTroca, iDip, iAnis, iZee
 write(110,*) '==================================================='
 write(110,*) 'tpc, pulso_width, Bmax, angle = ', tpc, pulso_width, Bmax, angle
 write(110,*) '==================================================='
 write(110,*) 'iSp, iSl = ', iSp, iSl
 write(110,*) ' '
 write(110,*) 'P, N_sp, N_sl, N_st, N_dip, N_pulso = ', P, N_sp, N_sl, N_st, N_dip, N_pulso
 write(110,*) ' '
 write(110,*) 'alpha, idamp', alpha, idamp
 write(110,*) '==================================================='
 close(110)
 
 !write(arqui,"('Anima.xyz')")
 !inquire(file = trim(dir3) // trim(arqui),exist = direx)
 !if (direx .eqv. .false.) then
 ! open(unit=iuni_xyz,file = trim(dir3) // trim(arqui),iostat=ierro)
 !end if
 
 !write(arqui,"('Pulso.xyz')")
 !inquire(file = trim(dir3) // trim(arqui),exist = direx)
 !if (direx .eqv. .false.) then
 ! open(unit=iuni_pul,file = trim(dir3) // trim(arqui),iostat=ierro)
 !end if

 call flush()

end subroutine inicial

subroutine rede1()
 use parametros
 use variaveis
 implicit none
 integer :: i,j,k,k1,k2,k3,kanis
 real(prec) :: delx,dely
 
 open(2,file='config_ini.xyz')
 read(2,*) N_st
 read(2,*)
 
 allocate(Sx(N_st),Sy(N_st),Sz(N_st),Snx(N_st),Sny(N_st),Snz(N_st),posixy(N_st,2))
 allocate(anisx(N_sp),anisy(N_sp))
 
 posixy = 0.0_prec
 Sx = 0.0_prec; Sy = 0.0_prec; Sz = 0.0_prec; Snx = 0.0_prec; Sny = 0.0_prec; Snz = 0.0_prec;  !====!
 anisx = 0.0_prec; anisy = 0.0_prec  !====!
 
 do i = 1,N_st
  read(2,*) k,posixy(i,1),posixy(i,2),Sx(i),Sy(i),Sz(i)
 end do
 
 close(2)

 k = -2
 kanis = 0
 do j = 1,P
  do i = 1,P
   if (mod(i+j,2).ne.0) then
    k = k + 3
    k1 = k + 1
    k2 = k + 2
    kanis = kanis + 1
    
    if (mod(j,2).ne.0) then

     anisx(kanis) = 1.0_prec

    else

     anisy(kanis) = 1.0_prec

    end if
   end if
  end do
 end do
 
 allocate(rijx(1:N_st*(N_st-3)),rijy(1:N_st*(N_st-3)),distij(1:N_st*(N_st-3)))

 distij = 0.0_prec; rijx = 0.0_prec; rijy = 0.0_prec

 k3 = 0

 do i = 1,N_sp
  k = 3*(i-1)+1
  do k1 = 1,3
   k2 = k + (k1 - 1)
   
   do j = 1,k-1
    k3 = k3 + 1
!----------------------------------------------------------------------!  
    rijx(k3) = posixy(k2,1)-posixy(j,1)
    rijy(k3) = posixy(k2,2)-posixy(j,2)
    distij(k3) = 1.0_prec/sqrt(rijx(k3)**2 + rijy(k3)**2)
   
    rijx(k3) = rijx(k3)*distij(k3)
    rijy(k3) = rijy(k3)*distij(k3)
    distij(k3) = distij(k3)**3
!----------------------------------------------------------------------!
   end do
   
   do j = k+3,N_st
    k3 = k3 + 1
!----------------------------------------------------------------------!  
    rijx(k3) = posixy(k2,1)-posixy(j,1)
    rijy(k3) = posixy(k2,2)-posixy(j,2)
    distij(k3) = 1.0_prec/sqrt(rijx(k3)**2 + rijy(k3)**2)
   
    rijx(k3) = rijx(k3)*distij(k3)
    rijy(k3) = rijy(k3)*distij(k3)
    distij(k3) = distij(k3)**3
!----------------------------------------------------------------------!
   end do
   
  end do
  
 end do
 
 print*, 'Inicio ok'
 
end subroutine rede1

subroutine rede()
 use parametros
 use variaveis
 implicit none
 integer :: i,j,k,k1,k2,k3,kanis
 real(prec) :: delx,dely
 
 allocate(Sx(N_st),Sy(N_st),Sz(N_st),Snx(N_st),Sny(N_st),Snz(N_st),posixy(N_st,2))
 allocate(anisx(N_sp),anisy(N_sp))
          
 posixy = 0.0_prec
 Sx = 0.0_prec; Sy = 0.0_prec; Sz = 0.0_prec; Snx = 0.0_prec; Sny = 0.0_prec; Snz = 0.0_prec;  !====!
 anisx = 0.0_prec; anisy = 0.0_prec  !====!

 k = -2
 kanis = 0
 do j = 1,P
  do i = 1,P
   if (mod(i+j,2).ne.0) then
    k = k + 3
    k1 = k + 1
    k2 = k + 2
    kanis = kanis + 1
    
    posixy(k,1) = 0.5_prec*iA*(1.0_prec*i - L2)
    posixy(k,2) = 0.5_prec*iA*(1.0_prec*j - L2)
    
    if (mod(j,2).ne.0) then
     posixy(k1,1) = posixy(k,1) - alpha
     posixy(k1,2) = posixy(k,2)
     posixy(k2,1) = posixy(k,1) + alpha
     posixy(k2,2) = posixy(k,2)
     anisx(kanis) = 1.0_prec
     Sx(k) = iSp
     Sx(k1) = iSl
     Sx(k2) = iSl
    else
     posixy(k1,1) = posixy(k,1) 
     posixy(k1,2) = posixy(k,2) - alpha
     posixy(k2,1) = posixy(k,1) 
     posixy(k2,2) = posixy(k,2) + alpha
     anisy(kanis) = 1.0_prec
     Sy(k) = iSp
     Sy(k1) = iSl
     Sy(k2) = iSl
    end if
   end if
  end do
 end do
 
 delx = maxval(posixy(:,1)) - minval(posixy(:,1))
 dely = maxval(posixy(:,2)) - minval(posixy(:,2))
 
 posixy(:,1) = posixy(:,1) - 0.5_prec*delx - minval(posixy(:,1))
 posixy(:,2) = posixy(:,2) - 0.5_prec*dely - minval(posixy(:,2))
 
 !write(*,'(1x,A,I5,A,I5)') 'Rede iniciada com ',N_sp,' ilhas e um total de ',N_st,' Spins'

 allocate(rijx(1:N_st*(N_st-3)),rijy(1:N_st*(N_st-3)),distij(1:N_st*(N_st-3)))

 distij = 0.0_prec; rijx = 0.0_prec; rijy = 0.0_prec

 k3 = 0

 do i = 1,N_sp
  k = 3*(i-1)+1
  do k1 = 1,3
   k2 = k + (k1 - 1)
   
   do j = 1,k-1
    k3 = k3 + 1
!----------------------------------------------------------------------!  
    rijx(k3) = posixy(k2,1)-posixy(j,1)
    rijy(k3) = posixy(k2,2)-posixy(j,2)
    distij(k3) = 1.0_prec/sqrt(rijx(k3)**2 + rijy(k3)**2)
   
    rijx(k3) = rijx(k3)*distij(k3)
    rijy(k3) = rijy(k3)*distij(k3)
    distij(k3) = distij(k3)**3
!----------------------------------------------------------------------!
   end do
   
   do j = k+3,N_st
    k3 = k3 + 1
!----------------------------------------------------------------------!  
    rijx(k3) = posixy(k2,1)-posixy(j,1)
    rijy(k3) = posixy(k2,2)-posixy(j,2)
    distij(k3) = 1.0_prec/sqrt(rijx(k3)**2 + rijy(k3)**2)
   
    rijx(k3) = rijx(k3)*distij(k3)
    rijy(k3) = rijy(k3)*distij(k3)
    distij(k3) = distij(k3)**3
!----------------------------------------------------------------------!
   end do
   
  end do
  
 end do

  
end subroutine rede 

subroutine rk4(Ndim,x,h,yx,yy,yz,dydx,dydy,dydz,youtx,youty,youtz,derivs)
 use parametros, only : prec
 implicit none
 integer(4),intent(in) :: Ndim
 real(prec),intent(in) :: x,h
 real(prec),dimension(Ndim),intent(in) :: yx,yy,yz,dydx,dydy,dydz
 real(prec),dimension(Ndim),intent(out) :: youtx,youty,youtz

 real(prec),dimension(Ndim) :: dymx,dytx,ytx
 real(prec),dimension(Ndim) :: dymy,dyty,yty
 real(prec),dimension(Ndim) :: dymz,dytz,ytz
 real(prec) :: h6,hh,xh
 
 hh = 0.5_prec*h
 h6 = h/6.0_prec
 xh = x + hh
 
 ytx(:) = yx(:) + hh*dydx(:)
 yty(:) = yy(:) + hh*dydy(:)
 ytz(:) = yz(:) + hh*dydz(:)
 call derivs(Ndim,xh,ytx,yty,ytz,dytx,dyty,dytz) 
 
 ytx(:) = yx(:) + hh*dytx(:)
 yty(:) = yy(:) + hh*dyty(:)
 ytz(:) = yz(:) + hh*dytz(:)
 call derivs(Ndim,xh,ytx,yty,ytz,dymx,dymy,dymz)
 
 ytx(:) = yx(:) + h*dymx(:)
 yty(:) = yy(:) + h*dymy(:)
 ytz(:) = yz(:) + h*dymz(:)
 dymx(:) = dytx(:) + dymx(:)
 dymy(:) = dyty(:) + dymy(:)
 dymz(:) = dytz(:) + dymz(:)
 call derivs(Ndim,x+h,ytx,yty,ytz,dytx,dyty,dytz)
 
 youtx(:) = yx(:) + h6*(dydx(:) + dytx(:) + 2.0_prec*dymx(:))
 youty(:) = yy(:) + h6*(dydy(:) + dyty(:) + 2.0_prec*dymy(:))
 youtz(:) = yz(:) + h6*(dydz(:) + dytz(:) + 2.0_prec*dymz(:))
 
 return
end subroutine rk4

subroutine Adams_Bashforth_Moulton(Ndim,x,h,yx,yy,yz,youtx,youty,youtz,derivs)
 use parametros
 use ABM
 implicit none
 integer(4),intent(in) :: Ndim
 real(prec),intent(in) :: x,h
 real(prec),dimension(Ndim),intent(in) :: yx,yy,yz
 real(prec),dimension(Ndim),intent(out) :: youtx,youty,youtz
 
 real(prec),dimension(Ndim) :: Predx,Corrx,pasx
 real(prec),dimension(Ndim) :: Predy,Corry,pasy
 real(prec),dimension(Ndim) :: Predz,Corrz,pasz
 !real(prec),dimension(1:3) :: tt
 
 !---------------------------------------------------------------------!
 pasx = 55.0_prec*F3x - 59.0_prec*F2x + 37.0_prec*F1x - 9.0_prec*F0x
 pasy = 55.0_prec*F3y - 59.0_prec*F2y + 37.0_prec*F1y - 9.0_prec*F0y
 pasz = 55.0_prec*F3z - 59.0_prec*F2z + 37.0_prec*F1z - 9.0_prec*F0z
 Predx = yx + h24*pasx
 Predy = yy + h24*pasy
 Predz = yz + h24*pasz
 !---------------------------------------------------------------------!
 
 call derivs(Ndim,x+h,Predx,Predy,Predz,F4x,F4y,F4z)
 
 !---------------------------------------------------------------------!
 pasx = 9.0_prec*F4x + 19.0_prec*F3x - 5.0_prec*F2x + F1x
 pasy = 9.0_prec*F4y + 19.0_prec*F3y - 5.0_prec*F2y + F1y
 pasz = 9.0_prec*F4z + 19.0_prec*F3z - 5.0_prec*F2z + F1z
 Corrx = yx + h24*pasx
 Corry = yy + h24*pasy
 Corrz = yz + h24*pasz
 !---------------------------------------------------------------------!
 
 F0x = F1x ; F0y = F1y ; F0z = F1z 
 F1x = F2x ; F1y = F2y ; F1z = F2z
 F2x = F3x ; F2y = F3y ; F2z = F3z
 Erro_passo = 0.0_prec

 !do i = 1,Ndim
 ! tt = (abs((Corrector(i,:) - Predictor(i,:))/Corrector(i,:))) 
 ! if (maxval(tt) .gt. Erro_passo) Erro_passo = maxval(tt)
 !end do
 
 !Erro_passo = 19.0_prec*Erro_passo/270.0_prec
 !yout = Corrector
 !Erro_Total = Erro_Total + Erro_passo
 
 youtx = Corrx
 youty = Corry
 youtz = Corrz
 
 !!!!call derivs(Ndim,x+h,Corrector,F3)!!!!
 
 return
 
end subroutine Adams_Bashforth_Moulton

subroutine inicia_rkABM()
 use parametros
 use variaveis
 use ABM
 implicit none

 h24 = h/24.0_prec 
 
 allocate(EffFieldx(N_st)); allocate(EffFieldy(N_st)); allocate(EffFieldz(N_st))
 allocate(F0x(N_st),F1x(N_st),F2x(N_st),F3x(N_st),F4x(N_st))
 allocate(F0y(N_st),F1y(N_st),F2y(N_st),F3y(N_st),F4y(N_st))
 allocate(F0z(N_st),F1z(N_st),F2z(N_st),F3z(N_st),F4z(N_st))
 
 EffFieldx = 0.0_prec; EffFieldy = 0.0_prec; EffFieldz = 0.0_prec
 F0x = 0.0_prec; F1x = 0.0_prec; F2x = 0.0_prec; F3x = 0.0_prec; F4x = 0.0_prec;
 F0y = 0.0_prec; F1y = 0.0_prec; F2y = 0.0_prec; F3y = 0.0_prec; F4y = 0.0_prec;
 F0z = 0.0_prec; F1z = 0.0_prec; F2z = 0.0_prec; F3z = 0.0_prec; F4z = 0.0_prec;

 !!-------------------------------------------------------------------!!
 !! Inicia os três primeiros passos utilizando o método de Runge-Kutta!!
 !!-------------------------------------------------------------------!!
 call derivs(N_st,time,Sx,Sy,Sz,F0x,F0y,F0z)
 call rk4(N_st,time,h,Sx,Sy,Sz,F0x,F0y,F0z,Snx,Sny,Snz,derivs)
 Sx = Snx; Sy = Sny; Sz = Snz
 
 call derivs(N_st,time,Sx,Sy,Sz,F1x,F1y,F1z)
 call rk4(N_st,time,h,Sx,Sy,Sz,F1x,F1y,F1z,Snx,Sny,Snz,derivs)
 Sx = Snx; Sy = Sny; Sz = Snz

 call derivs(N_st,time,Sx,Sy,Sz,F2x,F2y,F2z)
 call rk4(N_st,time,h,Sx,Sy,Sz,F2x,F2y,F2z,Snx,Sny,Snz,derivs)
 Sx = Snx; Sy = Sny; Sz = Snz
 
 !!-------------------------------------------------------------------!!
 !! Inicia os três primeiros passos utilizando o método de Runge-Kutta!!
 !!-------------------------------------------------------------------!!
 
end subroutine inicia_rkABM

subroutine derivs(N,x,yx,yy,yz,dydx,dydy,dydz)
 use parametros
 use variaveis, only : damp,idamp,iSp,iSl,N_sp
 use ABM
 implicit none
 integer,intent(in) :: N
 real(prec),intent(in) :: x
 real(prec),dimension(N),intent(in) :: yx,yy,yz
 real(prec),dimension(N),intent(out) :: dydx,dydy,dydz
 
 integer :: i,k,k1,k2
 real(prec),dimension(N) :: torquex,torquey,torquez,gammx,gammy,gammz

! dydx = 0.0_prec ; dydy = 0.0_prec ; dydz = 0.0_prec
! torquex = 0.0_prec ; torquey = 0.0_prec ; torquez = 0.0_prec
! gammx = 0.0_prec ; gammy = 0.0_prec ; gammz = 0.0_prec
 
 call Campo_Efetivo()
 
! do i = 1,N
  torquex = yy*EffFieldz - yz*EffFieldy
  torquey = yz*EffFieldx - yx*EffFieldz
  torquez = yx*EffFieldy - yy*EffFieldx
! end do
 
 !Max_Torque = MaxVal((torquex(:)**2 + torquey(:)**2 + torquez(:)**2))
 
! do i = 1,N
  gammx = yy*torquez - yz*torquey
  gammy = yz*torquex - yx*torquez
  gammz = yx*torquey - yy*torquex
! end do

 do i = 1,N_sp
  k = 3*(i-1) + 1
  k1 = k + 1
  k2 = k + 2
  
  dydx(k) = -idamp*(torquex(k) + (damp/iSp)*gammx(k))
  dydx(k1) = -idamp*(torquex(k1) + (damp/iSl)*gammx(k1))
  dydx(k2) = -idamp*(torquex(k2) + (damp/iSl)*gammx(k2))
  
  dydy(k) = -idamp*(torquey(k) + (damp/iSp)*gammy(k))
  dydy(k1) = -idamp*(torquey(k1) + (damp/iSl)*gammy(k1))
  dydy(k2) = -idamp*(torquey(k2) + (damp/iSl)*gammy(k2))
  
  dydz(k) = -idamp*(torquez(k) + (damp/iSp)*gammz(k))
  dydz(k1) = -idamp*(torquez(k1) + (damp/iSl)*gammz(k1))
  dydz(k2) = -idamp*(torquez(k2) + (damp/iSl)*gammz(k2))

 end do
 
 return

end subroutine derivs 

subroutine derivs1(N,x,yx,yy,yz,dydx,dydy,dydz)
 use parametros
 use variaveis, only : damp,idamp,iSp,iSl,N_sp
 use ABM
 implicit none
 integer,intent(in) :: N
 real(prec),intent(in) :: x
 real(prec),dimension(N),intent(in) :: yx,yy,yz
 real(prec),dimension(N),intent(out) :: dydx,dydy,dydz
 
 integer :: i,k,k1,k2
 real(prec),dimension(N) :: torquex,torquey,torquez,gammx,gammy,gammz

 call Campo_Efetivo()
 
 !EffFieldx = EffFieldx + Bmax*sin(omega*x)
 !EffFieldy = EffFieldy + Bmax*sin(omega*x)
 !EffFieldz = EffFieldz + Bmax*sin(omega*x)
 
 torquex = yy*EffFieldz - yz*EffFieldy
 torquey = yz*EffFieldx - yx*EffFieldz
 torquez = yx*EffFieldy - yy*EffFieldx

 gammx = yy*torquez - yz*torquey
 gammy = yz*torquex - yx*torquez
 gammz = yx*torquey - yy*torquex

 do i = 1,N_sp
  k = 3*(i-1) + 1
  k1 = k + 1
  k2 = k + 2
  
  dydx(k) = -idamp*(torquex(k) + (damp/iSp)*gammx(k))
  dydx(k1) = -idamp*(torquex(k1) + (damp/iSl)*gammx(k1))
  dydx(k2) = -idamp*(torquex(k2) + (damp/iSl)*gammx(k2))
  
  dydy(k) = -idamp*(torquey(k) + (damp/iSp)*gammy(k))
  dydy(k1) = -idamp*(torquey(k1) + (damp/iSl)*gammy(k1))
  dydy(k2) = -idamp*(torquey(k2) + (damp/iSl)*gammy(k2))
  
  dydz(k) = -idamp*(torquez(k) + (damp/iSp)*gammz(k))
  dydz(k1) = -idamp*(torquez(k1) + (damp/iSl)*gammz(k1))
  dydz(k2) = -idamp*(torquez(k2) + (damp/iSl)*gammz(k2))

 end do
 
 return

end subroutine derivs1 

subroutine pulso()
!======================================================================!
! Aplica um pulso de campo externo no sistema em um determinado angulo !
!																	   !
! Os dados nesta parte não são coletados							   !
!======================================================================!
 use parametros
 use variaveis!, only : h,N_pulso,angle,B_pulso,N_st,h,Spin,Spin_new,posixy,Bextx,Bexty,Bextz,iuni_pul
 use ABM!, only : F3
 integer :: i,j

 do i = 1,N_pulso
  Bextx = B_pulso(i)*cos(angle)
  Bexty = B_pulso(i)*sin(angle)
  Bextz = 0.0_prec
  call derivs(N_st,i*h+h,Sx,Sy,Sz,F3x,F3y,F3z)
  call Adams_Bashforth_Moulton(N_st,i*h,h,Sx,Sy,Sz,Snx,Sny,Snz,derivs)
  Sx = Snx
  Sy = Sny
  Sz = Snz
  !write(iuni_pul,'(1x,I4)') N_st
  !write(iuni_pul,*) ' '
  !do j = 1,N_st
  ! write(iuni_pul,'(I5,5(1x,f15.7))') j,posixy(j,:),Spin(j,:)
  !end do
 end do
 !flush(13)
 
 Bextx = 0.0_prec
 Bexty = 0.0_prec
 Bextz = 0.0_prec

 return
end subroutine pulso 

subroutine Campo_Efetivo()
 use parametros
 use variaveis
 use ABM
 implicit none
 integer :: i,j,k,k1,k2,k3
 real(prec) :: Dip1
 
 EffFieldx = 0.0_prec
 EffFieldy = 0.0_prec
 EffFieldz = 0.0_prec

 do i = 1,N_sp
  k = 3*(i-1) + 1
  k1 = k + 1
  k2 = k + 2
  
!======================================================================!
!==========================Campo Troca=================================!
!======================================================================!

  EffFieldx(k) = EffFieldx(k) + iTroca*(Sx(k1) + Sx(k2))
  EffFieldy(k) = EffFieldy(k) + iTroca*(Sy(k1) + Sy(k2))
  EffFieldz(k) = EffFieldz(k) + iTroca*(Sz(k1) + Sz(k2))
  
  EffFieldx(k1) = EffFieldx(k1) + iTroca*Sx(k)
  EffFieldy(k1) = EffFieldy(k1) + iTroca*Sy(k)
  EffFieldz(k1) = EffFieldz(k1) + iTroca*Sz(k)
  
  EffFieldx(k2) = EffFieldx(k2) + iTroca*Sx(k)
  EffFieldy(k2) = EffFieldy(k2) + iTroca*Sy(k)
  EffFieldz(k2) = EffFieldz(k2) + iTroca*Sz(k)
  
!======================================================================!
!==========================Campo Anisotropia===========================!
!======================================================================!

  EffFieldx(k) = EffFieldx(k) + 2.0_prec*iAnis*Sx(k)*anisx(i)
  EffFieldy(k) = EffFieldy(k) + 2.0_prec*iAnis*Sy(k)*anisy(i)
  
!======================================================================!
!==========================Campo Externo===============================!
!======================================================================!

  EffFieldx(k) = EffFieldx(k) + iZee*Bextx
  EffFieldy(k) = EffFieldy(k) + iZee*Bexty
  EffFieldz(k) = EffFieldz(k) + iZee*Bextz
  
  EffFieldx(k1) = EffFieldx(k1) + iZee*Bextx
  EffFieldy(k1) = EffFieldy(k1) + iZee*Bexty
  EffFieldz(k1) = EffFieldz(k1) + iZee*Bextz
  
  EffFieldx(k2) = EffFieldx(k2) + iZee*Bextx
  EffFieldy(k2) = EffFieldy(k2) + iZee*Bexty
  EffFieldz(k2) = EffFieldz(k2) + iZee*Bextz
  
 end do
 
!======================================================================!
!======================================================================!
!======================================================================!

!======================================================================!
!==========================Campo Dipolar===============================!
!======================================================================!
 
 !if (iDip .ne. 0.0_prec) then
  k3 = 0
  do i = 1,N_sp
   k = 3*(i-1) + 1
   do k1 = 1,3
    k2 = k + (k1-1)
    
    do j = 1,k-1
   
     k3 = k3 + 1
   
     Dip1 = 3.0_prec*(Sx(j)*rijx(k3) + Sy(j)*rijy(k3))

     EffFieldx(k2) = EffFieldx(k2) + iDip*(Dip1*rijx(k3) - Sx(j))*distij(k3)
     EffFieldy(k2) = EffFieldy(k2) + iDip*(Dip1*rijy(k3) - Sy(j))*distij(k3)
     EffFieldz(k2) = EffFieldz(k2) - iDip*Sz(j)*distij(k3)

    end do

    do j = k+3,N_st
   
     k3 = k3 + 1
   
     Dip1 = 3.0_prec*(Sx(j)*rijx(k3) + Sy(j)*rijy(k3))

     EffFieldx(k2) = EffFieldx(k2) + iDip*(Dip1*rijx(k3) - Sx(j))*distij(k3)
     EffFieldy(k2) = EffFieldy(k2) + iDip*(Dip1*rijy(k3) - Sy(j))*distij(k3)
     EffFieldz(k2) = EffFieldz(k2) - iDip*Sz(j)*distij(k3)

    end do
   
   end do
  
  end do
 !end if
!======================================================================!
!======================================================================!
!======================================================================!

end subroutine Campo_Efetivo

subroutine Energia()
 use parametros
 use variaveis
 implicit none
 integer :: i,j,k,k1,k2,k3
 real(prec) :: E_troca, E_Aniso, E_Dip, E_Zee, E_total, Dip1
 
 E_troca = 0.0_prec
 E_Aniso = 0.0_prec
 E_Zee = 0.0_prec
 E_Dip = 0.0_prec
 
 do i = 1,N_sp
 
  k = 3*(i-1) + 1
  k1 = k + 1
  k2 = k + 2
  
!======================================================================!
!==========================Energia Troca===============================!
!======================================================================!
  E_troca = E_troca - Sx(k)*(Sx(k1) + Sx(k2))
  E_troca = E_troca - Sy(k)*(Sy(k1) + Sy(k2))
  E_troca = E_troca - Sz(k)*(Sz(k1) + Sz(k2))
!======================================================================!
!==========================Energia Anisotropia=========================!
!======================================================================!
  E_Aniso = E_aniso - (Sx(k)*anisx(i))**2
  E_Aniso = E_aniso - (Sy(k)*anisy(i))**2 
!======================================================================!
!==========================Energia Zeeman==============================!
!======================================================================!
  E_Zee = E_Zee - ( Bextx*Sx(k) + Bexty*Sy(k) + Bextz*Sz(k) )
  E_Zee = E_Zee - ( Bextx*Sx(k1) + Bexty*Sy(k1) + Bextz*Sz(k1) )
  E_Zee = E_Zee - ( Bextx*Sx(k2) + Bexty*Sy(k2) + Bextz*Sz(k2) )
!======================================================================!
!==========================Energia Dipolar=============================!
!======================================================================!

 end do
 
! if (iDip .ne. 0.0_prec) then
  k3 = 0
  do i = 1,N_sp
   k = 3*(i-1) + 1
   do k1 = 1,3
    k2 = k + (k1-1)   
    do j = 1,k-1
     k3 = k3 + 1
     Dip1 = 3.0_prec*(Sx(j)*rijx(k3)  + Sy(j)*rijy(k3))*(Sx(k2)*rijx(k3) + Sy(k2)*rijy(k3))
     E_Dip = E_Dip - (Dip1 - (Sx(j)*Sx(k2) + Sy(j)*Sy(k2) + Sz(j)*Sz(k2)))*distij(k3)
    end do
  
    do j = k+3,N_st
     k3 = k3 + 1
     Dip1 = 3.0_prec*(Sx(j)*rijx(k3)  + Sy(j)*rijy(k3))*(Sx(k2)*rijx(k3) + Sy(k2)*rijy(k3))
     E_Dip = E_Dip - (Dip1 - (Sx(j)*Sx(k2) + Sy(j)*Sy(k2) + Sz(j)*Sz(k2)))*distij(k3)
    end do

   end do
 
  end do
! end if
 E_Troca = iTroca*E_troca
 E_Aniso = iAnis*E_Aniso
 E_Zee = iZee*E_Zee
 E_dip = 0.5_prec*iDip*E_Dip 
 
 E_total = E_Troca + E_Aniso + E_Zee + E_Dip
 
 write(iuni_e,'(7(1x,f15.7))') time, E_troca, E_Aniso, E_Zee, E_Dip, E_total
  
end subroutine Energia

subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8  ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

subroutine magnetizacao()
 use parametros
 use variaveis
 implicit none
 integer :: i,k,k1,k2
 real(prec) :: MxP, MyP, MzP, MxL, MyL, MzL, MxT, MyT, MzT
 
 MxP = 0.0_prec; MyP = 0.0_prec; MzP = 0.0_prec
 MxL = 0.0_prec; MyL = 0.0_prec; MzL = 0.0_prec
 MxT = 0.0_prec; MyT = 0.0_prec; MzT = 0.0_prec
 
 MxT = sum(Sx(:))
 MyT = sum(Sy(:))
 MzT = sum(Sz(:))
 
 do i = 1,N_sp
 
  k = 3*(i - 1) + 1
  k1 = k + 1
  k2 = k + 2
  
  MxP = MxP + Sx(k)
  MyP = MyP + Sy(k)
  MzP = MzP + Sz(k)
  
  MxL = MxL + Sx(k1) + Sx(k2)
  MyL = MyL + Sy(k1) + Sy(k2)
  MzL = MzL + Sz(k1) + Sz(k2)
  
 end do
 
 write(iuni_mag,'(4(1x,e15.7))')  time, Ms*MxT, Ms*MyT, Ms*MzT
 write(iuni_mgSp,'(4(1x,e15.7))') time, Ms*MxP, Ms*MyP, Ms*MzP
 write(iuni_mgSl,'(4(1x,e15.7))') time, Ms*MxL, Ms*MyL, Ms*MzL
 
end subroutine magnetizacao

subroutine histerese
 use variaveis
 use ABM
 implicit none
 integer :: i,j,k,N_hist,N_s
 real(prec),dimension(:),allocatable :: Bhx,Bhy,dataMagx,dataMagy
 real(prec) :: dBh,avex,varx,avey,vary
 character(60) :: arquihist
 
 N_hist = 50
 N_s = 1000
 dBh = 2.0_prec*Bmax/real(N_hist-1,prec)
 allocate(Bhx(2*N_hist),Bhy(2*N_hist)) 
 
 do i = 1,N_hist
  Bhx(i) = Bmax - dBh*(i-1)
  Bhy(i) = Bmax - dBh*(i-1)
 end do
 do i = N_hist,2*N_hist
  Bhx(i) = -Bmax + dBh*(i-N_hist)
  Bhy(i) = -Bmax + dBh*(i-N_hist)
 end do
 
 Bhx = Bhx*cos(angle)
 Bhy = Bhy*sin(angle)
 
 write(arquihist,'("histerese_rede_",I2.2,".dat")') Lxy
 open(887,file=arquihist)
 time = 0.0_prec
 !do k = 1,5
 
 allocate(dataMagx(N_s),dataMagy(N_s))
 
 print*,'Inicio da Histerese'
 
 do i = 1,2*N_hist
 
  Bextx = Bhx(i)
  Bexty = Bhy(i)
  Bextz = 0.0_prec
  
  print*, 'Campo:',i
  
  do j = 1,1000000
   call derivs(N_st,time+h,Sx,Sy,Sz,F3x,F3y,F3z)
   call Adams_Bashforth_Moulton(N_st,time,h,Sx,Sy,Sz,Snx,Sny,Snz,derivs)
   Sx = Snx
   Sy = Sny
   Sz = Snz
   time = time + h
  end do
  
  dataMagx = 0.0_prec
  dataMagy = 0.0_prec
   
  do j = 1,N_s
   call derivs(N_st,time+h,Sx,Sy,Sz,F3x,F3y,F3z)
   call Adams_Bashforth_Moulton(N_st,time,h,Sx,Sy,Sz,Snx,Sny,Snz,derivs)
   Sx = Snx
   Sy = Sny
   Sz = Snz
   time = time + h
   dataMagx(j) = sum(Sx)
   dataMagy(j) = sum(Sy)
  end do
   
  dataMagx = dataMagx*Ms
  dataMagy = dataMagy*Ms
 
  call avevar(N_s,dataMagx,avex,varx)
  call avevar(N_s,dataMagy,avey,vary)
  
  print*,'Parte 2:', i
 
  write(887,*) Bextx,Bexty,avex,varx,avey,vary
  call flush()
 end do

 !end do
 
 close(887)

end subroutine histerese

subroutine avevar(N,data1,ave,var)
 implicit none
 integer,intent(in) :: N
 real(8),dimension(N),intent(in) :: data1
 real(8),intent(out) :: ave,var
 real(8),dimension(N) :: s
 
 ave = sum(data1)/real(N,8)
 s = data1 - ave
 var =  dot_product(s,s)
 var = (var-sum(s)**2/real(N,8))/real(N-1,8)
 
 return
end subroutine avevar
 
program teste
 use variaveis
 use ABM
 implicit none
 integer :: i,j,k,ii,ierro
 real(prec) :: tini,tfin
 character(60) :: arqcomp
 logical :: direx
 
 call timestamp ( )
 
 call inicial
 call inicia_rkABM
 call histerese
 
 call timestamp ( )

end program teste

