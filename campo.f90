subroutine campo(N,arqui,h,x,y,Bx,By)
 implicit none
 !Input: arqui > Nome do arquivo contendo dados das nanoilhas
 !Input: N > Grade N x N para cálculo do campo
 !Input: h > Altura da 'ponta' do MFM
 !Output: x,y,Bx,By > Posição x,y e campos Bx,By na grade N x N
 integer,intent(in) :: N
 character(60),intent(in) :: arqui
 real(8),intent(in) :: h
 real(8),intent(out) :: x(N),y(N),Bx(N,N),By(N,N)
 !f2py intent(in) :: N, arqui,h
 !f2py intent(out) :: x,y,Bx,By
 !---------------------------------------------------------------------!
 integer :: Nt
 real(8),dimension(:),allocatable :: x0,y0,sx,sy,sz
 integer :: i,j,k
 real(8) :: xmax,xmin,dx
 real(8) :: drx,dry,d3
 !---------------------------------------------------------------------!
 open(10,file=arqui)
 read(10,*) Nt
 read(10,*) 
 allocate(x0(Nt),y0(Nt),sx(Nt),sy(Nt),sz(Nt))
 do i = 1,Nt
  read(10,*) k,x0(i),y0(i),sx(i),sy(i),sz(i)
 end do
 close(10)
 !---------------------------------------------------------------------!
 xmax = maxval(x0)
 xmin = minval(x0)
 dx = 1.5d0*(xmax - xmin)/real(N,8)
 do i = 1,N
  x(i) = 1.5d0*xmin + i*dx
  y(i) = 1.5d0*xmin + i*dx
 end do
 !---------------------------------------------------------------------! 
 Bx = 0.0d0
 By = 0.0d0
 open(12,file='teste.dat')
 do i = 1,N
  do j = 1,N
   do k = 1,Nt
    drx = x(i) - x0(k)
    dry = y(j) - y0(k)
    d3 = 1.0d0/sqrt(drx*drx + dry*dry + h*h)
    drx = drx*d3
    dry = dry*d3
    d3 = d3*d3*d3
    Bx(i,j) = Bx(i,j) + (3.0d0*drx*(drx*sx(k)+dry*sy(k)+h*sz(k)) - sx(k))*d3
    By(i,j) = By(i,j) + (3.0d0*dry*(drx*sx(k)+dry*sy(k)+h*sz(k)) - sy(k))*d3
   end do
   write(12,*) x(i),y(j),Bx(i,j),By(i,j)
  end do
 end do
 close(12)
 
 return
end subroutine campo

subroutine difcampo(N,Bx1,By1,Bx2,By2,Bx,By)
 implicit none
 !Intent(in): N,Bx1,By1,Bx2,By2
 !Intent(out): Bx,By
 integer,intent(in) :: N
 real(8),dimension(N,N),intent(in) :: Bx1,By1,Bx2,By2
 real(8),dimension(N,N),intent(out) :: Bx,By
 !f2py intent(in) :: N,Bx1,By1,Bx2,By2
 !f2py intent(out) :: Bx,By
 
 Bx = Bx2 - Bx1
 By = By2 - By1
 
 return
end subroutine difCampo

subroutine campo1(N,arqui,h,x,y,Bx,By)
 implicit none
 !Input: arqui > Nome do arquivo contendo dados das nanoilhas
 !Input: N > Grade N x N para cálculo do campo
 !Input: h > Altura da 'ponta' do MFM
 !Output: x,y,Bx,By > Posição x,y e campos Bx,By na grade N x N
 integer,intent(in) :: N
 character(60),intent(in) :: arqui
 real(8),intent(in) :: h
 real(8),intent(out) :: x(N),y(N),Bx(N,N),By(N,N)
 !f2py intent(in) :: N, arqui,h
 !f2py intent(out) :: x,y,Bx,By
 !---------------------------------------------------------------------!
 integer :: Nt
 real(8),dimension(:),allocatable :: x0,y0,sx,sy,sz
 integer :: i,j,k
 real(8),dimension(N*N,2) :: rxy
 real(8),allocatable :: dxy(:),dr(:,:),Bx0(:),By0(:)
 real(8) :: xmax,xmin,dx
 real(8) :: dip1
 !---------------------------------------------------------------------!
 open(10,file=arqui)
 read(10,*) Nt
 read(10,*) 
 allocate(x0(Nt),y0(Nt),sx(Nt),sy(Nt),sz(Nt))
 do i = 1,Nt
  read(10,*) k,x0(i),y0(i),sx(i),sy(i),sz(i)
 end do
 close(10)
 !---------------------------------------------------------------------!
 xmin = minval(x0)
 xmax = maxval(x0)
 dx = (xmax-xmin)/real(N,8)
 k = 0
 do i = 1,N
  do j = 1,N
   k = k + 1
   rxy(k,1) = xmin + i*dx
   rxy(k,2) = xmin + j*dx
  end do
  x(i) = xmin + i*dx
  y(i) = xmin + i*dx
 end do
 print*, 'ok'
 allocate(dr(Nt*N*N,2),dxy(Nt*N*N),Bx0(N*N),By0(N*N))
 k = 0
 do j = 1,N*N
  do i = 1,Nt
   k = k + 1
   dr(k,1) = rxy(j,1)-x0(i)
   dr(k,2) = rxy(j,2)-y0(i)
   dxy(k) = 1.0d0/sqrt(dr(k,1)**2+dr(k,2)**2+h**2)
   dr(k,:) = dr(k,:)*dxy(k)
   dxy(k) = dxy(k)**3
  end do
 end do
 
 Bx0 = 0.0d0
 By0 = 0.0d0
 
 print*, 'ok'
 k = 0
 do i = 1,N*N
  do j = 1,Nt
   k = k + 1
   dip1 = 3.0d0*(sx(j)*dr(k,1)+sy(j)*dr(k,2)+h*sz(j))
   Bx0(i) = Bx0(i) + (dip1*dr(k,1) - sx(j))
   By0(i) = By0(i) + (dip1*dr(k,2) - sy(j))
  end do
 end do
 k = 0
 print*, 'ok'
 open(12,file='testeca.dat')
 do i = 1,N
  do j = 1,N
   k = k + 1
   Bx(i,j) = Bx0(k)
   By(i,j) = By0(k)
   write(12,*) x(i),y(j),Bx0(k),By0(k)
  end do
 end do
 close(12)
 print*, 'ok'
 return

end subroutine campo1

program teste
 implicit none
 
 integer,parameter :: N = 21
 character(60):: arqui
 real(8) :: x(N),y(N),Bx(N,N),By(N,N)
 arqui = 'conf0_10.xyz'
 
 call campo1(N,arqui,0.0d0,x,y,Bx,By)
 print*, 'ok'
end program teste

