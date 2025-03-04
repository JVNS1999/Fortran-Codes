program SpinFlip
 implicit none
 integer::i,j,N
 real(8),dimension(3)::y,dydx,Bext,yout
 real(8)::x,h

 N = 1000000
 y = 1.0d0
 y = y/sqrt(y(1)**2 + y(2)**2 + y(3)**2)
 dydx = 0.0d0
 Bext = 0.0d0
 Bext(1) = 2.0d0
 h = 0.0001d0
 x = 0.0d0

 open(10,file='spin.xyz')
 do i = 1,N
  call campo(3,x,y,Bext,dydx)
  call rk4(3,y,dydx,Bext,x,h,yout,campo)
  y = yout
  if (mod(i,1000).eq.0.0d0) then
   write(10,*) 2
   write(10,*) ' '
   write(10,*) 0.0d0,0.0d0,y
   write(10,*) 0.0d0,0.0d0,0.5d0*Bext
  end if
  x = x + h
 end do
 Bext = -Bext
 N = 2*N
 do i = 1,N
  call campo(3,x,y,Bext,dydx)
  call rk4(3,y,dydx,Bext,x,h,yout,campo)
  y = yout
  if (mod(i,1000).eq.0.0d0) then
   write(10,*) 2
   write(10,*) ' '
   write(10,*) 0.0d0,0.0d0,y
   write(10,*) 0.0d0,0.0d0,0.5d0*Bext
  end if
  x = x + h
 end do
 
 close(10)

end program SpinFlip

subroutine campo(n,x,y,B,dydt)
 implicit none
 integer,intent(in)::n
 real(8),intent(in)::x
 real(8),dimension(n),intent(in)::y,B
 real(8),dimension(n),intent(out)::dydt
 real(8),dimension(n)::t,g
 real(8)::a,ia,modT

 t = 0.0d0
 g = 0.0d0
 a = 0.01d0
 ia = -1.0d0/(1.0d0 + a**2)

 t(1) = y(2)*B(3) - y(3)*B(2)
 t(2) = y(3)*B(1) - y(1)*B(3)
 t(3) = y(1)*B(2) - y(2)*B(1)

 modT = sqrt(t(1)**2 + t(2)**2 + t(3)**2)

 g(1) = y(2)*t(3) - y(3)*t(2)
 g(2) = y(3)*t(1) - y(1)*t(3)
 g(3) = y(1)*t(2) - y(2)*t(1)

 dydt = ia*(t + a*g)

 return
end subroutine campo 

subroutine rk4(n,y,dydx,B,x,h,yout,derivs)
 implicit none
 integer,intent(in)::n
 real(8),dimension(n),intent(in)::y,dydx,B
 real(8),intent(in)::x,h
 real(8),dimension(n),intent(out)::yout
 external::derivs
 real(8)::h6,hh,xh
 real(8),dimension(n)::dym,dyt,yt

 hh = 0.5d0*h
 h6 = h/6.0d0

 xh=x+hh
	yt=y+hh*dydx
	call derivs(n,xh,yt,B,dyt)
	
	yt=y+hh*dyt
	call derivs(n,xh,yt,B,dym)
	
	yt=y+h*dym
	dym=dyt+dym
	call derivs(n,x+h,yt,B,dyt)
	
	yout=y+h6*(dydx+dyt+2.0d0*dym)

 return
end subroutine rk4
 
