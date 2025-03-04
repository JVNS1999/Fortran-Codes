program man
   implicit none

   integer, parameter:: spint = selected_int_kind(8) !Suporta de -(2^15-1) ate (2^15-1)
   integer, parameter:: dpint = selected_int_kind(10) !Suporta de -2^63 ate 2^63
   integer, parameter:: qpint = selected_int_kind(20) !Suporta de -2^127 ate 2^127

   integer, parameter:: spreal = selected_real_kind(6, 37) !Single precision
   integer, parameter:: dpreal = selected_real_kind(15, 307) !Double precision
   integer, parameter:: qpreal = selected_real_kind(33, 4931) !Quadruple precision

   !Condicoes para integracao
   real (kind = dpreal), parameter:: to = 0.000d+00
   real (kind = dpreal), parameter:: tf = 1.000d+02
   real (kind = dpreal), parameter:: y1o = 1.000d+00
   real (kind = dpreal), parameter:: y2o = 0.000d+00
   real (kind = dpreal), parameter:: y3o = 1.000d+00
   real (kind = dpreal), parameter:: h = 1.000d-03
   
   !Parametros do sistema
   real (kind = dpreal), parameter:: p = 3.000d+00
   real (kind = dpreal), parameter:: b = 1.000d+00
   real (kind = dpreal), parameter:: r = 25.000d+00

   integer (kind = dpint), parameter:: npass = (tf-to)/h
   real (kind = dpreal), parameter:: h2 = 5.000d-01*h
   real (kind = dpreal), parameter:: h6 = h/6.000d+00

   integer (kind = dpint):: ipass
   integer (kind = dpint):: dstep
   real (kind = dpreal):: t
   real (kind = dpreal), dimension(3):: y
   real (kind = dpreal), dimension(3):: yt
   real (kind = dpreal), dimension(3):: k1
   real (kind = dpreal), dimension(3):: k2
   real (kind = dpreal), dimension(3):: k3
   real (kind = dpreal), dimension(3):: k4
   real (kind = dpreal), dimension(0:npass,4):: date

   k1 = 0.000d+00
   k2 = 0.000d+00
   k3 = 0.000d+00
   k4 = 0.000d+00
   date = 0.000d+00

   y = [y1o, y2o, y3o]

   date(0,:) = [to, y1o, y2o, y3o]

   do ipass = 1,npass
      yt = y
      k1(:) = [-p*(yt(1)-yt(2)), -yt(2)+(r-1.000D+00-yt(3))*yt(1), -b*yt(3)+yt(1)*yt(2)]
      yt = y+h2*k1
      k2(:) = [-p*(yt(1)-yt(2)), -yt(2)+(r-1.000D+00-yt(3))*yt(1), -b*yt(3)+yt(1)*yt(2)]
      yt = y+h2*k2
      k3(:) = [-p*(yt(1)-yt(2)), -yt(2)+(r-1.000D+00-yt(3))*yt(1), -b*yt(3)+yt(1)*yt(2)]
      yt = y+h*k3
      k4(:) = [-p*(yt(1)-yt(2)), -yt(2)+(r-1.000D+00-yt(3))*yt(1), -b*yt(3)+yt(1)*yt(2)]

      t = t+h
      y = y+h6*[k1+2.000d0*k2+2.000d0*k3+k4]

      date(ipass, :) = [t, y(1), y(2), y(3)] 

   end do
   
   open(unit = 1, file = './q3_resultados.dat', status = 'unknown', action = 'write')

   if (npass .le. 1000) then
      do ipass = 0,npass
         write(1,'(4ES23.15E2)') date(ipass,:)
      end do
   else
      dstep = npass/1000
      do ipass = 0,npass,dstep
         write(1,'(4ES23.15E2)') date(ipass,:)
      end do
   end if

   close(1)
end program man