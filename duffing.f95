MODULE duffing
    IMPLICIT NONE

    INTEGER                               ::  N, wno
    REAL(8)                               ::  h, Tf, a, b, k, f, w
    REAL(8), DIMENSION(:),   ALLOCATABLE  ::  t(:), wlist(:)
    REAL(8), DIMENSION(:,:), ALLOCATABLE  ::  y(:,:), yw(:,:)

    CONTAINS
    SUBROUTINE initw()
        ALLOCATE(wlist(wno))
    END SUBROUTINE initw

    SUBROUTINE init()
        N = CEILING(Tf/h)
        ALLOCATE (t(N),y(N,2), yw(N,wno))
    END SUBROUTINE init

    SUBROUTINE dealloc()
        t=0; y=0
        DEALLOCATE (t, y)
    END SUBROUTINE dealloc
    
    SUBROUTINE dealloc2()
        wlist=0; yw=0
        DEALLOCATE (wlist, yw)
    END SUBROUTINE dealloc2

    SUBROUTINE int_dffing()
        INTEGER                :: i
        REAL(8), DIMENSION(2)  :: k1, k2, k3, k4

        DO i=1,N-1
            k1 = h*f1(y(i,:))              
            k2 = h*f1(y(i,:) + k1/2) 
            k3 = h*f1(y(i,:) + k2/2) 
            k4 = h*f1(y(i,:) + k3)   

            y(i+1,:) = y(i,:) + (k1 + 2*k2 + 2*k3 + k4)/6
            t(i+1)   = t(i) + h          
        ENDDO

    END SUBROUTINE int_dffing

    SUBROUTINE int_forc_dffing()
        INTEGER                :: i, j
        REAL(8), DIMENSION(2)  :: k1, k2, k3, k4
        
        DO j = 1, wno
            w = wlist(j)
            DO i=1,N-1
                k1 = h*f2(t(i), y(i,:))              
                k2 = h*f2(t(i) + h/2, y(i,:) + k1/2) 
                k3 = h*f2(t(i) + h/2, y(i,:) + k2/2) 
                k4 = h*f2(t(i) + h  , y(i,:) + k3  )    

                y(i+1,:) = y(i,:) + (k1 + 2*k2 + 2*k3 + k4)/6
                t(i+1)   = t(i) + h            
            ENDDO
            yw(:,j) = y(:,1)
        ENDDO
    END SUBROUTINE int_forc_dffing

    FUNCTION f1(z)
        REAL(8), DIMENSION(2)  ::  f1, z
        
        f1(1) = z(2)
        f1(2) = -a*z(1) - b*(z(1))**3
        
        RETURN
    END FUNCTION

    FUNCTION f2(x,z)
        REAL(8)                ::  x
        REAL(8), DIMENSION(2)  ::  f2, z
        
        f2(1) = z(2)
        f2(2) = - k*z(2) - a*z(1) - b*(z(1))**3 + f*DCOS(w*x)

        RETURN
    END FUNCTION


END MODULE duffing