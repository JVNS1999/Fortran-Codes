Program Estudos
    Implicit None
    Real(8) :: A, B, C, root1, root2, discriminant

    Print*, 'Programa para resolver uma equação quadrática Ax^2 + Bx + C = 0.'
    Print*
    Write(*,'(a)', Advance = 'No') 'Enter A:' ! Advance não faz o cursor descer o resultado !
    Read(*,*) A
    Write(*,'(a)', Advance = 'No') 'Enter B:'
    Read(*,*) B
    Print*
    Write(*,'(a)', Advance = 'No') 'Enter C:'
    Read(*,*) C

    discriminant = B**2 - 4*A*C

    If(discriminant .ge. 0) Then ! Raízes Reais
        root1 = (-B + sqrt(discriminant))/(2.d0*A)
        root2 = (-B - sqrt(discriminant))/(2.d0*A)
        Print*, 'Root 1 =', root1
        Print*, 'Root 2 =', root2
    Else
    Print*, 'Discriminante =', discriminant
    Print*, 'Não é possível resolver'
    End If
    
End Program Estudos