module variaveis
    implicit none
    integer :: Nx,N,seed,Nterma,Nmc,Ntemp
    integer,dimension(:),allocatable :: Viz
    integer,dimension(:),allocatable :: S
    real(8),dimension(:),allocatable :: x,y
    real(8),dimension(:),allocatable :: E,M
    real(8),dimension(:),allocatable :: Temp,Beta0
    real(8) :: Etot,Mtot
    real(8) :: Beta
    real(8) :: E1,E2,Cv
    real(8) :: M1,M2,Sup

end module variaveis

subroutine inicial()
    use variaveis, only : seed,Nx,N,S,Viz,x,y,Etot,Mtot
    implicit none
    integer :: i,j,k,ni,nj
    real(8) :: dx,dy
    
    print*, 'Entre com Nx:'
    read(*,*) Nx
    print*, 'Entre com a seed:'
    read(*,*) seed

    N = Nx**2
    allocate(S(N),Viz(4*N),x(N),y(N))
    call srand(seed)

    S = 1
    do i = 1,N
        if (rand() < 0.5d0) then
            S(i) = -1
        end if
    end do

    k = 0
    do j = 1,Nx
        do i = 1,Nx
            k = k + 1
            x(k) = 1.0d0*(i-1)
            y(k) = 1.0d0*(j-1)
        end do
    end do

    k = 0
    do i = 1,N
        do j = 1,N
            do ni = -1,1
                do nj = -1,1
                    dx = x(i) - x(j) + real(ni*Nx)
                    dy = y(i) - y(j) + real(nj*Nx)
                    if (sqrt(dx**2 + dy**2) < 1.1d0 .and. sqrt(dx**2 + dy**2) > 0.d0) then
                        k = k + 1
                        Viz(k) = j
                    end if
                end do
            end do
        end do
    end do
    Etot = 0.0d0
    Mtot = 0.0d0
    Mtot = 1.0d0*sum(S(:))
    k = 0

    do i = 1,N
        do j = (i-1)*4 + 1,(i-1)*4 + 4
            Etot = Etot - 1.0d0*S(i)*S(Viz(j))
        end do
    end do
    Etot = 0.5d0*Etot

    return
end subroutine inicial

subroutine Metropolis()
    use variaveis, only : N,S,Viz,Etot,Mtot,Beta
    implicit none
    integer :: i,j,k
    real(8) :: dE

    do i = 1,N
        j = int(N*rand()) + 1
        dE = 0.0d0
        do k = (j-1)*4 + 1,(j-1)*4 + 4
            dE = dE + S(Viz(k))
        end do
        dE = 2.0d0*S(j)*dE
        if (rand() < exp(-dE*Beta)) then
            S(j) = -S(j)
            Etot = Etot + dE
            Mtot = Mtot - 2.0d0*S(j)
        end if
    end do
    return
end subroutine Metropolis

subroutine Termaliza()
    use variaveis, only : Nterma
    implicit none
    integer :: i
    do i = 1,Nterma
        call Metropolis
    end do
    return
end subroutine Termaliza

subroutine Passos()
    use variaveis, only : Nmc,Etot,Mtot,E,M
    implicit none
    integer :: i

    E = 0.0d0
    M = 0.0d0   ! M = np.null(Nmc)
    do i = 1,Nmc
        call Metropolis
        E(i) = Etot
        M(i) = Mtot
    end do

    return
end subroutine Passos

subroutine simulacao()
    use variaveis, only : N,Nmc,Nterma,Ntemp,Temp,beta0,beta,E,M,E1,E2,Cv,M1,M2,Sup
    implicit none
    integer :: i
    real(8) :: Ti,Tf,dT

    print*, 'Quantas Temperaturas:'
    read(*,*) Ntemp
    print*, 'Temperatura Inicial:'
    read(*,*) Ti
    print*, 'Temperaura Final:'
    read(*,*) Tf
    print*, 'Nmc?:'
    read(*,*) Nmc
    print*, 'Nterma?:'
    read(*,*) Nterma

    dT = (Tf - Ti)/Ntemp
    allocate(Temp(Ntemp),Beta0(Ntemp))
    do i = 1,Ntemp
        Temp(i) = Ti + (i-1)*dT       ! Temp = np.linspace(Ti,Tf,Ntemp) 
        Beta0(i) = 1.0d0/Temp(i)
    end do
    
    allocate(E(Nmc),M(Nmc))
    open(11,file='resultados.dat')
    do i = 1,Ntemp
        beta = beta0(i)
        call Termaliza
        call Passos

        call samples

        write(11,*) Temp(i),E1,sqrt(E2 - E1**2)/real(N,8),Cv,M1,sqrt(M2 - M1**2)/real(N,8),sup
    end do
    close(11)

    return
end subroutine simulacao

subroutine samples
    use variaveis, only : N,Nmc,E,M,E1,E2,M1,M2,Cv,Sup,beta
    implicit none
    integer :: i

    E1 = 0.0d0
    E2 = 0.0d0
    M1 = 0.0d0
    M2 = 0.0d0

    E1 = sum(1.0d0*E)/real(Nmc,8)
    M1 = sum(abs(1.0d0*M))/real(Nmc,8)
    do i = 1,Nmc
        E2 = E2 + E(i)**2
        M2 = M2 + M(i)**2
    end do
    E2 = E2/real(Nmc,8)
    M2 = M2/real(Nmc,8)

    Cv = (beta**2)*(E2 - E1**2)/real(N,8)
    Sup = beta*(M2 - M1**2)/real(N,8)

    E1 = E1/real(N,8)
    M1 = M1/real(N,8)
    E2 = E2/real(N,8)
    M2 = M2/real(N,8)

    return
end subroutine samples

program main
    use variaveis, only : N,S
    implicit none
    integer :: i
    call inicial

    open(10,file='config1.dat')
    do i = 1,N
        write(10,*) S(i)
    end do
    close(10)

    call simulacao

    open(10,file='config.dat')
    do i = 1,N
        write(10,*) S(i)
    end do
    close(10)

end program main