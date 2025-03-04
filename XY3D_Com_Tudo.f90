
!--------------------------!=========!----------------------------------
!--------------------------! Modules !----------------------------------
!--------------------------!=========!----------------------------------

!=======================================================================
!	Global Parameters
!=======================================================================
MODULE Constants

IMPLICIT NONE
REAL*8, PARAMETER:: Pi=4.0*ATAN(1.0),ZPi=2.0*Pi
REAL*8           :: delta
REAL*8           :: Temper,Beta
REAL*8           :: temp_i,dTempP,dTempS
REAL*8 			 :: neg_Zbeta
INTEGER          :: L, LOver2, LL, LLOver2, invLL
INTEGER          :: NEquil,NBoltz,NWolf,NOver,NHistogram
INTEGER			 :: sim_s, tot_data

END MODULE Constants
!=======================================================================
!	 Define spin arrays and Lattice Indexes 
!=======================================================================
MODULE Arrays

IMPLICIT NONE
REAL*8, ALLOCATABLE  :: Sx(:),Sy(:),Sz(:)
INTEGER, ALLOCATABLE :: Iep(:),Iem(:)
INTEGER, ALLOCATABLE :: Jep(:),Jem(:),IJplaquete(:)
INTEGER, ALLOCATABLE :: IBlack(:),IWhite(:)

END MODULE Arrays
!----------------------------------------------------------------------
module outputs
    implicit none  
    real*8,allocatable :: temp_vec(:)
	real*8,allocatable :: MagP(:),Mag2P(:),Mag4P(:)
	real*8,allocatable :: MagZ(:),Mag2Z(:),Mag4Z(:)
	real*8,allocatable :: MagZstag(:),Mag2Zstag(:),Mag4Zstag(:)
    real*8,allocatable :: neg_en(:),en2(:)
  end module outputs
!----------------------------------------------------------------------
!=======================================================================
!--------------------! Histogram Parameters !---------------------------
!=======================================================================
MODULE H_Parameters

IMPLICIT NONE
REAL*8               :: E_i,E_f,deltaE
INTEGER              :: NH
INTEGER, ALLOCATABLE :: H(:)

END MODULE H_Parameters
!=======================================================================
MODULE MPI_arrays

implicit none
INCLUDE 'mpif.h'
integer::ierr,tot_proc,proc_id,mpistatus(MPI_STATUS_SIZE)

END MODULE MPI_arrays
!=======================================================================

!--------------------------!=============!------------------------------
!--------------------------! Subroutines !------------------------------
!--------------------------!=============!------------------------------

!=======================================================================
!	Read several parameters
!=======================================================================
SUBROUTINE LerDados
	USE Arrays
	USE Constants
    USE H_Parameters
    USE outputs
    USE MPI_arrays
	IMPLICIT NONE
	INTEGER:: ih, Nlat

	OPEN(1,FILE="Dados.dat")

		READ(1,*)L! Lattice size (L)

        READ(1,*)NLat!how many times the lattice will be swept
		NEquil = NLat*L**2
		READ(1,*)NBoltz,NOver,NWolf! Metropolis, Wolf, Overrrelaxation
		READ(1,*)NHistogram! Points to build Histogram
								
		READ(1,*)delta ! Anisotropy at z

		READ(1,*)temp_i,dTempP,sim_s !temp step

		READ(1,*)E_i,E_f,deltaE

	CLOSE(1)

	dTempS = dTempP/(1.0d0*sim_s)
	LL=L*L;LOver2=L/2; LLOver2=LL/2; invLL= 1.0/LL

	tot_data = sim_s*NHistogram
	
	ALLOCATE(Sx(LL),Sy(LL),Sz(LL))
	ALLOCATE(Iep(LL),Iem(LL),Jep(LL),Jem(LL),IJplaquete(LL))
    ALLOCATE(IBlack(LLOver2),IWhite(LLOver2))
	ALLOCATE(MagP(tot_data),Mag2P(tot_data),Mag4P(tot_data))
	ALLOCATE(MagZ(tot_data),Mag2Z(tot_data),Mag4Z(tot_data))
	ALLOCATE(MagZstag(tot_data),Mag2Zstag(tot_data),Mag4Zstag(tot_data))
	ALLOCATE(neg_en(tot_data),en2(tot_data))
	ALLOCATE(temp_vec(0:sim_s-1))

	NH = (E_f - E_i)/deltaE + 1
	
	ALLOCATE(H(NH))

    Loop_H: DO ih = 1 , NH
        H(ih) = 0
	END DO Loop_H

END SUBROUTINE LerDados
!=======================================================================
!	Set Random Seed defined by computer clock
!=======================================================================
SUBROUTINE init_random_seed
    USE MPI_arrays
	integer :: i, n, clock
	integer,dimension(:),allocatable :: iseed
   
	call random_seed(size = n)
	allocate(iseed(n))
	call system_clock(count=clock)
	iseed = clock + 37*(/ (i - 1, i = 1, n) /)
	call random_seed(put = iseed*(proc_id+7))
	deallocate(iseed)

END SUBROUTINE init_random_seed
!=======================================================================
!	Set neighbour table
!=======================================================================
SUBROUTINE Table()
	USE Arrays
	USE Constants
	IMPLICIT NONE
	INTEGER :: ISitio,Ik,k,i,j,iy,ix,ixp,ixm,iyp,iym,is,ilin

	ISitio = 0
	Ik = 1
	k = 0

    Loop_InitVetors: DO i = 1 , LL ! Initialize index vectors
    
		Iep(i) = 0
		Iem(i) = 0
		Jep(i) = 0
		Jem(i) = 0
		IJplaquete(i) = 0
	END DO Loop_InitVetors

    Loop_yDirection: DO iy = 1 , L

		iyp = iy + 1
		iym = iy - 1
		IF(iyp .GT. L)iyp = iyp - L
		IF(iym .LT. 1)iym = iym + L

        Loop_xDirection: DO ix = 1 , L

			ixp = ix + 1
			ixm = ix - 1

			IF(ixp .GT. L)ixp = ixp - L
			IF(ixm .LT. 1)ixm = ixm + L

			ISitio = ISitio + 1
			Iep(ISitio) = ixp + (iy - 1)*L
			Iem(ISitio) = ixm + (iy - 1)*L
			Jep(ISitio) = ix + (iyp - 1)*L
			Jem(ISitio) = ix + (iym - 1)*L
            IJPlaquete(ISitio) = ixp + (iyp - 1)*L
            
		END DO Loop_xDirection

		Ik = -Ik
		is = MOD(iy,2)
		ilin = (iy -1)*L + is

        LOOP_IBlackIWhite: DO ix = 1 , LOver2

			k = k + 1
			IBlack(k) = ilin + 2 * ix-1
            IWhite(k) = IBlack(k) + ik
            
		END DO LOOP_IBlackIWhite

	END DO Loop_yDirection

END SUBROUTINE Table
!=======================================================================
!	Set Initial Configuration
!=======================================================================
SUBROUTINE InitCond()
	USE Arrays
	USE Constants
	IMPLICIT NONE
	REAL*8, DIMENSION(LL) :: Z,Phi,Sxy
    REAL*8, DIMENSION(LL) :: Rz,Rphi
    real*8                :: norm
	INTEGER               :: i,site

    CALL Random_Number(Rz)
	CALL Random_Number(Rphi)

    DO i = 1 , LL
        Sz(i)    = 1.0D0 - 2.0d0*Rz(i)
        Phi(i)   = ZPi*Rphi(i)
    END DO

    DO i = 1 , LL
        Sxy(i)   = SQRT(1.0D0 - Sz(i)**2)
    END DO

    DO i = 1 , LL
        Sx(i) = Sxy(i)*COS(Phi(i))
        Sy(i) = Sxy(i)*SIN(Phi(i))
    END DO

END SUBROUTINE InitCond
!=======================================================================
!	Calculate Energy
!=======================================================================
SUBROUTINE Energy(Neg_E)
	USE Arrays
	USE Constants
	IMPLICIT NONE
	INTEGER :: i,i1,i2
	REAL*8  :: Neg_E

	Neg_E = 0.0
    LOOP_Lattice: DO i = 1 , LL

		i1 = iep(i) ; i2 = jep(i)
			Neg_E = Neg_E + Sx(i)*(Sx(i1)+ Sx(i2)) + &
				Sy(i)*(Sy(i1)+ Sy(i2)) + &
                delta*Sz(i)*(Sz(i1)+ Sz(i2))
                
	END DO LOOP_Lattice

END SUBROUTINE Energy
!=======================================================================
!	Calculate Magnetization
!=======================================================================
SUBROUTINE Magnetization(Mx,My,Mz,MxStag,MyStag,MzStag)
	USE Arrays
	USE Constants
	IMPLICIT NONE
	INTEGER :: i,i1,i2
	REAL*8  :: Mx,My,Mz
	Real*8  :: Mxstag,Mystag,Mzstag
    
	Mx = 0.0 ; My = 0.0 ; Mz = 0.0 ; MxStag = 0.0 ; MyStag = 0.0 ; MzStag = 0.0

    LOOP_Lattice: DO i = 1 , LLover2

		i1 = IBlack(i)
		i2 = IWhite(i)

		Mx = Mx + Sx(i1) + Sx(i2)
		My = My + Sy(i1) + Sy(i2)
		Mz = Mz + Sz(i1) + Sz(i2)

		MxStag = MxStag + Sx(i1) - Sx(i2)
		MyStag = MyStag + Sy(i1) - Sy(i2)
		MzStag = MzStag + Sz(i1) - Sz(i2)	

	END DO LOOP_Lattice

END SUBROUTINE Magnetization
!=======================================================================
!	Do Metropolis
!=======================================================================
SUBROUTINE Metropolis(neg_Eold)
	USE Arrays
	USE Constants
	IMPLICIT NONE

	REAL*8		      :: Z,Phi
	REAL*8		      :: Vx,Vy,Vz,neg_Eold,p,neg_DE
    REAL*8, DIMENSION(LL) :: Rz,Rphi,R
	REAL*8, DIMENSION(LL) :: Sxaux,Syaux,Szaux,Sxy
	INTEGER               :: i,isite

    CALL Random_number(R)
    CALL Random_Number(Rz)
    CALL Random_Number(Rphi)
    
    DO i = 1 , LL
        
        Z    = 1.0D0 - 2.0d0*Rz(i)
        Phi   = ZPi*Rphi(i)
        Sxy(i)   = SQRT(1.0D0 - Z*Z)
        Sxaux(i) = sxy(i)*COS(phi)
        Syaux(i) = sxy(i)*SIN(phi)
        Szaux(i) = Z

	END DO

    LOOP_BOLTZ: DO isite = 1 , LL

		Vx =Sx(iep(isite)) + Sx(iem(isite)) + &
		Sx(jep(isite)) + Sx(jem(isite))

		Vy =Sy(iep(isite)) + Sy(iem(isite)) + &
		Sy(jep(isite)) + Sy(jem(isite))

		Vz =Sz(iep(isite)) + Sz(iem(isite)) + &
		Sz(jep(isite)) + Sz(jem(isite))

		neg_DE =Vx*(Sxaux(isite) - Sx(isite)) + &
		Vy*(Syaux(isite) - Sy(isite)) + &
		delta*Vz*(Szaux(isite) - Sz(isite))

        IF(neg_DE.GT.0.0) THEN

            Sx(isite) = Sxaux(isite)
            Sy(isite) = Syaux(isite)
            Sz(isite) = Szaux(isite)
            neg_Eold = neg_Eold + neg_DE

         ELSE

            p = EXP(beta*neg_DE)

            IF(R(isite).LE.p)THEN
                Sx(isite) = Sxaux(isite)
                Sy(isite) = Syaux(isite)
                Sz(isite) = Szaux(isite)
                neg_Eold = neg_Eold + neg_DE
            ENDIF

        ENDIF

    END DO	LOOP_BOLTZ

END SUBROUTINE Metropolis
!=======================================================================
!	Do Overrelaxation (Around z)
!=======================================================================
SUBROUTINE Overrelaxation()
	USE Arrays
	USE Constants
	IMPLICIT NONE
	INTEGER	:: i,isite,xplus,xminus,yplus,yminus
	REAL*8  :: Vx,Vy,Vz,V,Vx2,Vy2
	REAL*8  :: A1,A2,X,Y
  
	! Azimuthal Rotation (Energy is conserved)
    LOOP_White: DO i = 1 , LLover2

		isite = iwhite(i)

		xplus = iep(isite)
		xminus = iem(isite)
		yplus = jep(isite)
		yminus = jem(isite)

		Vx = Sx(xplus) + Sx(xminus) + Sx(yplus) + Sx(yminus)
		Vy = Sy(xplus) + Sy(xminus) + Sy(yplus) + Sy(yminus)

		Vx2 = Vx*Vx ; Vy2 = Vy*Vy
		V = 1.0/(Vx2 + vy2)
		A1 = V*Vx2 - V*Vy2
		A2 = 2.0*V*Vx*Vy

		X = Sx(isite)
		Y = Sy(isite)

		Sx(isite) = X*A1 + Y*A2
		Sy(isite) = X*A2 - Y*A1
									
    END DO LOOP_White
    
    LOOP_Black:	DO i = 1 , LLover2

		isite = iblack(i)

		xplus = iep(isite)
		xminus = iem(isite)
		yplus = jep(isite)
		yminus = jem(isite)

		Vx = Sx(xplus) + Sx(xminus) + Sx(yplus) + Sx(yminus)
		Vy = Sy(xplus) + Sy(xminus) + Sy(yplus) + Sy(yminus)

		Vx2 = Vx*Vx ; Vy2 = Vy*Vy
		V = 1.0/(Vx2 + vy2)
		A1 = V*Vx2 - V*Vy2
		A2 = 2.0*V*Vx*Vy

		X = Sx(isite)
		Y = Sy(isite)

		Sx(isite) = X*A1 + Y*A2
		Sy(isite) = X*A2 - Y*A1

    END DO LOOP_Black

END SUBROUTINE Overrelaxation
!=======================================================================
!	Do Wolff (In XY plane)
!=======================================================================
SUBROUTINE Wolff()
	USE Arrays
	USE Constants
	IMPLICIT NONE
	INTEGER 		:: i,n,k,j
	REAL*8  		:: Phi,X,Y,S1,S2,DE
	REAL*8, DIMENSION(2)   	:: Rinit
	INTEGER,DIMENSION(LL) 	:: Isite,Iac
	REAL*8,DIMENSION(LL)  	:: Rxplus,Rxminus,Ryplus,Ryminus
	REAL*8,DIMENSION(LL)  	:: Rw
	

	CALL Random_Number(Rxplus)
	CALL Random_Number(Rxminus)
	CALL Random_Number(Ryplus)
	CALL Random_Number(Ryminus)
	CALL Random_Number(Rinit)

    LOOP_Init:	DO i = 1 , LL
        Isite(i) = 0
        Iac(i)  = 0
    END DO LOOP_Init

	Phi = Zpi*Rinit(1)			! Define Random direction
	X = COS(Phi)
	Y = SIN(Phi)
				
	! Chose start point
	n = 1 
	Isite(n) = LL*Rinit(2)+1
	Iac(Isite(n)) = 1

    LOOP_Lattice:	DO i = 1 , LL

		! Check end cluster
		! If end is reached, flip cluster
		IF(i .GT. n)THEN		

			DO j = 1 , LL
				IF(IAC(j) .NE. 0)THEN		
					Sx(j) = Sx(j) - 2.0*Rw(j)*X
					Sy(j) = Sy(j) - 2.0*Rw(j)*Y
				ENDIF
			END DO
			RETURN

		 ELSE

			k = Isite(i)
			S1 = X*Sx(k) + Y*Sy(k)
			Rw(k) = S1

			! +X bond
			IF(Iac(iep(k)) .EQ. 0)THEN
				S2 = X*Sx(iep(k)) + Y*Sy(iep(k))
				DE = EXP(neg_Zbeta*S1*S2)
				IF(Rxplus(i) .GT. DE)THEN
					n = n + 1
					isite(n) = iep(k)
					iac(isite(n)) = 1
				ENDIF
			ENDIF

			! +Y bond
			IF(Iac(jep(k)) .EQ. 0)THEN
				S2 = X*Sx(jep(k)) + Y*Sy(jep(k))
				DE = EXP(neg_Zbeta*S1*S2)
				IF(Ryplus(i) .GT. DE)THEN
					n = n + 1
					isite(n) = jep(k)
					iac(isite(n)) = 1
				ENDIF
			ENDIF


			! -X bond
			IF(Iac(iem(k)) .EQ. 0)THEN
				S2 = X*Sx(iem(k)) + Y*Sy(iem(k))
				DE = EXP(neg_Zbeta*S1*S2)
				IF(Rxminus(i) .GT. DE)THEN
					n = n + 1
					isite(n) = iem(k)
					iac(isite(n)) = 1
				ENDIF
			ENDIF

			! -Y bond
			IF(Iac(jem(k)) .EQ. 0)THEN
				S2 = X*Sx(jem(k)) + Y*Sy(jem(k))
				DE = EXP(neg_Zbeta*S1*S2)
				IF(Ryminus(i) .GT. DE)THEN
					n = n + 1
					isite(n) = jem(k)
					iac(isite(n)) = 1
				ENDIF
			ENDIF

        ENDIF
        
	END DO Loop_Lattice

END SUBROUTINE Wolff
!=======================================================================
!	Build Histogram
!=======================================================================
SUBROUTINE Histogram(E)
	USE H_Parameters
	IMPLICIT NONE
	INTEGER		:: i
	REAL*8		:: E

	i = (E-E_i)/deltaE
	H(i) = H(i) + 1
    
END SUBROUTINE Histogram
!-----------------------------------------------------------------------
subroutine get_acumulators(i_data,neg_E,Mx,My,Mz,MxStag,MyStag,MzStag)
    use constants
    use arrays
    use outputs
    implicit none
    real*8 :: neg_E
	REAL*8  :: Mx,My,Mz,MxStag,MyStag,MzStag
	real*8  :: M2Z, M2ZStag, M2P
    integer :: i_data
    
    !-------------------------------------
    !--- Get M M**2 M**4 E E**2 E**4 -----
    !-------------------------------------
	
	M2P = Mx*Mx + My*My
	Mag2P(i_data) = M2P ; MagP(i_data) = sqrt(M2P)
	MagZ(i_data) = Mz; M2Z = Mz*Mz 
	Mag2Z(i_data) = M2Z; Mag4Z(i_data) = M2Z*M2Z
	neg_en(i_data) = neg_E; en2(i_data) = neg_E*neg_E; Mag4P(i_data) = M2P*M2P
	M2ZStag = MzStag*MzStag
	MagZstag(i_data) = MzStag; Mag2Zstag(i_data) = M2ZStag; Mag4Zstag(i_data) = M2ZStag*M2ZStag
    
end subroutine get_acumulators
!-------------------------------------------------------------------------------------------
subroutine data_archive()
    use constants
    use arrays
    use outputs
    use MPI_arrays
    implicit none
    character*10 :: file_id
    character*50 :: data_name
    integer ::j,k,kont,itemp
    real*8,dimension(13) :: vec
    
    !-----------------------------
    !--- Create data container ---
    !-----------------------------
    
    write(file_id, '(I3)') L
    data_name = 'data_L=' // trim(adjustl(file_id)) // '.dat'
    open(1,file = trim(data_name))
    
    kont = 0
    itemp = 0

    do j =1,tot_data
      
      kont = kont + 1
      vec(1) = temp_vec(itemp); vec(2) = kont
      vec(3) = magP(j); vec(4) = Mag2P(j); vec(5) = abs(MagZ(j))
	  vec(6) = Mag2Z(j); vec(7) = Mag4Z(j)
	  vec(8) = -1.0*neg_en(j); vec(9) = en2(j);  vec(10) = Mag4P(j)
	  vec(11) = abs(MagZstag(j)); vec(12) = Mag2Zstag(j); vec(13) = Mag4Zstag(j)
    
      ! write(1,*)vec(1),int(vec(2)),vec(3),vec(4),vec(5),vec(6) ! Use if process is serial

        if(proc_id .eq. 0)then

            write(1,*)vec(1),int(vec(2)),vec(3),vec(4),vec(5),vec(6),vec(7),vec(8),vec(9),vec(10),vec(11),vec(12),vec(13)

            do k=1,tot_proc-1

                call MPI_Recv(vec,13,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
                write(1,*)vec(1),int(vec(2)),vec(3),vec(4),vec(5),vec(6),vec(7),vec(8),vec(9),vec(10),vec(11),vec(12),vec(13)

            enddo

         else

            call MPI_Send(vec,13,MPI_DOUBLE_PRECISION,0,proc_id,MPI_COMM_WORLD,ierr)

        endif

        IF(mod(kont,NHistogram).EQ.0)THEN
			kont = 0; itemp = itemp + 1
        END IF

    enddo
    
    close(1)

end subroutine data_archive
!------------------------------------------------------------------------------------------


!=======================================================================
!==========!--------------!=============================================
!==========! MAIN PROGRAM !=============================================
!==========!--------------!=============================================
!=======================================================================
PROGRAM XY
    USE Arrays
    USE outputs
	USE Constants
	USE H_Parameters
    USE MPI_arrays
	IMPLICIT NONE
    INTEGER :: ih,iboltz,iover,iwolf,iequil
    INTEGER :: i,j,site,i_data,iserial
    REAL*8  :: norm
	REAL*8  :: neg_E
	REAL*8  :: Mx,My,Mz,MxStag,MyStag,MzStag
    REAL 	:: t_init,t_end
    
    !====================
    !-	Initialize mpi -=
    !====================
    Call MPI_INIT(ierr)
    !==============================
    !--	Get tot_proc and proc_id -=
    !==============================
    Call MPI_COMM_SIZE(MPI_COMM_WORLD, tot_proc, ierr)
    Call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierr)

	CALL LerDados()						! Read parameters
	CALL cpu_time(t_init)				! Start Counting the simulation time 
	CALL init_random_seed()				! Get Random Seed
	CALL Table()						! Set neighbour table
    CALL InitCond()						! Set initial conditions
    
    i_data = 0
	Temper = temp_i + (proc_id+1)*dTempP
	
    LOOP_TEMP: DO iserial = 0,sim_s-1

        Temper = Temper - iserial*dTempS
        temp_vec(iserial)=Temper
        Beta = 1.0/Temper
        neg_Zbeta = -2.0*beta

        CALL Energy(neg_E)					! Initialize Energy
            
        LOOP_EQUILIBRIUM: DO iequil = 1 , nequil
		
			LOOP_OVERRELAXATION: DO iover = 1 , nover
				CALL Overrelaxation()				
			END DO LOOP_OVERRELAXATION

            LOOP_METROPOLIS: DO iboltz = 1 , nboltz
                CALL Metropolis(neg_E)
            END DO LOOP_METROPOLIS

            LOOP_WOLF: DO iwolf = 1 , nwolf
            	CALL Wolff()
            	CALL Energy(neg_E)
			END DO LOOP_Wolf
			
			! LOOP_METROPOLIS2: DO iboltz = 1 , nboltz
			! 	CALL Metropolis(neg_E)
			! END DO LOOP_METROPOLIS2

			! if (iserial .eq. sim_s - 1 .and. proc_id .eq. 0) then
			! 	CALL Magnetization(Mx,My,Mz)		! <----- Used to check convergence only
			! 	WRITE(2,*)iequil,sqrt(Mx*Mx+My*My)*invLL,abs(Mz)*invLL,-neg_E*invLL
			! endif

        END DO LOOP_Equilibrium

        ! Used to photograph lattice after thermalization
        ! DO j = 1,L
        !     DO i = 1,L
        !         site = i+L*(j-1)
        !         norm = sqrt(Sx(site)**2+Sy(site)**2+Sz(site)**2)  
        !         WRITE(4,*)i,j,0,Sx(site),Sy(site),Sz(site),norm
        !     ENDDO
        ! ENDDO

        LOOP_Histogram: DO iequil = 1 , nhistogram

			LOOP_OVerrelaxation2: DO iover = 1 , nover
				CALL Overrelaxation()				
			END DO LOOP_OVerrelaxation2

            LOOP_Metropolis3: DO iboltz = 1 , nboltz
                CALL Metropolis(neg_E)
            END DO LOOP_Metropolis3

            LOOP_Wolf2: DO iwolf = 1 , nwolf
            	CALL Wolff()
            	CALL Energy(neg_E)			! I guess this command is out of place
			END DO LOOP_Wolf2
			
			! LOOP_Metropolis4: DO iboltz = 1 , nboltz
			! 	CALL Metropolis(neg_E)
			! END DO LOOP_Metropolis4

            !IF(E .GT. E_i .and. E .LT. E_f) CALL Histogram(E)
            Call Magnetization(Mx,My,Mz,MxStag,MyStag,MzStag)
            i_data = i_data + 1
            CALL get_acumulators(i_data,neg_E,Mx,My,Mz,MxStag,MyStag,MzStag)

        END DO LOOP_Histogram

    END DO LOOP_TEMP
    
    CALL data_archive()

	! DO ih = 1 , NH
	! 	WRITE(4,*)deltaE*ih,H(ih)
	! END DO
    
    Call MPI_Finalize (ierr)
	CALL cpu_time(t_end)
	! WRITE(*,*)'# tempo de maquina:',t_end - t_init
	
END PROGRAM XY
