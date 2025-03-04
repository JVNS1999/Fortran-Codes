MODULE PARAMETERS
IMPLICIT NONE

INTEGER(KIND=8),PARAMETER:: nstep = 1e6
INTEGER(KIND=8),PARAMETER:: nstepwl = 1e5
INTEGER(KIND=8),PARAMETER:: nstepheat_aleat_num = 531581		!Prime[Prime[10^3 + Prime[500]]]
INTEGER(KIND=8),PARAMETER:: POWER = 1e7
INTEGER(KIND=8),PARAMETER:: sobra = 1e3

REAL(KIND=16),PARAMETER:: LNF_INI = 1.q0, LNF_END = 1.q-8
REAL(KIND=16),PARAMETER:: expoente_funcao_lnf = 0.5q0
REAL(KIND=16),PARAMETER:: cond_flat_hist = 20.0q-2
REAL(KIND=16),PARAMETER:: cond_min_flat = QEXT(nstepwl)*1.q-2
REAL(KIND=16),PARAMETER:: DELTAE = 10.q0
REAL(KIND=16),PARAMETER:: alphae = 1.q0/DELTAE

REAL(KIND=16),PARAMETER:: DDIP = 1.q0
REAL(KIND=16),PARAMETER:: kb = 1.q0
REAL(KIND=16),PARAMETER:: TEMP_MAX = 1.5q1, TEMP_MIN = 1.0q-1, DTEMP = 1.q-2
REAL(KIND=16),PARAMETER:: rcut_frac_lado = QSQRT(1.q0/2.q0)

REAL(KIND=16),PARAMETER:: cte_ge = 1.q3

INTEGER(KIND=8),PARAMETER:: nstep_emin = 1e4
INTEGER(KIND=8),PARAMETER:: namostras_emin = 1e0
INTEGER(KIND=8),PARAMETER:: nstep_emax = 1e4
INTEGER(KIND=8),PARAMETER:: namostras_emax = 1e0
REAL(KIND=8),PARAMETER:: DTEMP_emin = 0.5q0, DTEMP_emax = 0.5q0

LOGICAL(KIND=4),PARAMETER:: INITIAL_ALEAT = .TRUE.
LOGICAL(KIND=4),PARAMETER:: periodic_condition = .TRUE.
LOGICAL(KIND=4),PARAMETER:: print_histograma = .FALSE.

END MODULE PARAMETERS

!===============================================================================

MODULE variaveis
USE PARAMETERS
IMPLICIT NONE

REAL(16),DIMENSION(:),ALLOCATABLE:: sx, sy, sz
REAL(16),DIMENSION(:),ALLOCATABLE:: sx_ini, sy_ini, sz_ini
REAL(16),DIMENSION(:),ALLOCATABLE:: sxmin, symin, szmin
REAL(16),DIMENSION(:),ALLOCATABLE:: sxmax, symax, szmax
REAL(16),DIMENSION(:),ALLOCATABLE:: dx, dy, dz
REAL(16),DIMENSION(:),ALLOCATABLE:: enerx, enery, enerz
REAL(16),DIMENSION(:),ALLOCATABLE:: enerx_ini, enery_ini, enerz_ini

REAL(16),DIMENSION(POWER):: r3, nrx, nry, nrz
INTEGER(8),DIMENSION(POWER):: viz
INTEGER(8),DIMENSION(:),ALLOCATABLE:: nviz

REAL(16),DIMENSION(:),ALLOCATABLE:: lnge, ge
INTEGER(8),DIMENSION(:),ALLOCATABLE:: h, hacc

REAL(16),DIMENSION(:),ALLOCATABLE:: Eviz

REAL(16):: x, y, z, p, r, s, T, beta, E, DE
REAL(16):: E_ini, Enew, Eold, hmed, hmini, hmax
REAL(16):: em, e2m, emin, emax, DELTA, ce
REAL(16):: PI, rcut
REAL(16):: esc1, esc2, esc3
REAL(16):: lnf, f

INTEGER(8):: i, j, k, m, n, kx, ky, kz, ind1, ind2
INTEGER(8):: sflip, cont, cont2, TAM
INTEGER(8):: eff1, eff2

 CHARACTER(80):: arq, dir, arq2, dir2
 CHARACTER(40):: str1, str2, str3
LOGICAL(4):: path, flat

!========== VARIAVEIS DE ENTRADA DE LEITURA DO PROGRAMA ==========
INTEGER(KIND=8):: DIMENSAO
INTEGER(KIND=8):: L
REAL(KIND=16):: LADO
INTEGER(KIND=8):: nspins

CONTAINS
	SUBROUTINE INICIALIZA_PROGRAMA_E_VETORES
	USE IFPORT
	USE PARAMETERS
	IMPLICIT NONE
	CALL RANDOM_SEED()
	
	DO i=1, nstepheat_aleat_num
		CALL RANDOM_NUMBER(s)
	ENDDO
	
	OPEN(10,FILE="position_spins.dat")
	
	READ(10,*)DIMENSAO
	READ(10,*)L
	READ(10,*)LADO
	READ(10,*)nspins
	
	WRITE(1,*)">>> DIMENSAO = ",DIMENSAO
	WRITE(1,*)">>> LADO = ", L
	WRITE(1,*)">>> NUM. DE SPINS = ", nspins
	
	ALLOCATE(sx(nspins) , sy(nspins) , sz(nspins))
	ALLOCATE(sx_ini(nspins) , sy_ini(nspins) , sz_ini(nspins))
	ALLOCATE(dx(nspins) , dy(nspins) , dz(nspins))
	ALLOCATE(sxmin(nspins) , symin(nspins) , szmin(nspins))
	ALLOCATE(sxmax(nspins) , symax(nspins) , szmax(nspins))
	
	DO i=1, nspins
		READ(10,*) dx(i), dy(i), dz(i), sx(i), sy(i), sz(i)
	ENDDO
	
	IF(EOF(10)) WRITE(1,*)">>> Leu todo o arquivo: position_spins.dat"
	CLOSE(10)
	
	IF(INITIAL_ALEAT)	THEN
		DO i=1, nspins
			CALL RANDOM_NUMBER(s)
			r = 2.q0*s-1.q0
			r = r/QABS(r)
			sx(i) = r*sx(i)
			sy(i) = r*sy(i)
			sz(i) = r*sz(i)
		ENDDO
	ENDIF
	
	FORALL(i=1:nspins)
		sx_ini(i) = sx(i)
		sy_ini(i) = sy(i)
		sz_ini(i) = sz(i)
	END FORALL
	
	WRITE(1,*)"			"
	WRITE(1,*)"			"
	WRITE(1,*)"			"
	
	RETURN
	END SUBROUTINE
	
!===============================================================================
	
	SUBROUTINE find_neighbors_2d
	USE IFPORT
	USE PARAMETERS
	IMPLICIT NONE
	INTEGER(8):: period_cond_factor
	CALL RANDOM_SEED()
	
	ALLOCATE(nviz(0:nspins))
	
	IF(periodic_condition)	THEN
		period_cond_factor = 1
	ELSE
		period_cond_factor = 0
	ENDIF
	
	rcut = rcut_frac_lado*LADO
	cont = 0
	DO i=1, nspins
		DO j=1, nspins
			IF(j == i)	CYCLE
			
			DO kx=-period_cond_factor, period_cond_factor
			DO ky=-period_cond_factor, period_cond_factor
				
				x = dx(j) - dx(i) + QEXT(kx)*LADO
				y = dy(j) - dy(i) + QEXT(ky)*LADO
				r = QSQRT(x*x + y*y)
				IF(r .LE. rcut)		THEN
					cont = cont + 1
					viz(cont) = j
					nrx(cont) = x/r
					nry(cont) = y/r
					r3(cont) = 1.q0/r/r/r
				ENDIF
				
			ENDDO
			ENDDO
		ENDDO
		nviz(i) = cont
	ENDDO
	nviz(0) = 0
	
	r = DBLE(nspins)/DBLE(L*L)
	x = PI*(rcut_frac_lado**2)*((r)**2)*(LADO**4)	
!	WRITE(*,*) cont, x
	WRITE(1,*) ">>> Contagem exata de vizinhos = ", cont
	WRITE(1,*) ">>> Previsao da total de vizinhos = ", x
	WRITE(1,*) ">>> Erro (prev/cont) = ", x/cont
	WRITE(1,*) "		"
	WRITE(1,*) "		"
	WRITE(1,*) "		"
	
	RETURN
	END SUBROUTINE
	
!===============================================================================
	
	SUBROUTINE energia_2d
	USE IFPORT
	USE PARAMETERS
	IMPLICIT NONE
	CALL RANDOM_SEED()
	
	ALLOCATE(enerx(nspins) , enery(nspins) , enerz(nspins))
	ALLOCATE(enerx_ini(nspins) , enery_ini(nspins) , enerz_ini(nspins))
	
	FORALL(i=1:nspins)
		enerx(i) = 0.q0
		enery(i) = 0.q0
		enerz(i) = 0.q0
	END FORALL
!	WRITE(*,*)enerx(17)
	
	e = 0.q0
	DO i=1, nspins
		enerx(i) = 0.q0
		enery(i) = 0.q0
		DO j=Nviz(i-1)+1, Nviz(i)
			esc1 = sx(viz(j))*r3(j)
			esc2 = sx(viz(j))*nrx(j) + sy(viz(j))*nry(j)
			esc3 = esc2*nrx(j)*r3(j)
			enerx(i) = enerx(i) + DDIP*(esc1 - 3.q0*esc3)
			
			esc1 = sy(viz(j))*r3(j)
			esc2 = sx(viz(j))*nrx(j) + sy(viz(j))*nry(j)
			esc3 = esc2*nry(j)*r3(j)
			enery(i) = enery(i) + DDIP*(esc1 - 3.q0*esc3)
		ENDDO
		E = E + sx(i)*enerx(i) + sy(i)*enery(i)
	ENDDO
	E = 0.5q0*E
	
	e_ini = e
	enerx_ini = enerx
	enery_ini = enery
	enerz_ini = enerz
	
	WRITE(1,*)">>> Energia inicial do sistema: ", E
	WRITE(1,*)">>> Energia inicial por spin  : ", E/QEXT(nspins)
	WRITE(1,*)"			"
	WRITE(1,*)"			"
	WRITE(1,*)"			"
	
	RETURN
	END SUBROUTINE
	
!===============================================================================
	
	SUBROUTINE WANG_LANDAU_single_spin_flip_2d
	USE IFPORT
	USE PARAMETERS
	IMPLICIT NONE
	CALL RANDOM_SEED()
	
	cont = 0
	DO i=1, nstepwl
		CALL RANDOM_NUMBER(r)
		r = r*QEXT(nspins)
		sflip = 1 + INT8(r)
		
		DE = -2.q0*( sx(sflip)*enerx(sflip) + sy(sflip)*enery(sflip) )
!		Eold = E
		Enew = E + DE
		
		m = INT8(alphae*(E-emin))
		n = INT8(alphae*(Enew-emin))
		
		IF((m<eff1).OR.(n<eff1).OR.(m>eff2).OR.(n>eff2))	THEN
			cont = cont + 1
			WRITE(1,*)"POSSIVEL ERRO: estamos encontrando novas" 
			WRITE(1,*)"energias maximas ou minimas", cont
			WRITE(1,*)"			"
			WRITE(1,*)"			"
!			PAUSE
		ENDIF
		
		IF(m<eff1)	eff1 = m
		IF(N<eff1)	eff1 = n
		IF(m>eff2)	eff2 = m
		IF(n>eff2)	eff2 = n
		
		p = QEXP(lnge(m) - lnge(n))
		CALL RANDOM_NUMBER(r)
		
		IF(r .LE. p)	THEN
			E = Enew
			sx(sflip) = -sx(sflip)
			sy(sflip) = -sy(sflip)
			lnge(n) = lnge(n) + lnf
			H(n) = H(n) + 1
			
			DO j=Nviz(sflip-1)+1, Nviz(sflip)
				esc1 = sx(sflip)*r3(j)
				esc2 = sx(sflip)*nrx(j) + sy(sflip)*nry(j)
				esc3 = esc2*nrx(j)*r3(j)
				enerx(viz(j)) = enerx(viz(j)) + 2.q0*DDIP*(esc1 - 3.q0*esc3)
				
				esc1 = sy(sflip)*r3(j)
				esc2 = sx(sflip)*nrx(j) + sy(sflip)*nry(j)
				esc3 = esc2*nry(j)*r3(j)
				enery(viz(j)) = enery(viz(j)) + 2.q0*DDIP*(esc1 - 3.q0*esc3)
			ENDDO
		ELSE
			lnge(m) = lnge(m) + lnf
			H(m) = H(m) + 1
		ENDIF
	ENDDO
	
	RETURN
	END SUBROUTINE
	
	
	
	
	
!===============================================================================
!================= AQUI ENTRA O WANG-LANDAU PARA O CASO EM 3D ==================
!===============================================================================
	
	
	
	
	
	SUBROUTINE find_neighbors_3d
	USE IFPORT
	USE PARAMETERS
	IMPLICIT NONE
	INTEGER(8):: period_cond_factor
	CALL RANDOM_SEED()
	
	ALLOCATE(nviz(0:nspins))
	
	IF(periodic_condition)	THEN
		period_cond_factor = 1
	ELSE
		period_cond_factor = 0
	ENDIF
	
	rcut = rcut_frac_lado*LADO
	cont = 0
	DO i=1, nspins
		DO j=1, nspins
			IF(j == i)	CYCLE
			
			DO kx=-period_cond_factor, period_cond_factor
			DO ky=-period_cond_factor, period_cond_factor
			DO kz=-period_cond_factor, period_cond_factor
				x = dx(j) - dx(i) + QEXT(kx)*LADO
				y = dy(j) - dy(i) + QEXT(ky)*LADO
				z = dz(j) - dz(i) + QEXT(kz)*LADO
				r = QSQRT(x*x + y*y + z*z)
				IF(r .LE. rcut)		THEN
					cont = cont + 1
					viz(cont) = j
					nrx(cont) = x/r
					nry(cont) = y/r
					nrz(cont) = z/r
					r3(cont) = 1.q0/r/r/r
				ENDIF
				
			ENDDO
			ENDDO
			ENDDO
		ENDDO
		nviz(i) = cont
	ENDDO
	nviz(0) = 0
	
	r = (DBLE(nspins)/DBLE(L*L*L))
	x = (4.d0/3.d0)*PI*(rcut_frac_lado**3)*((r)**2)*(LADO**6)	
!	WRITE(*,*) cont, x
	WRITE(1,*) ">>> Contagem exata de vizinhos = ", cont
	WRITE(1,*) ">>> Previsao da total de vizinhos = ", x
	WRITE(1,*) ">>> Erro (prev/cont) = ", x/cont
	WRITE(1,*) "		"
	WRITE(1,*) "		"
	WRITE(1,*) "		"
	
	
	
	RETURN
	END SUBROUTINE
	
!===============================================================================
	
	SUBROUTINE energia_3d
	USE IFPORT
	USE PARAMETERS
	IMPLICIT NONE
	CALL RANDOM_SEED()
	
	ALLOCATE(enerx(nspins) , enery(nspins) , enerz(nspins))
	ALLOCATE(enerx_ini(nspins) , enery_ini(nspins) , enerz_ini(nspins))
	
	FORALL(i=1:nspins)
		enerx(i) = 0.q0
		enery(i) = 0.q0
		enerz(i) = 0.q0
	END FORALL
!	WRITE(*,*)enerx(17)
	
	e = 0.q0
	DO i=1, nspins
		enerx(i) = 0.q0
		enery(i) = 0.q0
		enerz(i) = 0.q0
		DO j=Nviz(i-1)+1, Nviz(i)
			esc1 = sx(viz(j))*r3(j)
			esc2 = sx(viz(j))*nrx(j) + sy(viz(j))*nry(j) + sz(viz(j))*nrz(j)
			esc3 = esc2*nrx(j)*r3(j)
			enerx(i) = enerx(i) + DDIP*(esc1 - 3.q0*esc3)
			
			esc1 = sy(viz(j))*r3(j)
			esc2 = sx(viz(j))*nrx(j) + sy(viz(j))*nry(j) + sz(viz(j))*nrz(j)
			esc3 = esc2*nry(j)*r3(j)
			enery(i) = enery(i) + DDIP*(esc1 - 3.q0*esc3)
			
			esc1 = sz(viz(j))*r3(j)
			esc2 = sx(viz(j))*nrx(j) + sy(viz(j))*nry(j) + sz(viz(j))*nrz(j)
			esc3 = esc2*nrz(j)*r3(j)
			enerz(i) = enerz(i) + DDIP*(esc1 - 3.q0*esc3)
		ENDDO
		E = E + sx(i)*enerx(i) + sy(i)*enery(i) + sz(i)*enerz(i)
	ENDDO
	E = 0.5q0*E
	
	e_ini = e
	enerx_ini = enerx
	enery_ini = enery
	enerz_ini = enerz
	
	WRITE(1,*)">>> Energia inicial do sistema: ", E
	WRITE(1,*)">>> Energia inicial por spin  : ", E/QEXT(nspins)
	WRITE(1,*)"			"
	WRITE(1,*)"			"
	WRITE(1,*)"			"
	
	RETURN
	END SUBROUTINE
	
!===============================================================================
	
	SUBROUTINE WANG_LANDAU_single_spin_flip_3d
	USE IFPORT
	USE PARAMETERS
	IMPLICIT NONE
	CALL RANDOM_SEED()
	
	cont = 0
	DO i=1, nstepwl
		CALL RANDOM_NUMBER(r)
		r = r*QEXT(nspins)
		sflip = 1 + INT8(r)
		
!		DE = -2.q0*( sx(sflip)*enerx(sflip) + sy(sflip)*enery(sflip) + sz(sflip)*enerz(sflip))
		DE = -2.q0*sx(sflip)*enerx(sflip)
		DE = DE -2.q0*sy(sflip)*enery(sflip)
		DE = DE -2.q0*sz(sflip)*enerz(sflip)
!		Eold = E
		Enew = E + DE
		
		m = INT8(alphae*(E-emin))
		n = INT8(alphae*(Enew-emin))
		
		IF((m<eff1).OR.(n<eff1).OR.(m>eff2).OR.(n>eff2))	THEN
			cont = cont + 1
			WRITE(1,*)"POSSIVEL ERRO: estamos encontrando novas" 
			WRITE(1,*)"energias maximas ou minimas", cont
			WRITE(1,*)"			"
			WRITE(1,*)"			"
!			PAUSE
		ENDIF
		
		IF(m<eff1)	eff1 = m
		IF(N<eff1)	eff1 = n
		IF(m>eff2)	eff2 = m
		IF(n>eff2)	eff2 = n
		
		p = QEXP(lnge(m) - lnge(n))
		CALL RANDOM_NUMBER(r)
		
		IF(r .LE. p)	THEN
			E = Enew
			sx(sflip) = -sx(sflip)
			sy(sflip) = -sy(sflip)
			sz(sflip) = -sz(sflip)
			lnge(n) = lnge(n) + lnf
			H(n) = H(n) + 1
			
			DO j=Nviz(sflip-1)+1, Nviz(sflip)
				esc1 = sx(sflip)*r3(j)
				esc2 = sx(sflip)*nrx(j) + sy(sflip)*nry(j) + sz(sflip)*nrz(j)
				esc3 = esc2*nrx(j)*r3(j)
				enerx(viz(j)) = enerx(viz(j)) + 2.q0*DDIP*(esc1 - 3.q0*esc3)
				
				esc1 = sy(sflip)*r3(j)
				esc2 = sx(sflip)*nrx(j) + sy(sflip)*nry(j) + sz(sflip)*nrz(j)
				esc3 = esc2*nry(j)*r3(j)
				enery(viz(j)) = enery(viz(j)) + 2.q0*DDIP*(esc1 - 3.q0*esc3)
				
				esc1 = sz(sflip)*r3(j)
				esc2 = sx(sflip)*nrx(j) + sy(sflip)*nry(j) + sz(sflip)*nrz(j)
				esc3 = esc2*nrz(j)*r3(j)
				enerz(viz(j)) = enerz(viz(j)) + 2.q0*DDIP*(esc1 - 3.q0*esc3)
			ENDDO
		ELSE
			lnge(m) = lnge(m) + lnf
			H(m) = H(m) + 1
		ENDIF
	ENDDO
	
	
	RETURN
	END SUBROUTINE
	
!===============================================================================
	
	
	
	
!===============================================================================
!=========== SUBROTINAS FLATLIST E STATISTICAL INDEPENDEM DA DIMENSAO ==========
!===============================================================================
	
	SUBROUTINE flat_list
	USE IFPORT
	USE PARAMETERS
	IMPLICIT NONE
	CALL RANDOM_SEED()
	
	z = QEXT(eff2-eff1)
	z = 1.q0/z
	
	hmed = 0.q0
	DO i=eff1, eff2
		hmed = hmed + QEXT(H(i))
	ENDDO
	hmed = z*hmed
	
	hmini = 1.q100
	DO i=ind1,ind2
		IF(H(i) == 0)					CYCLE
		IF(H(I) .LE. cond_min_flat)		CYCLE
		IF(QEXT(H(i)) .LE. hmini)		hmini = QEXT(H(i))
	ENDDO
	
	hmax = QEXT(MAXVAL(H))
	
	r = (hmax-hmed )/hmed
	s = (hmed-hmini)/hmed
	
	IF((r < cond_flat_hist).AND.(s < cond_flat_hist))	THEN
		flat = .TRUE.
	ENDIF
	
	IF(print_histograma)	THEN
		OPEN(123,FILE="histograma.dat")
		
		z = QEXT(SUM(H))
		Z = 1.q0/Z
		DO i=eff1, eff2
			x = DELTAE*QEXT(i) + emin
			WRITE(123,'(F16.8,F16.8)') x, Z*H(i)
		ENDDO
		PAUSE
		
		CLOSE(123)
	ENDIF
	
	RETURN
	END SUBROUTINE
	
!===============================================================================
	
	SUBROUTINE wang_landau_statistical
	USE IFPORT
	USE PARAMETERS
	IMPLICIT NONE
	CALL RANDOM_SEED()
	
	!========== PRIMEIRO, NORMALIZANDO g(E) ==========
	z = MAXVAL(lnge)
	x = cte_ge - z
	DO i=eff1, eff2
		lnge(i) = lnge(i) + x
	ENDDO
	WRITE(1,*) z, MAXVAL(lnge), MINVAL(lnge)
	
	
	!========== VAMOS AGORA, A ESTATISTICA DO PROBLEMA ==========
	OPEN(10,FILE="energia_wl/temp_ener.dat")
	OPEN(11,FILE="sheat_wl/temp_sheat.dat")

	s = 1.q0/QEXT(nspins)
	T = TEMP_MAX
	DO WHILE(T .GE. TEMP_MIN)
		z = 0.q0
		em = 0.q0
		e2m = 0.q0
		beta = 1.q0/T/kb
		DO i=eff1, eff2
			DELTA = QEXT(i)*DELTAE
			E = QEXT(i)*DELTAE + EMIN
			
			z = z + QEXP(lnge(i))*QEXP(-beta*DELTA)
			em = em + E*QEXP(lnge(i))*QEXP(-beta*DELTA)
			e2m = e2m + E*E*QEXP(lnge(i))*QEXP(-beta*DELTA)
		ENDDO
		em = em/z
		e2m = e2m/z
		ce = (e2m-em*em)*kb*beta*beta
		
		WRITE(10,110)T, s*em
		WRITE(11,111)T, s*ce
		
		T = T - DTEMP
	ENDDO
	
	CLOSE(10)
	CLOSE(11)
	
110	FORMAT(2F20.10)
111	FORMAT(2F20.10)	
	
	RETURN
	END SUBROUTINE
	
!===============================================================================

END MODULE variaveis

!===============================================================================

PROGRAM wang_landau_spin_ice_dipolar_1D_2D_3D
USE IFPORT
USE PARAMETERS
USE variaveis
IMPLICIT NONE
CALL RANDOM_SEED()

OPEN(1,FILE="status_program_wang_ladau.txt")

path = MAKEDIRQQ("energia_wl")
path = MAKEDIRQQ("sheat_wl")
path = MAKEDIRQQ("magnetizacao_wl")
path = MAKEDIRQQ("suscetibilidade_wl")

PI = 4.q0*QATAN(1.q0)
!WRITE(*,*) PI

CALL INICIALIZA_PROGRAMA_E_VETORES

IF(DIMENSAO == 1)			THEN
	
	
	
	
	
ELSEIF(DIMENSAO == 2)		THEN
	CALL find_neighbors_2d
	CALL energia_2d
	
	CALL find_indice_emin_2d
	CALL find_indice_emax_2d
	
	eff1 = 0
	TAM = INT8((emax-emin)/DELTAE)
	eff2 = TAM
	ind1 = eff1 - sobra
	ind2 = eff2 + sobra
	
	ALLOCATE(lnge(ind1:ind2) , ge(ind1:ind2))
	ALLOCATE(H(ind1:ind2) , hacc(ind1:ind2))
	ALLOCATE(Eviz(ind1:ind2))
	
	DO i=ind1, ind2
		Eviz(i) = QEXT(i)*DELTAE + emin !+ 0.5q0*DELTAE
	ENDDO
	
	sx = sx_ini
	sy = sy_ini
	e = e_ini
	enerx = enerx_ini
	enery = enery_ini
	
	DO i=ind1, ind2
		lnge(i) = 0.q0
		ge(i) = 1.q0
		H(i) = 0
		Hacc(i) = 0
	ENDDO
	
	lnf = LNF_INI
	flat = .FALSE.
	DO WHILE(lnf .GE. LNF_END)
		CALL WANG_LANDAU_single_spin_flip_2d
		CALL flat_list
	
		IF(flat)	THEN
			h = 0
			lnf = lnf*expoente_funcao_lnf
			WRITE(1,*)"LNF = ", lnf
			flat = .FALSE.
		ENDIF
	ENDDO
	
	CALL wang_landau_statistical
	
	
	
	
	
ELSEIF(DIMENSAO == 3)		THEN
	CALL find_neighbors_3d
	CALL energia_3d
	
	CALL find_indice_emin_3d
	CALL find_indice_emax_3d
	
	eff1 = 0
	TAM = INT8((emax-emin)/DELTAE)
	eff2 = TAM
	ind1 = eff1 - sobra
	ind2 = eff2 + sobra
	
	ALLOCATE(lnge(ind1:ind2) , ge(ind1:ind2))
	ALLOCATE(H(ind1:ind2) , hacc(ind1:ind2))
	ALLOCATE(Eviz(ind1:ind2))
	
	DO i=ind1, ind2
		Eviz(i) = QEXT(i)*DELTAE + emin !+ 0.5q0*DELTAE
	ENDDO
	
	sx = sx_ini
	sy = sy_ini
	sz = sz_ini
	e = e_ini
	enerx = enerx_ini
	enery = enery_ini
	enerz = enerz_ini
	
	DO i=ind1, ind2
		lnge(i) = 0.q0
		ge(i) = 1.q0
		H(i) = 0
		Hacc(i) = 0
	ENDDO
	
	lnf = LNF_INI
	flat = .FALSE.
	DO WHILE(lnf .GE. LNF_END)
		CALL WANG_LANDAU_single_spin_flip_3d
		CALL flat_list
	
		IF(flat)	THEN
			h = 0
			lnf = lnf*expoente_funcao_lnf
			WRITE(1,*)"LNF = ", lnf
			flat = .FALSE.
		ENDIF
	ENDDO
	
	CALL wang_landau_statistical
	
	
ELSE
	
	WRITE(1,*)">>> A algo de exagerado na dimensao!" 
	
ENDIF







 CLOSE(1)

STOP
END PROGRAM wang_landau_spin_ice_dipolar_1D_2D_3D



!===============================================================================
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!===============================================================================

	SUBROUTINE find_indice_emin_2d
	USE IFPORT
	USE PARAMETERS
	USE variaveis
	IMPLICIT NONE
	REAL(KIND=16),PARAMETER:: DDIP_eff = DDIP
	REAL(16),DIMENSION(nspins):: sx1, sy1
	REAL(16),DIMENSION(nspins):: sx1min, sy1min
	REAL(16),DIMENSION(nspins):: enerx1, enery1
	REAL(16):: e1
	CALL RANDOM_SEED()
	
	sx1 = sx_ini
	sy1 = sy_ini
	e1 = e_ini
	enerx1 = enerx_ini
	enery1 = enery_ini
	
	emin = e1
	DO k=1, namostras_emin
		T = TEMP_MAX
		DO i=1, nstepheat_aleat_num
			CALL RANDOM_NUMBER(r)
		ENDDO
		DO WHILE(T .GE. TEMP_MIN)
			beta = 1.q0/T/kb
			DO i=1, nstep_emin
				CALL RANDOM_NUMBER(r)
				r = r*QEXT(nspins)
!				sflip = 1 + KIQINT(r)			!KIQINT() LEVA REAL*16 P/ INTEIRO*8
				sflip = 1 + INT8(r)				!INT8() LEVA QUALQUER NUMERO PARA INTEIRO*8
				DE = -2.q0*( sx1(sflip)*enerx1(sflip) + sy1(sflip)*enery1(sflip) )
				IF(DE .LE. 0.q0)	THEN
					E1 = E1 + DE
					sx1(sflip) = -sx1(sflip)
					sy1(Sflip) = -sy1(Sflip)
					
					DO j=nviz(sflip-1)+1, nviz(sflip)
						esc1 = sx1(sflip)*r3(j)
						esc2 = sx1(sflip)*nrx(j) + sy1(sflip)*nry(j)
						esc3 = esc2*nrx(j)*r3(j)
						enerx1(viz(j)) = enerx1(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
						
						esc1 = sy1(sflip)*r3(j)
						esc2 = sx1(sflip)*nrx(j) + sy1(sflip)*nry(j)
						esc3 = esc2*nry(j)*r3(j)
						enery1(viz(j)) = enery1(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
					ENDDO
				ELSE
					p = QEXP(-beta*DE)
					CALL RANDOM_NUMBER(r)
					IF(r .LT. p)		THEN
						E1 = E1 + DE
						sx1(sflip) = -sx1(sflip)
						sy1(Sflip) = -sy1(Sflip)
					
						DO j=nviz(sflip-1)+1, nviz(sflip)
							esc1 = sx1(sflip)*r3(j)
							esc2 = sx1(sflip)*nrx(j) + sy1(sflip)*nry(j)
							esc3 = esc2*nrx(j)*r3(j)
							enerx1(viz(j)) = enerx1(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
							
							esc1 = sy1(sflip)*r3(j)
							esc2 = sx1(sflip)*nrx(j) + sy1(sflip)*nry(j)
							esc3 = esc2*nry(j)*r3(j)
							enery1(viz(j)) = enery1(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
						ENDDO
				!	ELSE
				!		Mantain the same	
					ENDIF
				ENDIF
				
				IF(e1 .LT. emin)	THEN
					emin = e1
					sx1min = sx1
					sy1min = sy1
				ENDIF
			ENDDO
			T = T - DTEMP_emin
		ENDDO
		WRITE(*,*) emin
	ENDDO
	
	sxmin = sx1min
	symin = sy1min
	
	WRITE(1,*)">>> Estimativa para energia minima: ", emin
	WRITE(1,*)"			"
	WRITE(1,*)"			"
	
	OPEN(123,FILE="config_emin.xyz")
	
	WRITE(123,*)nspins
	WRITE(123,*)"T = ", T
	DO i=1, nspins
		WRITE(123,555) dx(i), dy(i), dz(i), sx1min(i), sy1min(i), sz(i)
	ENDDO
	
	CLOSE(123)
555	FORMAT(1x,' C ',3F13.7,' atom_vector ',3F13.7)	
	
	RETURN
	END SUBROUTINE
	
!===============================================================================

	SUBROUTINE find_indice_emax_2d
	USE IFPORT
	USE PARAMETERS
	USE variaveis
	IMPLICIT NONE
	REAL(KIND=16),PARAMETER:: DDIP_eff = -DDIP
	REAL(16),DIMENSION(nspins):: sx2, sy2
	REAL(16),DIMENSION(nspins):: sx2max, sy2max
	REAL(16),DIMENSION(nspins):: enerx2, enery2
	REAL(16):: e2
	CALL RANDOM_SEED()
	
	sx2 = sx_ini
	sy2 = sy_ini
	e2 = e_ini
	enerx2 = -enerx_ini
	enery2 = -enery_ini
	
	emax = e2
	DO k=1, namostras_emax
		T = TEMP_MAX
		DO i=1, nstepheat_aleat_num
			CALL RANDOM_NUMBER(r)
		ENDDO
		DO WHILE(T .GE. TEMP_MIN)
			beta = 1.q0/T/kb
			DO i=1, nstep_emax
				CALL RANDOM_NUMBER(r)
				r = r*QEXT(nspins)
!				sflip = 1 + KIQINT(r)			!KIQINT() LEVA REAL*16 P/ INTEIRO*8
				sflip = 1 + INT8(r)				!INT8() LEVA QUALQUER NUMERO PARA INTEIRO*8
				DE = -2.q0*( sx2(sflip)*enerx2(sflip) + sy2(sflip)*enery2(sflip) )
				IF(DE .LE. 0.q0)	THEN
					E2 = E2 + DE
					sx2(sflip) = -sx2(sflip)
					sy2(Sflip) = -sy2(Sflip)
					
					DO j=nviz(sflip-1)+1, nviz(sflip)
						esc1 = sx2(sflip)*r3(j)
						esc2 = sx2(sflip)*nrx(j) + sy2(sflip)*nry(j)
						esc3 = esc2*nrx(j)*r3(j)
						enerx2(viz(j)) = enerx2(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
						
						esc1 = sy2(sflip)*r3(j)
						esc2 = sx2(sflip)*nrx(j) + sy2(sflip)*nry(j)
						esc3 = esc2*nry(j)*r3(j)
						enery2(viz(j)) = enery2(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
					ENDDO
				ELSE
					p = QEXP(-beta*DE)
					CALL RANDOM_NUMBER(r)
					IF(r .LT. p)		THEN
						E2 = E2 + DE
						sx2(sflip) = -sx2(sflip)
						sy2(Sflip) = -sy2(Sflip)
					
						DO j=nviz(sflip-1)+1, nviz(sflip)
							esc1 = sx2(sflip)*r3(j)
							esc2 = sx2(sflip)*nrx(j) + sy2(sflip)*nry(j)
							esc3 = esc2*nrx(j)*r3(j)
							enerx2(viz(j)) = enerx2(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
							
							esc1 = sy2(sflip)*r3(j)
							esc2 = sx2(sflip)*nrx(j) + sy2(sflip)*nry(j)
							esc3 = esc2*nry(j)*r3(j)
							enery2(viz(j)) = enery2(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
						ENDDO
				!	ELSE
				!		Mantain the same	
					ENDIF
				ENDIF
				
				IF(e2 .LT. emax)	THEN
					emax = e2
					sx2max = sx2
					sy2max = sy2
				ENDIF
			ENDDO
			T = T - DTEMP_emax
!			WRITE(*,*) T
!			PAUSE
		ENDDO
		WRITE(*,*) -emax
	ENDDO
	
	emax = -emax
	sxmax = sx2max
	symax = sy2max
	
	WRITE(1,*)">>> Estimativa para energia maxima: ", emax
	WRITE(1,*)"			"
	WRITE(1,*)"			"
	
	OPEN(123,FILE="config_emax.xyz")
	
	WRITE(123,*)nspins
	WRITE(123,*)"T = ", T
	DO i=1, nspins
		WRITE(123,555) dx(i), dy(i), dz(i), sx2max(i), sy2max(i), sz(i)
	ENDDO
	
	CLOSE(123)
555	FORMAT(1x,' C ',3F13.7,' atom_vector ',3F13.7)	
	
	RETURN
	END SUBROUTINE
	
	
	
	
!===============================================================================
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!===============================================================================

	SUBROUTINE find_indice_emin_3d
	USE IFPORT
	USE PARAMETERS
	USE variaveis
	IMPLICIT NONE
	REAL(KIND=16),PARAMETER:: DDIP_eff = DDIP
	REAL(16),DIMENSION(nspins):: sx1, sy1, sz1
	REAL(16),DIMENSION(nspins):: sx1min, sy1min, sz1min
	REAL(16),DIMENSION(nspins):: enerx1, enery1, enerz1
	REAL(16):: e1
	CALL RANDOM_SEED()
	
	sx1 = sx_ini
	sy1 = sy_ini
	sz1 = sz_ini
	e1 = e_ini
	enerx1 = enerx_ini
	enery1 = enery_ini
	enerz1 = enerz_ini
	
	emin = e1
	DO k=1, namostras_emin
		T = TEMP_MAX
		DO i=1, nstepheat_aleat_num
			CALL RANDOM_NUMBER(r)
		ENDDO
		DO WHILE(T .GE. TEMP_MIN)
			beta = 1.q0/T/kb
			DO i=1, nstep_emin
				CALL RANDOM_NUMBER(r)
				r = r*QEXT(nspins)
!				sflip = 1 + KIQINT(r)			!KIQINT() LEVA REAL*16 P/ INTEIRO*8
				sflip = 1 + INT8(r)				!INT8() LEVA QUALQUER NUMERO PARA INTEIRO*8
				DE = -2.q0*sx1(sflip)*enerx1(sflip) 
				DE = DE -2.q0*sy1(sflip)*enery1(sflip)
				DE = DE -2.q0*sz1(sflip)*enerz1(sflip)
				IF(DE .LE. 0.q0)	THEN
					E1 = E1 + DE
					sx1(sflip) = -sx1(sflip)
					sy1(Sflip) = -sy1(Sflip)
					sz1(sflip) = -sz1(sflip)
					
					DO j=nviz(sflip-1)+1, nviz(sflip)
						esc1 = sx1(sflip)*r3(j)
						esc2 = sx1(sflip)*nrx(j) + sy1(sflip)*nry(j) + sz1(sflip)*nrz(j)
						esc3 = esc2*nrx(j)*r3(j)
						enerx1(viz(j)) = enerx1(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
						
						esc1 = sy1(sflip)*r3(j)
						esc2 = sx1(sflip)*nrx(j) + sy1(sflip)*nry(j) + sz1(sflip)*nrz(j)
						esc3 = esc2*nry(j)*r3(j)
						enery1(viz(j)) = enery1(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
						
						esc1 = sz1(sflip)*r3(j)
						esc2 = sx1(sflip)*nrx(j) + sy1(sflip)*nry(j) + sz1(sflip)*nrz(j)
						esc3 = esc2*nrz(j)*r3(j)
						enerz1(viz(j)) = enerz1(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
					ENDDO
				ELSE
					p = QEXP(-beta*DE)
					CALL RANDOM_NUMBER(r)
					IF(r .LT. p)		THEN
						E1 = E1 + DE
						sx1(sflip) = -sx1(sflip)
						sy1(Sflip) = -sy1(Sflip)
						sz1(sflip) = -sz1(sflip)
					
						DO j=nviz(sflip-1)+1, nviz(sflip)
							esc1 = sx1(sflip)*r3(j)
							esc2 = sx1(sflip)*nrx(j) + sy1(sflip)*nry(j) + sz1(sflip)*nrz(j)
							esc3 = esc2*nrx(j)*r3(j)
							enerx1(viz(j)) = enerx1(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
							
							esc1 = sy1(sflip)*r3(j)
							esc2 = sx1(sflip)*nrx(j) + sy1(sflip)*nry(j) + sz1(sflip)*nrz(j)
							esc3 = esc2*nry(j)*r3(j)
							enery1(viz(j)) = enery1(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
							
							esc1 = sz1(sflip)*r3(j)
							esc2 = sx1(sflip)*nrx(j) + sy1(sflip)*nry(j) + sz1(sflip)*nrz(j)
							esc3 = esc2*nrz(j)*r3(j)
							enerz1(viz(j)) = enerz1(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
						ENDDO
				!	ELSE
				!		Mantain the same	
					ENDIF
				ENDIF
				
				IF(e1 .LT. emin)	THEN
					emin = e1
					sx1min = sx1
					sy1min = sy1
					sz1min = sz1
				ENDIF
			ENDDO
			T = T - DTEMP_emin
		ENDDO
		WRITE(*,*) emin
	ENDDO
	
	sxmin = sx1min
	symin = sy1min
	szmin = sz1min
	
	WRITE(1,*)">>> Estimativa para energia minima: ", emin
	WRITE(1,*)"			"
	WRITE(1,*)"			"
	
	OPEN(123,FILE="config_emin.xyz")
	
	WRITE(123,*)nspins
	WRITE(123,*)"T = ", T
	DO i=1, nspins
		WRITE(123,555) dx(i), dy(i), dz(i), sx1min(i), sy1min(i), sz1min(i)
	ENDDO
	
	CLOSE(123)
555	FORMAT(1x,' C ',3F13.7,' atom_vector ',3F13.7)	
	
	RETURN
	END SUBROUTINE
	
!===============================================================================

	SUBROUTINE find_indice_emax_3d
	USE IFPORT
	USE PARAMETERS
	USE variaveis
	IMPLICIT NONE
	REAL(KIND=16),PARAMETER:: DDIP_eff = -DDIP
	REAL(16),DIMENSION(nspins):: sx2, sy2, sz2
	REAL(16),DIMENSION(nspins):: sx2max, sy2max, sz2max
	REAL(16),DIMENSION(nspins):: enerx2, enery2, enerz2
	REAL(16):: e2
	CALL RANDOM_SEED()
	
	sx2 = sx_ini
	sy2 = sy_ini
	sz2 = sz_ini
	e2 = e_ini
	enerx2 = -enerx_ini
	enery2 = -enery_ini
	enerz2 = -enerz_ini
	
	emax = e2
	DO k=1, namostras_emax
		T = TEMP_MAX
		DO i=1, nstepheat_aleat_num
			CALL RANDOM_NUMBER(r)
		ENDDO
		DO WHILE(T .GE. TEMP_MIN)
			beta = 1.q0/T/kb
			DO i=1, nstep_emax
				CALL RANDOM_NUMBER(r)
				r = r*QEXT(nspins)
!				sflip = 1 + KIQINT(r)			!KIQINT() LEVA REAL*16 P/ INTEIRO*8
				sflip = 1 + INT8(r)				!INT8() LEVA QUALQUER NUMERO PARA INTEIRO*8
				DE = -2.q0*sx2(sflip)*enerx2(sflip)
				DE = DE -2.q0*sy2(sflip)*enery2(sflip)
				DE = DE -2.q0*sz2(sflip)*enerz2(sflip)
				IF(DE .LE. 0.q0)	THEN
					E2 = E2 + DE
					sx2(sflip) = -sx2(sflip)
					sy2(sflip) = -sy2(sflip)
					sz2(sflip) = -sz2(sflip)
					
					DO j=nviz(sflip-1)+1, nviz(sflip)
						esc1 = sx2(sflip)*r3(j)
						esc2 = sx2(sflip)*nrx(j) + sy2(sflip)*nry(j) + sz2(sflip)*nrz(j)
						esc3 = esc2*nrx(j)*r3(j)
						enerx2(viz(j)) = enerx2(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
						
						esc1 = sy2(sflip)*r3(j)
						esc2 = sx2(sflip)*nrx(j) + sy2(sflip)*nry(j) + sz2(sflip)*nrz(j)
						esc3 = esc2*nry(j)*r3(j)
						enery2(viz(j)) = enery2(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
						
						esc1 = sz2(sflip)*r3(j)
						esc2 = sx2(sflip)*nrx(j) + sy2(sflip)*nry(j) + sz2(sflip)*nrz(j)
						esc3 = esc2*nrz(j)*r3(j)
						enerz2(viz(j)) = enerz2(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
					ENDDO
				ELSE
					p = QEXP(-beta*DE)
					CALL RANDOM_NUMBER(r)
					IF(r .LT. p)		THEN
						E2 = E2 + DE
						sx2(sflip) = -sx2(sflip)
						sy2(sflip) = -sy2(sflip)
						sz2(sflip) = -sz2(sflip)
					
						DO j=nviz(sflip-1)+1, nviz(sflip)
							esc1 = sx2(sflip)*r3(j)
							esc2 = sx2(sflip)*nrx(j) + sy2(sflip)*nry(j) + sz2(sflip)*nrz(j)
							esc3 = esc2*nrx(j)*r3(j)
							enerx2(viz(j)) = enerx2(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
							
							esc1 = sy2(sflip)*r3(j)
							esc2 = sx2(sflip)*nrx(j) + sy2(sflip)*nry(j) + sz2(sflip)*nrz(j)
							esc3 = esc2*nry(j)*r3(j)
							enery2(viz(j)) = enery2(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
							
							esc1 = sz2(sflip)*r3(j)
							esc2 = sx2(sflip)*nrx(j) + sy2(sflip)*nry(j) + sz2(sflip)*nrz(j)
							esc3 = esc2*nrz(j)*r3(j)
							enerz2(viz(j)) = enerz2(viz(j)) + 2.q0*DDIP_eff*(esc1 - 3.q0*esc3)
						ENDDO
				!	ELSE
				!		Mantain the same	
					ENDIF
				ENDIF
				
				IF(e2 .LT. emax)	THEN
					emax = e2
					sx2max = sx2
					sy2max = sy2
					sz2max = sz2
				ENDIF
			ENDDO
			T = T - DTEMP_emax
!			WRITE(*,*) T
!			PAUSE
		ENDDO
		WRITE(*,*) -emax
	ENDDO
	
	emax = -emax
	sxmax = sx2max
	symax = sy2max
	szmax = sz2max
	WRITE(1,*)">>> Estimativa para energia maxima: ", emax
	WRITE(1,*)"			"
	WRITE(1,*)"			"
	
	OPEN(123,FILE="config_emax.xyz")
	
	WRITE(123,*)nspins
	WRITE(123,*)"T = ", T
	DO i=1, nspins
		WRITE(123,555) dx(i), dy(i), dz(i), sx2max(i), sy2max(i), sz2max(i)
	ENDDO
	
	CLOSE(123)
555	FORMAT(1x,' C ',3F13.7,' atom_vector ',3F13.7)	
	
	RETURN
	END SUBROUTINE
