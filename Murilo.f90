Module Parametros
    Integer,Parameter :: L = 8
    Integer :: IE,ID,JC,JB,W,MCs,E,M,dE,Soma,Ind
    Integer :: K,L2,I,J,B,s(L,L),IP(L),IM(L), II
    Real :: Aleatorio1,Aleatorio2,Aleatorio3,PA,T0 = 2.39
    Real :: EAux,MAux
    Real :: Energia(0:15) 
End Module Parametros

Subroutine Teste()
    Use Parametros
    Implicit None

    Open(25, file ='EMT22L8106.dat', status = 'Unknown')
    MCs = 2000
    L2 = L*L
    Do I = 1, L
        Do J = 1, L
            S(i,j) = 1
        End Do
    End Do

    M = L2
    E = -2*L2
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !              Lista de Vizinhos             !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Do I = 1, L
        IP(I) = I + 1
        IM(I) = I - 1
    End Do

    IP(L) = 1
    IM(1) = L

    Do I = 1,9,2
        Energia(I) = Exp(-2.d0*(I-5.d0)/T0)
    End Do

    Do 100 II = 1, MCs
        Do K = 1,L2
            Aleatorio1 = Ran()
            Aleatorio2 = Ran()
            I = Int(Aleatorio1*(L)) + 1 
            J = Int(Aleatorio2*(L)) + 1
            B = -S(i,j)
            IE = IM(I)
            ID = IP(I)
            JC = IM(J)
            JB = IP(J)
            Soma = S(I,JC) + S(I,JB) + S(IE,J) + S(ID,J)
            Ind = S(i,j)*Soma
            dE = 2*Ind
            Ind = 5 + Ind
            PA = Energia(Ind)

            Aleatorio3 = Ran()
            If (Aleatorio3 <= PA) Then
                S(i,j) = B
                M = M + 2*S(i,j)
                E = E + dE
            End If
        End Do

    Write(2,3) E,M
    100 Continue

    3 Format(2x,I4,2x,I4)
End Subroutine Teste


Subroutine Propriedades()
    Use Parametros
    Implicit None
    Integer :: IT, MCsD
    Real(8) :: T,B0
    Real(8) :: Den, EPP, MPP,Ex
    Real(8) :: Emed,Mmed,M2med,E2med, M4med,C,SU, U

    Open(33, file='ECMSUT22L8106.dat', status = 'Unknown')
    T0 = 2.39
    B0 = 1.d0/T0
    T = 2.38
    MCs = 2000
    MCsD = 500 !! Termalização
    L2 = L*L
    Do 105 IT = 1,25
        Open(25, file ='EMT22L8106.dat', status = 'Unknown')
        B = 1.d0/T
        DEN = 0.d0
        Emed = 0.d0
        E2med = 0.d0
        Mmed = 0.d0
        M2med = 0.d0
        M4med = 0.d0
        Do 101 i = 1, MCsD
            Read(2,*)E,M
        101 Continue 

        Do 100 i = MCsD + 1, MCs
            Read(2,*) E,M
            EPP = E
            MPP = Abs(M)
            Ex = Exp((B0 - B)*E)
            DEN = DEN + Ex
            Emed = Emed + (EPP*Ex)
            E2med = E2med + (EPP*EPP*Ex)
            Mmed = Mmed + (MPP*Ex)
            M2med = M2med + (MPP*MPP*Ex)
            M4med = M4med + (MPP*MPP*MPP*MPP*Ex)
        100 Continue

        Emed = Emed/DEN
        E2med = E2med/DEN
        C = (E2med - Emed*Emed)/(1.d0*T*T)
        Mmed = Mmed/DEN
        M2med = M2med/DEN
        SU = (M2med - Mmed*Mmed)/(1.d0*T)
        M4med = M4med/DEN
        U = (1 - (M4med/(3.d0*M2med*M2med)))
        Write(33,3) T, Emed, C, Mmed, SU, U

        3 Format (F18.8, 2x, F18.8, 2x, F18.8, 2x, F18.8, 2x, F18.8, 2x, F18.8)
        Close(25)

        !! Variação da Temperatura !!

        105 T = T + 0.001
        Close(33)
End Subroutine Propriedades        

Subroutine Erros()
    Use Parametros
    Implicit None

    Integer :: IB,Box,NB,IBox,MCsD
    Real(8) :: T,B0
    Real(8) :: Den,EPP,MPP,Ex,EE,EC,ESU,EU,EM,DenB 
    Real(8) :: Emed,Mmed,M2med,E2med,M4med,C,SU,U,C2,SU2,U2
    Real(8) :: EmedB,MmedB,M2medB,E2medB,M4medB,CB,SUB,UB

    Open(33, file = 'ECMSUT236L845A6er.dat', status = 'Unknown')
    T0 = 2.36
    B0 = 1.d0/T0
    T = 2.36
    MCs = 2000
    MCsD = 500
    L2 = L*L
    Box = 100000
    NB = MCs/Box
    DEN = 0.d0
    Emed = 0.d0
    E2med = 0.d0
    Mmed = 0.d0
    M2med = 0.d0
    M4med = 0.d0
    C = 0.d0
    SU = 0.d0
    U = 0.d0
    C2 = 0.d0
    SU2 = 0.d0
    U2 = 0.d0

    Open(25, file = 'EMT236L845a6.dat', status= 'Old')
    Do 101 I = 1, MCsD
        Read(25,*)E,M
    101 Continue
    IBox = MCsD
    Do 105 IB = 1, NB
        B = 1.d0/T
        DenB = 0.d0
        EmedB = 0.d0
        E2medB = 0.d0
        MmedB = 0.d0
        M2medB = 0.d0
        M4medB = 0.d0

        Do 100 I = IBox + 1, IBox + Box
            Read(2,*) E,M
            EPP = E
            MPP = Abs(M)
            EPP = EPP/L2

            MPP = MPP/L2
            EX = Exp((B0-B)*E)
            DenB = DenB + EX
            EmedB = EmedB + (EPP*EX)
            E2medB = E2medB + (EPP*EPP*EX)
            MmedB = MmedB + (MPP*EX)
            M2medB = M2medB + (MPP*MPP*EX)
            M4medB = M4medB + (MPP*MPP*MPP*MPP*EX)

        100 Continue

            EmedB = EmedB/DenB 
            E2medB = E2medB/DenB
            CB = (E2medB - EmedB*EmedB)/(1.0d0*T*T)
            MmedB = MmedB/DenB
            M2medB = M2medB/DenB
            SUB = (M2medB - MmedB*MmedB)/(1.0d0*T)
            M4medB = M4medB/DenB
            UB = (1 - (M4medB/(3.0d0*M2medB*M2medB)))
            Emed = Emed + EmedB
            E2med = E2med + E2medB
            Mmed = Mmed + MmedB
            M2med = M2med + M2medB
            M4med = M4med + M4medB
            C = C + CB
            SU = SU + SUB
            U = U + UB
            C2 = C2 + (CB*CB)
            SU2 = SU2 + (SUB*SUB)
            U2 = U2 + (UB*UB)
            Den = Den + DenB

    105 IBox = IBox + Box
    Emed = Emed/NB
    E2med = E2med/NB
    Mmed = Mmed/NB
    M2med = M2med/NB
    M4med = M4med/NB
    C = C/NB
    SU = SU/NB
    U = U/NB
    C2 = C2/NB
    SU2 = SU2/NB
    U2 = U2/NB

    EE = Sqrt((E2med - Emed*Emed)/(NB - 1.0d0))
    EM = Sqrt((M2med - Mmed*Mmed)/(NB - 1.0d0))
    EC = Sqrt((C2 - C*C)/(NB - 1.0d0))
    ESU = Sqrt((SU2 - SU*SU)/(NB - 1.0d0))
    EU = Sqrt((U2 - U*U)/(NB - 1.0d0))
    write(33,2) T,MCs,Box,NB
    write(33,1) Den
    1 format(2x,'DEN', F10.0)
    2 FORMAT('T =', F7.4, 2x,'MCS =', I8, 2x,'box =', I8, 3x,'box =', I6)
    write(33,3) Emed, C, Mmed, SU, U
    3 FORMAT(F10.6, 2x, F10.6, 2x, F10.6, 2x, F10.6, 2x, F10.6, 2x, F10.6)
    write(33,3) EE,EC,EM,ESU,EU

    close(25)
    close(33)
End Subroutine Erros


Program Main 
    Use Parametros
    Implicit None

    Call Teste
    Call Propriedades
    Call Erros

End Program Main