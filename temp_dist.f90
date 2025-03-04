Program Paint 
       Implicit None
       Logical :: Ok_para_pintar = .False.
       Real(8) :: Velocida_Vento, Temp, Umida_relat

       Print*, 'Entre com o valor da Velocidade do Vento:'
       Print*
       Print*, 'Entre com a Umidade Relativa:'
       Print*
       Print*, 'Entre com a Temperatura:'
       Read(*,*)Velocida_Vento, Umida_relat, Temp

       If(Umida_relat > 85.d0) then
              Ok_para_pintar = .False.
       Else If (Velocida_Vento >= 3.d0) Then
              Ok_para_pintar = .False.
       Else If (Temp <10.0 .OR. Temp > 30.0) Then
              Ok_para_pintar = .False.
       Else
              Ok_para_pintar = .True.
       End if

       Print*; Print*
       Print*,'OK TO PAINT = : ', Ok_para_pintar

End Program 

