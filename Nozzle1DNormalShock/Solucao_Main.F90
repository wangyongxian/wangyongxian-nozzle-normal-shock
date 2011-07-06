module Main

use Initialization
use output
use variaveis
use Solucao_CDS
use Solucao_CDS_UDS
use Solucao_UDS
use solvers_1D
use Solucao_Analitica

implicit none

contains

subroutine solucao_numerica

	integer :: it
    real*8 :: k    ! auxiliar
    real*8 :: M_in ! número de Mach na entrada
    
	tcpu = timef() ! zera cronômetro
	
    ue = 0.0d0
    !p     = p_cam/((1.0d0+(gama-1.0d0)*(Mach**2)/2.0d0)**(gama/(gama-1.0d0)))
    ro    = p / ( R * T )
    ue    = u(1:n-1)
    pl    = 0.0d0
    roe   = ro(1:n-1)
    pl_in = 0.0d0
    bc    = 0.0d0
    
    ! inicialização na entrada da tubeira
    u_in  = ue(1)
    u_out  = ue(N-1)
    T_in  = T_cam - 0.5d0*(gama-1.0d0)*(u_in**2)/(gama*R)
    T_out  = T_cam - 0.5d0*(gama-1.0d0)*(u_out**2)/(gama*R)
    
    M_in = u_in / dsqrt ( gama * R * T_in )
    M_out = u_out / dsqrt ( gama * R * T_out )
    
    k = 1.0d0 + 0.5d0 * ( gama - 1.0d0 ) * ( M_in ** 2 )
    p_in = p_cam / ( k ** ( gama / ( gama - 1.0d0 ) ) )
    k = 1.0d0 + 0.5d0 * ( gama - 1.0d0 ) * ( M_out ** 2 )
    p_outN = p_out
    
    p(1)  = 2.0d0*p_in - p(2)
    u(1)  = 2*u(2) - u(3)
    T(1)  = -T(2) + 2*T_in
    ro(1) = p(1) / ( R * T(1) )
    
    ! inicialização na saída da tubeira
    p(n)  = 2.0d0*p_out - p(n-1)  !2.0d0*p(n-1) - p(n-2)
    u(n)  = 2.0d0*u(n-1) - u(n-2)
   ! T(n)  = 2.0d0*T(n-1) - T(n-2)
    T(N) = -T(N-1) + 2*T_out
    ro(n) = p(n) / ( R * T(n) )
    de = 0.0d0
    ds = 0.0d0
    
    !residuo
    call coeficientes_e_fontes_qml_cds_uds
    call Norma_L1( n, aw, ap, ae, bp, u, Residuo_U )
	Residuo_U_o = Residuo_U
	
	call coeficientes_fontes_massa_cds_uds
	call Norma_L1( n, aw, ap, ae, bp, pl, Residuo_P )
	Residuo_P_o = Residuo_P
	
	call coeficientes_e_fontes_energia_cds_uds
	call Norma_L1( n, aw, ap, ae, bp, T, Residuo_T )
	Residuo_T_o = Residuo_T
	
	
	open(8, file='norma_l1.txt')
	write(8,4)
	4 format(t1,'iteracao',t18 ,'Residuo T',t38 ,'Residuo U',t58,'Residuo Pl')
	
    do it = 1, iteracao
    
      !if ((it > 5000).and.(it < 50000)) beta = 0.9d0
      !if (it >= 50000 ) beta = 1.0d0
	
	  ! Atualizando campos para novo avanço
	  u_o  = u
	  ue_o = ue
	  p_o = p
	  T_o = T
	  ro_o = ro
	   
	  ! atualização da pressão na entrada da tubeira
      p_ia = p_in
      M_in = u_in / dsqrt ( gama * R * T_in )
      k = 1.0d0 + 0.5d0 * ( gama - 1.0d0 ) * ( M_in ** 2 )
      p_in = p_cam / ( k ** ( gama / ( gama - 1.0d0 ) ) )
      p(1) = 2.0d0*p_in - p(2)
      pl_in = p_in - p_ia
       
      p_outa = p_outN
      M_out = u_out / dsqrt ( gama * R * T_out )
      k = 1.0d0 + 0.5d0 * ( gama - 1.0d0 ) * ( M_out ** 2 )
      p_outN = p_out / ( k ** ( gama / ( gama - 1.0d0 ) ) )
      p(N) = 2.0d0*p_outN - p(N-1)
      pl_out = p_outN - p_outa
	

	  ! cálculo dos coeficientes e termos fontes
	  call coeficientes_e_fontes_qml_cds_uds
	   
	  ! solução do sistema de equações
	  call tdma (N,ap,aw,ae,bp,u)	
	  
	  call Norma_L1 (N,ap,aw,ae,bp,u,Residuo_U)
	  Residuo_U = Residuo_U/Residuo_U_o
	    
      ! cálculo dos coeficientes do método SIMPLEC
	  call coeficientes_simplec
	   
      ! cálculos das velocidades na face leste
	  call calculo_velocidades_face_cds_uds
!-----------------------------------------------------		  
      ! inicialização na entrada da tubeira
      u_in  = ue(1)
      u_out  = ue(N-1)

      T_in = T_cam - 0.5d0*(gama-1.0d0)*(u_in**2)/(gama*R)
      T_out = T_cam - 0.5d0*(gama-1.0d0)*(u_out**2)/(gama*R)
    
	  ! cálculo dos coef e fontes da energia
	  call coeficientes_e_fontes_energia_cds_uds
	  
	  ! solução do sistema de equações
	  call tdma (N,ap,aw,ae,bp,T)
	  
	  call Norma_L1 (N,ap,aw,ae,bp,T,Residuo_T)
	  Residuo_T = Residuo_T/Residuo_T_o
	   
	  ! cálculo da massa específica
      ro = p / ( R * T )
       
	  call calculo_massa_especifica_nas_faces
	  
	  ! cálculo dos coef e fontes da massa
      call coeficientes_fontes_massa_cds_uds
	  
      ! solução do sistema de equações
      call tdma (N,ap,aw,ae,bp,pl)
      
      call Norma_L1 (N,ap,aw,ae,bp,pl,Residuo_P)
	  Residuo_P = Residuo_P/Residuo_P_o
	   
      call correcoes_com_plinha
      call calculo_massa_especifica_nas_faces
      
      !T_ex   = 0.5d0 * ( T(n-1) + T(n) )
      
      write(8,16) it, Residuo_T, Residuo_U, Residuo_P
	  16 format (i11,5x,3(1pe20.13))
	   
!-----------------------------------------------------	      

	end do
    close(8)
    
    ! coeficiente da pressao
    call coeficientes_fontes_massa_cds_uds
    call gera_arq_coef(2)
    !coeficiente da Energia
    !call coeficientes_e_fontes_energia_cds
    !call gera_arq_coef(3)
    !coeficiente da QML
    !call coeficientes_e_fontes_qml_cds
    !call gera_arq_coef(4)
    
    ! entrada da tubeira
    u(1)  = u_in
    p(1)  = p_in
    T(1)  = T_in

    ! saída da tubeira
    T_ex  = 0.5d0 * ( T(n-1) + T(n) )
    p_ex  = 0.5d0 * ( p(n-1) + p(n) )
    u(n)  = ue(n-1)
    p(n)  = p_ex
    T(n)  = T_ex
    
    ro_in = p_in / ( R * T_in )
    ro(1) = ro_in
    ro_ex = p_ex / ( R * T_ex )
    ro(n) = ro_ex

    
    call calcula_empuxo
    call calcula_coeficiente_descarga
    call calcula_fluxo_massa
    call calcula_coeficiente_descarga
    
	tcpu = timef()

  end subroutine solucao_numerica

 subroutine coeficientes_simplec
    integer ::i
    
    ds(1) = 0.0d0
    ds(n) = 0.0d0

    do i = 2, n-1
       ds(i) = sp(i) / ( ap(i) + aw(i) + ae(i) )
    end do

    de(1)   = ds(2)
    de(n-1) = ds(n-1)

    do i = 2,n-2
       de(i) = 0.5d0 * ( ds(i) + ds(i+1) )
    end do

	

  end subroutine coeficientes_simplec
  
  
subroutine calcula_fluxo_massa
        
        m = roe*ue*Se
        
end subroutine calcula_fluxo_massa

!-------------------------------------------------

subroutine calcula_coeficiente_descarga
    
    call calcula_fluxo_massa
    Cd = m/(dsqrt(gama*R*T))
    
end subroutine calcula_coeficiente_descarga

!-------------------------------------------------

subroutine calcula_empuxo
        
    Empuxo = ro*u*Sp*u
        
end subroutine calcula_empuxo

!-------------------------------------------------

subroutine correcoes_com_plinha

    real*8 :: fator ! auxiliar
    integer ::i
    
    ! pressão
    p = p + pl
    !p(N) = P_out

    ! massa específica nodal
    ro = ro + pl / ( R * T )

    ! velocidade média
    do i = 1, n-1
       ue(i) = ue(i) - de(i) * ( pl(i+1) - pl(i) )
    end do

    ! velocidade na entrada da tubeira
    u_in = ue(1)

    ! velocidade nodal
    do i = 2,n-1
       u(i) = u(i) - 0.5d0 * ds(i) * ( pl(i+1) - pl(i-1) )
    end do

    fator = 2.0d0 * xp(2) / ( xp(3) - xp(2) )
	u(1) = u(2) - fator * ( u(3) - u(2) )

    fator = 2.0d0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
	u(n) = u(n-1) + fator * ( u(n-1) - u(n-2) )
	
end subroutine correcoes_com_plinha



!-------------------------------------------------

  subroutine atualizar_ficticios_massa
	! Calculando fictício P = 1
	pl(1) = 2.0d0*pl(1+1) - pl(1+2)

	! Calculando fictício P = N 
	pl(N) = 2.0d0*pl(N-1) - pl(N-2)    
  
  end subroutine atualizar_ficticios_massa 

!-------------------------------------------------

  subroutine corrigir_velocidades
    integer ::i
     !arrumar
	u(1) = - u(2) + 2.0d0*u(1)
   
	! Calculando para volumes internos
	do i = 2, N-1
	   u(i) = u(i) - ds(i)*(pl(i+1)-pl(i-1))/2.0d0
	end do

	! Calculando fictício P = N
	u(N) = u(N-1) + (u(N-1)-u(N-2))

  end subroutine corrigir_velocidades

!-------------------------------------------------

  subroutine corrigir_velocidades_faces
    integer ::i
    
    do i = 2, N-2
	   ue(i) = ue(i) - de(i)*(pl(i+1)-pl(i))
	end do 

  end subroutine corrigir_velocidades_faces

!-------------------------------------------------

subroutine calculo_massa_especifica
    
    ro = p/(R*T)
    
end subroutine calculo_massa_especifica

!-------------------------------------------------

subroutine calculo_massa_especifica_nas_faces
    integer ::i
    do i = 1, n-1

       roe(i) = ro(i) + beta * 0.5d0 * ( ro(i+1) - ro(i) )

    end do
    
end subroutine calculo_massa_especifica_nas_faces

end module Main