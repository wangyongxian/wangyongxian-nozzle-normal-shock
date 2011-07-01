module coeficientes
  
use variaveis  
use dados
  
implicit none
  
contains

!-------------------------------------------------


  subroutine coeficientes_e_fontes_qml
	real*8  :: fator, dx ! auxiliar
	real*8 :: oeste, leste ! auxiliares
	integer ::i
	! volume 1 (fictício)

    fator = 2.0d0 * xp(2) / ( xp(3) - xp(2) )
    aw(1) =  0.0d0
    ap(1) =  1.0d0
    ae(1) = -1.0d0
    bp(1) = -fator * ( u(3) - u(2) )
    bf(1) = 0.0d0
    bc(1) = 0.0d0
    ! volumes internos

    do i = 2, n-1

       dx = xe(i) - xe(i-1)

       aw(i) = - roe(i-1) * ue(i-1) * se(i-1)

       ae(i) = 0.0d0

       ap(i) = ro_o(i) * sp(i) * dx / dt - ( aw(i) + ae(i) )

       bf(i) = 0.0d0 !- pi * f(i) * ro(i) * (u(i)**2) * raio(i) * dx / 4.0d0   

       bp(i) = ro_o(i) * sp(i) * dx * u_o(i) / dt + bf(i)   &
             + 0.5d0  * sp(i) * ( p(i-1) - p(i+1) )
             
       oeste = roe(i-1) * ue(i-1) * se(i-1) * ( u(i) - u(i-1) )

       leste = roe(i) * ue(i) * se(i) * ( u(i+1) - u(i) )

       bc(i) = 0.5d0 * beta * ( oeste - leste )

    end do

    ! volume n (fictício)

    fator = 2.0d0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) = -1.0d0
    ap(n) =  1.0d0
    ae(n) =  0.0d0
    bp(n) =  fator * ( u(n-1) - u(n-2) )
    bf(n) =  0.0d0
    bc(n) = 0.0d0
    bp = bp + bc
    
  end subroutine coeficientes_e_fontes_qml

!-------------------------------------------------
  
  subroutine calculo_velocidades_face

	real*8 :: massa_e ! massa do volume E
    real*8 :: massa_p ! massa do volume P
    real*8 :: somap, somae
    integer ::i
    ! média antiga

    do i = 2, n-2

       massa_p = ro_o(i) * sp(i) * (xe(i) - xe(i-1))

       massa_e = ro_o(i+1) * sp(i+1) * (xe(i+1) - xe(i))

       somap = aw(i)*u(i-1) + ae(i)*u(i+1)
 
       somae = aw(i+1)*u(i) + ae(i+1)*u(i+2)

       ue(i) = (-somap - somae + bc(i) + bc(i+1) + bf(i) + bf(i+1)       &
             + (massa_p+massa_e)*ue_o(i)/dt - 2.0d0*se(i)*(p(i+1)-p(i)))  &
             / (ap(i)+ap(i+1))

    end do

    ue(1) = 0.5d0 * ( u(1) + u(2) )

    ue(n-1) = 0.5d0 * ( u(n) + u(n-1) )
 
  end subroutine calculo_velocidades_face

!-------------------------------------------------


  subroutine coeficientes_e_fontes_energia
	real*8  :: fator, dx ! auxiliar
	real*8 :: oeste, leste ! auxiliares
	integer ::i
	! volume 1 (fictício)

    aw(1) = 0.0d0
    ap(1) = 1.0d0
    ae(1) = 1.0d0
    bp(1) = 2.0d0 * T_in
    bc(1) = 0.0d0
    ! volumes internos

    do i = 2, n-1

       dx = xe(i) - xe(i-1)

       aw(i) = - cp * roe(i-1) * ue(i-1) * se(i-1)

       ae(i) = 0.0d0

       ap(i) = cp * ro_o(i) * sp(i) * dx / dt - ( aw(i) + ae(i) )

       bp(i) = cp * ro_o(i) * sp(i) * dx * T_o(i) / dt        &
             + sp(i) * dx * ( p(i) - p_o(i) ) / dt              &
             + sp(i) * u(i) * ( p(i+1) - p(i-1) ) * 0.5d0     ! &
             !+ pi * f(i) * ro(i) * (u(i)**3) * raio(i) * dx / 4  
       oeste = roe(i-1) * ue(i-1) * se(i-1) * ( T(i) - T(i-1) )

       leste = roe(i) * ue(i) * se(i) * ( T(i+1) - T(i) )

       bc(i) = 0.5d0 *beta* cp * ( oeste - leste )
    end do

    ! volume n (fictício)

    fator = 2.0d0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) =  -1.0d0
    ap(n) =  1.0d0
    ae(n) =  0.0d0
    bc(n) = 0.0d0
    bp(n) =  fator * ( T(n-1) - T(n-2) )
    !bp(n) =  2.0d0*T_out
! volume n (fictício)

    ! atualiza o termo independente

    bp = bp + bc
  end subroutine coeficientes_e_fontes_energia
  
!-------------------------------------------------

  subroutine coeficientes_fontes_massa

    real*8 :: fator, dx ! auxiliar
    real*8 :: oeste, leste ! auxiliares
    integer ::i 
    ! volume 1 (fictício)

    bc(1) = 0.0d0
    aw(1) = 0.0d0
    ap(1) = 1.0d0
    ae(1) = 1.0d0
    bp(1) = 2.0d0 * pl_in

    ! volumes internos

    do i = 2, n-1

       dx = xe(i) - xe(i-1)

       aw(i) = - roe(i-1) * de(i-1) * se(i-1)          &
               - ue(i-1)  * se(i-1) / ( R * T(i-1) )

       ae(i) = - roe(i) * de(i) * se(i)

       ap(i) = ( sp(i) * dx / dt + ue(i) * se(i) ) / ( R * T(i) )   &
             + roe(i-1) * de(i-1) * se(i-1)                              &
             + roe(i)   * de(i)   * se(i)

       bp(i) = - ( ( ro(i) - ro_o(i) ) * sp(i) * dx / dt              &
             + ro(i) * ue(i) * se(i) - ro(i-1) * ue(i-1) * se(i-1) )
       
       oeste = se(i-1) * ( ro(i) - ro(i-1) ) * ue(i-1)

       leste = se(i)   * ( ro(i+1) - ro(i) ) * ue(i)

       bc(i) = 0.5d0 * beta * ( oeste - leste )
       
    end do
 ! volume n (fictício)

    bc(n) = 0.0d0
    !fator = 2.0d0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) = 1.0d0
    ap(n) =  1.0d0
    ae(n) =  0.0d0
    !bp(n) =  fator * ( pl(n-1) - pl(n-2) )
    bp(n) =  2.0d0*pl_out
    ! atualiza o termo independente

    bp = bp + bc
  end subroutine coeficientes_fontes_massa

!-------------------------------------------------

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

!-------------------------------------------------

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

end module coeficientes

