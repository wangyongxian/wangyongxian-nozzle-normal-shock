module Solucao_CDS
  
use variaveis  
!use dados
  
implicit none
  
contains

!-------------------------------------------------


  subroutine coeficientes_e_fontes_qml_cds
	real*8  :: fator, dx ! auxiliar
	real*8 :: oeste, leste ! auxiliares
	integer ::i
	! volume 1 (fict�cio)

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

    ! volume n (fict�cio)

    fator = 2.0d0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) = -1.0d0
    ap(n) =  1.0d0
    ae(n) =  0.0d0
    bp(n) =  fator * ( u(n-1) - u(n-2) )
    bf(n) =  0.0d0
    bc(n) = 0.0d0
    bp = bp + bc
    
  end subroutine coeficientes_e_fontes_qml_cds

!-------------------------------------------------
  
  subroutine calculo_velocidades_face_cds

	real*8 :: massa_e ! massa do volume E
    real*8 :: massa_p ! massa do volume P
    real*8 :: somap, somae
    integer ::i
    ! m�dia antiga

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
 
  end subroutine calculo_velocidades_face_cds

!-------------------------------------------------


  subroutine coeficientes_e_fontes_energia_cds
	real*8  :: fator, dx ! auxiliar
	real*8 :: oeste, leste ! auxiliares
	integer ::i
	! volume 1 (fict�cio)

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

    ! volume n (fict�cio)

    fator = 2.0d0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) =  -1.0d0
    ap(n) =  1.0d0
    ae(n) =  0.0d0
    bc(n) = 0.0d0
    bp(n) =  fator * ( T(n-1) - T(n-2) )
    !bp(n) =  2.0d0*T_out
! volume n (fict�cio)

    ! atualiza o termo independente

    bp = bp + bc
  end subroutine coeficientes_e_fontes_energia_cds
  
!-------------------------------------------------

  subroutine coeficientes_fontes_massa_cds

    real*8 :: dx, fator ! auxiliar
    real*8 :: oeste, leste ! auxiliares
    integer ::i 
    ! volume 1 (fict�cio)

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
 ! volume n (fict�cio)

    bc(n) = 0.0d0
    fator = 2.0d0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) = -1.0d0
    ap(n) =  1.0d0
    ae(n) =  0.0d0
    bp(n) =  fator * ( pl(n-1) - pl(n-2) )
    !bp(n) =  2.0d0*pl_out
    ! atualiza o termo independente

    bp = bp + bc
  end subroutine coeficientes_fontes_massa_cds

!-------------------------------------------------

end module Solucao_CDS

