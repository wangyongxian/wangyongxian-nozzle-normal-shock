module Solucao_CDS_UDS
  
use variaveis  
  
implicit none
  
contains

!-------------------------------------------------


  subroutine coeficientes_e_fontes_qml_cds_uds
	real*8  :: fator, dx ! auxiliar
	integer ::i
	
	! volume 1 (fictício)
    fator = 2.0d0 * xp(2) / ( xp(3) - xp(2) )
    aw(1) =  0.0d0
    ap(1) =  1.0d0
    ae(1) = -1.0d0
    bp(1) = -fator * ( u(3) - u(2) )
   ! bc = 0.0d0
    ! volumes internos
    do i = 2, n-1

       dx = xe(i) - xe(i-1)

       aw(i) = -roe(i-1) * ue(i-1) * se(i-1)

       ae(i) = 0.0d0

       ap(i) = ro(i) * sp(i) * dx / dt + ( roe(i) * ue(i) * se(i) )

       bp(i) = ro_o(i) * sp(i) * dx * u_o(i) / dt   &
             - 0.5d0 * sp(i) * ( p(i+1) - p(i-1) ) &
             + 0.5d0 * beta * roe(i-1) * Se(i-1) * ue(i-1) * ( u_o(i) - u_o(i-1) )  &
             - 0.5d0 * beta * roe(i) * Se(i) * ue(i) * ( u_o(i+1) - u_o(i) )
       !bc(i) = 0.5d0 * beta * ( 0.5d0 * beta * roe(i-1) * Se(i-1) * ue(i-1) * ( u_o(i) - u_o(i-1) ) - 0.5d0 * beta * roe(i) * Se(i) * ue(i) * ( u_o(i+1) - u_o(i) ) )
    end do
    
    ! volume n (fictício)
    fator = 2.0d0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) = -1.0d0
    ap(n) =  1.0d0
    ae(n) =  0.0d0
    bp(n) =  fator * ( u(n-1) - u(n-2) )
    
  end subroutine coeficientes_e_fontes_qml_cds_uds

!-------------------------------------------------
  
  subroutine calculo_velocidades_face_cds_uds

	real*8 :: massa_e ! massa do volume E
    real*8 :: massa_p ! massa do volume P
    real*8 :: somap, somae
    integer ::i

    do i = 2, n-2

       massa_p = ro_o(i) * sp(i) * (xe(i) - xe(i-1))

       massa_e = ro_o(i+1) * sp(i+1) * (xe(i+1) - xe(i))

       somap = aw(i)*u_o(i-1) + ae(i)*u_o(i+1)
 
       somae = aw(i+1)*u_o(i) + ae(i+1)*u_o(i+2)

       ue(i) = ( -somap -somae + (massa_p+massa_e)*ue_o(i)/dt - 2.0d0*se(i)*(p(i+1)-p(i)) &
                -0.5d0 * beta * roe(i) * ue(i) * se(i) * ( u_o(i+1) - u_o(i) ) &
                -0.5d0 * beta * roe(i+1) * ue(i+1) * se(i+1) * ( u_o(i+2) - u_o(i+1) ) &
                +0.5d0 * beta * roe(i-1) * ue(i-1) * se(i-1) * ( u_o(i) - u_o(i-1) )   &
                +0.5d0 * beta * roe(i) * ue(i) * se(i) * ( u_o(i+1) - u_o(i) ) ) &
               / (ap(i)+ap(i+1))

    end do

    ue(1) = 0.5d0 * ( u(1) + u(2) )

    ue(n-1) = 0.5d0 * ( u(n) + u(n-1) )
 
  end subroutine calculo_velocidades_face_cds_uds

!-------------------------------------------------


  subroutine coeficientes_e_fontes_energia_cds_uds
	real*8  :: fator, dx ! auxiliar
	integer ::i
	
	! volume 1 (fictício)
    aw(1) = 0.0d0
    ap(1) = 1.0d0
    ae(1) = 1.0d0
    bp(1) = 2.0d0 * T_in
    
    ! volumes internos
    do i = 2, n-1

       dx = xe(i) - xe(i-1)

       aw(i) = -cp * roe(i-1) * ue(i-1) * se(i-1)

       ae(i) = 0.0d0

       ap(i) = cp * ro(i) * sp(i) * dx / dt + cp * roe(i) * ue(i) * se(i)

       bp(i) = cp * ro_o(i) * sp(i) * dx * T_o(i) / dt        &
             + sp(i) * dx * ( p(i) - p_o(i) ) / dt              &
             + sp(i) * u(i) * ( p(i+1) - p(i-1) ) * 0.5d0 &
             + 0.5d0 * beta * cp * roe(i-1) * ue(i-1) * se(i-1) * ( T(i) - T(i-1) ) &
             - 0.5d0 * beta * cp * roe(i) * ue(i) * se(i) * ( T(i+1) - T(i) )

    end do

    ! volume n (fictício)
    fator = 2.0d0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) =  -1.0d0
    ap(n) =  1.0d0
    ae(n) =  0.0d0
    bp(n) =  fator * ( T(n-1) - T(n-2) )
    !bp(n) =  2.0d0*T_out

  end subroutine coeficientes_e_fontes_energia_cds_uds
  
!-------------------------------------------------

  subroutine coeficientes_fontes_massa_cds_uds

    real*8 :: dx, fator ! auxiliar
    integer ::i 
    
    ! volume 1 (fictício)
    aw(1) = 0.0d0
    ap(1) = 1.0d0
    ae(1) = 1.0d0
    bp(1) = 2.0d0 * pl_in

    ! volumes internos
    do i = 2, n-1

       dx = xe(i) - xe(i-1)

       aw(i) = - roe(i-1) * de(i-1) * se(i-1)          &
               - ue_o(i-1)  * se(i-1) / ( R * T(i-1) )

       ae(i) = - roe(i) * de(i) * se(i)       

       ap(i) = ( sp(i) * dx / dt + ue_o(i) * se(i) ) / ( R * T(i) )   &
             + roe(i-1) * de(i-1) * se(i-1)                              &
             + roe(i) * de(i)   * se(i)

       bp(i) = ro_o(i) * sp(i) * dx / dt              &
               +roe(i)*ue_o(i)*se(i) &
               -roe(i-1)*ue_o(i-1)*se(i-1) &
               -0.5d0*beta*(ro_o(i+1)-ro_o(i))*ue_o(i)*se(i) &
               +0.5d0*beta*(ro_o(i)-ro_o(i-1))*ue_o(i-1)*se(i-1) &
               -(sp(i)*dx/dt+ue_o(i)*Se(i))*ro_o(i) &
               +Se(i-1)*ue_o(i-1)*ro_o(i-1) & ! - daqui pra baixo
               -Se(i)*roe(i)*ue_o(i) &
               +Se(i-1)*roe(i-1)*ue_o(i-1)
       
    end do

    ! volume n (fictício)
    fator = 2.0d0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) = -1.0d0
    ap(n) =  1.0d0
    ae(n) =  0.0d0
    bp(n) =  fator * ( pl(n-1) - pl(n-2) )
    !bp(n) =  2.0d0*pl_out

  end subroutine coeficientes_fontes_massa_cds_uds

!-------------------------------------------------

end module Solucao_CDS_UDS

