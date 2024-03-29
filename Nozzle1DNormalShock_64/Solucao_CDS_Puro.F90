module Solucao_CDS_Puro
  
use variaveis  
!use dados
  
implicit none
  
contains

!-------------------------------------------------


  subroutine coeficientes_e_fontes_qml_cds_p
	real*16  :: fator, dx ! auxiliar
	integer ::i
	! volume 1 (fict�cio)

    fator = 2.0q0 * xp(2) / ( xp(3) - xp(2) )
    aw(1) =  0.0q0
    ap(1) =  1.0q0
    ae(1) = -1.0q0
    bp(1) = -fator * ( u(3) - u(2) )
    bf(1) = 0.0q0
    bc(1) = 0.0q0
    ! volumes internos

    do i = 2, n-1

       dx = xe(i) - xe(i-1)

       aw(i) = -0.5q0*roe(i-1) * ue(i-1) * se(i-1)

       ae(i) = 0.5q0*roe(i) * ue(i) * se(i)

       ap(i) = ro(i) * sp(i) * dx / dt + ( aw(i) + ae(i) )

       bp(i) = ro_o(i) * sp(i) * dx * u_o(i) / dt  &
             - 0.5q0  * sp(i) * ( p(i+1) - p(i-1) )
             

    end do

    ! volume n (fict�cio)

    fator = 2.0q0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) = -1.0q0
    ap(n) =  1.0q0
    ae(n) =  0.0q0
    bp(n) =  fator * ( u(n-1) - u(n-2) )
    bf(n) =  0.0q0
    bc(n) = 0.0q0
    !bp = bp + bc
    
  end subroutine coeficientes_e_fontes_qml_cds_p

!-------------------------------------------------
  
  subroutine calculo_velocidades_face_cds_p

	real*16 :: massa_e ! massa do volume E
    real*16 :: massa_p ! massa do volume P
    real*16 :: somap, somae
    integer ::i
    ! m�dia antiga

    do i = 2, n-2

       massa_p = ro_o(i) * sp(i) * (xe(i) - xe(i-1))

       massa_e = ro_o(i+1) * sp(i+1) * (xe(i+1) - xe(i))

       somap = aw(i)*u(i-1) + ae(i)*u(i+1)
 
       somae = aw(i+1)*u(i) + ae(i+1)*u(i+2)

       ue(i) = (-somap - somae + (massa_p+massa_e)*ue_o(i)/dt - 2.0q0*se(i)*(p(i+1)-p(i)) &
                -0.25q0 * beta * roe(i) * ue(i) * se(i) * ( u_o(i+1) - u_o(i) ) &
                -0.25q0*beta*roe(i+1) * ue(i+1) * se(i+1) * ( u_o(i+2) - u_o(i+1) ) &
                +0.25q0*beta*roe(i-1) * ue(i-1) * se(i-1) * ( u_o(i) - u_o(i-1) )   &
                +0.25q0*beta*roe(i) * ue(i) * se(i) * ( u_o(i) - u_o(i-1) ) ) &
               / (ap(i)+ap(i+1))

    end do

    ue(1) = 0.5q0 * ( u(1) + u(2) )

    ue(n-1) = 0.5q0 * ( u(n) + u(n-1) )
 
  end subroutine calculo_velocidades_face_cds_p

!-------------------------------------------------


  subroutine coeficientes_e_fontes_energia_cds_p
	real*16  :: fator, dx ! auxiliar
	integer ::i

	! volume 1 (fict�cio)
    aw(1) = 0.0q0
    ap(1) = 1.0q0
    ae(1) = 1.0q0
    bp(1) = 2.0q0 * T_in
    bc(1) = 0.0q0
    ! volumes internos

    do i = 2, n-1

       dx = xe(i) - xe(i-1)

       aw(i) = -0.5q0*cp * roe(i-1) * ue(i-1) * se(i-1)

       ae(i) = 0.5q0*cp * roe(i) * ue(i) * se(i)

       ap(i) = cp * ro(i) * sp(i) * dx / dt + cp * roe(i) * ue(i) * se(i) + ( aw(i) + ae(i) )

       bp(i) = cp * ro_o(i) * sp(i) * dx * T_o(i) / dt        &
             + sp(i) * dx * ( p(i) - p_o(i) ) / dt              &
             + sp(i) * u(i) * ( p(i+1) - p(i-1) ) * 0.5q0 

    end do

    ! volume n (fict�cio)
    fator = 2.0q0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) =  -1.0q0
    ap(n) =  1.0q0
    ae(n) =  0.0q0
    bc(n) = 0.0q0
    bp(n) =  fator * ( T(n-1) - T(n-2) )

  end subroutine coeficientes_e_fontes_energia_cds_p
  
!-------------------------------------------------

  subroutine coeficientes_fontes_massa_cds_p

    real*16 :: dx ! auxiliar
    integer ::i 
    ! volume 1 (fict�cio)

    bc(1) = 0.0q0
    aw(1) = 0.0q0
    ap(1) = 1.0q0
    ae(1) = 1.0q0
    bp(1) = 2.0q0 * pl_in

    ! volumes internos
    do i = 2, n-1

       dx = xe(i) - xe(i-1)

       aw(i) = - 0.5q0*(ro(i)+ro(i-1)) * de(i-1) * se(i-1)          &
               - ue(i-1)  * se(i-1) / ( R * T(i-1) )

       ae(i) = - 0.5q0*(ro(i)+ro(i+1)) * de(i) * se(i)      &
               - ue(i-1)  * se(i-1) / ( R * T(i-1) )  

       ap(i) = ( sp(i) * dx / dt + ue(i) * se(i) ) / ( R * T(i) )   &
             - 0.5q0*(ro(i)+ro(i-1)) * de(i-1) * se(i-1)                              &
             + 0.5q0*(ro(i)+ro(i+1))   * de(i)   * se(i)

       bp(i) = ro_o(i) * sp(i) * dx / dt              &
               +0.5q0*(ro_o(i)+ro_o(i+1))*ue_o(i)*se(i) &
               -0.5q0*(ro_o(i)+ro_o(i-1))*ue_o(i-1)*se(i-1) &
               -(sp(i)*dx/dt+ue_o(i)*Se(i))*ro_o(i) &
               +Se(i-1)*ue_o(i-1)*ro_o(i-1) & ! - daqui pra baixo
               -0.5q0*Se(i)*(ro_o(i)+ro_o(i+1))*ue_o(i) &
               +0.5q0*Se(i-1)*(ro_o(i)+ro_o(i-1))*ue_o(i-1)
       
    end do

 ! volume n (fict�cio)
    bc(n) = 0.0q0
    !fator = 2.0q0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) = 1.0q0
    ap(n) =  1.0q0
    ae(n) =  0.0q0
    !bp(n) =  fator * ( pl(n-1) - pl(n-2) )
    bp(n) =  2.0q0*pl_out

  end subroutine coeficientes_fontes_massa_cds_p

!-------------------------------------------------

end module Solucao_CDS_Puro

