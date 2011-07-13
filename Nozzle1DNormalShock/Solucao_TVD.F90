module Solucao_TVD
  
use variaveis  
!use dados
  
implicit none
  
contains

!-------------------------------------------------
function PSI(R)
real*8 PSI
real*8 R
!Van Leer
PSI = (R+ABS(R))/(1.0d0+R)
!Van Albada
!PSI = (R+R**2.0d0)/(1.0d0+R**2.0d0)
!Min-Mod
!if (R <= 0.0d0) then
!    PSI = 0.0d0
!else
!    PSI = MIN(R,1.0d0)
!end if
!Superbee de Roe
!    PSI = MAX(0.0d0,MIN(2.0d0*R,1.0d0),MIN(R,2.0d0))
!Sweby
!    PSI = MAX(0.0d0,MIN(BETA_TVD*R,1.0d0),MIN(R,BETA_TVD))
!QUICK
!    PSI = MAX(0.0d0,MIN(2.0d0*R,(3.0d0+R)/4.0d0,2.0d0))
!UMIST
!    PSI = MAX(0.0d0,MIN(2.0d0*R,(1.0d0+3.0d0*R)/4.0d0,(3.0d0+R)/4.0d0,2))
end function PSI

function re(i,F)
real*8 ::re
integer ::i
real*8, DIMENSION(N) ::F
    re = (F(i)-F(i-1))/(F(i+1)-F(i))
end function re

function rw(i, F)
real*8 ::rw
integer ::i
real*8, DIMENSION(N) ::F
    rw = (F(i-1)-F(i-2))/(F(i)-F(i-1))
end function rw

  subroutine coeficientes_e_fontes_qml_tvd
	real*8  :: fator, dx ! auxiliar
	real*8 :: oeste, leste ! auxiliares
	integer ::i
	! volume 1 (fictício)

    fator = 2.0d0 * xp(2) / ( xp(3) - xp(2) )
    aw(1) =  0.0d0
    ap(1) =  1.0d0
    ae(1) = -1.0d0
    bp(1) = -fator * ( u(3) - u(2) )
    
    ! volumes internos
    do i = 2, n-1

       dx = xe(i) - xe(i-1)

       aw(i) = -roe(i-1) * ue(i-1) * se(i-1)

       ae(i) = 0.0d0

       ap(i) = ro(i) * sp(i) * dx / dt + ( roe(i) * ue(i) * se(i) )

       bp(i) = ro_o(i) * sp(i) * dx * u_o(i) / dt  &
             - 0.5d0   * sp(i) * ( p(i+1) - p(i-1) ) &
             + 0.5d0   * roe(i-1) * ue(i-1) * se(i-1) * PSI(rw(i,u)) * ( u_o(i) - u_o(i-1) ) &
             - 0.5d0   * roe(i) * ue(i) * se(i) * PSI(re(i,u)) * ( u_o(i) - u_o(i-1) )
    end do

    ! volume n (fictício)
    fator = 2.0d0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) = -1.0d0
    ap(n) =  1.0d0
    ae(n) =  0.0d0
    bp(n) =  fator * ( u(n-1) - u(n-2) )
    
  end subroutine coeficientes_e_fontes_qml_tvd

!-------------------------------------------------
  
  subroutine calculo_velocidades_face_tvd

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

       ue(i) = (-somap - somae + (massa_p+massa_e)*ue_o(i)/dt - 2.0d0*se(i)*(p(i+1)-p(i)) &
                -0.25d0 * beta * roe(i) * ue(i) * se(i) * ( u_o(i+1) - u_o(i) ) &
                -0.25d0*beta*roe(i+1) * ue(i+1) * se(i+1) * ( u_o(i+2) - u_o(i+1) ) &
                +0.25d0*beta*roe(i-1) * ue(i-1) * se(i-1) * ( u_o(i) - u_o(i-1) )   &
                +0.25d0*beta*roe(i) * ue(i) * se(i) * ( u_o(i) - u_o(i-1) ) ) &
               / (ap(i)+ap(i+1))

    end do

    ue(1) = 0.5d0 * ( u(1) + u(2) )

    ue(n-1) = 0.5d0 * ( u(n) + u(n-1) )
    
  end subroutine calculo_velocidades_face_tvd

!-------------------------------------------------


  subroutine coeficientes_e_fontes_energia_tvd
	real*8  :: fator, dx ! auxiliar
	real*8 :: oeste, leste ! auxiliares
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

       bp(i) = cp * ro_o(i) * sp(i) * dx * T_o(i) / dt           &
             + sp(i) * dx * ( p(i) - p_o(i) ) / dt               &
             + sp(i) * u(i) * ( p(i+1) - p(i-1) ) * 0.5d0        &
             + 0.5d0 * cp * roe(i-1) * ue(i-1) * se(i-1) * PSI(rw(i,T)) * ( T(i) - T(i-1) )  &
             - 0.5d0 * cp * roe(i) * ue(i) * se(i) * PSI(re(i,T)) * ( T(i+1) - T(i) )

    end do

    ! volume n (fictício)
    fator = 2.0d0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) =  -1.0d0
    ap(n) =  1.0d0
    ae(n) =  0.0d0
    bp(n) =  fator * ( T(n-1) - T(n-2) )
    !bp(n) =  2.0d0*T_out
! volume n (fictício)

  end subroutine coeficientes_e_fontes_energia_tvd
  
!-------------------------------------------------

  subroutine coeficientes_fontes_massa_tvd

    real*8 :: fator, dx ! auxiliar
    real*8 :: oeste, leste ! auxiliares
    integer ::i 
    
    ! volume 1 (fictício)
    aw(1) = 0.0d0
    ap(1) = 1.0d0
    ae(1) = 1.0d0
    bp(1) = 2.0d0 * pl_in

    ! volumes internos
    do i = 2, n-1

       dx = xe(i) - xe(i-1)

       aw(i) = - 0.5d0*(ro(i)+ro(i-1)) * de(i-1) * se(i-1)          &
               - ue(i-1)  * se(i-1) / ( R * T(i-1) )

       ae(i) = - 0.5d0*(ro(i)+ro(i+1)) * de(i) * se(i)       

       ap(i) = ( sp(i) * dx / dt + ue(i) * se(i) ) / ( R * T(i) )   &
             - 0.5d0*(ro(i)+ro(i-1)) * de(i-1) * se(i-1)                              &
             + 0.5d0*(ro(i)+ro(i+1))   * de(i)   * se(i)

       bp(i) = ro_o(i) * sp(i) * dx / dt              &
               +0.5d0*(ro_o(i)+ro_o(i+1))*ue_o(i)*se(i) &
               -0.5d0*(ro_o(i)+ro_o(i-1))*ue_o(i-1)*se(i-1) &
               -0.5d0*Se(i)*PSI(re(i,ro))*(ro_o(i)+ro_o(i+1))*ue_o(i) &
               +0.5d0*Se(i-1)*PSI(rw(i,ro))*(ro_o(i)+ro_o(i-1))*ue_o(i-1) &
               -(Sp(i)*dx/dt+Se(i)*ue_o(i))*ro_o(i) & ! - daqui pra baixo
               -Se(i-1)*ue_o(i-1)*ro_o(i-1) & 
               -0.5d0*(ro_o(i)+ro_o(i+1))*ue_o(i)*se(i) &
               +0.5d0*(ro_o(i)+ro_o(i-1))*ue_o(i-1)*se(i-1)
    end do

 ! volume n (fictício)
    !fator = 2.0d0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) = 1.0d0
    ap(n) =  1.0d0
    ae(n) =  0.0d0
    !bp(n) =  fator * ( pl(n-1) - pl(n-2) )
    bp(n) =  2.0d0*pl_out

  end subroutine coeficientes_fontes_massa_tvd

!-------------------------------------------------

end module Solucao_TVD

