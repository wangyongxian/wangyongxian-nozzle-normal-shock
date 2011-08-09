module Solucao_TVD
  
use variaveis  
!use dados
  
implicit none
  
contains

!-------------------------------------------------
function PSI(R)
real*16 PSI
real*16 R
real*16 BETA_TVD
real*16 ALPHA_TVD
real*16 GAMA_TVD
!Van Leer
!PSI = (R+ABS(R))/(1.0q0+R)
!Van Albada 1
!PSI = (R+R**2.0q0)/(1.0q0+R**2.0q0)
!Van Albada 2
!PSI = 2.0q0*R/(R**2.0q0+1.0q0)
!Min-Mod
!if (R <= 0.0q0) then
!    PSI = 0.0q0
!else
!    PSI = MIN(R,1.0q0)
!end if
!Superbee de Roe
    !PSI = MAX(0.0q0,MIN(2.0q0*R,1.0q0),MIN(R,2.0q0))
!Sweby muito promissor
!1.0q0<=BETA_TVD <=2.0q0
!BETA_TVD = 1.0q0
!    PSI = MAX(0.0q0,MIN(BETA_TVD*R,1.0q0),MIN(R,BETA_TVD))
!QUICK
!    PSI = MAX(0.0q0,MIN(2.0q0*R,(3.0q0+R)/4.0q0,2.0q0))
!UMIST
    PSI = MAX(0.0q0,MIN(2.0q0*R,(1.0q0+3.0q0*R)/4.0q0,(3.0q0+R)/4.0q0,2.0q0))
!http://en.wikipedia.org/wiki/Flux_limiter
!Koren ruim
!PSI = max(0.0q0,min(2.0q0*R,(1.0q0+2.0q0*R)/3.0q0,2.0q0))
!HQUICK nao funcionou
!PSI = 2.0q0*(R+ABS(R))/(R+3.0q0)
!HCUS nao funcionou
!PSI = 1.5q0*(R+ABS(R))/(R+2.0q0) 
!CHARM parece promissor
!if(R>0.0q0)then
!PSI = R*(3.0q0*R+1.0q0)/(R+1)**2.0q0
!else
!PSI = 0.0q0
!end if
!monotonized central (MC) +-
!PSI = max(0.0q0,min(2.0q0*R,0.5*(1.0q0+R),2.0q0))
!OSHER promissor
!1.0q0<=BETA_TVD <=2.0q0
BETA_TVD=1.9q0
!PSI = max(0.0q0,min(R,BETA_TVD))
!ospre promissor
!PSI = 1.5q0*(R**2.0q0+R)/(R**2.0q0+R+1.0q0)
!smart muito promissor
!PSI = max(0.0q0,min(2.0q0*R,(0.25q0+0.75q0*R),4.0q0))
!TOPUS parece promissor
!ALPHA_TVD = -2.0q0, 0.0q0, 2.0q0
ALPHA_TVD = -2.0q0
!PSI = max(0.0q0,0.5q0*(ABS(R)+R)*((-0.5q0+1.0q0)*R**2.0q0+(ALPHA_TVD+4.0q0)*R+(-0.5q0*ALPHA_TVD+3.0q0))/(1.0q0+ABS(R))**3.0q0)
!FDPUS-C1 parece promissor
!PSI = max(0.0q0,0.5q0*(ABS(R)+R)*(4.0q0*R**2.0q0+12.0*R)/(1.0q0+ABS(R))**4.0q0)
!SDPUS-C1 parece promissor
!GAMA - 4.0q0, 6.0q0,8.0q0,10.0q0,12.0q0
GAMA_TVD = 12.0q0
!PSI = max(0.0q0,0.5q0*(ABS(R)+R)*((-8.0q0+2.0q0*GAMA_TVD)*R**3.0q0+(40.0q0-4.0q0*GAMA_TVD)*R**2.0q0+2.0q0*GAMA_TVD*R)/(1.0q0+ABS(R))**5.0q0)
end function PSI

function re(i,F)
real*16 ::re
logical a
integer ::i
real*16, DIMENSION(N) ::F
if (F(i+1)==F(i)) then
    re = 1.0q0
else
    re = (F(i)-F(i-1))/(F(i+1)-F(i))
end if
a = isnan(re)
if (a == .true.) then
re = 1.0q0
end if
!soma = F(i+1)-F(i)
end function re

function rw(i, F)
real*16 ::rw
integer ::i
real*16, DIMENSION(N) ::F
if (F(i) == F(i-1)) then
    rw=1.0q0
else
    rw = (F(i-1)-F(i-2))/(F(i)-F(i-1))
end if

end function rw

  subroutine coeficientes_e_fontes_qml_tvd
	real*16  :: fator, dx ! auxiliar
	integer ::i
	! volume 1 (fictício)

    fator = 2.0q0 * xp(2) / ( xp(3) - xp(2) )
    aw(1) =  0.0q0
    ap(1) =  1.0q0
    ae(1) = -1.0q0
    bp(1) = -fator * ( u(3) - u(2) )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    i = 2
    dx = xe(i) - xe(i-1)

       aw(i) = -roe(i-1) * ue(i-1) * se(i-1)

       ae(i) = 0.0q0

       ap(i) = ro(i) * sp(i) * dx / dt + ( roe(i) * ue(i) * se(i) )
        
       bp(i) = ro_o(i) * sp(i) * dx * u_o(i) / dt  &
             - 0.5q0   * sp(i) * ( p(i+1) - p(i-1) ) &
             + 0.5q0   * roe(i-1) * ue(i-1) * se(i-1) * PSI(1.0q0) * ( u(i) - u(i-1) ) &
             - 0.5q0   * roe(i) * ue(i) * se(i) * PSI(re(i,u)) * ( -u(i) + u(i-1) )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! volumes internos
    do i = 3, n-1

       dx = xe(i) - xe(i-1)

       aw(i) = -roe(i-1) * ue(i-1) * se(i-1)

       ae(i) = 0.0q0

       ap(i) = ro(i) * sp(i) * dx / dt + ( roe(i) * ue(i) * se(i) )

       bp(i) = ro_o(i) * sp(i) * dx * u_o(i) / dt  &
             - 0.5q0   * sp(i) * ( p(i+1) - p(i-1) ) &
             + 0.5q0   * roe(i-1) * ue(i-1) * se(i-1) * PSI(rw(i,u)) * ( u(i) - u(i-1) ) &
             - 0.5q0   * roe(i) * ue(i) * se(i) * PSI(re(i,u)) * ( -u(i) + u(i+1) )
    end do

    ! volume n (fictício)
    fator = 2.0q0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) = -1.0q0
    ap(n) =  1.0q0
    ae(n) =  0.0q0
    bp(n) =  fator * ( u(n-1) - u(n-2) )
    
  end subroutine coeficientes_e_fontes_qml_tvd

!-------------------------------------------------
  
  subroutine calculo_velocidades_face_tvd

	real*16 :: massa_e ! massa do volume E
    real*16 :: massa_p ! massa do volume P
    real*16 :: somap, somae
    integer ::i
    ! média antiga
    i=2
       massa_p = ro_o(i) * sp(i) * (xe(i) - xe(i-1))

       massa_e = ro_o(i+1) * sp(i+1) * (xe(i+1) - xe(i))

       somap = aw(i)*u(i-1) + ae(i)*u(i+1)
 
       somae = aw(i+1)*u(i) + ae(i+1)*u(i+2)

       ue(i) = (-somap - somae + (massa_p+massa_e)*ue_o(i)/dt - 2.0q0*se(i)*(p(i+1)-p(i)) &
                -0.5q0 * PSI(re(i,u_o)) * roe(i) * ue(i) * se(i) * ( u(i+1) - u(i) ) &
                -0.5q0 * PSI(re(i,u_o)) * roe(i+1) * ue(i+1) * se(i+1) * ( u(i+2) - u(i+1) ) &
                +0.5q0 * PSI(1.0q0) * roe(i-1) * ue(i-1) * se(i-1) * ( u(i) - u(i-1) )   &
                +0.5q0 * PSI(1.0q0) * roe(i) * ue(i) * se(i) * ( u(i) - u(i-1) ) ) &
               / (ap(i)+ap(i+1))
               

    do i = 3, n-2

       massa_p = ro_o(i) * sp(i) * (xe(i) - xe(i-1))

       massa_e = ro_o(i+1) * sp(i+1) * (xe(i+1) - xe(i))

       somap = aw(i)*u(i-1) + ae(i)*u(i+1)
 
       somae = aw(i+1)*u(i) + ae(i+1)*u(i+2)

       ue(i) = (-somap - somae + (massa_p+massa_e)*ue_o(i)/dt - 2.0q0*se(i)*(p(i+1)-p(i)) &
                -0.5q0 * PSI(re(i,u_o)) * roe(i) * ue(i) * se(i) * ( u(i+1) - u(i) ) &
                -0.5q0 * PSI(re(i,u_o)) * roe(i+1) * ue(i+1) * se(i+1) * ( u(i+2) - u(i+1) ) &
                +0.5q0 * PSI(rw(i,u_o)) * roe(i-1) * ue(i-1) * se(i-1) * ( u(i) - u(i-1) )   &
                +0.5q0 * PSI(rw(i,u_o)) * roe(i) * ue(i) * se(i) * ( u(i) - u(i-1) ) ) &
               / (ap(i)+ap(i+1))

    end do

    ue(1) = 0.5q0 * ( u(1) + u(2) )

    ue(n-1) = 0.5q0 * ( u(n) + u(n-1) )
    
  end subroutine calculo_velocidades_face_tvd

!-------------------------------------------------


  subroutine coeficientes_e_fontes_energia_tvd
	real*16  :: fator, dx ! auxiliar
	integer ::i
	! volume 1 (fictício)

    aw(1) = 0.0q0
    ap(1) = 1.0q0
    ae(1) = 1.0q0
    bp(1) = 2.0q0 * T_in
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    i=2
    dx = xe(i) - xe(i-1)

       aw(i) = -cp * roe(i-1) * ue(i-1) * se(i-1)

       ae(i) = 0.0q0

       ap(i) = cp * ro(i) * sp(i) * dx / dt + cp * roe(i) * ue(i) * se(i)

       bp(i) = cp * ro_o(i) * sp(i) * dx * T_o(i) / dt           &
             + sp(i) * dx * ( p(i) - p_o(i) ) / dt               &
             + 0.5q0 * sp(i) * u(i) * ( p(i+1) - p(i-1) )        &
             + 0.5q0 * cp * roe(i-1) * ue(i-1) * se(i-1) * PSI(1.0q0) * ( T(i) - T(i-1) )  &
             - 0.5q0 * cp * roe(i) * ue(i) * se(i) * PSI(re(i,T)) * ( T(i+1) - T(i) )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! volumes internos
    do i = 3, n-1

       dx = xe(i) - xe(i-1)

       aw(i) = -cp * roe(i-1) * ue(i-1) * se(i-1)

       ae(i) = 0.0q0

       ap(i) = cp * ro(i) * sp(i) * dx / dt + cp * roe(i) * ue(i) * se(i)

       bp(i) = cp * ro_o(i) * sp(i) * dx * T_o(i) / dt           &
             + sp(i) * dx * ( p(i) - p_o(i) ) / dt               &
             + 0.5q0 * sp(i) * u(i) * ( p(i+1) - p(i-1) )         &
             + 0.5q0 * cp * roe(i-1) * ue(i-1) * se(i-1) * PSI(rw(i,T)) * ( T(i) - T(i-1) )  &
             - 0.5q0 * cp * roe(i) * ue(i) * se(i) * PSI(re(i,T)) * ( T(i+1) - T(i) )

    end do

    ! volume n (fictício)
    fator = 2.0q0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) =  -1.0q0
    ap(n) =  1.0q0
    ae(n) =  0.0q0
    bp(n) =  fator * ( T(n-1) - T(n-2) )
    !bp(n) =  2.0q0*T_out
! volume n (fictício)

  end subroutine coeficientes_e_fontes_energia_tvd
  
!-------------------------------------------------

  subroutine coeficientes_fontes_massa_tvd

    real*16 :: dx, fator ! auxiliar
    integer ::i 
    real*16 ::Mp_ro,Me_ro,Mw_ro,Me_u,Mw_u,bpn
    
    ! volume 1 (fictício)
    aw(1) = 0.0q0
    ap(1) = 1.0q0
    ae(1) = 1.0q0
    bp(1) = 2.0q0 * pl_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
i=2
dx = xe(i) - xe(i-1)
    Mp_ro = sp(i)*dx/dt + ue(i)*se(i)
    Me_ro = 0.0q0
    Mw_ro = -ue(i-1)*se(i-1)
    Me_u = se(i)*(ro(i) + 0.5q0*PSI(1.0q0)*(ro(i+1)-ro(i)))
    Mw_u = -se(i-1)*(ro(i-1)+0.5q0*PSI(1.0q0)*(ro(i)-ro(i-1)))
    bpn = ro_o(i)*sp(i)*dx/dt + se(i)*ue(i)*(ro(i)+0.5q0*PSI(1.0q0)*(ro(i+1)-ro(i))) &
          -se(i-1)*ue(i-1)*(ro(i-1)+0.5q0*PSI(1.0q0)*(ro(i)-ro(i-1)))
         
       aw(i) = Mw_u * de(i-1) &
               + Mw_ro / ( R * T(i-1) )

       ae(i) = Me_ro / (R* T(i+1)) &
               - Me_u * de(i) 

       ap(i) = Mp_ro / ( R * T(i) )   &
              + Me_u * de(i)    &
              - Mw_u * de(i-1)

       bp(i) = bpn  &
               -(Mp_ro*ro(i)  &
               +Me_ro*ro(i+1) &
               +Mw_ro*ro(i-1) &
               +Me_u*ue(i)    &
               +Mw_u*ue(i-1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! volumes internos
    do i = 3, n-1
        dx = xe(i) - xe(i-1)

    Mp_ro = sp(i)*dx/dt + ue(i)*se(i)
    Me_ro = 0.0q0
    Mw_ro = -ue(i-1)*se(i-1)
    Me_u = se(i)*(ro(i) + 0.5q0*PSI(re(i,ro))*(ro(i+1)-ro(i)))
    Mw_u = -se(i-1)*(ro(i-1)+0.5q0*PSI(rw(i,ro))*(ro(i)-ro(i-1)))
    bpn = ro_o(i)*sp(i)*dx/dt + se(i)*ue(i)*(ro(i)+0.5q0*PSI(re(i,ro))*(ro(i+1)-ro(i))) &
          -se(i-1)*ue(i-1)*(ro(i-1)+0.5q0*PSI(rw(i,ro))*(ro(i)-ro(i-1)))
         
       aw(i) = Mw_u * de(i-1) &
               + Mw_ro / ( R * T(i-1) )

       ae(i) = Me_ro / (R* T(i+1)) &
               - Me_u * de(i) 

       ap(i) = Mp_ro / ( R * T(i) )   &
              + Me_u * de(i)    &
              - Mw_u * de(i-1)

       bp(i) = bpn  &
               -(Mp_ro*ro(i)  &
               +Me_ro*ro(i+1) &
               +Mw_ro*ro(i-1) &
               +Me_u*ue(i)    &
               +Mw_u*ue(i-1))
 
    end do

 ! volume n (fictício)
    fator = 2.0q0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    aw(n) = -1.0q0
    ap(n) =  1.0q0
    ae(n) =  0.0q0
    bp(n) =  fator * ( pl(n-1) - pl(n-2) )
    !bp(n) =  2.0q0*pl_out

  end subroutine coeficientes_fontes_massa_tvd

!-------------------------------------------------

subroutine calculo_massa_especifica_nas_faces_tvd
    integer ::i
    
    roe(1) = ro(1) + 0.5q0*PSI(1.0q0)*( ro(1+1) - ro(1) )
    
    do i = 2, n-1
       roe(i) = ro(i) + 0.5q0*PSI(re(i,ro))*( ro(i+1) - ro(i) )
    end do
    
end subroutine calculo_massa_especifica_nas_faces_tvd

end module Solucao_TVD

