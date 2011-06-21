Module NormalShock1D

use variaveis
implicit none

contains

subroutine PrandtRelation(Mach1, Mach2)
    real*8, intent(in) ::Mach1
    real*8, intent(out) ::Mach2
    
    Mach2 = 1.0d0/Mach1
    
end subroutine PrandtRelation

subroutine Mach2Calc(gama, Mach1, Mach2)
    
    real*8, intent(in)  ::gama
    real*8, intent(in)  ::Mach1
    real*8, intent(out)::Mach2
    Mach2 = dsqrt((1.0d0+((gama-1.0d0)/2.0d0)*Mach1**2.0d0)/(gama*Mach1**2.0d0-(gama-1.0d0)/2.0d0))

end subroutine Mach2Calc
    
subroutine P2Calc(gama, Mach1, P1, P2)
    real*8, intent(in) ::gama
    real*8, intent(in) ::Mach1
    real*8, intent(in) ::P1
    real*8, intent(out) ::P2
    
    P2 = (1.0d0+2.0d0*gama*(Mach1**2.0d0-1.0d0)/(gama+1.0d0))*P1

end subroutine

subroutine RO2Calc(gama, Mach1, RO1, RO2)
    real*8, intent(in) ::gama
    real*8, intent(in) ::Mach1
    real*8, intent(in) ::RO1
    real*8, intent(out) ::RO2
    
    RO2 = (((gama+1.0d0)*Mach1**2.0d0)/(2.0d0+(gama-1.0d0)*Mach1**2.0d0))*RO1

end subroutine

subroutine T2Calc(gama, Mach1, T1, T2)
    real*8, intent(in) ::gama
    real*8, intent(in) ::Mach1
    real*8, intent(in) ::T1
    real*8, intent(out) ::T2
    
    T2 = ((1.0d0+2.0d0*gama*(Mach1**2.0d0-1.0d0)/(gama+1.0d0))*1.0d0/(((gama+1.0d0)*Mach1**2.0d0)/(2.0d0+(gama-1.0d0)*Mach1**2.0d0)))*T1

end subroutine

subroutine U2Calc(gama, Mach1, U1, U2)
    real*8, intent(in) ::gama
    real*8, intent(in) ::Mach1
    real*8, intent(in) ::U1
    real*8, intent(out) ::U2

    U2 = 1.0d0/((((gama+1.0d0)*Mach1**2.0d0)/(2.0d0+(gama-1.0d0)*Mach1**2.0d0))*U1)

end subroutine

subroutine ShockLocationCalc(pe, p01, Ae, At, ratio)
    real*8, intent(in) ::pe
    real*8, intent(in) ::p01
    real*8, intent(in) ::Ae
    real*8, intent(in) ::At
    real*8, intent(out)::ratio
    
    real*8 ::factor
    factor = pe*Ae/(p01*At)
    
    Me = dsqrt(-1.0d0/(gama-1.0d0)+dsqrt(1.0d0/(gama-1.0d0)**2.0d0 + &
                                    (2.0d0/(gama-1.0d0))*((2.0d0/(gama+1.0d0))**((gama+1.0d0)/(gama-1.0d0)))* & 
                                    ((factor)**2.0d0)) )
    

end subroutine

subroutine P0Calc(gama, P, Mach, P0)
    real*8, intent(in) ::gama
    real*8, intent(in) ::P
    real*8, intent(in) ::Mach
    real*8, intent(out) ::P0
    
    P0 = ((1.0d0+((gama-1.0d0)/2.0d0)*Mach**2.0d0)**(gama/(gama-1.0d0)))*P

end subroutine P0Calc

subroutine T0Calc(gama, T, Mach, T0)
    real*8, intent(in) ::gama
    real*8, intent(in) ::T
    real*8, intent(in) ::Mach
    real*8, intent(out) ::T0
    
    T0 = (1.0d0+((gama-1.0d0)/2.0d0)*Mach**2.0d0)*T

end subroutine T0Calc

subroutine Ro0Calc(gama, Ro, Mach, Ro0)
    real*8, intent(in) ::gama
    real*8, intent(in) ::Ro
    real*8, intent(in) ::Mach
    real*8, intent(out) ::Ro0
    
    Ro0 = ((1.0d0+((gama-1.0d0)/2.0d0)*Mach**2.0d0)**(gama/(gama-1.0d0)))*Ro

end subroutine Ro0Calc


end module NormalShock1D

  !subroutine GIBBS_MACH1D_envia6 ( cpt, gamat, Rt, Yit )
   ! real*8,              intent(out) :: cpt, gamat, Rt
    !real*8,dimension(Ne),intent(out) :: Yit
    !real*8,dimension(Ne),intent(inout) :: Yit
    !real*8,dimension(Ne),intent(in) :: Yit
!  end subroutine GIBBS_MACH1D_envia6