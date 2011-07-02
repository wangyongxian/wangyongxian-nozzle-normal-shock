Module NormalShock1D

use variaveis
implicit none

contains

subroutine PrandtRelation(Mach1, Mach2)
    real*8, intent(in)  ::Mach1
    real*8, intent(out) ::Mach2
    
    Mach2 = 1.0d0/Mach1
    
end subroutine PrandtRelation

subroutine Mach2Calc(gama, Mach1, Mach2)
    
    real*8, intent(in)  ::gama
    real*8, intent(in)  ::Mach1
    real*8, intent(out) ::Mach2
    Mach2 = dsqrt((1.0d0+((gama-1.0d0)/2.0d0)*Mach1**2.0d0)/(gama*Mach1**2.0d0-(gama-1.0d0)/2.0d0))

end subroutine Mach2Calc

subroutine Mach1Calc(gama, Mach1, Mach2)
    
    real*8, intent(in)  ::gama
    real*8, intent(out) ::Mach1
    real*8, intent(in)  ::Mach2
    Mach1 = dsqrt((1.0d0+((gama-1.0d0)/2.0d0)*Mach2**2.0d0)/(gama*Mach2**2.0d0-(gama-1.0d0)/2.0d0))

end subroutine Mach1Calc
    
subroutine P2Calc(gama, Mach1, P1, P2)
    real*8, intent(in)  ::gama
    real*8, intent(in)  ::Mach1
    real*8, intent(in)  ::P1
    real*8, intent(out) ::P2
    
    P2 = (1.0d0+2.0d0*gama*(Mach1**2.0d0-1.0d0)/(gama+1.0d0))*P1

end subroutine

subroutine RO2Calc(gama, Mach1, Ro1, Ro2)
    real*8, intent(in) ::gama
    real*8, intent(in) ::Mach1
    real*8, intent(in) ::Ro1
    real*8, intent(out) ::Ro2
    
    Ro2 = (((gama+1.0d0)*Mach1**2.0d0)/(2.0d0+(gama-1.0d0)*Mach1**2.0d0))*Ro1

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
    real*8, intent( out) ::U2

    U2 = 1.0d0/((((gama+1.0d0)*Mach1**2.0d0)/(2.0d0+(gama-1.0d0)*Mach1**2.0d0))*U1)

end subroutine

subroutine ShockLocationCalc(gama, Pe, P01, P02, Ae, A, AA, At)
    real*8, intent(in) ::gama
    real*8, intent(in) ::Pe
    real*8, intent(in) ::P01
    real*8, intent(inout) ::P02
    real*8, intent(in) ::Ae
    real*8, intent(in) ::A
    real*8, intent(out)::AA
    real*8, intent(out)::At
    real*8 ::Me
    real*8 ::factor
    real*8 ::Mach
    
    factor = pe*A/(p01)
    !Calcula o mach de saida
    Me = dsqrt(-1.0d0/(gama-1.0d0)+dsqrt(1.0d0/(gama-1.0d0)**2.0d0 + &
                                    (2.0d0/(gama-1.0d0))*((2.0d0/(gama+1.0d0))**((gama+1.0d0)/(gama-1.0d0)))* & 
                                    ((1.0d0/factor)**2.0d0)) )
    
    !calculo da pressao de estagnacao apos o choque
    call P0Calc(gama, Pe, Me, P02)                           
    
    Mach = ROOT(gama, P01, P02, 10.0d0)
    
    call AARatioCalc(gama, Mach, AA)
    
    call NewACalc(gama, Me, Ae, At)
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

subroutine AARatioCalc(gama, Mach, AA)
    real*8, intent(in)  ::gama
    real*8, intent(in)  ::Mach
    real*8, intent(out) ::AA
    
    AA = dsqrt((1.0d0/Mach**2.0d0)*((2.0d0/(gama+1.0d0))*(1.0d0+((gama-1.0d0)/2.0d0)*Mach**2.0d0))**((gama+1.0d0)/(gama-1.0d0)))

end subroutine AARatioCalc

subroutine NewACalc(gama, Mach, A, At)
    real*8, intent(in)  ::gama
    real*8, intent(in)  ::Mach
    real*8, intent(in)  ::A
    real*8, intent(out) ::At
    
    At = A/dsqrt((1.0d0/Mach**2.0d0)*((2.0d0/(gama+1.0d0))*(1.0d0+((gama-1.0d0)/2.0d0)*Mach**2.0d0))**((gama+1.0d0)/(gama-1.0d0)))

end subroutine NewACalc

subroutine ACalc(gama, R, T, a)
    real*8, intent(in)  ::gama
    real*8, intent(in)  ::R
    real*8, intent(in)  ::T
    real*8, intent(out) ::a
    
    a = dsqrt(gama*R*T)

end subroutine ACalc

subroutine UCalc(Mach, a, u)
    real*8, intent(in)  ::Mach
    real*8, intent(in)  ::a
    real*8, intent(out) ::u
    
    u = Mach*a
    
end subroutine UCalc

REAL*8 FUNCTION ROOT(gama, P01, P02, CHUTE)
    REAL*8 ::F
    REAL*8 ::CHUTE
    REAL*8 ::FL
    REAL*8 ::M
    REAL*8 ::P01
    REAL*8 ::P02
    REAL*8 ::gama
    REAL*8 ::D      !denominador
    REAL*8 ::DL     !derivada do denominador
    REAL*8 ::N      !numerador
    REAL*8 ::NL     !derivada do numerador
    REAL*8 ::M2     !mach2
    REAL*8 ::M2_N   !mach2
    REAL*8 ::M2_NL  !mach2
    REAL*8 ::M2_D   !mach2
    REAL*8 ::M2_DL  !mach2
    REAL*8 ::P      !potencia
  M  = CHUTE
  DO
    N = 1.0d0+0.5d0*(gama-1.0d0)*M**2.0d0
    NL = M*(gama-1.0d0)
    M2_N = (1.0d0+0.5d0*(gama-1.0d0)*M**2.0d0)
    M2_NL = NL
    M2_D = (gama*M**2.0d0-0.5d0*(gama-1.0d0)) 
    M2_DL = 2.0d0*gama*M 
    M2 = M2_N/M2_D
    D = 1.0d0+0.5d0*(gama-1.0d0)*M2
    DL = 0.5d0*(gama-1.0d0)*((M2_NL*M2_D - M2_DL*M2_N)/(M2_D**2.0d0))
    P = gama/(gama-1.0d0)
    F = 1.0d0+2.0d0*gama*(M**2.0d0-1.0d0)/(gama+1.0d0) &
        -(P02/P01)*(N/D)**P
    FL = 4.0d0*gama*M/(gama+1.0d0)-(P02/P01)*(gama/(gama-1))* &
        ((N/D)**(P-1.0d0))*((NL*D - DL*N)/D**2.0d0)
    ROOT  = M - F/FL
    IF (ABS((M-ROOT)/ROOT) .LT. 1E-06) RETURN
    M = ROOT
  END DO
END FUNCTION

subroutine ShockFinder(AA, N, pos)
integer, intent(inout) ::pos
real*8, intent(in)  ::AA
!real*8 ::AAp
integer ::N
integer ::i
do i=N-1, 1, -1
    if (aap(i-1) < AA .and. AA < aap(i+1) ) then
       pos = i 
       exit
    end if
end do
end subroutine ShockFinder

function ThroatFinder()
    integer :: ThroatFinder 
    integer ::i
    ThroatFinder = -1
    do i=1, N
        if ( Se(i) == Pi*rg**2.0d0 ) then
            ThroatFinder = i
            exit
        end if
    end do
end function



end module NormalShock1D

