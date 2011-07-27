Module NormalShock1D

use variaveis
implicit none

contains

subroutine PrandtRelation(Mach1, Mach2)
    real*16, intent(in)  ::Mach1
    real*16, intent(out) ::Mach2
    
    Mach2 = 1.0q0/Mach1
    
end subroutine PrandtRelation

subroutine Mach2Calc(gama, Mach1, Mach2)
    
    real*16, intent(in)  ::gama
    real*16, intent(in)  ::Mach1
    real*16, intent(out) ::Mach2
    Mach2 = qsqrt((1.0q0+((gama-1.0q0)/2.0q0)*Mach1**2.0q0)/(gama*Mach1**2.0q0-(gama-1.0q0)/2.0q0))

end subroutine Mach2Calc

subroutine Mach1Calc(gama, Mach1, Mach2)
    
    real*16, intent(in)  ::gama
    real*16, intent(out) ::Mach1
    real*16, intent(in)  ::Mach2
    Mach1 = qsqrt((1.0q0+((gama-1.0q0)/2.0q0)*Mach2**2.0q0)/(gama*Mach2**2.0q0-(gama-1.0q0)/2.0q0))

end subroutine Mach1Calc
    
subroutine P2Calc(gama, Mach1, P1, P2)
    real*16, intent(in)  ::gama
    real*16, intent(in)  ::Mach1
    real*16, intent(in)  ::P1
    real*16, intent(out) ::P2
    
    P2 = (1.0q0+2.0q0*gama*(Mach1**2.0q0-1.0q0)/(gama+1.0q0))*P1

end subroutine

subroutine RO2Calc(gama, Mach1, Ro1, Ro2)
    real*16, intent(in) ::gama
    real*16, intent(in) ::Mach1
    real*16, intent(in) ::Ro1
    real*16, intent(out) ::Ro2
    
    Ro2 = (((gama+1.0q0)*Mach1**2.0q0)/(2.0q0+(gama-1.0q0)*Mach1**2.0q0))*Ro1

end subroutine

subroutine T2Calc(gama, Mach1, T1, T2)
    real*16, intent(in) ::gama
    real*16, intent(in) ::Mach1
    real*16, intent(in) ::T1
    real*16, intent(out) ::T2
    
    T2 = ((1.0q0+2.0q0*gama*(Mach1**2.0q0-1.0q0)/(gama+1.0q0))*1.0q0/(((gama+1.0q0)*Mach1**2.0q0)/(2.0q0+(gama-1.0q0)*Mach1**2.0q0)))*T1

end subroutine

subroutine U2Calc(gama, Mach1, U1, U2)
    real*16, intent(in) ::gama
    real*16, intent(in) ::Mach1
    real*16, intent(in) ::U1
    real*16, intent( out) ::U2

    U2 = 1.0q0/((((gama+1.0q0)*Mach1**2.0q0)/(2.0q0+(gama-1.0q0)*Mach1**2.0q0))*U1)

end subroutine

subroutine ShockLocationCalc(gama, Pe, P01, P02, Ae, A, AA, At)
    real*16, intent(in) ::gama
    real*16, intent(in) ::Pe
    real*16, intent(in) ::P01
    real*16, intent(inout) ::P02
    real*16, intent(in) ::Ae
    real*16, intent(in) ::A
    real*16, intent(out)::AA
    real*16, intent(out)::At
    real*16 ::Me
    real*16 ::factor
    real*16 ::Mach
    
    factor = pe*A/(p01)
    !Calcula o mach de saida
    Me = qsqrt(-1.0q0/(gama-1.0q0)+qsqrt(1.0q0/(gama-1.0q0)**2.0q0 + &
                                    (2.0q0/(gama-1.0q0))*((2.0q0/(gama+1.0q0))**((gama+1.0q0)/(gama-1.0q0)))* & 
                                    ((1.0q0/factor)**2.0q0)) )
    
    !calculo da pressao de estagnacao apos o choque
    call P0Calc(gama, Pe, Me, P02)                           
    
    Mach = ROOT(gama, P01, P02, 10.0q0)
    
    call AARatioCalc(gama, Mach, AA)
    
    call NewACalc(gama, Me, Ae, At)
    
end subroutine

subroutine P0Calc(gama, P, Mach, P0)
    real*16, intent(in) ::gama
    real*16, intent(in) ::P
    real*16, intent(in) ::Mach
    real*16, intent(out) ::P0
    
    P0 = ((1.0q0+((gama-1.0q0)/2.0q0)*Mach**2.0q0)**(gama/(gama-1.0q0)))*P

end subroutine P0Calc

subroutine T0Calc(gama, T, Mach, T0)
    real*16, intent(in) ::gama
    real*16, intent(in) ::T
    real*16, intent(in) ::Mach
    real*16, intent(out) ::T0
    
    T0 = (1.0q0+((gama-1.0q0)/2.0q0)*Mach**2.0q0)*T

end subroutine T0Calc

subroutine Ro0Calc(gama, Ro, Mach, Ro0)
    real*16, intent(in) ::gama
    real*16, intent(in) ::Ro
    real*16, intent(in) ::Mach
    real*16, intent(out) ::Ro0
    
    Ro0 = ((1.0q0+((gama-1.0q0)/2.0q0)*Mach**2.0q0)**(gama/(gama-1.0q0)))*Ro

end subroutine Ro0Calc

subroutine AARatioCalc(gama, Mach, AA)
    real*16, intent(in)  ::gama
    real*16, intent(in)  ::Mach
    real*16, intent(out) ::AA
    
    AA = qsqrt((1.0q0/Mach**2.0q0)*((2.0q0/(gama+1.0q0))*(1.0q0+((gama-1.0q0)/2.0q0)*Mach**2.0q0))**((gama+1.0q0)/(gama-1.0q0)))

end subroutine AARatioCalc

subroutine NewACalc(gama, Mach, A, At)
    real*16, intent(in)  ::gama
    real*16, intent(in)  ::Mach
    real*16, intent(in)  ::A
    real*16, intent(out) ::At
    
    At = A/qsqrt((1.0q0/Mach**2.0q0)*((2.0q0/(gama+1.0q0))*(1.0q0+((gama-1.0q0)/2.0q0)*Mach**2.0q0))**((gama+1.0q0)/(gama-1.0q0)))

end subroutine NewACalc

subroutine ACalc(gama, R, T, a)
    real*16, intent(in)  ::gama
    real*16, intent(in)  ::R
    real*16, intent(in)  ::T
    real*16, intent(out) ::a
    
    a = qsqrt(gama*R*T)

end subroutine ACalc

subroutine UCalc(Mach, a, u)
    real*16, intent(in)  ::Mach
    real*16, intent(in)  ::a
    real*16, intent(out) ::u
    
    u = Mach*a
    
end subroutine UCalc

real*16 FUNCTION ROOT(gama, P01, P02, CHUTE)
    real*16 ::F
    real*16 ::CHUTE
    real*16 ::FL
    real*16 ::M
    real*16 ::P01
    real*16 ::P02
    real*16 ::gama
    real*16 ::D      !denominador
    real*16 ::DL     !derivada do denominador
    real*16 ::N      !numerador
    real*16 ::NL     !derivada do numerador
    real*16 ::M2     !mach2
    real*16 ::M2_N   !mach2
    real*16 ::M2_NL  !mach2
    real*16 ::M2_D   !mach2
    real*16 ::M2_DL  !mach2
    real*16 ::P      !potencia
  M  = CHUTE
  DO
    N = 1.0q0+0.5q0*(gama-1.0q0)*M**2.0q0
    NL = M*(gama-1.0q0)
    M2_N = (1.0q0+0.5q0*(gama-1.0q0)*M**2.0q0)
    M2_NL = NL
    M2_D = (gama*M**2.0q0-0.5q0*(gama-1.0q0)) 
    M2_DL = 2.0q0*gama*M 
    M2 = M2_N/M2_D
    D = 1.0q0+0.5q0*(gama-1.0q0)*M2
    DL = 0.5q0*(gama-1.0q0)*((M2_NL*M2_D - M2_DL*M2_N)/(M2_D**2.0q0))
    P = gama/(gama-1.0q0)
    F = 1.0q0+2.0q0*gama*(M**2.0q0-1.0q0)/(gama+1.0q0) &
        -(P02/P01)*(N/D)**P
    FL = 4.0q0*gama*M/(gama+1.0q0)-(P02/P01)*(gama/(gama-1))* &
        ((N/D)**(P-1.0q0))*((NL*D - DL*N)/D**2.0q0)
    ROOT  = M - F/FL
    IF (ABS((M-ROOT)/ROOT) .LT. 1E-06) RETURN
    M = ROOT
  END DO
END FUNCTION

subroutine ShockFinder(AA, N, pos)
integer, intent(inout) ::pos
real*16, intent(in)  ::AA
!real*16 ::AAp
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
    !real*16 ::teste1,teste2, teste3
    real*16 ::loc
    
    ThroatFinder = -1
    !teste3 = Se(1)
    loc = Se(1)
    do i=1, N-1
    !teste1 = Pi*rg**2.0q0 
    !teste2 = Se(i)
        if ( Se(i) == Pi*rg**2.0q0 ) then
            ThroatFinder = i
            exit
        end if
        if (loc > Se(i)) then
        ThroatFinder = i
     !   teste3 = Se(i)
        loc = Se(i)
        end if
    end do
    
end function

function MachAACalc(razao, Mach)
    real*16 ::razao
    real*16 ::MachN
    real*16 ::Mach
    real*16 ::Machlx
    real*16 ::Machx
    real*16 ::MachAACalc
    integer ::i
    
    do i=1, 50
        Machlx = (-1.0q0/(Mach*Mach))*(((2.0q0/(gama+1.0q0))*(1.0q0+(gama-1.0q0)*Mach*Mach/2.0q0))**((gama+1.0q0)/(2.0q0*gama-2.0q0)))+((2.0q0/(gama+1.0q0))*(1.0q0+(gama-1.0q0)*Mach*Mach/2.0q0))**((gama+1.0q0)/(2.0q0*gama-2.0q0)-1.0q0);
        Machx = -razao + (1.0q0/Mach)*(((2.0q0/(gama+1.0q0))*(1.0q0+(gama-1.0q0)*Mach*Mach/2.0q0))**((gama+1.0q0)/(2.0q0*gama-2.0q0)))
        MachN = Mach - Machx/Machlx
        Mach = MachN
    end do
    
    MachAACalc = Mach        
    
end function

subroutine ShockNumFinder(p,local)
real*16, dimension(:) , intent(in) ::p
integer, intent(out) ::local

real*16 ::diff
integer ::i

diff = 0.0q0
do i=2,n-1
if((p(i+1)-p(i)) > diff) then
local = i
diff = p(i+1)-p(i)
end if
end do

end subroutine ShockNumFinder

end module NormalShock1D

