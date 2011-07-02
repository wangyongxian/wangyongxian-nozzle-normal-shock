module Solucao_Analitica
use variaveis
use NormalShock1D

implicit none

contains

subroutine solucao_analitica_init

    real*8 :: MachN, Machx, Machlx, razao, AA, P02, At
    integer ::i, j, pos 
    
    cp = gama*R/(gama-1.0d0)
    
    !localiza o choque
    call ShockLocationCalc(gama, P_out, P_cam, P02, Se(N-1),(rin/rg)**2.0d0, AA, At)
    
    call ShockFinder(AA, N, pos)
    !call ShockLocationCalc(1.4d0,0.5d0, 1.0d0, 0.0d0, 3.0d0, AA)
    
    do j=1, N
        if (xp(j) < (lc+ln/2.0d0) ) then
            Mach(j) = 0.1d0
            razao = Sp(j)/(Pi*(rg**2.0d0))
        elseif ( j < pos) then
            razao = Sp(j)/(Pi*(rg**2.0d0))
            Mach(j) = 2.0d0
        else
            razao = Sp(j)/At
            Mach(j) = 0.1d0
        end if
        
        Mach(j) = MachAACalc(razao,Mach(j))
        
        !do i=1, 50
        !    Machlx = (-1.0d0/(Mach(j)*Mach(j)))*(((2.0d0/(gama+1.0d0))*(1.0d0+(gama-1.0d0)*Mach(j)*Mach(j)/2.0d0))**((gama+1.0d0)/(2.0d0*gama-2.0d0)))+((2.0d0/(gama+1.0d0))*(1.0d0+(gama-1.0d0)*Mach(j)*Mach(j)/2.0d0))**((gama+1.0d0)/(2.0d0*gama-2.0d0)-1.0d0);
        !    Machx = -razao + (1.0d0/Mach(j))*(((2.0d0/(gama+1.0d0))*(1.0d0+(gama-1.0d0)*Mach(j)*Mach(j)/2.0d0))**((gama+1.0d0)/(2.0d0*gama-2.0d0)))
        !    MachN = Mach(j) - Machx/Machlx
        !    Mach(j) = MachN
        !end do
    end do
    
    do i=1, N
    if (i < pos) then
        p(i) = P_cam/(1.0d0+(gama-1.0d0)*(Mach(i)**2.0d0)/2.0d0)**(gama/(gama-1.0d0))
        T(i) = T_cam/(1.0d0+(gama-1.0d0)*(Mach(i)**2.0d0)*0.5d0)
        ro(i) = (P_cam/(T_cam*R))*(1.0d0+(gama-1.0d0)*(Mach(i)**2.0d0)/2.0d0)**(-1.0d0/(gama-1.0d0))
        u(i) = Mach(i)*(gama*R*T_cam*(1.0d0+(gama-1.0d0)*(Mach(i)**2.0d0)/2.0d0)**(-1.0d0))**(0.5d0)
    else
        p(i) = P02/(1.0d0+(gama-1.0d0)*(Mach(i)**2.0d0)/2.0d0)**(gama/(gama-1.0d0))
        T(i) = T_cam/(1.0d0+(gama-1.0d0)*(Mach(i)**2.0d0)*0.5d0)
        ro(i) = (P02/(T_cam*R))*(1.0d0+(gama-1.0d0)*(Mach(i)**2.0d0)/2.0d0)**(-1.0d0/(gama-1.0d0))
        u(i) = Mach(i)*(gama*R*T_cam*(1.0d0+(gama-1.0d0)*(Mach(i)**2.0d0)/2.0d0)**(-1.0d0))**(0.5d0)
    end if
    end do
    Ma = ro*Sp*u
    
    Ua = u
    Ta = T    
    Pa = p
    roa = ro
    ro_o = ro
    roe = ro
    p_o = p
    u_o = u
    ue_o = u
    T_o = T
    
    
end subroutine solucao_analitica_init


end module Solucao_Analitica