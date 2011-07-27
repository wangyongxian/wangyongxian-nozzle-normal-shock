module Solucao_Analitica
use variaveis
use Initialization
use NormalShock1D

implicit none

contains

subroutine GenerateSolutionFile(numeroNos)
    integer, intent(in) ::numeroNos
    integer ::i
    call init_alloc(numeroNos)
    call solucao_analitica_init(numeroNos)
    12 format(11(1pe21.13))
    open(12,file='solucao_analitica.txt')
    do i=1, numeroNos
        write(12,12) xp(i), xe(i), u(i), ue(i), t(i), p(i), ro(i), Mach(i), m(i), cd(i), empuxo(i)
    end do
    13 format(t3, 'Xp' ,t24 ,'Xe' ,t45 ,'Velocidade nó' ,t66 ,'Velocidade face' ,t87 ,'Temperatura' ,t108 ,&
            'Pressão',t129 ,'Ro' ,t150 ,'Mach', t171, 'Fluxo de Massa', t192, 'Coeficiente de Descarga', t213, 'Empuxo')
    write(12,13)
    close(12)
    call dealloc()
end subroutine

subroutine solucao_analitica_init(N)
    integer, intent(in) ::N
    
    real*16 :: razao, AA, P02, At
    integer ::i, j, pos 
    
    cp = gama*R/(gama-1.0q0)
    if (P_OUT /= -1.0q0) then
        !localiza o choque
        call ShockLocationCalc(gama, P_out, P_cam, P02, Se(N-1),(rin/rg)**2.0q0, AA, At)
        
        call ShockFinder(AA, N, pos)
    else
        !para a analitica ignorar o choque
        pos = N + 1
    end if
    
    do j=1, N
    
        if (xp(j) < (lc+ln/2.0q0) ) then
            Mach(j) = 0.1q0
            razao = Sp(j)/(Pi*(rg**2.0q0))
        elseif ( j < pos) then
            razao = Sp(j)/(Pi*(rg**2.0q0))
            Mach(j) = 2.0q0
        else
            razao = Sp(j)/At
            Mach(j) = 0.1q0
        end if
        
        Mach(j) = MachAACalc(razao,Mach(j))
        
    end do
    
    do i=1, N
    
        if (i < pos) then
            p(i) = P_cam/(1.0q0+(gama-1.0q0)*(Mach(i)**2.0q0)/2.0q0)**(gama/(gama-1.0q0))
            T(i) = T_cam/(1.0q0+(gama-1.0q0)*(Mach(i)**2.0q0)*0.5q0)
            ro(i) = (P_cam/(T_cam*R))*(1.0q0+(gama-1.0q0)*(Mach(i)**2.0q0)/2.0q0)**(-1.0q0/(gama-1.0q0))
            u(i) = Mach(i)*(gama*R*T_cam*(1.0q0+(gama-1.0q0)*(Mach(i)**2.0q0)/2.0q0)**(-1.0q0))**(0.5q0)
        else
            p(i) = P02/(1.0q0+(gama-1.0q0)*(Mach(i)**2.0q0)/2.0q0)**(gama/(gama-1.0q0))
            T(i) = T_cam/(1.0q0+(gama-1.0q0)*(Mach(i)**2.0q0)*0.5q0)
            ro(i) = (P02/(T_cam*R))*(1.0q0+(gama-1.0q0)*(Mach(i)**2.0q0)/2.0q0)**(-1.0q0/(gama-1.0q0))
            u(i) = Mach(i)*(gama*R*T_cam*(1.0q0+(gama-1.0q0)*(Mach(i)**2.0q0)/2.0q0)**(-1.0q0))**(0.5q0)
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
    ue = u
    ue_o = u
    T_o = T
    
    !calculo do fluxo de massa
    m = roe*ue*Se
    !calculo do coeficiente de descarga
    Cd = m/(qsqrt(gama*R*T))
    !calculo do empuxo
    Empuxo = ro*u*Sp*u
    
    
end subroutine solucao_analitica_init


end module Solucao_Analitica