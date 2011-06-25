module dados


use variaveis

!-------------------------------------------------

implicit none

contains


!-------------------------------------------------

  subroutine le_dados

   ! ver = system('notepad dados_entrada.ent') ! lista dados

    open(7,file='entrada_dados.txt')

	read(7,*) ! +Propriedades
	read(7,*) P_cam
	read(7,*) T_cam
    read(7,*) N
    read(7,*) dt
	read(7,*) iteracao
	read(7,*) Lt
	read(7,*) fDarcy
	read(7,*) Rgases
	read(7,*) gama
	read(7,*) rin
	read(7,*) rg
	read(7,*) Lc
	read(7,*) Ln
	read(7,*) Beta
	read(7,*) ! +Saida
    read(7,*) title
    read(7,*) caso
    read(7,*) ! +Saida para o programa Richardson
    read(7,*) richardson_1
    read(7,*) richardson_2
    read(7,*) richardson_3
    read(7,*) richardson_4
    read(7,*) ! 012345   1         2         3         4         5         6         7         8         9
    read(7,*) ! +Graficos 1 - para ligar 0 para desligar
    read(7,*) graf_m ! fluxo de massa
    read(7,*) graf_t ! temperatura
    read(7,*) graf_v ! velocidade
    read(7,*) graf_ro ! massa especifica
    read(7,*) graf_p ! pressao
    read(7,*) graf_e ! empuxo
    read(7,*) graf_cdesc ! coeficiente de descarga
    read(7,*) graf_dom ! coeficiente de descarga
    
    close(7)

  end subroutine le_dados

!-------------------------------------------------

  subroutine inicializacao
    real*8 :: MachN, Machx, Machlx, razao, dx
    integer ::i, j
    ! alocação de memória
    allocate (xp(N),xe(N),Sp(N),Se(N),M(N),Me(N),Raio(N),u(N),p(N), T(N), ro(N))
    allocate (pl(N),u_o(N),ue(N),ue_o(N),ds(N),de(N), p_o(N), ro_o(N), roe(N), ropA(N))
    allocate (afu(N),atu(N),btu(N),bpru(N))
	allocate (aw(N),ap(N),ae(N),bp(N), bf(N), bc(N))
  	allocate (Empuxo(N), Mach(n), Mache(N), Ma(N), Ua(N), Cd(N), Ta(N), T_o(N), Pa(N), roa(N), f(N))
  	f = 0.0d0
  	T_o=0.0d0
  	Ta = 0.0d0
  	ropA = 0.0d0
  	Cd = 0.0d0
  	Ua = 0.0d0
  	Empuxo = 0.0d0
    afu = 0.0d0
    atu = 0.0d0
    btu = 0.0d0
    bpru = 0.0d0
    aw = 0.0d0
    ap = 0.0d0
    ae = 0.0d0
    bp = 0.0d0
    p = 0.0d0
    pa = 0.0d0
    p_o = 0.0d0
    pl = 0.0d0
    u = 0.0d0 
    u_o = 0.0d0
    ue = 0.0d0
    ue_o = 0.0d0
    roe = 0.0d0
    ro = 0.0d0
    ro_o = 0.0d0
    roa = 0.0d0
    T = 0.0d0
    Sp = 0.0d0
    Se = 0.0d0
    xp = 0.0d0
    xe = 0.0d0
    Mach = 0.0d0
    Mache = 0.0d0
    Ma = 0.0d0
    bc = 0.0d0
    bf = 0.0d0
   
    ! Cálculo do dx 2 ficticios
    dx = Lt/(N-2)
	
	! Cálculo do Pi
	Pi = dacos (-1.0d0)

    ! Cálculo xp internos e nos contornos
    xp(1) = 0.0d0
    do i = 2, N-1
       xp(i) = (i-2.0d0)*dx + (dx/2.0d0)
    end do
    xp(N) = Lt

    ! Calculando xe da face leste
    do i = 1, N-1 
       xe(i) = (i-(1.0d0))*dx
	end do
    xe(N) = 0.0d0
    ! Cálculo do raio
    do i = 1, N
        if(xp(i) >= Lc) then
           Raio(i) = rg + ((rin-rg)/2.0d0)*(1.0d0+cos(2.0d0*PI*(xp(i)-Lc)/Ln))
        else
           Raio(i) = rin
        end if
    end do
  


	! Calculando a área no ponto p
	do i = 1, N
	   Sp(i) = Pi*(Raio(i)**2)
	end do

	! Calculando a área na face leste
	do i = 1, N-1
	   if(xe(i) >= Lc) then
	        Se(i) = Pi*((rg + ((rin-rg)/2.0d0)*(1.0d0+cos(2.0d0*PI*(xe(i)-Lc)/Ln)))**2)
	   else
            Se(i) = Pi*(rin**2)
       end if
	end do  

    cp = gama*Rgases/(gama-1.0d0)
    
    do j=1, N
        if (xp(j) < (lc+ln/2.0d0)) then
            Mach(j) = 0.1d0
        else
            Mach(j) = 2.0d0
        end if
        razao = Sp(j)/(Pi*(rg**2))
        do i=1, 50
            Machlx = (-1.0/(Mach(j)*Mach(j)))*(((2.0/(gama+1.0))*(1.0+(gama-1.0)*Mach(j)*Mach(j)/2.0))**((gama+1.0)/(2.0*gama-2.0)))+((2.0/(gama+1.0))*(1.0+(gama-1.0)*Mach(j)*Mach(j)/2.0))**((gama+1.0)/(2.0*gama-2.0)-1.0);
            Machx = -razao + (1.0d0/Mach(j))*(((2.0d0/(gama+1.0d0))*(1.0d0+(gama-1.0d0)*Mach(j)*Mach(j)/2.0d0))**((gama+1.0d0)/(2.0d0*gama-2.0d0)))
            MachN = Mach(j) - Machx/Machlx
            Mach(j) = MachN
        end do
    end do
    
    p = P_cam/(1.0d0+(gama-1.0d0)*(Mach**2)/2.0d0)**(gama/(gama-1.0d0))
    T = T_cam/(1.0d0+(gama-1.0d0)*(Mach**2)*0.5d0)
    ro = (P_cam/(T_cam*rgases))*(1.0d0+(gama-1.0d0)*(Mach**2)/2.0d0)**(-1.0d0/(gama-1.0d0))
    u = Mach*(gama*rgases*T_cam*(1.0d0+(gama-1.0d0)*(Mach**2)/2.0d0)**(-1))**(0.5d0)
    Ma = ro*Sp*u
    !solucao analitica para as faces
    !do j=1, N-1
    !    if (xe(j) < (lc+ln/2.0d0)) then
    !        Mache(j) = 0.1d0
    !    else
    !        Mache(j) = 2.0d0
    !    end if
    !    razao = Se(j)/(Pi*(rg**2))
    !    do i=1, 50
    !        Machlx = (-1.0/(Mache(j)*Mache(j)))*(((2.0/(gama+1.0))*(1.0+(gama-1.0)*Mache(j)*Mache(j)/2.0))**((gama+1.0)/(2.0*gama-2.0)))+((2.0/(gama+1.0))*(1.0+(gama-1.0)*Mache(j)*Mache(j)/2.0))**((gama+1.0)/(2.0*gama-2.0)-1.0);
    !        Machx = -razao + (1.0d0/Mache(j))*(((2.0d0/(gama+1.0d0))*(1.0d0+(gama-1.0d0)*Mache(j)*Mache(j)/2.0d0))**((gama+1.0d0)/(2.0d0*gama-2.0d0)))
    !        MachN = Mache(j) - Machx/Machlx
    !        Mache(j) = MachN
    !    end do
    !    roe(j) = (P_cam/(T_cam*rgases))*(1.0d0+(gama-1.0d0)*(Mache(j)**2)/2.0d0)**(-1.0d0/(gama-1.0d0))
    !    ue(j) = Mache(j)*(gama*rgases*T_cam*(1.0d0+(gama-1.0d0)*(Mache(j)**2)/2.0d0)**(-1))**(0.5d0)
    !end do
    Ua = u
    roa = ro
    ro_o = ro
    roe = ro
    p_o = p
    u_o = u
    ue_o = u
    ropA = ro
    Ta = T    
    Pa = p
    T_o=T
    
  end subroutine inicializacao

!-------------------------------------------------

end module dados
