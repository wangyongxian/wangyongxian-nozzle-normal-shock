module dados


use variaveis

!-------------------------------------------------

implicit none

contains


!-------------------------------------------------

  subroutine le_dados

    ver = system('notepad dados_entrada.ent') ! lista dados

    open(7,file='dados_entrada.ent')

	read(7,*) caso
	read(7,*) P0
	read(7,*) T0
    read(7,*) N
    read(7,*) deltat
	read(7,*) iteracao
	read(7,*) Lt
	read(7,*) fator
	read(7,*) Rgases
	read(7,*) gama
	read(7,*) rin
	read(7,*) rg
	read(7,*) Lc
	read(7,*) Ln
	read(7,*) Beta
    read(7,*) title
    close(7)

  end subroutine le_dados

!-------------------------------------------------

  subroutine mostra_dados

    integer :: comp

    comp = len(trim(adjustl(caso)))

    write(10,1) trim(adjustl(caso)),N,deltat,iteracao,Lt,fator

    1 format(/,2x,'DADOS',//,  &
                 a<comp>,' = caso',/, &	 
				 8x,i8,  ' = número de volumes de controle',/, &
				 1pe16.8,' = número de avanços no tempo',/,    &
				 8x,i8,  ' = número de iterações externas',/,  & 
				 1pe16.8,' = comprimento domínio de cálculo',/, &
				 1pe16.8,' = viscosidade absoluta do fluido',/, &
				 1pe16.8,' = massa específica do fluido',/,     &
				 1pe16.8,' = fator de atrito de Dárcy',/,   &
				 1pe16.8,' = velocidade inicial no duto',/, &  
				 1pe16.8,' = coeficiente angular do diâmetro ',/)
				 

  end subroutine mostra_dados

!-------------------------------------------------

  subroutine inicializacao
    real*8 :: MachN, Machx, Machlx, razao
    ! alocação de memória
    allocate (x(N),xe(N),A(N),Ae(N),M(N),Me(N),Raio(N),u(N),p(N), T(N), rop(N))
    allocate (plinha(N),u_o(N),ue(N),ue_o(N),ds(N),de(N), p_o(N), rop_o(N), roe(N))
    allocate (afu(N),atu(N),btu(N),bpru(N))
	allocate (awu(N),aPu(N),aeu(N),bPu(N))
	allocate (awt(N),aPt(N),aet(N),bPt(N))
  	allocate (awplinha(N),aPplinha(N),aeplinha(N),bPplinha(N), Empuxo(N), Mach(n), Mache(N), Ma(N))
  	Empuxo = 0.0d0
    afu = 0.0d0
    atu = 0.0d0
    btu = 0.0d0
    bpru = 0.0d0
    awu = 0.0d0
    aPu = 0.0d0
    aeu = 0.0d0
    bPu = 0.0d0
    awt = 0.0d0
    aPt = 0.0d0
    aet = 0.0d0
    bPt = 0.0d0
    awplinha = 0.0d0
    aPplinha = 0.0d0
    aeplinha = 0.0d0
    bPplinha = 0.0d0
    p = 0.0d0
    p_o = 0.0d0
    plinha = 0.0d0
    u = 0.0d0 
    u_o = 0.0d0
    ue = 0.0d0
    ue_o = 0.0d0
    roe = 0.0d0
    rop = 0.0d0
    rop_o = 0.0d0
    T = 0.0d0
    A = 0.0d0
    Ae = 0.0d0
    x = 0.0d0
    xe = 0.0d0
    Mach = 0.0d0
    Mache = 0.0d0
    Ma = 0.0d0
   
    ! Cálculo do deltax 2 ficticios
    deltax = Lt/(N-2.0d0)
	
	! Cálculo do Pi
	Pi = dacos (-1.0d0)

    ! Cálculo xp internos e nos contornos
    x(1) = 0.0d0
    do i = 2, N-1
       x(i) = (i-2.0d0)*deltax + (deltax/2.0d0)
    end do
    x(N) = Lt

    ! Calculando xe da face leste
    do i = 1, N-1 
       xe(i) = (i-(1.0d0))*deltax
	end do
    xe(N) = 0.0d0
    ! Cálculo do raio
    do i = 1, N
        if(x(i) >= Lc) then
           Raio(i) = rg + ((rin-rg)/2.0d0)*(1.0d0+cos(2.0d0*PI*(x(i)-Lc)/Ln))
        else
           Raio(i) = rin
        end if
    end do
  


	! Calculando a área no ponto p
	do i = 1, N
	   A(i) = Pi*(Raio(i)**2)
	end do

	! Calculando a área na face leste
	do i = 1, N-1
	   if(xe(i) >= Lc) then
	        Ae(i) = Pi*((rg + ((rin-rg)/2.0d0)*(1.0d0+cos(2.0d0*PI*(xe(i)-Lc)/Ln)))**2)
	   else
            Ae(i) = Pi*(rin**2)
       end if
	end do  

    cp = gama*Rgases/(gama-1.0d0)
    
    do j=1, N
        if (x(j) < (lc+ln/2.0d0)) then
            Mach(j) = 0.1d0
        else
            Mach(j) = 2.0d0
        end if
        razao = A(j)/(Pi*(rg**2))
        do i=1, 50
            Machlx = (-1.0/(Mach(j)*Mach(j)))*(((2.0/(gama+1.0))*(1.0+(gama-1.0)*Mach(j)*Mach(j)/2.0))**((gama+1.0)/(2.0*gama-2.0)))+((2.0/(gama+1.0))*(1.0+(gama-1.0)*Mach(j)*Mach(j)/2.0))**((gama+1.0)/(2.0*gama-2.0)-1.0);
            Machx = -razao + (1.0d0/Mach(j))*(((2.0d0/(gama+1.0d0))*(1.0d0+(gama-1.0d0)*Mach(j)*Mach(j)/2.0d0))**((gama+1.0d0)/(2.0d0*gama-2.0d0)))
            MachN = Mach(j) - Machx/Machlx
            Mach(j) = MachN
        end do
    
        p(j) = P0*(1.0d0+(gama-1.0d0)*(Mach(j)**2)/2.0d0)**(-gama/(gama-1.0d0))
        T(j) = T0*(1.0d0+(gama-1.0d0)*(Mach(j)**2)/2.0d0)**(-1)
        rop(j) = (P0/(T0*rgases))*(1.0d0+(gama-1.0d0)*(Mach(j)**2)/2.0d0)**(-1.0d0/(gama-1.0d0))
        u(j) = Mach(j)*(gama*rgases*T0*(1.0d0+(gama-1.0d0)*(Mach(j)**2)/2.0d0)**(-1))**(0.5d0)
        Ma(j) = rop(j)*A(j)*u(j)
    end do
    
    do j=1, N-1
        if (xe(j) < (lc+ln/2.0d0)) then
            Mache(j) = 0.1d0
        else
            Mache(j) = 2.0d0
        end if
        razao = Ae(j)/(Pi*(rg**2))
        do i=1, 50
            Machlx = (-1.0/(Mache(j)*Mache(j)))*(((2.0/(gama+1.0))*(1.0+(gama-1.0)*Mache(j)*Mache(j)/2.0))**((gama+1.0)/(2.0*gama-2.0)))+((2.0/(gama+1.0))*(1.0+(gama-1.0)*Mache(j)*Mache(j)/2.0))**((gama+1.0)/(2.0*gama-2.0)-1.0);
            Machx = -razao + (1.0d0/Mache(j))*(((2.0d0/(gama+1.0d0))*(1.0d0+(gama-1.0d0)*Mache(j)*Mache(j)/2.0d0))**((gama+1.0d0)/(2.0d0*gama-2.0d0)))
            MachN = Mache(j) - Machx/Machlx
            Mache(j) = MachN
        end do
    
        roe(j) = (P0/(T0*rgases))*(1.0d0+(gama-1.0d0)*(Mache(j)**2)/2.0d0)**(-1.0d0/(gama-1.0d0))
        ue(j) = Mache(j)*(gama*rgases*T0*(1.0d0+(gama-1.0d0)*(Mache(j)**2)/2.0d0)**(-1))**(0.5d0)
    end do

    rop_o = rop
    p_o = p
    u_o = u
    ue_o = ue
    
    
  end subroutine inicializacao

!-------------------------------------------------

end module dados
