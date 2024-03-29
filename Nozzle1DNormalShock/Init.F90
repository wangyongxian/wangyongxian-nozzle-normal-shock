module Initialization

use variaveis
use NormalShock1D

implicit none

contains

!-------------------------------------------------

  subroutine conf_file_read

    ver = system('notepad entrada_dados.txt') ! lista dados

    open(7,file='entrada_dados.txt')

	read(7,*) ! +Propriedades
	read(7,*) P_cam
	read(7,*) P_out
	read(7,*) T_cam
    read(7,*) N
    read(7,*) dt
	read(7,*) iteracao
	read(7,*) freq
	read(7,*) Lt
	read(7,*) fDarcy
	read(7,*) R
	read(7,*) gama
	read(7,*) rin
	read(7,*) rg
	read(7,*) Lc
	read(7,*) Ln
	read(7,*) Beta
	read(7,*) tipo
	read(7,*) ! +Razao de refino
    read(7,*) RazaoRef
    read(7,*) Niveis
	read(7,*) ! +Saida
    read(7,*) title
    read(7,*) caso
    read(7,*) ! +caminho para o aplicativo de Richardson
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
    read(7,*) graf_mach ! mach
    read(7,*) graf_e ! empuxo
    read(7,*) graf_cdesc ! coeficiente de descarga
    read(7,*) graf_dom ! coeficiente de descarga
    read(7,*) ! +Resultados (1 - para mostrar 0 - para nao mostrar)
    read(7,*) res_result ! coeficiente de descarga
    read(7,*) res_iter ! coeficiente de descarga
    read(7,*) res_coef ! abrir arquivo de coeficientes
    read(7,*) gerar_analitico ! gerar solucao analitica
    read(7,*) N_Sol ! numero de n�s para a malha da solucao analitica
    
    close(7)

    allocate(files(Niveis))
    !USE IFPORT
    !$MAXPATH
    !richardson_path = trim(adjustl(richardson_path))
    !richardson_exe =  trim(adjustl(richardson_exe))
    !richardson_exe =  richardson_path // richardson_exe
    !richardson_1 = trim(adjustl(richardson_1))
    !richardson_2 = trim(adjustl(richardson_2))
    !richardson_3 = trim(adjustl(richardson_3))
    !richardson_4 = trim(adjustl(richardson_4))

  end subroutine conf_file_read

!-------------------------------------------------

  subroutine init_alloc(N)
    integer, intent(in)  ::N
    real*8 ::dx
    integer ::i
    ! aloca��o de mem�ria
    allocate (xp(N),xe(N),Sp(N),Se(N),M(N),Me(N),Raio(N),u(N),p(N), T(N), ro(N))
    allocate (pl(N),u_o(N),ue(N),ue_o(N),ds(N),de(N), p_o(N), ro_o(N), roe(N))
    allocate (afu(N),atu(N),btu(N),bpru(N))
	allocate (aw(N),ap(N),ae(N),bp(N), bf(N), bc(N), aww(N))
  	allocate (Empuxo(N), Mach(n), Mache(N), Ma(N), Ua(N), Cd(N), Ta(N), T_o(N), Pa(N), roa(N), f(N), AAp(N))
  	
  	
  	aww = 0.0d0
  	AAp = 0.0d0
  	f = 0.0d0
  	T_o=0.0d0
  	Ta = 0.0d0
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
   
    ! C�lculo do dx 2 ficticios
    dx = Lt/(N-2)
	
	! C�lculo do Pi
	Pi = dacos (-1.0d0)

    ! C�lculo xp internos e nos contornos
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
    ! C�lculo do raio
    do i = 1, N
        if(xp(i) >= Lc) then
           Raio(i) = rg + ((rin-rg)/2.0d0)*(1.0d0+cos(2.0d0*PI*(xp(i)-Lc)/Ln))
        else
           Raio(i) = rin
        end if
    end do

	! Calculando a �rea no ponto p
	do i = 1, N
	   Sp(i) = Pi*(Raio(i)**2.0d0)
	end do
	   AAp = (Raio/rg)**2.0d0

	! Calculando a �rea na face leste
	do i = 1, N-1
	   if(xe(i) >= Lc) then
	        Se(i) = Pi*((rg + ((rin-rg)/2.0d0)*(1.0d0+cos(2.0d0*PI*(xe(i)-Lc)/Ln)))**2.0d0)
	   else
            Se(i) = Pi*(rin**2.0d0)
       end if
	end do  
!Paper: "Weighted Essentially Non-oscillatory Schemes on Triangular Meshes"

 !Paper: "Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic Conservation Laws"

 !Paper: "Efficient Implementation of Weighted ENO Schemes"

 !Also GoogleScholar search WENO and look into the Journal of Physics, Journal of Fluid Mechanics, and AIAA websites.

end subroutine init_alloc

!-------------------------------------------------

subroutine dealloc()
    deallocate (xp,xe,Sp,Se,M,Me,Raio,u,p,T,ro)
    deallocate (pl,u_o,ue,ue_o,ds,de, p_o, ro_o, roe)
    deallocate (afu,atu,btu,bpru)
	deallocate (aw,ap,ae,bp, bf, bc, aww)
  	deallocate (Empuxo, Mach, Mache, Ma, Ua, Cd, Ta, T_o, Pa, roa, f, AAp)
end subroutine dealloc


end module Initialization