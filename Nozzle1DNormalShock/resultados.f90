module resultados

use variaveis
use coeficientes
use solvers_1D
use NormalShock1D


implicit none

contains

! -----------------------------------------------

  subroutine solucao_numerica

	integer :: it
    real*8 :: k    ! auxiliar
    real*8 :: M_in ! número de Mach na entrada
    
	tcpu = timef() ! zera cronômetro
	
    ue = 0.0d0
    !p     = p_cam/((1.0d0+(gama-1.0d0)*(Mach**2)/2.0d0)**(gama/(gama-1.0d0)))
    ro    = p / ( R * T )
    ue    = u(1:n-1)
    pl    = 0.0d0
    roe   = ro(1:n-1)
    pl_in = 0.0d0
    bc    = 0.0d0
    
    ! inicialização na entrada da tubeira
    u_in  = ue(1)
    u_out  = ue(N-1)
    T_in  = T_cam - 0.5d0*(gama-1.0d0)*(u_in**2)/(gama*R)
    T_out  = T_cam - 0.5d0*(gama-1.0d0)*(u_out**2)/(gama*R)
    
    M_in = u_in / dsqrt ( gama * R * T_in )
    M_out = u_out / dsqrt ( gama * R * T_out )
    
    k = 1.0d0 + 0.5d0 * ( gama - 1.0d0 ) * ( M_in ** 2 )
    p_in = p_cam / ( k ** ( gama / ( gama - 1.0d0 ) ) )
    k = 1.0d0 + 0.5d0 * ( gama - 1.0d0 ) * ( M_out ** 2 )
    p_outN = p_out
    
    p(1)  = 2.0d0*p_in - p(2)
    u(1)  = 2*u(2) - u(3)
    T(1)  = -T(2) + 2*T_in
    ro(1) = p(1) / ( R * T(1) )
    
    ! inicialização na saída da tubeira
    p(n)  = 2.0d0*p_out - p(n-1)  !2.0d0*p(n-1) - p(n-2)
    u(n)  = 2.0d0*u(n-1) - u(n-2)
   ! T(n)  = 2.0d0*T(n-1) - T(n-2)
    T(N) = -T(N-1) + 2*T_out
    ro(n) = p(n) / ( R * T(n) )
    de = 0.0d0
    ds = 0.0d0
    
    !residuo
    call coeficientes_e_fontes_qml
    call Norma_L1( n, aw, ap, ae, bp, u, Residuo_U )
	Residuo_U_o = Residuo_U
	
	call coeficientes_fontes_massa
	call Norma_L1( n, aw, ap, ae, bp, pl, Residuo_P )
	Residuo_P_o = Residuo_P
	
	call coeficientes_e_fontes_energia
	call Norma_L1( n, aw, ap, ae, bp, T, Residuo_T )
	Residuo_T_o = Residuo_T
	
	
	open(8, file='norma_l1.txt')
	write(8,4)
	4 format(t1,'iteracao',t18 ,'Residuo T',t38 ,'Residuo U',t58,'Residuo Pl')
	
    do it = 1, iteracao
    
        !if ((it > 5000).and.(it < 50000)) beta = 0.9d0
        !if (it >= 50000 ) beta = 1.0d0
	
	   ! Atualizando campos para novo avanço
	   u_o  = u
	   ue_o = ue
	   p_o = p
	   T_o = T
	   ro_o = ro
	   
	   ! atualização da pressão na entrada da tubeira
       p_ia = p_in
       M_in = u_in / dsqrt ( gama * R * T_in )
       k = 1.0d0 + 0.5d0 * ( gama - 1.0d0 ) * ( M_in ** 2 )
       p_in = p_cam / ( k ** ( gama / ( gama - 1.0d0 ) ) )
       p(1) = 2.0d0*p_in - p(2)
       pl_in = p_in - p_ia
       
       p_outa = p_outN
       M_out = u_out / dsqrt ( gama * R * T_out )
       k = 1.0d0 + 0.5d0 * ( gama - 1.0d0 ) * ( M_out ** 2 )
       p_outN = p_out / ( k ** ( gama / ( gama - 1.0d0 ) ) )
       p(N) = 2.0d0*p_outN - p(N-1)
       pl_out = p_outN - p_outa
	

	   ! cálculo dos coeficientes e termos fontes
	   call coeficientes_e_fontes_qml
	   
	   ! solução do sistema de equações
	   call tdma (N,ap,aw,ae,bp,u)	
	  
	   call Norma_L1 (N,ap,aw,ae,bp,u,Residuo_U)
	   Residuo_U = Residuo_U/Residuo_U_o
	    
       ! cálculo dos coeficientes do método SIMPLEC
	   call coeficientes_simplec
	   
       ! cálculos das velocidades na face leste
	   call calculo_velocidades_face
!-----------------------------------------------------		  
     ! inicialização na entrada da tubeira
      u_in  = ue(1)
      u_out  = ue(N-1)

      T_in = T_cam - 0.5d0*(gama-1.0d0)*(u_in**2)/(gama*R)
      T_out = T_cam - 0.5d0*(gama-1.0d0)*(u_out**2)/(gama*R)
    
	  ! cálculo dos coef e fontes da energia
	  call coeficientes_e_fontes_energia
	  
	  ! solução do sistema de equações
	  call tdma (N,ap,aw,ae,bp,T)
	  
	  call Norma_L1 (N,ap,aw,ae,bp,T,Residuo_T)
	  Residuo_T = Residuo_T/Residuo_T_o
	   
	  ! cálculo da massa específica
      ro = p / ( R * T )
       
	  call calculo_massa_especifica_nas_faces
	  
	  ! cálculo dos coef e fontes da massa
      call coeficientes_fontes_massa
	  
      ! solução do sistema de equações
      call tdma (N,ap,aw,ae,bp,pl)
      
      call Norma_L1 (N,ap,aw,ae,bp,pl,Residuo_P)
	  Residuo_P = Residuo_P/Residuo_P_o
	   
      call correcoes_com_plinha
      call calculo_massa_especifica_nas_faces
      
      !T_ex   = 0.5d0 * ( T(n-1) + T(n) )
      
      write(8,16) it, Residuo_T, Residuo_U, Residuo_P
	  16 format (i11,5x,3(1pe20.13))
	   
!-----------------------------------------------------	      

	end do
    close(8)
    
    ! coeficiente da pressao
    call coeficientes_fontes_massa
    call gera_arq_coef(2)
    !coeficiente da Energia
    !call coeficientes_e_fontes_energia
    !call gera_arq_coef(3)
    !coeficiente da QML
    !call coeficientes_e_fontes_qml
    !call gera_arq_coef(4)
    
    ! entrada da tubeira
    u(1)  = u_in
    p(1)  = p_in
    T(1)  = T_in

    ! saída da tubeira
    T_ex  = 0.5d0 * ( T(n-1) + T(n) )
    p_ex  = 0.5d0 * ( p(n-1) + p(n) )
    u(n)  = ue(n-1)
    p(n)  = p_ex
    T(n)  = T_ex
    
    ro_in = p_in / ( R * T_in )
    ro(1) = ro_in
    ro_ex = p_ex / ( R * T_ex )
    ro(n) = ro_ex

	tcpu = timef()

  end subroutine solucao_numerica

!-------------------------------------------------
subroutine gera_txt
    integer ::i
    integer*4    :: var(8)    ! data e hora
    character*20 :: vardate   ! data e hora
    character*20 :: vartime   ! data e hora
    character*20 :: varzone   ! data e hora
    character*2  :: aux1,aux2
    character*4  :: aux3
    character*50 :: aux
    character*12 :: dia       ! data da simulação
    character*8  :: hora      ! horário da simulação
    character*62 :: head      ! título do gráfico + dia

    open(10, file="resultados.txt")
        
    call date_and_time(vardate,vartime,varzone,var)
    
    write(aux,*) var(5)
    aux1 = trim(adjustl(aux))
    write(aux,*) var(6)
    aux2 = trim(adjustl(aux))
    write(aux,*) var(7)
    aux3 = trim(adjustl(aux))
    hora = trim(aux1)//':'//trim(aux2)//':'//aux3
    
    head = trim(title)//" "//trim(dia)//"'"
    
    !write(10,2) trim(adjustl(title)), dia, hora
    !2 format(/,'Título = ', a<comp1>, &
    !        //,5x,'Dia = ',a12,5x,'Hora = ',a8)
    !comp = len(trim(adjustl(caso)))

    write(10,1)  P_cam, P_out, T_cam, N, dt, iteracao, Lt, fDarcy, R, Gama, rin, rg, Lc, Ln, Beta

    1 format(/,2x,'DADOS',//,  &
                ! a<comp>,  ' = caso',/, &	 
				 1pe16.8,  ' = pressao de estagnacao',/, &
				 1pe16.8,  ' = pressao da saida',/, &
				 1pe16.8,  ' = temperatura de estagnacao',/, &
				 8x,i8,    ' = número de volumes de controle',/, &
				 1pe16.8,  ' = número de avanços no tempo',/, &
				 8x,i8,    ' = número de iterações (ciclo total)',/, &
				 1pe16.8,  ' = comprimento do domínio de cálculo',/, &
				 1pe16.8,  ' = fator de atrito de Darcy',/, &
				 1pe16.8,  ' = R constante do gas',/,    &
				 1pe16.8,  ' = Gama',/,  & 
				 1pe16.8,  ' = Rin',/,  & 
				 1pe16.8,  ' = Rg',/,  & 
				 1pe16.8,  ' = Lc',/,  & 
				 1pe16.8,  ' = Ln',/,  & 
				 1pe16.8,  ' = Beta',/)

    9 format(i4,4x,4(1pe21.11))
	
	write(10,13)
    13 format(//,t1,'volume',t13,'U(Analitica)',t34,'U(Numerica)',t55,'Erro')

	do i = 1, N
	  write(10,9) i, ua(i), u(i), (ua(i) - u(i))
	end do
	
	!Temperarura
	write(10,14)
    14 format(//,t1,'volume',t13,'T(Analitica)',t34,'T(Numerica)',t55,'Erro')
	do i = 1, N
	    
	  write(10,9) i, Ta(i), T(i), (Ta(i) - T(i))
	end do
	
	!pressao
	write(10,15)
    15 format(//,t1,'volume',t13,'P(Analitica)',t34,'P(Numerica)',t55,'Erro')
	do i = 1, N
	  write(10,9) i, Pa(i), P(i), (Pa(i) - p(i))
	end do
	
	!Mach
	write(10,16)
    16 format(//,t1,'volume',t13,'Mach(Analitico)',t34,'Mach(Numerico)',t55,'Erro')
	do i = 1, N 
	  write(10,9) i, Mach(i), (u(i)/dsqrt(gama*R*T(i))), (Mach(i) - (u(i)/dsqrt(gama*R*T(i))))
	end do
	
	!massa especifica
	write(10,17)
    17 format(//,t1,'volume',t13,'RO(Analitico)',t34,'RO(Numerico)',t55,'Erro')
	do i = 1, N
	  write(10,9) i, roa(i), ro(i), (roa(i) - ro(i))
	end do

    write(10,10) tcpu
    10 format(/, f14.3, ' = tempo de processamento (segundos)')
    
	close(10)
    
end subroutine gera_txt

subroutine gera_graficos
    integer ::i
    real*8 ::T0
    
    call calcula_empuxo
    call calcula_coeficiente_descarga
    call calcula_fluxo_massa
    
    open(23, file='fm.dat')
    do i = 1, N
	  write(23,*) xp(i), Ma(i), M(i)
	end do
    close(23)
    
    open(23, file='u.dat')
    do i = 1, N
	  write(23,48) xp(i), Ua(i), U(i), raio(i)*10000.d0
	  48 format(4(1pe27.18))
	end do
	
    close(23)
    
    open(23, file='Mach.dat')
    do i = 1, N
	  write(23,48) xp(i), Mach(i), (u(i)/dsqrt(gama*R*T(i))), raio(i)
	end do
	
	close(23)
    open(23, file='ro.dat')
    do i = 1, N
	  write(23,48) xp(i), roa(i), ro(i), raio(i)
	end do
	
	open(23, file='p.dat')
    do i = 1, N
	  write(23,48) xp(i), pa(i), p(i), raio(i)
	end do
	
    close(23)
    open(22,file='dominio.dat')
	do i = 1, N
	  write(22,*) xp(i), raio(i)
	end do
	close(22)
	
    open(22,file='T.dat')
	do i = 1, N
	    call T0Calc(gama,T(i), (u(i)/dsqrt(gama*R*T(i))), T0)
	    !call T0Calc(gama,Ta(i), Mach(i), T0)
	  write(22,48) xp(i), Ta(i), T(i), raio(i)
	end do
	close(22)
    
    if (graf_p) then
        ver = system('wgnuplot p.gnu')
    end if
    if (graf_cdesc) then
        ver = system('wgnuplot coef_descarga.gnu')
    end if
    if (graf_e) then
        ver = system('wgnuplot empuxo.gnu')
    end if
    if (graf_dom) then
        ver = system('wgnuplot dominio.gnu')
    end if
    if (graf_ro) then
        ver = system('wgnuplot ro.gnu')
    end if
    if (graf_m) then
        ver = system('wgnuplot fm.gnu')
    end if
    if (graf_v) then
        ver = system('wgnuplot U.gnu')
    end if
    if (graf_t) then
       ver = system('wgnuplot T.gnu')
    end if
    if(graf_mach)then
       ver = system('wgnuplot Mach.gnu')
    end if
    if(res_coef) then
    ver = system('coeficientes.txt')
    end if 
    if (res_result) then
        ver = system('resultados.txt')
    end if
    
    if (res_iter) then
        ver = system('norma_l1.txt')
    end if
    
end subroutine gera_graficos

subroutine gera_arq_coef(prop)
    integer, intent(in) ::prop
    !prop 2 - PL
    !prop 3 - Energia
    !prop 4 - QML
    integer ::i
    open(14,file='coeficientes.txt')
    2 format('Coeficientes de Pl', /,t4, 'volume', t13, 'x', t34,'oeste',t55,'central',t76,'leste',t97,'fonte',/)
    3 format('Coeficientes Energia', /,t4, 'volume', t13, 'x', t34,'oeste',t55,'central',t76,'leste',t97,'fonte',/)
    4 format('Coeficientes QML', /,t4, 'volume', t13, 'x', t34,'oeste',t55,'central',t76,'leste',t97,'fonte',/)
    5 format(i4, 4x, 5(1pe21.11))
    select case(prop)
        case(2)
            write(14, 2)
        case(3)
            write(14, 3)
        case(4)
            write(14, 4)
    end select
    
    do i=1,N
        write(14, 5) i, xp(i), aw(i), ap(i), ae(i), bp(i)
    end do
    
    close(14)
end subroutine gera_arq_coef
subroutine solucao_analitica()

    
end subroutine solucao_analitica


!-------------------------------------------------

end module resultados
