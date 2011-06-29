module resultados

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
    T_in  = T_cam - 0.5d0*(gama-1.0d0)*(u_in**2)/(gama*R)
    
    M_in = u_in / dsqrt ( gama * R * T_in )
    k = 1.0d0 + 0.5d0 * ( gama - 1.0d0 ) * ( M_in ** 2 )
    p_in = p_cam / ( k ** ( gama / ( gama - 1.0d0 ) ) )
    
    p(1)  = 2.0d0*p_in - p(2)
    u(1)  = 2*u(2) - u(3)
    T(1)  = -T(2) + 2*T_in
    ro(1) = p(1) / ( R * T(1) )
    
    ! inicialização na saída da tubeira
    p(n)  = 2.0d0*p(n-1) - p(n-2)
    u(n)  = 2.0d0*u(n-1) - u(n-2)
    T(n)  = 2.0d0*T(n-1) - T(n-2)
    ro(n) = p(n) / ( R * T(n) )
    de = 0.0d0
    ds = 0.0d0
    
    do it = 1, iteracao
	
	   ! atualização da pressão na entrada da tubeira
       p_ia = p_in
       M_in = u_in / dsqrt ( gama * R * T_in )
       k = 1.0d0 + 0.5d0 * ( gama - 1.0d0 ) * ( M_in ** 2 )
       p_in = p_cam / ( k ** ( gama / ( gama - 1.0d0 ) ) )
       p(1) = 2.0d0*p_in - p(2)
       pl_in = p_in - p_ia
	
       ! Atualizando campos para novo avanço
	   u_o  = u
	   ue_o = ue
	   p_o = p
	   T_o = T
	   ro_o = ro
       
	   ! cálculo dos coeficientes e termos fontes
	   call coeficientes_e_fontes_qml
	   
	   ! solução do sistema de equações
	   call tdma (N,ap,aw,ae,bp,u)	
	   
       ! cálculo dos coeficientes do método SIMPLEC
	   call coeficientes_simplec
	   
       ! cálculos das velocidades na face leste
	   call calculo_velocidades_face
!-----------------------------------------------------		  
      T_in = T_cam - 0.5d0*(gama-1.0d0)*(u_in**2)/(gama*R)
    
	  ! cálculo dos coef e fontes da energia
	  call coeficientes_e_fontes_energia
	  
	  ! solução do sistema de equações
	  call tdma (N,ap,aw,ae,bp,T)
	  
	  ! cálculo da massa específica
      ro = p / ( R * T )
       
	  call calculo_massa_especifica_nas_faces
	  
	  ! cálculo dos coef e fontes da massa
      call coeficientes_fontes_massa
	  
      ! solução do sistema de equações
      call tdma (N,ap,aw,ae,bp,pl)
      
      call correcoes_com_plinha
      call calculo_massa_especifica_nas_faces
      
      T_ex   = 0.5d0 * ( T(n-1) + T(n) )
!-----------------------------------------------------	      

	end do
    
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
subroutine escreve_dados
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
    integer :: comp, comp1

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
	  write(10,9) i, Mach(i), M(i), (Mach(i) - M(i))
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
    ver = system('resultados.txt')
    
end subroutine escreve_dados

  subroutine escreve
    integer ::i
    
    !call calcula_empuxo
    !call calcula_coeficiente_descarga
    !call calcula_fluxo_massa

    
    
    open(23, file='fm.dat')
    do i = 1, N
	  write(23,*) xp(i), Ma(i), M(i)
	end do
    close(23)
   ! ver = system('wgnuplot fm.gnu')
    
    open(23, file='u.dat')
    do i = 1, N
	  write(23,48) xp(i), Ua(i), U(i), raio(i)*10000.d0
	  48 format(4(1pe27.18))
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
	  write(22,48) xp(i), Ta(i), T(i), raio(i)*10000
	end do
	close(22)
    
    if (graf_p) then
        ver = system('wgnuplot p.gnu')
    end if
    if (graf_cdesc) then
    
    end if
    if (graf_e) then
    
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
    
  end subroutine escreve

subroutine solucao_analitica()


!campo de temperatura
!campo de velocidade
!campo de pressao
!campo de ro
!posicao do choque

 !   call AARatioCalc(gama, Mach, AA)
!    call Mach2Calc(gama,mach,mach2)
    
    
end subroutine solucao_analitica


!-------------------------------------------------

end module resultados
