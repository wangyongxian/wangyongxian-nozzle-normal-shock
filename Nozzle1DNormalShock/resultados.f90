module resultados

use coeficientes
use solvers_1D

implicit none

contains

! -----------------------------------------------

  subroutine solucao_numerica

	integer :: it
    real*8 :: k    ! auxiliar
    real*8 :: M_in ! número de Mach na entrada
    
	tcpu = timef() ! zera cronômetro
ue=0.0d0
    p     = p_cam/((1.0d0+(gama-1.0d0)*(Mach**2)/2.0d0)**(gama/(gama-1.0d0)))
    ro    = p / ( Rgases * T )
    ue    = u(1:n-1)
    pl    = 0.0d0
    roe   = ro(1:n-1)
    pl_in = 0.0d0
    bc    = 0.0d0
    
    ! inicialização na entrada da tubeira
    u_in  = ue(1)
    T_in  = T_cam - 0.5d0*(gama-1.0d0)*(u_in**2)/(gama*Rgases)
    
    M_in = u_in / dsqrt ( gama * Rgases * T_in )
    k = 1.0d0 + 0.5d0 * ( gama - 1.0d0 ) * ( M_in ** 2 )
    p_in = p_cam / ( k ** ( gama / ( gama - 1.0d0 ) ) )
    
    p(1)  = 2.0d0*p_in - p(2)
    u(1)  = 2*u(2) - u(3)
    T(1)  = -T(2) + 2*T_in
    ro(1) = p(1) / ( Rgases * T(1) )
    
    ! inicialização na saída da tubeira
    p(n)  = 2.0d0*p(n-1) - p(n-2)
    u(n)  = 2.0d0*u(n-1) - u(n-2)
    T(n)  = 2.0d0*T(n-1) - T(n-2)
    ro(n) = p(n) / ( Rgases * T(n) )
    de = 0.0d0
    ds = 0.0d0
    do it = 1, iteracao
	
	   ! atualização da pressão na entrada da tubeira
       p_ia = p_in
       M_in = u_in / dsqrt ( gama * Rgases * T_in )
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
      T_in = T_cam - 0.5d0*(gama-1.0d0)*(u_in**2)/(gama*Rgases)
    
	  ! cálculo dos coef e fontes da energia
	  call coeficientes_e_fontes_energia
	  
	  ! solução do sistema de equações
	  call tdma (N,ap,aw,ae,bp,T)
	  
	  ! cálculo da massa específica
      ro = p / ( Rgases * T )
       
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
    
    ro_in = p_in / ( Rgases * T_in )
    ro(1) = ro_in
    ro_ex = p_ex / ( Rgases * T_ex )
    ro(n) = ro_ex


	tcpu = timef()
	
    ! escrita da variável primária e sua visualização
    call escreve
    call escreve_dados

  end subroutine solucao_numerica

!-------------------------------------------------
subroutine escreve_dados
    integer ::i
    
    open(10, file="resultados.txt")

	write(10,1)
    1 format(/,t4,'volume',t13,'x versus velocidades u (nodais)',/)

	do i = 1, N
	  write(10,2) i, xp(i), u(i)
      2 format(i4,4x,2(1pe21.11))
	end do
	
	write(10,3)
    3 format(//,t4,'volume',t13,'x versus velocidades u (faces)',/)

	do i = 1, N
	  if (i == N) then
	     xe(N) = xe(N-1)
		 ue(N) = ue(N-1)
	  end if
	  write(10,4) i, xe(i), ue(i)
      4 format(i4,4x,2(1pe21.11))
	end do

	write(10,5)
    5 format(//,t4,'volume',t13,'x versus fluxo de massa nas faces',/)

	do i = 1, N 
	  if (i == N) then
	     xe(N) = xe(N-1)
		 Me(N) = Me(N-1)
	  end if
	  write(10,7) i, xe(i), Me(i)
      7 format(i4,4x,2(1pe21.11))
	end do

	close(9)

	write(10,8)
    8 format(//,t4,'volume',t13,'x versus pressões (p e plinha)',/)


	do i = 1, N
	  write(10,9) i, xp(i), p(i), pl(i)
      9 format(i4,4x,3(1pe21.11))
	end do
    !velocudade nodal
	write(10,13)
    13 format(//,t4,'volume',t13,'xp',t13,'u',t13,'ua',t13,'erro',/)

	do i = 1, N
	  write(10,9) i, u(i), ua(i), (ua(i) - u(i))
	end do
	
	!Temperarura
	write(10,14)
    14 format(//,t4,'volume',t13,'xp',t13,'u',t13,'ua',t13,'temperatura',/)
	do i = 1, N
	  write(10,9) i, T(i), Ta(i), (Ta(i) - T(i))
	end do
	
	!pressao
	write(10,15)
    15 format(//,t4,'volume',t13,'xp',t13,'u',t13,'ua',t13,'pressao',/)
	do i = 1, N
	  write(10,9) i, P(i), Pa(i), (Pa(i) - p(i))
	end do
	
	!Mach
	!write(10,16)
    !16 format(//,t4,'volume',t13,'xp',t13,'u',t13,'ua',t13,'mach',/)
	!do i = 1, N
	 ! write(10,9) i, M(i), Ma(i), (Ma(i) - M(i))
!	end do
	
	!massa especifica
	write(10,17)
    17 format(//,t4,'volume',t13,'xp',t13,'u',t13,'ua',t13,'ro',/)
	do i = 1, N
	  write(10,9) i, ro(i), roa(i), (roa(i) - ro(i))
	end do

	

	close(10)
    ver = system('resultados.txt')
    
end subroutine escreve_dados

  subroutine escreve
    integer ::i
    
    call calcula_empuxo
    call calcula_coeficiente_descarga
    call calcula_fluxo_massa

    
    
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
	  write(23,48) xp(i), ropA(i), ro(i), raio(i)
	end do
	
	open(23, file='p.dat')
    do i = 1, N
	  write(23,48) xp(i), p(i), p(i), raio(i)
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
    
  end subroutine escreve

!-------------------------------------------------

end module resultados
