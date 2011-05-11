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
    
   ! write(10,15) 
	!15 format (/,t4,'Iteração',6x,'Norma L1(n)/L1(0)',/)

    !call calcula_fluxo_massa
	
	! cálculo dos coeficientes e termos fontes
	!call coeficientes_e_fontes_qml
	
    
	!call norma (N,apu,-awu,-aeu,bpu,u,R)
	!R_o = R
	
	!open(8,file='Norma.dat')

	tcpu = timef() ! zera cronômetro

    p     = p_cam/((1.0d0+(gama-1.0d0)*(M**2)/2.0d0)**(gama/(gama-1.0d0)))
    ro    = p / ( Rgases * T )
    ue    = u(1:n-1)
    pl    = 0.0d0
    roe   = ro(1:n-1)
    pl_in = 0.0d0
    !bc    = 0.0d0
    
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
    
    
    do it = 1, iteracao
	
	   ! atualização da pressão na entrada da tubeira
       p_ia = p_in
       M_in = u_in / dsqrt ( gama * Rgases * T_in )
       k = 1.0d0 + 0.5d0 * ( gama - 1.0d0 ) * ( M_in ** 2 )
       p_in = p_cam / ( k ** ( gama / ( gama - 1.0d0 ) ) )
       p(1) = 2.0d0*p_in - p(2)
       pl_in = p_in - p_ia
	
	
	   ! cálculo dos coeficientes e termos fontes
	   call coeficientes_e_fontes_qml
	   
	   ! solução do sistema de equações
	   call tdma (N,aPu,-awu,-aeu,bPu,u)	
	   
	!   	write(8,16) it, R
	 !  16 format (i11,5x,1pe20.13)
	   
	!   call norma (N,apu,awu,aeu,bpu,u,R)
	!   R = R/R_o	   
	   
	   !write(8,16) it, R
     
       ! cálculo dos coeficientes do método SIMPLEC
	   call coeficientes_simplec
	   
       ! cálculos das velocidades na face leste
	   call calculo_velocidades_face
!-----------------------------------------------------		  
      T_in = T_cam - 0.5d0*(gama-1.0d0)*(u_in**2)/(gama*Rgases)
    
	  ! cálculo dos coef e fontes da energia
	  call coeficientes_e_fontes_energia
	  
	  ! solução do sistema de equações
	  call tdma (N,apT,-awT,-aeT,bpT,T)
	  
	  ! cálculo da massa específica
      ro = p / ( Rg * T )
       
	  !call calculo_massa_especifica
	  
	  call calculo_massa_especifica_nas_faces
	  
	  ! cálculo dos coef e fontes da massa
      call coeficientes_fontes_massa
	  
      ! solução do sistema de equações
      call tdma (N,aPplinha,-awplinha,-aeplinha,bPplinha,pl)
      
      ! atualizando pl fictícios da massa
      !call atualizar_ficticios_massa
	  
      ! corrigir a pressao e obter p(p)
      !call corrigir_pressao
	  
      !call corrigir_massa_especifica
      
      ! corrigir velocidades e obter u(p)
      !call corrigir_velocidades
	  
      ! corrigir velocidades das faces
     !call corrigir_velocidades_faces
      call correcoes_com_plinha
      call calculo_massa_especifica_nas_faces
      !call calculo_massa_especifica_nas_faces
!-----------------------------------------------------	      
	   ! Atualizando campos para novo avanço
	   u_o  = u
	   ue_o = ue
	   p_o = p
	   T_o = T
	   ro_o = ro

	end do

	!close(8)

	!write(10,16) it, R
	
	tcpu = timef()
	
	! escrita dos coeficientes e fontes velocidade 
	!call lista_coeficientes_velocidade

	! escrita dos coeficientes e fontes pressão 
	!call lista_coeficientes_pressao
   
    ! escrita da variável primária e sua visualização
    call escreve

   ! write(10,1) tcpu
    !1 format(/, f14.3, ' = tempo de processamento (segundos)')

  end subroutine solucao_numerica

!-------------------------------------------------


  subroutine escreve

    !integer :: j
	!real*8  :: Pref        
	
    ! Antes das tabelas um Pós-processamento
    !arrumar uin
    !u(1) = 0
	!u(N) = (u(N-1)+u(n))/2.0d0
	!p(1) = (p(1)+p(2))/2.0d0
	!p(N) = (p(N-1)+p(N))/2.0d0
	
	!Pref = p(1)
	!do i = 1, N
	!   p(i) = p(i) - Pref
	!end do
	
   !! call calcula_empuxo
   ! !call calcula_coeficiente_descarga
  !  call calcula_fluxo_massa
 !   
!	write(10,14) 
!	14 format(//,4x,'OUTROS RESULTADOS RELEVANTES',/)
!
!	write(10,1)
 !   1 format(/,t4,'volume',t13,'xp versus velocidades u (nodais)',/)
!
!	! abertura de arquivo para gravar resultados de u (numérico)
 !   open(7,file='U.dat')
!
!	do i = 1, N
!	  write( 7,2) i, xp(i), u(i)
!	  write(10,2) i, xp(i), u(i)
 !     2 format(i4,4x,2(1pe21.11))
!	end do
	
!	close(7)

!	write(10,3)
 !   3 format(//,t4,'volume',t13,'xp versus velocidades u (faces)',/)

!	do i = 1, N
!	  if (i == N) then
!	     xe(N) = xe(N-1)
!		 ue(N) = ue(N-1)
!	  end if
!	  write(10,4) i, xe(i), ue(i)
 !     4 format(i4,4x,2(1pe21.11))
!	end do

!	write(10,5)
  !  5 format(//,t4,'volume',t13,'xp versus fluxo de massa',/)

	! abertura de arquivo para gravar resultados de u (numérico)
!	do i = 1, N 
!	  write(10,7) i, xp(i), M(i)
   !   7 format(i4,4x,2(1pe21.11))
!	end do

	!write(10,8)
  !  8 format(//,t4,'volume',t13,'xp versus pressões (p e pl)',/)

	! abertura de arquivo para gravar resultados de u (numérico)
    !open(11,file='p.dat')

	!!do i = 1, N
	!  write(11,9) i, xp(i), p(i), pl(i)
	!  write(10,9) i, xp(i), p(i), pl(i)
    !  9 format(i4,4x,3(1pe21.11))
	!end do

	!close(11)
	
    ! adapta arquivo de comandos para fazer gráfico
    !open(7,file='U.gnu')
    !  do j = 1, 5
    !     read(7,*)
    !  end do
    !  write(7,10) head
    !  10 format("set title '",a62,/,"replot")
    !close(7)

    ! mostra o gráfico de U
    !ver = system('wgnuplot U.gnu')	

	! adapta arquivo de comandos para fazer gráfico
    !open(11,file='p.gnu')
    !  do j = 1, 5
    !     read(11,*)
    !  end do
    !  write(11,13) head
    !  13 format("set title '",a62,/,"replot")
    !close(11)

    ! mostra o gráfico de p
    !ver = system('wgnuplot p.gnu')

    ! mostra o dominio de calculo


    !ver = system('wgnuplot dominio.gnu')
    
    
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
