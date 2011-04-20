module resultados

use coeficientes
use solvers_1D

implicit none

contains

! -----------------------------------------------

  subroutine solucao_numerica

	integer :: it
    
   ! write(10,15) 
	!15 format (/,t4,'Iteração',6x,'Norma L1(n)/L1(0)',/)

    !call calcula_fluxo_massa
	
	! cálculo dos coeficientes e termos fontes
	!call coeficientes_e_fontes_qml
	
    
	!call norma (N,apu,-awu,-aeu,bpu,u,R)
	!R_o = R
	
	!open(8,file='Norma.dat')

	tcpu = timef() ! zera cronômetro
    
    do it = 1, iteracao
	
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
	  ! cálculo dos coef e fontes da energia
	  call coeficientes_e_fontes_energia
	  
	  ! solução do sistema de equações
	  call tdma (N,apT,-awT,-aeT,bpT,T)
	  
	  call calculo_massa_especifica
	  
	  call calculo_massa_especifica_nas_faces
	  
	  ! cálculo dos coef e fontes da massa
      call coeficientes_fontes_massa
	  
      ! solução do sistema de equações
      call tdma (N,aPplinha,-awplinha,-aeplinha,bPplinha,plinha)
      
      ! atualizando plinha fictícios da massa
      call atualizar_ficticios_massa
	  
      ! corrigir a pressao e obter p(p)
      call corrigir_pressao
	  
      call corrigir_massa_especifica
      
      ! corrigir velocidades e obter u(p)
      call corrigir_velocidades
	  
      ! corrigir velocidades das faces
     call corrigir_velocidades_faces
      
      call calculo_massa_especifica_nas_faces
!-----------------------------------------------------	      
	   ! Atualizando campos para novo avanço
	   u_o  = u
	   ue_o = ue
	   p_o = p
	   T_o = T
	   rop_o = rop

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
 !   1 format(/,t4,'volume',t13,'x versus velocidades u (nodais)',/)
!
!	! abertura de arquivo para gravar resultados de u (numérico)
 !   open(7,file='U.dat')
!
!	do i = 1, N
!	  write( 7,2) i, x(i), u(i)
!	  write(10,2) i, x(i), u(i)
 !     2 format(i4,4x,2(1pe21.11))
!	end do
	
!	close(7)

!	write(10,3)
 !   3 format(//,t4,'volume',t13,'x versus velocidades u (faces)',/)

!	do i = 1, N
!	  if (i == N) then
!	     xe(N) = xe(N-1)
!		 ue(N) = ue(N-1)
!	  end if
!	  write(10,4) i, xe(i), ue(i)
 !     4 format(i4,4x,2(1pe21.11))
!	end do

!	write(10,5)
  !  5 format(//,t4,'volume',t13,'x versus fluxo de massa',/)

	! abertura de arquivo para gravar resultados de u (numérico)
!	do i = 1, N 
!	  write(10,7) i, x(i), M(i)
   !   7 format(i4,4x,2(1pe21.11))
!	end do

	!write(10,8)
  !  8 format(//,t4,'volume',t13,'x versus pressões (p e plinha)',/)

	! abertura de arquivo para gravar resultados de u (numérico)
    !open(11,file='p.dat')

	!!do i = 1, N
	!  write(11,9) i, x(i), p(i), plinha(i)
	!  write(10,9) i, x(i), p(i), plinha(i)
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
	  write(23,*) x(i), Ma(i), M(i)
	end do
	
    close(23)
   ! ver = system('wgnuplot fm.gnu')
    
    open(23, file='u.dat')
    do i = 1, N
	  write(23,48) x(i), Ua(i), U(i), raio(i)*10000.d0
	  48 format(4(1pe27.18))
	end do
	
    close(23)
    
    open(23, file='ro.dat')
    do i = 1, N
	  write(23,48) x(i), ropA(i), rop(i), raio(i)
	end do
	
	open(23, file='p.dat')
    do i = 1, N
	  write(23,48) x(i), p(i), p(i), raio(i)
	end do
	
    close(23)
    open(22,file='dominio.dat')
	do i = 1, N
	  write(22,*) x(i), raio(i)
	end do
	close(22)
	
    open(22,file='T.dat')
	do i = 1, N
	  write(22,48) x(i), Ta(i), T(i), raio(i)*10000
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
