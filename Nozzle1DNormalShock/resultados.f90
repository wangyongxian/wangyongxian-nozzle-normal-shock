module resultados

use coeficientes
use solvers_1D

implicit none

contains

! -----------------------------------------------

  subroutine solucao_numerica

	integer :: it, ite  ! auxiliar 

    write(10,15) 
	15 format (/,t4,'Itera��o',6x,'Norma L1(n)/L1(0)',/)

	! c�lculo dos coeficientes e termos fontes
	call coeficientes_e_fontes_qml
	
	call norma (N,apu,awu,aeu,bpu,u,R)
	R_o = R
	
	open(8,file='Norma.dat')

	tcpu = timef() ! zera cron�metro
   
    do it = 1, 10!iteracao
	
	   ! c�lculo dos coeficientes e termos fontes
	   call coeficientes_e_fontes_qml
	   
	   ! solu��o do sistema de equa��es
	   call tdma (N,aPu,awu,aeu,bPu,u)	
	   
	   	write(8,16) it, R
	   16 format (i11,5x,1pe20.13)
	   
	   call norma (N,apu,awu,aeu,bpu,u,R)
	   R = R/R_o	   
	   
	   !write(8,16) it, R
     
       ! c�lculos das velocidades na face leste
	   call calculo_velocidades_face
	   
	   ! c�lculo dos coeficientes do m�todo SIMPLEC
	   call coeficientes_simplec
	   
	   do ite = 1, 2
		  
		  ! c�lculo dos coef e fontes da energia
		  call coeficientes_e_fontes_energia
		  
		  ! solu��o do sistema de equa��es
		  call tdma (N,apT,awT,aeT,bpT,T)
		  
		  call calculo_massa_especifica
		  
		  call calculo_massa_especifica_nas_faces
		  
		  ! c�lculo dos coef e fontes da massa
	      call coeficientes_fontes_massa
		  
	      ! solu��o do sistema de equa��es
	      call tdma (N,aPplinha,awplinha,aeplinha,bPplinha,plinha)
	      
	      ! atualizando plinha fict�cios da massa
	      call atualizar_ficticios_massa
		  
	      ! corrigir a pressao e obter p(p)
	      call corrigir_pressao
		  
	      call corrigir_massa_especifica
	      
	      call corrigir_massa_especifica_faces
	      
	      ! corrigir velocidades e obter u(p)
	      call corrigir_velocidades
		  
	      ! corrigir velocidades das faces
	      call corrigir_velocidades_faces
	      
	   end do

	   ! Atualizando campos para novo avan�o
	   u_o  = u
	   ue_o = ue

	end do

	close(8)

	write(10,16) it, R
	
	tcpu = timef()
	
	! escrita dos coeficientes e fontes velocidade 
	call lista_coeficientes_velocidade

	! escrita dos coeficientes e fontes press�o 
	call lista_coeficientes_pressao
   
    ! escrita da vari�vel prim�ria e sua visualiza��o
    call escreve_T

    write(10,1) tcpu
    1 format(/, f14.3, ' = tempo de processamento (segundos)')

  end subroutine solucao_numerica

!-------------------------------------------------

  ! Solu��o Anal�tica T_exato

  subroutine escreve_T

    integer :: j
	real*8  :: Pref        
	
    ! Antes das tabelas um P�s-processamento
    u(1) = Uin
	u(N) = (u(N-1)+u(n))/2.0d0
	p(1) = (p(1)+p(2))/2.0d0
	p(N) = (p(N-1)+p(N))/2.0d0
	
	Pref = p(1)
	do i = 1, N
	   p(i) = p(i) - Pref
	end do
	
	write(10,14) 
	14 format(//,4x,'OUTROS RESULTADOS RELEVANTES',/)

	write(10,1)
    1 format(/,t4,'volume',t13,'x versus velocidades u (nodais)',/)

	! abertura de arquivo para gravar resultados de u (num�rico)
    open(7,file='U.dat')

	do i = 1, N
	  write( 7,2) i, x(i), u(i)
	  write(10,2) i, x(i), u(i)
      2 format(i4,4x,2(1pe21.11))
	end do
	
	close(7)

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

	! abertura de arquivo para gravar resultados de u (num�rico)
    open(9,file='fm.dat')

	do i = 1, N 
	  if (i == N) then
	     xe(N) = xe(N-1)
		 Me(N) = Me(N-1)
	  end if
	  write( 9,7) i, xe(i), Me(i)
	  write(10,7) i, xe(i), Me(i)
      7 format(i4,4x,2(1pe21.11))
	end do

	close(9)

	write(10,8)
    8 format(//,t4,'volume',t13,'x versus press�es (p e plinha)',/)

	! abertura de arquivo para gravar resultados de u (num�rico)
    open(11,file='p.dat')

	do i = 1, N
	  write(11,9) i, x(i), p(i), plinha(i)
	  write(10,9) i, x(i), p(i), plinha(i)
      9 format(i4,4x,3(1pe21.11))
	end do

	close(11)
	
	! adapta arquivo de comandos para fazer gr�fico
    !open(8,file='Norma.txt')
    !  do j = 1, 6
    !     read(8,*)
    !  end do
    !  write(8,17) head
    !  17 format("set title '",a62,/,"replot")
    !close(8)

	! mostra o gr�fico de Ue
    !ver = system('wgnuplot Norma.txt')

    ! adapta arquivo de comandos para fazer gr�fico
    open(7,file='U.gnu')
      do j = 1, 5
         read(7,*)
      end do
      write(7,10) head
      10 format("set title '",a62,/,"replot")
    close(7)

    ! mostra o gr�fico de U
    ver = system('wgnuplot U.gnu')	

	! adapta arquivo de comandos para fazer gr�fico
    open(9,file='fm.gnu')
      do j = 1, 5
         read(9,*)
      end do
      write(9,12) head
      12 format("set title '",a62,/,"replot")
    close(9)

    ! mostra o gr�fico de fm
    ver = system('wgnuplot fm.gnu')

	! adapta arquivo de comandos para fazer gr�fico
    open(11,file='p.gnu')
      do j = 1, 5
         read(11,*)
      end do
      write(11,13) head
      13 format("set title '",a62,/,"replot")
    close(11)

    ! mostra o gr�fico de p
    ver = system('wgnuplot p.gnu')

  end subroutine escreve_T

!-------------------------------------------------

end module resultados
