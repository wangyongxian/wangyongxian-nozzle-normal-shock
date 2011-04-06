module coeficientes
  
! objetivo: calcular os coeficientes e termos fontes
!           das equações discretizadas

use variaveis  
use dados
  
implicit none
  
contains

!-------------------------------------------------

  subroutine lista_coeficientes

    write(10,4)
    4 format(//,4x,'COEFICIENTES DA VELOCIDADE',//,  &
             t4,'volume',t13,'X',t34,'oeste',t55,'central', &
                     t76,'leste',t97,'fonte',/)
    do i = 1, N
       write(10,2) i, X(i), awu(i), aPu(i), aeu(i), bPu(i)
    end do

    2 format(i4,4x,5(1pe21.11))

  end subroutine lista_coeficientes

!-------------------------------------------------

  subroutine lista_coeficientes_pressao

    write(10,10)
    10 format(//,4x,'COEFICIENTES DA PRESSÃO',//,  &
             t4,'volume',t13,'X',t34,'oeste',t55,'central', &
                     t76,'leste',t97,'fonte',/)
    do i = 1, N
       write(10,2) i, X(i), awplinha(i), aPplinha(i), aeplinha(i), bPplinha(i)
    end do

    2 format(i4,4x,5(1pe21.11))

  end subroutine lista_coeficientes_pressao

!-------------------------------------------------


  subroutine coeficientes_e_fontes_qml
	
	! Calculando o fluxo de massa na face leste
	do i = 1, N
	   Me(i) = ro*ue(i)*Ae(i) 
	end do

	! Fictício esquerdo P = 1
    awu(1) = 0.0d0
    aeu(1) = -1.0d0
    aPu(1) = 1.0d0
    bPu(1) = 2.0d0*Uin
	
	! volumes internos
    do i = 2, N-1
       afu(i)  = 0 
	   atu(i)  = 0 
	   btu(i)  = 0 
	   bpru(i) = - A(i)*(p(i+1)-p(i-1))/2.0d0

	   awu(i)  = Me(i-1)
	   aeu(i)  = 0 
       aPu(i)  = Me(i)
	   bpu(i)  = bpru(i)
	end do
      
    ! Fictício direito P = N
    awu(N) = 1.0d0
    aeu(N) = 0.0d0
    aPu(N) = 1.0d0
    bPu(N) = u(N-1)-u(N-2)

  end subroutine coeficientes_e_fontes_qml

!-------------------------------------------------


  subroutine coeficientes_e_fontes_energia
	
	! Calculando o fluxo de massa na face leste
	do i = 1, N
	   Me(i) = ro*ue(i)*Ae(i) 
	end do

	! Fictício esquerdo P = 1
    awu(1) = 0.0d0
    aeu(1) = -1.0d0
    aPu(1) = 1.0d0
    bPu(1) = T0-((gama-1.0d0)/2.0d0)*(Uin*Uin)/(Rgases*gama)
	
	! volumes internos
    do i = 2, N-1
	   awu(i)  = Me(i-1)
	   aeu(i)  = 0 
       aPu(i)  = Me(i)
	   bpu(i)  = cp*A(i)*(p(i+1)-p(i-1))/2.0d0
	end do
      
    ! Fictício direito P = N
    awu(N) = 1.0d0
    aeu(N) = 0.0d0
    aPu(N) = 1.0d0
    bPu(N) = u(N-1)-u(N-2)

  end subroutine coeficientes_e_fontes_energia
  
!-------------------------------------------------
  
  subroutine velocidades_ue

    real*8  :: sigmap, sigmaE ! auxiliar

	! Calculando ue no passo deltat+1
	ue(1) = Uin
	
	do i = 2, N-2
	   sigmap = awu(i)*u(i-1) + aeu(i)*u(i+1)
	   sigmaE = awu(i+1)*u(i) + aeu(i+1)*u(i+2)
	   ue(i) = (sigmap+sigmaE - 2.0d0*Ae(i)*(p(i+1)-p(i)))/(apu(i)+apu(i+1))
	end do
	
	ue(N-1) = (u(N-1)+u(N))/2.0d0
	
	Fat = (Uin*Ae(1))/(ue(N-1)*Ae(N-1))
	
  	ue(N-1) = Fat*ue(N-1)

  end subroutine velocidades_ue

!-------------------------------------------------

  subroutine coeficientes_simplec

    do i = 2, N-1
	   ds(i) = A(i)/(apu(i)-awu(i)-aeu(i))
    end do

	do i = 2, N-2
	   de(i) = (ds(i) + ds(i+1))/2.0d0
	end do

	! Calculando nos contornos
	de(1) = ds(1+1)
	de(N-1) = ds(N-1)

  end subroutine coeficientes_simplec

!-------------------------------------------------

  subroutine coeficientes_fontes_massa

    ! Calculando fictício P = 1	
    awplinha(1) = 0.0d0
    aeplinha(1) = 0.0d0
    aPplinha(1) = 1.0d0
    bpplinha(1) = 0.0d0

    ! Calculando os volumes internos
    do i = 2, N-1
       awplinha(i) = de(i-1)*Ae(i-1)
	   aeplinha(i) = de(i)*Ae(i)
	   aPplinha(i) = awplinha(i) + aeplinha(i)
	   bpplinha(i) = ue(i-1)*Ae(i-1) - ue(i)*Ae(i)
    end do
  
    ! Calculando fictício P = N 
    awplinha(N) = 0.0d0
    aeplinha(N) = 0.0d0
    aPplinha(N) = 1.0d0
    bpplinha(N) = 0.0d0

  end subroutine coeficientes_fontes_massa

!-------------------------------------------------

  subroutine atualizar_ficticios_massa

	! Calculando fictício P = 1
	plinha(1) = 2.0d0*plinha(1+1) - plinha(1+2)

	! Calculando fictício P = N 
	plinha(N) = 2.0d0*plinha(N-1) - plinha(N-2)    

  end subroutine atualizar_ficticios_massa 

!-------------------------------------------------

  subroutine corrigir_pressao
    
	do i = 1, N
	   p(i) = p(i) + plinha(i)
	end do 

  end subroutine corrigir_pressao

!-------------------------------------------------
  
  subroutine corrigir_velocidades
     
	! Calculando fictício P = 1
	u(1) = - u(1+1) + 2.0d0*Uin

	! Calculando para volumes internos
	do i = 2, N-1
	   u(i) = u(i) - ds(i)*(plinha(i+1)-plinha(i-1))/2.0d0
	end do

	! Calculando fictício P = N
	u(N) = u(N-1) + (u(N-1)-u(N-2))

  end subroutine corrigir_velocidades

!-------------------------------------------------

  subroutine corrigir_velocidades_faces

	! Calculando ue somente para volumes internos
    do i = 2, N-2
	   ue(i) = ue(i) - de(i)*(plinha(i+1)-plinha(i))
	end do 

  end subroutine corrigir_velocidades_faces

!-------------------------------------------------
subroutine calculo_massa_especifica

    do i=2, N-1
        rop(i) = plinha(i)/(Rgases*Temperatura(i))
    end do    
    
end subroutine calculo_massa_especifica

!-------------------------------------------------

subroutine calculo_massa_especifica_nas_faces
    
    do i=2, N-1
        roe(i) = rop(i)+Beta*(rop(i+1)-rop(i))/2.0d0
    end do
    
end subroutine calculo_massa_especifica_nas_faces

subroutine corrigir_massa_especifica
    do i=2, N-1
    rop(i) = rop(i)+P(i)/(Rgases*Temperatura(i))
    end do
end subroutine corrigir_massa_especifica

subroutine corrigir_massa_especifica_faces
    do i=2, N-1
    roe(i) = roe(i) + P(i)/(Rgases*Temperatura(i))
    end do
end subroutine corrigir_massa_especifica_faces


end module coeficientes

