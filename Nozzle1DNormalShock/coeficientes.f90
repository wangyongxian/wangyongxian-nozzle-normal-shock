module coeficientes
  
use variaveis  
use dados
  
implicit none
  
contains

!-------------------------------------------------

  !subroutine lista_coeficientes_velocidade

   ! write(10,4)
   ! 4 format(//,4x,'COEFICIENTES DA VELOCIDADE',//,  &
    !         t4,'volume',t13,'X',t34,'oeste',t55,'central', &
     !                t76,'leste',t97,'fonte',/)
    !do i = 1, N
     !  write(10,2) i, X(i), awu(i), aPu(i), aeu(i), bPu(i)
    !end do

    !2 format(i4,4x,5(1pe21.11))

  !end subroutine lista_coeficientes_velocidade

!-------------------------------------------------

  !subroutine lista_coeficientes_pressao

   ! write(10,10)
    !10 format(//,4x,'COEFICIENTES DA PRESSÃO',//,  &
     !        t4,'volume',t13,'X',t34,'oeste',t55,'central', &
      !               t76,'leste',t97,'fonte',/)
    !do i = 1, N
     !  write(10,2) i, X(i), awplinha(i), aPplinha(i), aeplinha(i), bPplinha(i)
   ! end do

    !2 format(i4,4x,5(1pe21.11))

  !end subroutine lista_coeficientes_pressao

!-------------------------------------------------


  subroutine coeficientes_e_fontes_qml
	
	real*8 :: bpUDS, bpB
	
	! Fictício esquerdo P = 1
    awu(1) = 0.0d0
    aeu(1) = -1.0d0
    aPu(1) = 1.0d0
    bPu(1) = -(2.0d0*x(2)*(u(3)-u(2))/(x(3)-x(2)))
	
	! volumes internos
    do i = 2, N-1
    !aw(i) = - rom(i-1) * um(i-1) * se(i-1)
	   awu(i)  = -rop(i-1)*u(i-1)*Ae(i-1)
	   aeu(i)  = 0.0d0
	   !ap(i) = roa(i) * sp(i) * dx / dt - ( aw(i) + ae(i) )
       aPu(i)  = rop_o(i)*A(i)*dx/dt - (awu(i) + aeu(i))
	   !bf(i) = - pi * f(i) * ro(i) * (u(i)**2) * rp(i) * dx / 4   
	   !bp(i) = roa(i) * sp(i) * dx * ua(i) / dt + bf(i)   &
       !      + 0.5d0  * sp(i) * ( p(i-1) - p(i+1) )
	   bpUDS = rop_o(i)*A(i)*u_o(i)*dx/dt &
	                - 0.5d0*A(i)*(p(i+1)-p(i-1))  &
	                - 0.25d0*Pi*fator*rop(i)*(u(i)**2)*raio(i)*dx
        
       !oeste = rom(i-1) * um(i-1) * se(i-1) * ( u(i) - u(i-1) )
       !leste = rom(i) * um(i) * se(i) * ( u(i+1) - u(i) )
       !bc(i) = 0.5d0 * beta * ( oeste - leste )
       	                
	   bpB = 0.5d0*beta*(rop(i-1)*u(i-1)*Ae(i-1)*(u(i)-u(i-1)) &
	            -rop(i)*u(i)*Ae(i)*(u(i+1)-u(i)))
	   bpu(i) = bpUDS + bpB
	end do
      
    ! Fictício direito P = N
    awu(N) = -1.0d0
    aeu(N) = 0.0d0
    aPu(N) = 1.0d0
    bPu(N) = (2.0d0*(x(N)-x(N-1))*(u(N-1)-u(N-2))/(x(N-1)-x(N-2)))
    
  end subroutine coeficientes_e_fontes_qml

!-------------------------------------------------
  
  subroutine calculo_velocidades_face

    real*8  :: sigmap, sigmaE, bcP, bcE, bfP, bfE ! auxiliar
   
	! Calculando ue no passo dt+1
	ue(1) = (u(1)+u(2))/2.0d0
	
	do i = 2, N-2
	   bcP = beta*((roe(i-1)*Ae(i-1)*ue(i-1))*(u(i)-u(i-1))-(roe(i)*Ae(i)*ue(i))*(u(i+1)-u(i)))/2.0d0
	   bcE = beta*((roe(i)*Ae(i)*ue(i))*(u(i+1)-u(i))-(roe(i+1)*Ae(i+1)*ue(i+1))*(u(i+2)-u(i+1)))/2.0d0
	   ! arrumar o raio hidraulico do volume de controle
	   bfP = 0 !-Pi*fator*rop(i)*u(i)*raio(i)*dx/4.0d0
	   bfE = 0 !-Pi*fator*rop(i+1)*u(i+1)*raio(i+1)*dx/4.0d0
	   sigmap = awu(i)*u(i-1) + aeu(i)*u(i+1)
	   sigmaE = awu(i+1)*u(i) + aeu(i+1)*u(i+2)
	   !ue(i) = (-sigmap-sigmaE+bcP+bcE+bfP+bfE+((rop_o(i)*A(i)*dx+rop_o(i+1)*A(i+1)*dx)/dt)*ue_o(i) - 2.0d0*Ae(i)*(p(i+1)-p(i)))/(apu(i)+apu(i+1))
	   ue(i) = (-sigmap-sigmaE+(rop_o(i)*A(i)*dx+rop_o(i+1)*A(i+1)*dx)*ue_o(i)/(dt*2.0d0) - 2.0d0*Ae(i)*(p(i+1)-p(i)))/(apu(i)+apu(i+1))
	end do
	
	ue(N-1) = (u(N-1)+u(N))/2.0d0
	!arrumar o uin u(1)>uin
	!Fat = (u(1)*Ae(1))/(ue(N-1)*Ae(N-1))
  	!ue(N-1) = Fat*ue(N-1)
 
  end subroutine calculo_velocidades_face

!-------------------------------------------------


  subroutine coeficientes_e_fontes_energia
	
	real*8 :: bUDS, bBeta

	! Fictício esquerdo P = 1
    awT(1) = 0.0d0
    aeT(1) = 1.0d0
    apT(1) = 1.0d0
    !arrumar o u(1) q eh uin
    bpT(1) = 2.0d0*(T0-(gama-1.0d0)*(u(1)**2)/(2.0d0*gama*Rgases))
	
	! volumes internos
    do i = 2, N-1
	   awT(i)  = -cp*rop(i-1)*u(i-1)*Ae(i-1)
	   aeT(i)  = 0 
       apT(i)  = cp*A(i)*rop_o(i)*dx/dt - (awT(i)+aeT(i))
!       aw(i) = - cp(i) * rom(i-1) * um(i-1) * se(i-1)       
!       ap(i) = cp(i) * roa(i) * sp(i) * dx / dt - ( aw(i) + ae(i) ) &
!             + h(i) * C_rec(i) * sn(i)                              &
!             + epsilon_cte * sigma * sn(i) * ( T(i) ** 3 )    

!       bp(i) = cp(i) * roa(i) * sp(i) * dx * Ta(i) / dt        &
!             + sp(i) * dx * ( p(i) - pa(i) ) / dt              &
!             + sp(i) * u(i) * ( p(i+1) - p(i-1) ) * 0.5d0      &
!             + pi * f(i) * ro(i) * (u(i)**3) * rp(i) * dx / 4  &
!             + h(i) * T_wall(i) * sn(i)                        &
!             + epsilon_cte * sigma * sn(i) * ( T_wall(i) ** 4 )
       
       bUDS = cp*A(i)*rop_o(i)*T_o(i)*dx/dt &
                +A(i)*(p(i)-p_o(i))*dx/dt &
                +A(i)*u(i)*(p(i+1)-p(i-1))*0.5d0 &
                +0.25d0*Pi*fator*rop(i)*raio(i)*dx*u(i)**3
       !oeste = rom(i-1) * um(i-1) * se(i-1) * ( T(i) - T(i-1) )
       !leste = rom(i) * um(i) * se(i) * ( T(i+1) - T(i) )
       !bc(i) = 0.5d0 * beta * cp(i) * ( oeste - leste )
       	   
	   bBeta = 0.5d0*beta*cp*(rop(i-1)*u(i-1)*Ae(i-1)*(T(i)-T(i-1)) &
	                -rop(i)*u(i)*Ae(i)*(T(i+1)-T(i)))
	                
	   bpT(i)  = bUDS + bBeta
	end do
      
    ! Fictício direito P = N
    awT(N) = -1.0d0
    aeT(N) = 0.0d0
    apT(N) = 1.0d0
    bpT(N) = (2.0d0*(x(N)-x(N-1))/(x(N-1)-x(N-2)))*(T(N-1)-T(N-2))

  end subroutine coeficientes_e_fontes_energia
  
!-------------------------------------------------

  subroutine coeficientes_fontes_massa

    real*8 :: bUDS, bBeta
    ! Calculando fictício P = 1	
    awplinha(1) = 0.0d0
    aeplinha(1) = 1.0d0
    aPplinha(1) = 1.0d0
    !arrumar
    bpplinha(1) = 2.0d0*plinha(1)

    ! Calculando os volumes internos
    do i = 2, N-1
    
      ! aw(i) = - rom(i-1) * de(i-1) * se(i-1)          &
      !        - um(i-1)  * se(i-1) / ( Rg(i-1) * T(i-1) )

      ! ae(i) = - rom(i) * de(i) * se(i)

      ! ap(i) = ( sp(i) * dx / dt + um(i) * se(i) ) / ( Rg(i) * T(i) )   &
      !       + rom(i-1) * de(i-1) * se(i-1)                              &
      !       + rom(i)   * de(i)   * se(i)

      ! bp(i) = - ( ( ro(i) - roa(i) ) * sp(i) * dx / dt              &
      !       + ro(i) * um(i) * se(i) - ro(i-1) * um(i-1) * se(i-1) )
             
             
       awplinha(i) = -rop(i-1)*Ae(i-1)*de(i-1) &
                        -u(i-1)*Ae(i-1)/(rgases*T(i-1)) 
	   aeplinha(i) = -rop(i)*Ae(i)*de(i)
	   aPplinha(i) = (A(i)*dx/dt+u(i)*Ae(i))/(rgases*T(i)) &
	                    +rop(i-1)*de(i-1)*Ae(i-1) &
	                    +rop(i)*de(i)*Ae(i) 
	                    
	   bUDS = -(A(i)*dx*(rop(i)-rop_o(i))/dt &
	            +rop(i)*Ae(i)*u(i) &
	            -rop(i-1)*Ae(i-1)*u(i-1))
       !oeste = se(i-1) * ( ro(i) - ro(i-1) ) * um(i-1)
       !leste = se(i)   * ( ro(i+1) - ro(i) ) * um(i)
       !bc(i) = 0.5d0 * beta * ( oeste - leste )
	   bBeta = 0.5d0*beta*(Ae(i-1)*(rop(i)-rop(i-1))*u(i-1)-Ae(i)*(rop(i+1)-rop(i))*u(i))
	   
	   bpplinha(i) = bUDS+bBeta
    end do
  
    ! Calculando fictício P = N 
    awplinha(N) = -1.0d0
    aeplinha(N) = 0.0d0
    aPplinha(N) = 1.0d0
    bpplinha(N) = (2.0d0*(x(N)-x(N-1))/(x(N-1)-x(N-2)))*(plinha(N-1)-plinha(N-2))

  end subroutine coeficientes_fontes_massa

!-------------------------------------------------

  subroutine coeficientes_simplec

    ds(1) = 0.0d0
    ds(N) = 0.0d0
    
    do i = 2, N-1
	   ds(i) = A(i)/(apu(i)+awu(i)+aeu(i))
    end do

    ! Calculando nos contornos
	de(1) = ds(2)
	de(N-1) = ds(N-1)
	
	do i = 2, N-2
	   de(i) = 0.5d0*(ds(i) + ds(i+1))
	end do

	

  end subroutine coeficientes_simplec


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
     !arrumar
	u(1) = - u(2) + 2.0d0*u(1)
   
	! Calculando para volumes internos
	do i = 2, N-1
	   u(i) = u(i) - ds(i)*(plinha(i+1)-plinha(i-1))/2.0d0
	end do

	! Calculando fictício P = N
	u(N) = u(N-1) + (u(N-1)-u(N-2))

  end subroutine corrigir_velocidades

!-------------------------------------------------

  subroutine corrigir_velocidades_faces

    do i = 2, N-2
	   ue(i) = ue(i) - de(i)*(plinha(i+1)-plinha(i))
	end do 

  end subroutine corrigir_velocidades_faces

!-------------------------------------------------

  subroutine corrigir_massa_especifica

    do i = 1, N
	   rop(i) = rop(i) + plinha(i)/(Rgases*T(i))
	end do 

  end subroutine corrigir_massa_especifica
  
!-------------------------------------------------

subroutine calculo_massa_especifica

    do i=1, N
        rop(i) = p(i)/(Rgases*T(i))
    end do    
    
end subroutine calculo_massa_especifica

!-------------------------------------------------

subroutine calculo_massa_especifica_nas_faces
    
    do i=1, N-1
        roe(i) = rop(i) + beta*(rop(i+1)-rop(i))/2.0d0
    end do
    
end subroutine calculo_massa_especifica_nas_faces

!-------------------------------------------------

subroutine calcula_fluxo_massa
    do i=1,N
        m(i) = rop(i)*u(i)*A(i)
    end do
end subroutine calcula_fluxo_massa

!-------------------------------------------------

subroutine calcula_coeficiente_descarga
    call calcula_fluxo_massa
    do i=1,N
        Cd(i) = m(i)/Ma(i)
    end do
end subroutine calcula_coeficiente_descarga

!-------------------------------------------------

subroutine calcula_empuxo
    do i=1, N
        Empuxo(i) = rop(i)*u(i)*A(i)*u(i)
    end do
end subroutine calcula_empuxo

!-------------------------------------------------

end module coeficientes

